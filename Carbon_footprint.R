#----------------------------------------------------------------
# This script attempts to estimate energy consumption and carbon emission of R script
# - Auto CPU TDP best-effort detection
# - Free carbon intensity via UK National Grid API
# - Sampler runs in background Rscript, samples process cpu & rss, optional nvidia-smi GPU power
#
# Defaults to review/override:
#  - cpu_watts_default = 65 (W)
#  - ram_watts_per_GiB = 0.3725 (W per GiB)
#  - carbon_fallback_gpkwh = 475 (g CO2 / kWh) -global fallback if regional API fails
#
# NOTE: This script uses only free endpoints (UK National Grid) and local heuristics.
#----------------------------------------------------------------


#----------------------------------------------------------------
#### Configurations ####

## Installing packages if necessary
required_pkgs <- c("peakRAM", "ps", "processx", "digest", "jsonlite", "httr")
for (p in required_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cloud.r-project.org")
}

## Loading packages
library(peakRAM)
library(ps)
library(processx)
library(digest)
library(jsonlite)
library(httr)
library(tidyverse)

#----------------------------------------------------------------
#### Creating a table with default CPU TDP values (to adapt) ####

## Small CPU TDP lookup table
.cpu_tdp_table <- list("Intel(R) Core(TM) i7-9700K" = 95,
                       "Intel(R) Core(TM) i7-8700K" = 95,
                       "Intel(R) Core(TM) i9-9900K" = 95,
                       "Intel(R) Core(TM) i7-10700K" = 125,
                       "Intel(R) Core(TM) i7-11700K" = 125,
                       "Intel(R) Xeon" = 85,
                       "AMD Ryzen 7 3700X" = 65,
                       "AMD Ryzen 9 3900X" = 105,
                       "AMD EPYC" = 180)

#----------------------------------------------------------------
#### Creating functions ####

## CPU TDP detection
detect_cpu_tdp_best_effort <- function() {
  model <- NA_character_
  
  # Linux / Unix: try lscpu then /proc/cpuinfo
  if (.Platform$OS.type == "unix") {
    try({
      out <- tryCatch(system("lscpu | grep 'Model name' | sed -E 's/Model name:\\s*//'", intern = TRUE), error = function(e) NA_character_)
      if (!is.na(out) && length(out) > 0 && nzchar(out[1])) model <- out[1]
    }, silent = TRUE)
    if (is.na(model) || !nzchar(model)) {
      try({
        info <- readLines("/proc/cpuinfo", warn = FALSE)
        mn <- grep("^model name\\s*:", info, value = TRUE)
        if (length(mn) >= 1) model <- sub("^model name\\s*:\\s*", "", mn[1])
      }, silent = TRUE)
    }
  }
  
  # Windows: using wmic
  if (.Platform$OS.type == "windows" && (is.na(model) || !nzchar(model))) {
    try({
      out <- suppressWarnings(system2("wmic", args = c("cpu", "get", "Name"), stdout = TRUE, stderr = FALSE))
      out <- out[nzchar(out)]
      if (length(out) >= 2) model <- trimws(out[2])
    }, silent = TRUE)
  }
  
  # Fallback to Sys.info
  if (is.na(model) || !nzchar(model)) {
    model <- paste(Sys.info()[c("sysname", "release", "machine")], collapse = " ")
  }
  
  tdp <- NA_real_; source <- "none"
  
  for (k in names(.cpu_tdp_table)) {
    if (grepl(k, model, ignore.case = TRUE)) {
      tdp <- .cpu_tdp_table[[k]]; source <- paste0("table_match:", k)
      break
    }
  }
  
  # Generic heuristics
  if (is.na(tdp)) {
    for (pat in c("Xeon", "EPYC", "Ryzen", "Intel", "Core")) {
      if (grepl(pat, model, ignore.case = TRUE)) {
        tdp <- switch(tolower(pat),
                      "xeon" = 85, "epyc" = 180, "ryzen" = 95,
                      "intel" = 65, "core" = 65, NA)
        source <- paste0("heuristic:", pat)
        break
      }
    }
  }
  list(model = model, tdp = tdp, source = source)
}

## Carbon intensity fetcher (free API): returns numeric gCO2/kWh

get_carbon_intensity <- function(fallback_g_per_kwh = 475, verbose = TRUE) {
  url <- "https://api.carbonintensity.org.uk/intensity"
  
  # trying GET with timeout
  res <- tryCatch(httr::GET(url, httr::timeout(5)), error = function(e) e)
  if (inherits(res, "error") || !httr::status_code(res) %in% c(200)) {
    if (verbose) warning("GB CI API request failed or returned non-200 --> fallback used: ", fallback_g_per_kwh, " gCO2/kWh")
    return(as.numeric(fallback_g_per_kwh))
  }
  
  # parsing safely as text then json
  txt <- tryCatch(httr::content(res, as = "text", encoding = "UTF-8"), error = function(e) NA_character_)
  if (is.na(txt) || !nzchar(txt)) {
    if (verbose) warning("GB CI API returned empty body --> fallback used: ", fallback_g_per_kwh)
    return(as.numeric(fallback_g_per_kwh))
  }
  parsed <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE), error = function(e) e)
  if (inherits(parsed, "error") || !is.list(parsed)) {
    if (verbose) warning("GB CI JSON parse failed or wrong type --> fallback used: ", fallback_g_per_kwh)
    return(as.numeric(fallback_g_per_kwh))
  }

  if (!"data" %in% names(parsed) || length(parsed$data) < 1 || !is.list(parsed$data[[1]])) {
    if (verbose) warning("GB CI payload unexpected structure --> fallback used: ", fallback_g_per_kwh)
    return(as.numeric(fallback_g_per_kwh))
  }
  first <- parsed$data[[1]]
  if (!is.list(first) || !"intensity" %in% names(first) || !is.list(first$intensity)) {
    if (verbose) warning("GB CI intensity missing or malformed --> fallback used: ", fallback_g_per_kwh)
    return(as.numeric(fallback_g_per_kwh))
  }
  intensity <- first$intensity

  if (!is.null(intensity$actual) && is.numeric(intensity$actual)) return(as.numeric(intensity$actual))
  if (!is.null(intensity$forecast) && is.numeric(intensity$forecast)) return(as.numeric(intensity$forecast))
  if (verbose) warning("GB CI actual/forecast missing --> fallback used: ", fallback_g_per_kwh)
  return(as.numeric(fallback_g_per_kwh))
}

## Checking if nvidia-smi available
has_nvidia_smi <- function() nzchar(Sys.which("nvidia-smi"))

## Sampler script writer
write_sampler_script <- function(path) {
  code <- '
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("sampler needs 5 args")
sentinel <- args[1]; out_csv <- args[2]; target_pid <- as.integer(args[3]); interval_ms <- as.integer(args[4]); use_gpu <- as.logical(args[5])
if (!requireNamespace("ps", quietly = TRUE)) install.packages("ps", repos="https://cloud.r-project.org")
library(ps)
writeLines("time_iso,cpu_user,cpu_system,cpu_total,rss_MiB,gpu_power_W", con = out_csv)
while (file.exists(sentinel)) {
  t <- format(Sys.time(), "%Y-%m-%dT%H:%M:%OS3")
  cpu_user <- cpu_system <- cpu_total <- rss_MiB <- gpu_power <- NA
  try({
    p <- ps_handle(target_pid)
    ct <- ps_cpu_times(p)
    cpu_user <- as.numeric(ct["user"]); cpu_system <- as.numeric(ct["system"])
    cpu_total <- sum(cpu_user, cpu_system, na.rm = TRUE)
    rss_MiB <- as.numeric(ps_memory_info(p)["rss"]) / (1024^2)
  }, silent = TRUE)
  if (use_gpu) {
    try({
      out <- system2("nvidia-smi", args = c("--query-gpu=power.draw","--format=csv,noheader,nounits"), stdout = TRUE, stderr = FALSE)
      if (length(out) >= 1) gpu_power <- as.numeric(out[1])
    }, silent = TRUE)
  }
  line <- paste(t, ifelse(is.na(cpu_user),"",cpu_user), ifelse(is.na(cpu_system),"",cpu_system), ifelse(is.na(cpu_total),"",cpu_total), ifelse(is.na(rss_MiB),"",rss_MiB), ifelse(is.na(gpu_power),"",gpu_power), sep = ",")
  write(line, file = out_csv, append = TRUE)
  Sys.sleep(interval_ms / 1000)
}
'
writeLines(code, con = path)
invisible(path)
}

## Main function to track carbon emission
monitor_emissions <- function(code,
    auto_detect_cpu = TRUE,
    cpu_watts = NULL,
    cpu_watts_default = 65,
    ram_watts_per_GiB = 0.3725,
    gpu_watts_fallback = NULL,
    carbon_region = "GB",
    carbon_fallback_gpkwh = 475,
    use_gpu_sampling = TRUE,
    sample_interval_ms = 200L,
    monitor_dir = tempdir(),
    verbose = TRUE) {
  
  # Preparing
  code_expr <- substitute(code)
  timestamp <- Sys.time()
  code_hash <- digest::digest(code_expr)
  
  # CPU watt selection
  cpu_detect <- NULL
  if (auto_detect_cpu && is.null(cpu_watts)) {
    cpu_detect <- detect_cpu_tdp_best_effort()
    cpu_watts_est <- if (!is.na(cpu_detect$tdp)) cpu_detect$tdp else cpu_watts_default
    cpu_watts_source <- if (!is.na(cpu_detect$tdp)) cpu_detect$source else "fallback_default"
  } else if (!is.null(cpu_watts)) {
    cpu_watts_est <- as.numeric(cpu_watts)
    cpu_watts_source <- "user_provided"
  } else {
    cpu_watts_est <- cpu_watts_default
    cpu_watts_source <- "fallback_default"
  }
  
  # Carbon intensity (free region only currently: GB)
  carbon_gpkwh <- carbon_fallback_gpkwh
  if (toupper(carbon_region) %in% c("GB", "UK")) {
    carbon_gpkwh <- tryCatch(get_carbon_intensity(fallback_g_per_kwh = carbon_fallback_gpkwh, verbose = verbose),
                             error = function(e) { if (verbose) warning("Carbon fetch error: ", conditionMessage(e)); carbon_fallback_gpkwh })
  } else {
    if (verbose) message("Region not supported for free endpoint; using fallback gCO2/kWh = ", carbon_fallback_gpkwh)
    carbon_gpkwh <- carbon_fallback_gpkwh
  }
  
  # Sampler setup
  sampler_script <- tempfile(fileext = ".R")
  sampler_csv <- tempfile(fileext = ".csv")
  sentinel_file <- tempfile()
  write_sampler_script(sampler_script)
  file.create(sentinel_file)
  rscript_bin <- file.path(R.home("bin"), "Rscript")
  use_gpu_flag <- use_gpu_sampling && has_nvidia_smi()
  sampler_proc <- processx::process$new(rscript_bin,
                                        args = c(sampler_script, sentinel_file, sampler_csv, as.character(Sys.getpid()), as.character(sample_interval_ms), as.character(use_gpu_flag)),
                                        stdout = "|", stderr = "|")
  if (verbose) message("Sampler started (PID): ", sampler_proc$get_pid())
  
  # Running target code with peakRAM
  if (verbose) message("Running target code (peakRAM)...")
  peak_res <- peakRAM::peakRAM(eval(code_expr))
  if (!is.data.frame(peak_res) || nrow(peak_res) < 1) stop("peakRAM failed to produce results")
  elapsed_sec_peak <- as.numeric(peak_res$Elapsed_Time_sec[1])
  peak_RAM_MiB_peakRAM <- as.numeric(peak_res$Peak_RAM_Used_MiB[1])
  
  # Stop sampler
  if (file.exists(sentinel_file)) file.remove(sentinel_file)
  # wait for process to stop and then clean up
  sampler_proc$wait(2000)
  if (sampler_proc$is_alive()) { sampler_proc$kill(); sampler_proc$wait(1000) }
  
  # Reading sampler csv safely
  sampler_df <- data.frame()
  if (file.exists(sampler_csv)) {
    raw <- readLines(sampler_csv, warn = FALSE)
    if (length(raw) > 1) {
      tmp <- utils::read.csv(text = raw, stringsAsFactors = FALSE)
      
      # coercing numeric columns
      suppressWarnings({
        tmp$cpu_user <- as.numeric(tmp$cpu_user)
        tmp$cpu_system <- as.numeric(tmp$cpu_system)
        tmp$cpu_total <- as.numeric(tmp$cpu_total)
        tmp$rss_MiB <- as.numeric(tmp$rss_MiB)
        tmp$gpu_power_W <- as.numeric(tmp$gpu_power_W)
      })
      sampler_df <- tmp
    }
  }
  
  # computing delta cpu sec: prefering sampler delta if available, else peakRAM/ps diff
  delta_proc_cpu_sec <- NA_real_
  if (nrow(sampler_df) >= 2 && !all(is.na(sampler_df$cpu_total))) {
    cpu_first <- sampler_df$cpu_total[1]; cpu_last <- sampler_df$cpu_total[nrow(sampler_df)]
    delta_sampler_cpu <- cpu_last - cpu_first
    if (!is.na(delta_sampler_cpu) && delta_sampler_cpu >= 0) delta_proc_cpu_sec <- delta_sampler_cpu
  }
  
  # fallback: trying ps times around evaluation (less precise)
  if (is.na(delta_proc_cpu_sec)) {
    try({
      ph <- ps_handle(Sys.getpid())
      ct_before <- ps_cpu_times(ph)

      ct_after <- ct_before
      delta_proc_cpu_sec <- as.numeric(ct_after["user"] + ct_after["system"])

      if (!is.finite(delta_proc_cpu_sec)) delta_proc_cpu_sec <- NA_real_
    }, silent = TRUE)
  }
  
  # elapsed time: prefer peakRAM value; if zero, deriving from sampler timestamps
  elapsed_sec <- elapsed_sec_peak
  if (is.na(elapsed_sec) || elapsed_sec <= 0) {
    if (nrow(sampler_df) >= 2) {
      t0 <- try(as.POSIXct(sampler_df$time_iso[1], format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"), silent = TRUE)
      t1 <- try(as.POSIXct(sampler_df$time_iso[nrow(sampler_df)], format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC"), silent = TRUE)
      if (!inherits(t0, "try-error") && !inherits(t1, "try-error") && !is.na(t0) && !is.na(t1)) {
        elapsed_sec <- as.numeric(difftime(t1, t0, units = "secs"))
      }
    }
  }
  if (is.na(elapsed_sec) || elapsed_sec <= 0) elapsed_sec <- elapsed_sec_peak
  
  # peak RAM: prefering sampler max if present
  peak_RAM_MiB_sampler <- if (nrow(sampler_df) >= 1 && any(!is.na(sampler_df$rss_MiB))) max(sampler_df$rss_MiB, na.rm = TRUE) else NA_real_
  peak_RAM_MiB <- if (!is.na(peak_RAM_MiB_sampler) && peak_RAM_MiB_sampler > 0) peak_RAM_MiB_sampler else peak_RAM_MiB_peakRAM
  if (is.na(peak_RAM_MiB)) peak_RAM_MiB <- 0
  
  # cores & avg CPU utilization fraction
  n_cores <- parallel::detectCores(logical = TRUE)
  avg_cpu_util_frac <- NA_real_
  if (!is.na(delta_proc_cpu_sec) && elapsed_sec > 0 && n_cores > 0) {
    avg_cpu_util_frac <- min(max(delta_proc_cpu_sec / (elapsed_sec * n_cores), 0), 1)
  }
  
  # GPU mean power if available
  gpu_power_mean_W <- NA_real_
  if (nrow(sampler_df) >= 1 && any(!is.na(sampler_df$gpu_power_W))) {
    gpu_power_mean_W <- mean(sampler_df$gpu_power_W, na.rm = TRUE)
  } else if (has_nvidia_smi() && use_gpu_sampling) {
    out <- tryCatch(system2("nvidia-smi", args = c("--query-gpu=power.draw", "--format=csv,noheader,nounits"), stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
    if (!is.na(out[1])) gpu_power_mean_W <- suppressWarnings(as.numeric(out[1]))
  } else if (!is.null(gpu_watts_fallback)) {
    gpu_power_mean_W <- as.numeric(gpu_watts_fallback)
  }
  
  # Energies (Joules)
  cpu_energy_J <- if (!is.na(avg_cpu_util_frac)) cpu_watts_est * avg_cpu_util_frac * elapsed_sec else NA_real_
  ram_GiB <- peak_RAM_MiB / 1024
  ram_power_W <- ram_GiB * ram_watts_per_GiB
  ram_energy_J <- ram_power_W * elapsed_sec
  gpu_energy_J <- if (!is.na(gpu_power_mean_W)) gpu_power_mean_W * elapsed_sec else NA_real_
  
  total_energy_J <- sum(c(ifelse(is.na(cpu_energy_J), 0, cpu_energy_J),
                          ifelse(is.na(ram_energy_J), 0, ram_energy_J),
                          ifelse(is.na(gpu_energy_J), 0, gpu_energy_J)), na.rm = TRUE)
  total_energy_kWh <- total_energy_J / 3.6e6
  carbon_g <- total_energy_kWh * carbon_gpkwh
  
  # Building result object
  result <- list(
    metadata = list(
      timestamp = as.character(timestamp),
      code_label = if (!is.null(code_expr)) deparse(code_expr)[1] else "",
      code_hash = code_hash,
      cpu_detect = cpu_detect,
      cpu_watts_used = cpu_watts_est,
      cpu_watts_source = cpu_watts_source,
      ram_watts_per_GiB = ram_watts_per_GiB,
      gpu_power_mean_W = gpu_power_mean_W,
      carbon_gpkwh = carbon_gpkwh
    ),
    timing = list(elapsed_sec = elapsed_sec, delta_proc_cpu_sec = delta_proc_cpu_sec, avg_cpu_util_frac = avg_cpu_util_frac, n_cores = n_cores),
    memory = list(peak_RAM_MiB = peak_RAM_MiB, peak_RAM_GiB = peak_RAM_MiB / 1024),
    energy = list(cpu_energy_J = cpu_energy_J, ram_energy_J = ram_energy_J, gpu_energy_J = gpu_energy_J, total_energy_J = total_energy_J, total_energy_kWh = total_energy_kWh, carbon_g = carbon_g),
    sampler_df = sampler_df,
    sampler_files = list(sampler_script = sampler_script, sampler_csv = sampler_csv)
  )
  class(result) <- "monitor_emissions_result"
  
  # Print summary
  if (verbose) {
    message("---- monitor_emissions summary ----")
    message("CPU model detected: ", result$metadata$cpu_detect$model)
    message("CPU watts used: ", result$metadata$cpu_watts_used, " (source: ", result$metadata$cpu_watts_source, ")")
    message("RAM watts per GiB: ", result$metadata$ram_watts_per_GiB)
    message("Carbon intensity (gCO2/kWh) used: ", result$metadata$carbon_gpkwh, " (region: ", carbon_region, ")")
    message(sprintf("Elapsed (s): %.4f | delta_cpu_sec: %s | avg_cpu_util_frac: %s", result$timing$elapsed_sec,
                    ifelse(is.na(result$timing$delta_proc_cpu_sec), "NA", format(result$timing$delta_proc_cpu_sec, digits = 4)),
                    ifelse(is.na(result$timing$avg_cpu_util_frac), "NA", format(result$timing$avg_cpu_util_frac, digits = 4))))
    message(sprintf("Peak RAM (MiB): %.1f | CPU energy (J): %s | RAM energy (J): %.2f | GPU mean (W): %s",
                    result$memory$peak_RAM_MiB,
                    ifelse(is.na(result$energy$cpu_energy_J), "NA", format(result$energy$cpu_energy_J, digits = 6)),
                    ifelse(is.na(result$energy$ram_energy_J), 0, result$energy$ram_energy_J),
                    ifelse(is.na(result$metadata$gpu_power_mean_W), "NA", result$metadata$gpu_power_mean_W)))
    message(sprintf("Total energy (kWh): %.9f | Estimated CO2 (g): %.6f", result$energy$total_energy_kWh, result$energy$carbon_g))
    message("-------------------------------------------")
  }
  
  invisible(result)
}

## Function to create a clean result table
as_emissions_long <- function(res) {

  data_list <- list(
    "Timestamp"                = as.character(res$metadata$timestamp),
    "Code hash"                = as.character(res$metadata$code_hash),
    "Carbon intensity"         = paste0(res$metadata$carbon_gpkwh, " gCO2/kWh"),
    
    # CPU
    "CPU model"                = as.character(res$metadata$cpu_detect$model %||% "NA"),
    "CPU power used"           = paste0(res$metadata$cpu_watts_used %||% "NA", " W"),
    "CPU power source"         = as.character(res$metadata$cpu_watts_source %||% "NA"),
    
    # RAM
    "RAM power per GiB"        = paste0(signif(res$metadata$ram_watts_per_GiB %||% NA, 3), " W/GiB"),
    
    # GPU
    "Mean GPU power"           = paste0(res$metadata$gpu_power_mean_W %||% "NA", " W"),
    
    # Timing
    "Elapsed time"             = paste0(signif(res$timing$elapsed_sec %||% NA,3), " s"),
    "Delta CPU time"           = paste0(signif(res$timing$delta_proc_cpu_sec %||% NA,3), " s"),
    "Average CPU utilization"  = paste0(round((res$timing$avg_cpu_util_frac %||% NA)*100,2), " %"),
    "Number of cores"          = res$timing$n_cores %||% "NA",
    
    # Memory
    "Peak RAM"                 = paste0(round(res$memory$peak_RAM_MiB %||% NA,1), " MiB"),
    "Peak RAM (GiB)"           = paste0(round(res$memory$peak_RAM_GiB %||% NA,3), " GiB"),
    
    # Energy
    "CPU energy"               = paste0(round(res$energy$cpu_energy_J %||% NA,3), " J"),
    "RAM energy"               = paste0(round(res$energy$ram_energy_J %||% NA,3), " J"),
    "GPU energy"               = paste0(round(res$energy$gpu_energy_J %||% NA,3), " J"),
    "Total energy"             = paste0(round(res$energy$total_energy_J %||% NA,3), " J"),
    "Total energy (kWh)"       = paste0(round(res$energy$total_energy_kWh %||% NA,9), " kWh"),
    "Estimated carbon"         = paste0(round(res$energy$carbon_g %||% NA,9), " g CO2"))
  
  tibble(
    Metric = names(data_list),
    Value  = unlist(data_list, use.names = FALSE))
}

#----------------------------------------------------------------
#### Examples of usage ####

## Example 1: auto-detecting CPU and using UK carbon intensity
res1 <- monitor_emissions(
  quote({
    Sys.sleep(1)
    x <- rnorm(2e6)
    sum(x)
  }),
  verbose = TRUE
)
as_emissions_long(res1)

## Example 2: providing known CPU watts (overriding TDP detection)
res2 <- monitor_emissions(
  quote({
    mat <- matrix(rnorm(1e7), ncol = 1000)
    mean(mat)
  }),
  cpu_watts = 95,
  verbose = TRUE)
as_emissions_long(res2)

## Example 3: disabling GPU sampling (for machines without nvidia-smi)
res3 <- monitor_emissions(
  quote({
    Sys.sleep(0.5)
    sum(rnorm(5e6))
  }),
  use_gpu_sampling = FALSE,
  verbose = TRUE)

as_emissions_long(res3)

## Example 4: auto-detecting CPU and using carbon intensity specified by user

res4 <- monitor_emissions(
  code = quote({
    x <- rnorm(1e6)
    mean(x)
  }),
  carbon_region="FR",
  carbon_fallback_gpkwh = 60,
  verbose = TRUE)

as_emissions_long(res4)
