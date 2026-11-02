#-------------------------------------------------------------------------------
## Reproducible and generalizable Cox regression script with time-varying covariates (counting-process format: one row per interval (start, stop])
#
# covers:
# - data import:
#   ¤ one table of baseline data (one row per individual)
#   ¤ one table of repeated measurements (one row per individual and measurement time)
# - construction of the counting-process dataset: value carried forward, optional lag, optional expiry of old values
# - checks of the constructed dataset, description of the measurement process
# - Cox models (cause-specific, one per cause when there are competing events), cluster-robust standard errors:
#   ¤ current value (main)
#   ¤ baseline value only
#   ¤ baseline value + change
#   ¤ lagged value
#   ¤ old values excluded
# - assumptions: independence, proportional hazards, linearity, multicollinearity (overall and by period), effect modification
# - remediation of non-proportional hazards (covariate x log(time) interaction)
# - visualization: forest plots, comparison of the models, trajectories of the covariates
#
# Note:
# - parameters to adapt are in "User settings"
# - the example uses survival::pbc (baseline) and survival::pbcseq (repeated measurements)
#-------------------------------------------------------------------------------

#----------------------------------------------------------------
#### Configurations ####

## Disabling memory torture
gctorture(FALSE)

## Installing and loading packages
pack_needed<-c("data.table","tidyverse","survival","survminer","conflicted","broom","here","car","qs2","patchwork","tidycmprsk","ggsurvfit")

is_installed<-pack_needed %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed==FALSE)){
  install.packages(pack_needed[!is_installed],repos="http://cran.us.r-project.org")
}
invisible(lapply(pack_needed, library, character.only=TRUE))

## Preventing package conflicts
conflict_prefer("select","dplyr")
conflict_prefer("filter","dplyr")
conflict_prefer("slice","dplyr")
conflict_prefer("summarize","dplyr")
conflict_prefer("lag","dplyr")
conflict_prefer("first","dplyr")

## Setting the working directory and the output folder
out_dir<-here::here("Cox regression with time-varying covariates")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
setwd(out_dir)
output_dir<-"Results"
if(!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
safe_path<-function(filename) file.path(output_dir, filename)

## Setting seed
seed<-123

#----------------------------------------------------------------
#### Creating functions ####

#- - - - - -
## Formatting numbers for display

test_format<-function(x){
  x<-as.numeric(x)
  sign_x<-if_else(x<0,"neg","pos")
  x<-abs(x)
  x_raw<-x
  if(is.na(x)|is.infinite(x)) x_raw<-0
  if(x_raw>=100){
    virg_pos<-str_locate(as.character(x_raw),"[.]")[1]
    if(!is.na(virg_pos)&as.numeric(substr(x_raw,virg_pos+1,virg_pos+1))>=5){
      x<-x+1
      x<-as.numeric(substr(x,1,virg_pos-1))
    }
  }
  if(is.na(x)|is.infinite(x)){
    x<-""
  }else{
    if(x<0.05|x>=10000){
      x<-format(signif(x,3), scientific=TRUE)
    }else{
      x_save<-x
      x<-signif(x,3)
      if(nchar(x)==6) x<-as.numeric(substr(x,1,5))
      if(nchar(x)==5){
        if(as.numeric(substr(x,5,5))>=5){
          x<-x+0.01
          x<-substr(x,1,4)
        }else{
          x<-substr(x,1,4)
        }
      }else{
        if(x>=1000) x<-as.character(signif(x_save,4)) else x<-as.character(x)
      }
    }
  }
  if(sign_x=="neg"&x_raw!=0) x<-paste("-",x,sep="")
  if(x=="0e+00") x<-"0"
  return(x)
}

vf<-function(x) vapply(x, test_format, character(1)) # test_format applied number by number

#- - - - - -
## Custom ggplot theme

theme_Gaia<-function(){
  theme_bw() +
    theme(strip.text=element_text(size=12, colour="black", face="bold"),
          strip.background=element_rect(fill="#CAE1FF",colour="black"),
          axis.text=element_text(size=12, color="black"),
          axis.title=element_text(size=13, face="bold", color="black"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=13, face="bold"),
          legend.key.size=unit(0.5, "cm"),
          axis.line=element_line(color="black", linewidth=0.1))
}

#- - - - - -
## Saving figures (any ggplot2 object, or a list of ggplots) with automatic size
#
# - ext: one or several formats, e.g. c("png","pdf","tiff"); one file per format (name_<date>.<ext>)
# - compression: TIFF compression (e.g., "lzw"); only used for tiff files
# - discrete axes: the size is computed from the number and length of the labels and from the text written inside the tiles (geom_text/geom_label)
# - continuous axes: panel_width x panel_height (cm)
# - facets with free space keep their relative sizes; fixed-ratio coordinates (coord_fixed, coord_sf) keep their aspect ratio
# - one plot: the panel size is set directly in the plot layout and the final size is measured on this layout (panels + axes + strips + legend), widened if the title/subtitle/caption is wider
# - several plots: one single legend, placed at legend_position ("top" = above all plots)
# - a patchwork object cannot be measured: width/height are used (default 25 x 18 cm)
# - width, height: optional (cm); override the automatic size
# - ...: other arguments passed to ggplot2::ggsave

save_plot<-function(plot, name, ncol=1, dpi=300, ext="png", panel_width=9, panel_height=7, legend_position=NULL, margin=0.5,
                    width=NULL, height=NULL, compression="lzw", ...){
  plots<-if(inherits(plot, "ggplot")) list(plot) else plot
  ncol<-min(ncol, length(plots))
  nrow_p<-ceiling(length(plots)/ncol)
  ext<-unique(tolower(ext))
  files<-paste0(name, "_", Sys.Date(), ".", ext) # one file per format
  
  # null device: needed by grid to measure text, without creating a file
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add=TRUE)
  
  pt_cm<-0.03528 # 1 point in cm
  
  # width (cm) of the widest text of a vector
  text_cm<-function(x, fontsize){
    x<-as.character(x[!is.na(x)])
    if(length(x)==0) return(0)
    grid::convertWidth(grid::grobWidth(grid::textGrob(x, gp=grid::gpar(fontsize=fontsize))), "cm", valueOnly=TRUE)
  }
  # value of a theme element of the plot (e.g., font size of the axis text)
  theme_el<-function(p, el, field, default){
    th<-tryCatch(ggplot2::theme_get() + p$theme, error=function(e) ggplot2::theme_get())
    v<-tryCatch(ggplot2::calc_element(el, th)[[field]], error=function(e) NULL)
    if(is.null(v) || length(v)==0 || is.na(v[1])) default else v
  }
  
  # panel area of one plot: returns its size (cm) and its layout (gtable) with the panels set to this size
  panel_dims<-function(p){
    b<-ggplot2::ggplot_build(p)
    g<-ggplot2::ggplotGrob(p)
    lay<-b$layout$layout
    pan<-g$layout[grepl("^panel", g$layout$name), ]
    cols<-sort(unique(pan$l)) # columns and rows of the layout occupied by panels
    rows<-sort(unique(pan$t))
    
    # with coord_flip, the x scale is drawn vertically
    flipped<-inherits(p$coordinates, "CoordFlip")
    h_ax<-if(flipped) "y" else "x"
    v_ax<-if(flipped) "x" else "y"
    sc_h<-b$layout[[paste0("panel_scales_", h_ax)]]
    sc_v<-b$layout[[paste0("panel_scales_", v_ax)]]
    idx_h<-lay[[paste0("SCALE_", toupper(h_ax))]]
    idx_v<-lay[[paste0("SCALE_", toupper(v_ax))]]
    
    # text layers (all labels need room, including repelled ones)
    cls<-sapply(p$layers, function(l) class(l$geom)[1])
    is_txt<-grepl("Text|Label", cls)
    txt<-if(any(is_txt)) as.character(unlist(lapply(b$data[which(is_txt)], function(d) d$label))) else character(0)
    has_txt<-length(txt)>0
    
    # text written inside the tiles (plain geom_text / geom_label, not repelled): the slots of discrete axes must be large enough to contain it
    is_tile<-cls %in% c("GeomText","GeomLabel")
    tile_txt<-if(any(is_tile)) as.character(unlist(lapply(b$data[which(is_tile)], function(d) d$label))) else character(0)
    tile_txt<-tile_txt[!is.na(tile_txt) & tile_txt!=""]
    tile_fs<-if(any(is_tile)) max(unlist(lapply(b$data[which(is_tile)], function(d) if("size" %in% names(d)) d$size else 3.88)), na.rm=TRUE) * ggplot2::.pt else 11
    tile_w<-if(length(tile_txt)>0) 1.15 * text_cm(unlist(strsplit(tile_txt, "\n", fixed=TRUE)), tile_fs) + 0.6 else 0 # widest line (cm): +15% and padding
    
    # number of labels stacked in the same slot of the discrete axis (e.g., dodged labels of several groups): the slot must be tall enough for all of them
    n_stack<-1
    if(any(is_tile)){
      v_discrete<-tryCatch(sc_v[[1]]$is_discrete(), error=function(e) FALSE)
      ax<-if(flipped) "x" else "y"
      n_stack<-max(1, unlist(lapply(b$data[which(is_tile)], function(dl){
        if(!v_discrete || !all(c(ax, "PANEL") %in% names(dl))) return(1)
        max(table(dl$PANEL, round(as.numeric(dl[[ax]]))))
      })))
    }
    tile_h<-if(length(tile_txt)>0) n_stack * (max(stringr::str_count(tile_txt, "\n")) + 1) * tile_fs * pt_cm * 1.4 + 0.5 else 0 # height of the stacked labels (cm)
    
    fs_x<-theme_el(p, "axis.text.x", "size", 11)
    fs_y<-theme_el(p, "axis.text.y", "size", 11)
    ang_x<-theme_el(p, "axis.text.x", "angle", 0)
    lab_of<-function(s) tryCatch(as.character(s$get_labels()), error=function(e) as.character(s$get_limits()))
    
    # horizontal size of one panel: discrete = one slot per label (widest of: label, text inside the tile); continuous = panel_width
    size_h<-function(s){
      if(s$is_discrete()){
        lab<-lab_of(s)
        fs<-if(flipped) fs_y else fs_x
        per<-if(!flipped && ang_x>=45) max(0.9, fs * pt_cm * 1.8) else max(0.9, text_cm(lab, fs) + 0.5)
        max(4, length(lab) * max(per, tile_w))
      }else{
        max(panel_width, if(has_txt) 3 * text_cm(txt, 11) else 0)
      }
    }
    # vertical size of one panel: discrete = one slot per label (height of its lines, or of the text inside the tile); continuous = panel_height
    size_v<-function(s){
      if(s$is_discrete()){
        lab<-lab_of(s)
        fs<-if(flipped) fs_x else fs_y
        n_l<-stringr::str_count(lab, "\n") + 1
        per<-pmax(if(has_txt) 1.1 else 0.8, n_l * fs * pt_cm * 1.25 + 0.35, tile_h)
        max(3, sum(per))
      }else{
        panel_height
      }
    }
    w<-sapply(seq_len(max(lay$COL)), function(i) max(sapply(unique(idx_h[lay$COL==i]), function(k) size_h(sc_h[[k]]))))
    h<-sapply(seq_len(max(lay$ROW)), function(i) max(sapply(unique(idx_v[lay$ROW==i]), function(k) size_v(sc_v[[k]]))))
    if(length(cols)!=length(w)) w<-rep(max(w), length(cols))
    if(length(rows)!=length(h)) h<-rep(max(h), length(rows))
    
    # facets with free space (unequal null units) keep their relative sizes, otherwise all panels are equal
    keep_prop<-function(sizes, units){
      nul<-tryCatch(as.numeric(units), error=function(e) rep(1, length(units)))
      if(length(nul)>1 && diff(range(nul))>1e-8) nul/max(nul) * max(sizes) else rep(max(sizes), length(sizes))
    }
    w<-keep_prop(w, g$widths[cols])
    h<-keep_prop(h, g$heights[rows])
    
    # fixed aspect ratio (coord_fixed, coord_sf...): the height follows the width
    free<-tryCatch(p$coordinates$is_free(), error=function(e) TRUE)
    if(!free){
      asp<-tryCatch(sum(as.numeric(g$heights[rows]))/sum(as.numeric(g$widths[cols])), error=function(e) 1)
      w<-rep(panel_width/length(w), length(w))
      h<-rep(sum(w) * asp/length(h), length(h))
      if(sum(h)>2 * panel_height){
        h<-h * 2 * panel_height/sum(h)
        w<-w * (sum(h)/asp)/sum(w)
      }
    }
    
    # panels set to the computed size in the layout
    g$widths[cols]<-grid::unit(w, "cm")
    g$heights[rows]<-grid::unit(h, "cm")
    list(w=sum(grid::convertWidth(g$widths[min(cols):max(cols)], "cm", valueOnly=TRUE)),
         h=sum(grid::convertHeight(g$heights[min(rows):max(rows)], "cm", valueOnly=TRUE)),
         grob=g)
  }
  
  if(any(sapply(plots, inherits, what="patchwork"))){
    # already assembled patchwork: cannot be measured
    message("A patchwork object cannot be measured: width/height used (default 25 x 18 cm); use width= and height= (cm) to change them")
    fig<-if(length(plots)==1) plots[[1]] else patchwork::wrap_plots(plots, ncol=ncol, guides="collect")
    if(!is.null(legend_position)) fig<-fig & ggplot2::theme(legend.position=legend_position)
    if(is.null(width)) width<-25
    if(is.null(height)) height<-18
    
  }else if(length(plots)==1){
    # one plot: panels set directly in the layout of the plot, size measured on this layout (panels + axes + strips + legend)
    p<-plots[[1]]
    if(!is.null(legend_position)) p<-p + ggplot2::theme(legend.position=legend_position)
    d<-panel_dims(p)
    message("Panel area computed: ", round(d$w, 1), " cm wide x ", round(d$h, 1), " cm high")
    fig<-d$grob
    w_layout<-sum(grid::convertWidth(fig$widths, "cm", valueOnly=TRUE))
    h_layout<-sum(grid::convertHeight(fig$heights, "cm", valueOnly=TRUE))
    # title, subtitle and caption are not part of the column widths: the figure must be at least as wide as the longest one
    w_text<-max(c(0, sapply(c("title","subtitle","caption"), function(l){
      s<-tryCatch(p$labels[[l]], error=function(e) NULL)
      if(is.null(s) || length(s)==0 || is.na(s[1])) return(0)
      1.1 * text_cm(strsplit(as.character(s), "\n", fixed=TRUE)[[1]], theme_el(p, paste0("plot.", l), "size", 11))
    })))
    if(is.null(width)) width<-max(w_layout * 1.02, w_text) + margin
    if(is.null(height)) height<-h_layout * 1.02 + margin
    
  }else{
    # several plots: panel area of each column and each row of plots = largest of its plots, assembled with patchwork
    sizes<-lapply(plots, panel_dims)
    idx<-matrix(seq_len(nrow_p * ncol), nrow=nrow_p, ncol=ncol, byrow=TRUE)
    w_col<-sapply(seq_len(ncol), function(j) max(sapply(idx[, j][idx[, j]<=length(plots)], function(k) sizes[[k]]$w)))
    h_row<-sapply(seq_len(nrow_p), function(i) max(sapply(idx[i, ][idx[i, ]<=length(plots)], function(k) sizes[[k]]$h)))
    message("Panel area computed: ", paste0(round(w_col, 1), collapse=" / "), " cm wide x ", paste0(round(h_row, 1), collapse=" / "), " cm high (per column / per row)")
    # one single legend: kept on the first plot only (all plots of a figure must share the same color scale), removed from the others
    if(!is.null(legend_position)){
      plots<-lapply(seq_along(plots), function(k) plots[[k]] + ggplot2::theme(legend.position=if(k==1) legend_position else "none"))
    }
    # legend on top of everything (above the titles of the first row): a dedicated first row is reserved for the collected legend
    top_legend<-identical(legend_position, "top")
    if(top_legend){
      # height of the legend (cm), measured on the first plot
      leg_h<-tryCatch({
        gl<-ggplot2::ggplotGrob(plots[[1]] + ggplot2::theme(legend.position="top"))
        top_panel<-min(gl$layout$t[grepl("^panel", gl$layout$name)])
        i<-which(grepl("guide-box", gl$layout$name) & gl$layout$t<top_panel)
        if(length(i)==0) 1.5 else sum(grid::convertHeight(gl$heights[unique(gl$layout$t[i])], "cm", valueOnly=TRUE))
      }, error=function(e) 1.5)
      # row 1 = area of the collected legend (all columns), then one area per plot (filled row by row)
      areas<-Reduce(c, c(list(patchwork::area(1, 1, 1, ncol)),
                         lapply(seq_along(plots), function(k) patchwork::area((k - 1) %/% ncol + 2, (k - 1) %% ncol + 1))))
      fig<-patchwork::wrap_plots(c(list(patchwork::guide_area()), plots), design=areas, widths=grid::unit(w_col, "cm"),
                                 heights=grid::unit(c(leg_h, h_row), "cm"), guides="collect")
    }else{
      fig<-patchwork::wrap_plots(plots, ncol=ncol, widths=grid::unit(w_col, "cm"), heights=grid::unit(h_row, "cm"), guides="collect")
    }
    # total size = size of the assembled figure (+2% and a margin)
    gt<-patchwork::patchworkGrob(fig)
    if(is.null(width)) width<-sum(grid::convertWidth(gt$widths, "cm", valueOnly=TRUE)) * 1.02 + margin
    if(is.null(height)) height<-sum(grid::convertHeight(gt$heights, "cm", valueOnly=TRUE)) * 1.02 + margin
  }
  
  # one file per format; compression is only passed to tiff devices
  for(k in seq_along(files)){
    args<-c(list(plot=fig, filename=files[k], width=width, height=height, units="cm", dpi=dpi, bg="white", limitsize=FALSE), list(...))
    if(ext[k] %in% c("tiff","tif")){
      args$compression<-compression
      args$device<-if(requireNamespace("ragg", quietly=TRUE)) ragg::agg_tiff else "tiff"
    }
    do.call(ggplot2::ggsave, args)
    message("Figure saved: ", files[k], " (", round(width, 1), " x ", round(height, 1), " cm)")
  }
  invisible(c(width=width, height=height))
}

#- - - - - -
## Safe file name from a label

to_filename<-function(label){
  label %>% str_replace_all("[^[:alnum:] ]", "") %>% str_replace_all("\\s+", "_") %>% str_trunc(60, ellipsis="")
}

#- - - - - -
## Tidy table of a Cox model: hazard ratio (HR) with 95% CI, direction of the association, formatted text
#
# - exponentiate=TRUE: HR per 1 unit of the covariate (cluster-robust standard errors when the model has cluster(id))
# - direction: positive / negative if the 95% CI excludes 1, uncertain otherwise

tidy_cox<-function(mod, label){
  broom::tidy(mod, exponentiate=TRUE, conf.int=TRUE) %>%
    mutate(direction=case_when(conf.low>1 ~ "positive association", conf.high<1 ~ "negative association", TRUE ~ "uncertain association"),
           HR=paste0(vf(estimate), " [", vf(conf.low), "; ", vf(conf.high), "]"),
           model=label)
}

#- - - - - -
## Tidy table of cox.zph (Schoenfeld test of the proportional hazards assumption)

tidy_ph<-function(z){
  as_tibble(data.frame(term=rownames(z$table), z$table)) %>%
    rename_with(~ c("term","chisq","df","p_value"), everything()) %>%
    mutate(PH_violated=p_value<alpha_ph)
}

#- - - - - -
## Tidy table of car::vif (VIF, or GVIF for factors with more than 2 levels)

tidy_vif<-function(model, thr){
  v<-car::vif(model)
  if(is.matrix(v)){
    as_tibble(v, rownames="term") %>% rename(GVIF=1, Df=2, `GVIF^(1/(2*Df))`=3) %>% mutate(VIF_flag=`GVIF^(1/(2*Df))`>sqrt(thr))
  }else tibble(term=names(v), VIF=v) %>% mutate(VIF_flag=VIF>thr)
}

#- - - - - -
## Building the counting-process dataset (one row per interval (tstart, tstop]) with time-varying covariates
#
# Arguments:
# - base: one row per individual
#   ¤ id, time (end of follow-up)
#   ¤ ev_code (0 = censored, 1, 2, ... = causes)
#   ¤ time-fixed covariates
# - meas: one row per measurement: id, day (measurement time, same unit as time), one column per time-varying covariate (names = tv_names)
# - lag: the measurement is only usable `lag` time units after it was made (0 = immediately) -> limits reverse causation within the interval
# - max_carry: a value is set to missing when it is older than max_carry time units (Inf = carried forward until the next measurement)
#
# Method: 
# - for each covariate, the last observed (non-missing) value is carried forward from its measurement time (tdc)
# - the time of this measurement is kept to compute the age of the value (since_<covariate>)
# - intervals are split when a value expires
#
# Returns: tmerge dataset with columns id, tstart, tstop, endpt (0 = no event in the interval, else the cause), covariates and since_<covariate>

build_long<-function(base, meas, tv_names, lag=0, max_carry=Inf){
  long<-tmerge(base, base, id=id, endpt=event(time, ev_code))
  for(v in tv_names){
    m<-meas[!is.na(meas[[v]]), c("id","day",v)]
    names(m)[3]<-"val"
    m$day<-m$day + lag
    m$last<-m$day
    long<-tmerge(long, m, id=id, tv_value=tdc(day, val), tv_last=tdc(day, last))
    if(is.finite(max_carry)){
      ex<-data.frame(id=m$id, day=m$day + max_carry, one=1)
      long<-tmerge(long, ex, id=id, tv_split=tdc(day, one)) # only used to split the intervals at the expiry time
      long$tv_split<-NULL
    }
    age_value<-long$tstart - long$tv_last
    long[[paste0("since_", v)]]<-age_value
    long[[v]]<-ifelse(!is.na(age_value) & age_value<max_carry, long$tv_value, NA_real_)
    long$tv_value<-NULL
    long$tv_last<-NULL
  }
  as_tibble(long)
}

#- - - - - -
## Checking the counting-process dataset: intervals ordered and non-empty, at most one event per individual, same person-time as the baseline data

check_long<-function(long, base){
  if(!all(long$tstop>long$tstart)) stop("Some intervals are empty or reversed (tstop <= tstart)")
  n_ev<-long %>% group_by(id) %>% summarise(n_ev=sum(endpt!=0), .groups="drop")
  if(any(n_ev$n_ev>1)) stop("Some individuals have more than one event")
  pt_long<-sum(long$tstop - long$tstart)
  pt_base<-sum(base$time)
  if(abs(pt_long - pt_base)>1e-6 * pt_base) warning("Person-time differs between the counting-process data and the baseline data")
  invisible(TRUE)
}

#- - - - - -
## Time scale: "followup" (start and stop = time since inclusion; age can then be a covariate) or "age" (start and stop = age, delayed entry at the age of inclusion)

add_time_scale<-function(long){
  if(time_scale=="age"){
    long$start<-long[[age_var]] + long$tstart/time_per_year
    long$stop<-long[[age_var]] + long$tstop/time_per_year
  }else{
    long$start<-long$tstart
    long$stop<-long$tstop
  }
  long
}

#- - - - - -
## Fitting a cause-specific Cox model on the counting-process data
#
# - events of the other causes are censored; cluster(id): robust standard errors, as an individual contributes several rows
# - data rows with a missing covariate are dropped by coxph (complete-case analysis)

fit_cox<-function(data, rhs, cause, cluster=TRUE){
  data$event_cs<-as.integer(data$endpt==cause)
  f<-as.formula(paste0("Surv(start, stop, event_cs) ~ ", paste(rhs, collapse=" + "), if(cluster) " + cluster(id)" else ""))
  coxph(f, data=data, ties="efron")
}

#- - - - - -
## VIF overall and by period of time (quartiles of the stop time): collinearity can change over follow-up

vif_by_period<-function(long, rhs){
  d<-long[, c("stop", rhs)] %>% drop_na() %>% mutate(across(where(is.character), factor))
  bins<-cut(d$stop, unique(quantile(d$stop, 0:4/4)), include.lowest=TRUE, labels=FALSE)
  bind_rows(lapply(c("All", sort(unique(bins))), function(b){
    dd<-if(b=="All") d else d[bins==as.integer(b), ]
    tryCatch(tidy_vif(lm(as.formula(paste("stop ~", paste(rhs, collapse=" + "))), data=dd), VIF_set) %>% mutate(period=ifelse(b=="All", "All", paste0("Quarter ", b))),
             error=function(e) NULL)
  }))
}

#- - - - - -
## Follow-up, event proportion, cumulative risk and incidence
#
# - ev_code: 0 = censored, otherwise the cause (codes)
# - time in the unit of the data
# - tpy = number of time units per year
# - horizon_years = horizon of the cumulative risk
# - median follow-up: reverse Kaplan-Meier (censoring is treated as the event), with quartiles and 95% CI (Brookmeyer-Crowley) in years
# - events_per_1000_persons: share of the individuals who had the event during the whole follow-up x 1000, exact binomial 95% CI
# - cumulative_risk_horizon_pct: Kaplan-Meier (1 - survival) at horizon_years, cause-specific (events of other causes are censored)
# - incidence_per_1000_py: events / person-years x 1000, exact Poisson 95% CI
# - person-years = sum of the follow-up times / tpy

rate_table<-function(time, ev_code, codes, labels, tpy, horizon_years=5){
  py<-sum(time)/tpy
  n<-length(time)
  fu<-survfit(Surv(time, as.integer(ev_code==0)) ~ 1)
  q<-quantile(fu, probs=c(0.25, 0.5, 0.75))
  fq<-as.numeric(q$quantile)/tpy
  fl<-as.numeric(q$lower)/tpy
  fh<-as.numeric(q$upper)/tpy
  method<-"reverse Kaplan-Meier"
  if(is.na(fq[2])){
    fq<-as.numeric(quantile(time, c(0.25, 0.5, 0.75)))/tpy
    fl<-rep(NA_real_, 3)
    fh<-rep(NA_real_, 3)
    method<-"median of the follow-up times (no censoring)"
  }
  bind_rows(lapply(codes, function(k){
    x<-sum(ev_code==k)
    b<-binom.test(x, n)$conf.int
    p<-poisson.test(x, py)$conf.int
    km<-survfit(Surv(time/tpy, as.integer(ev_code==k)) ~ 1)
    s<-summary(km, times=horizon_years, extend=TRUE)
    tibble(cause=labels[[as.character(k)]], persons=n, events=x,
           events_per_1000_persons=1000 * x/n, events_per_1000_persons_low=1000 * b[1], events_per_1000_persons_high=1000 * b[2],
           cumulative_risk_horizon_pct=100 * (1 - s$surv), cumulative_risk_low=100 * (1 - s$upper), cumulative_risk_high=100 * (1 - s$lower), horizon_years=horizon_years,
           person_years=py, incidence_per_1000_py=1000 * x/py, incidence_low=1000 * p[1], incidence_high=1000 * p[2],
           median_follow_up_years=fq[2], median_follow_up_low=fl[2], median_follow_up_high=fh[2], follow_up_q1_years=fq[1], follow_up_q3_years=fq[3],
           follow_up_method=method)
  }))
}

#- - - - - -
## Prevalence at inclusion

prevalence_table<-function(prevalent, label="Disease at inclusion"){
  x<-sum(prevalent==1, na.rm=TRUE)
  n<-sum(!is.na(prevalent))
  b<-binom.test(x, n)$conf.int
  tibble(condition=label, persons=n, prevalent_cases=x, prevalence_per_1000=1000 * x/n, prevalence_low=1000 * b[1], prevalence_high=1000 * b[2])
}

#----------------------------------------------------------------
#### User settings - edit with your own data ####

## Variables of the baseline table (one row per individual)
id_var<-"id" # individual identifier
time_var<-"time" # end of follow-up (event or censoring), from inclusion
event_var<-"status" # 0 = censored, 1, 2, ... = causes of event (with a single cause: 0/1)

## Causes: names = codes of event_var, values = labels; the FIRST one is the cause of interest, the others are competing events (cause-specific hazards)
cause_labels<-c("2"="Death", "1"="Transplant")

## Time-fixed covariates (baseline table)
fixed_covs<-c("age","sex")

## Time scale: "followup" (time since inclusion; age is a covariate) or "age" (age as the time scale; age_var = age at inclusion, in years)
time_scale<-"followup"
age_var<-"age"
time_per_year<-365.25 # number of time units (time_var) per year (used for the age scale)

## Time-varying covariates (measurement table): name, label, transform ("none", "log", "sqrt"); numeric or binary 0/1 variables
tv_dict<-tibble(name=c("bili","albumin","protime","ascites"),
                label=c("Bilirubin","Albumin","Prothrombin time","Ascites"),
                transform=c("log","none","none","none"))
meas_id_var<-"id"
meas_time_var<-"day" # time of the measurement, same unit and same origin as time_var

## Measurement process
lag_main<-0 # main analysis: value usable immediately after its measurement
lag_sens<-30 # sensitivity analysis: value usable 30 time units after its measurement
max_carry_main<-Inf # main analysis: value carried forward until the next measurement
max_carry_sens<-365 # sensitivity analysis: values older than 365 time units are set to missing

## Effect modification: variable (among fixed_covs) whose interaction with the other covariates is tested (NULL to skip)
interaction_var<-"sex"

## Thresholds
alpha_ph<-0.05 # alpha risk of the proportional hazards test
VIF_set<-5 # VIF above which collinearity is flagged
spline_df<-3 # degrees of freedom of the natural splines (linearity check)

#----------------------------------------------------------------
#### Data import and preprocessing (replace with your own data) ####

# - baseline table: id, follow-up time, event code, fixed covariates
# - measurement table: id, measurement time, time-varying covariates
# - example: survival::pbc (312 randomized patients) and survival::pbcseq (repeated measurements)

base_raw<-as_tibble(survival::pbc) %>% filter(id<=312)
meas_raw<-as_tibble(survival::pbcseq)

## Standardizing names: base = id, time, ev_code, covariates; meas = id, day, covariates
base<-base_raw %>%
  select(all_of(c(id_var, time_var, event_var, fixed_covs))) %>%
  rename(id=all_of(id_var), time=all_of(time_var), ev_code=all_of(event_var)) %>%
  mutate(across(where(is.character), factor))

tv_dict<-tv_dict %>% mutate(model_name=ifelse(transform=="none", name, paste0(transform, "_", name)))
meas<-meas_raw %>% select(all_of(c(meas_id_var, meas_time_var, tv_dict$name))) %>% rename(id=all_of(meas_id_var), day=all_of(meas_time_var))

## Transforming the time-varying covariates (non-finite values after transformation are set to missing)
for(i in seq_len(nrow(tv_dict))){
  x<-meas[[tv_dict$name[i]]]
  x<-switch(tv_dict$transform[i], none=x, log=log(x), sqrt=sqrt(x), stop("Unknown transform: ", tv_dict$transform[i]))
  x[!is.finite(x)]<-NA
  meas[[tv_dict$model_name[i]]]<-x
}
tv_names<-tv_dict$model_name
meas<-meas %>% filter(id %in% base$id) %>% select(id, day, all_of(tv_names))

if(any(base$time<=0)) stop("Follow-up times must be strictly positive")
cause_codes<-as.numeric(names(cause_labels))
cause_main<-cause_codes[1]
fixed_terms<-if(time_scale=="age") setdiff(fixed_covs, age_var) else fixed_covs # age is the time scale: not a covariate
main_terms<-c(fixed_terms, tv_names)

#----------------------------------------------------------------
#### Counting-process dataset ####

long_main<-build_long(base, meas, tv_names, lag_main, max_carry_main) %>% add_time_scale()
check_long(long_main, base)

# Example of the structure for one individual (one row per interval; the value of the covariates changes at each new measurement)
long_main %>% filter(id==unique(id)[2]) %>% select(id, tstart, tstop, endpt, all_of(tv_names)) %>% print(n=Inf)

#----------------------------------------------------------------
#### Descriptive statistics ####

## Individuals, events, person-time
desc_overview<-tibble(n_individuals=n_distinct(base$id),
                      n_rows=nrow(long_main),
                      person_time=sum(long_main$tstop - long_main$tstart),
                      events=paste(sapply(cause_codes, function(k) paste0(cause_labels[as.character(k)], ": ", sum(base$ev_code==k))), collapse="; "))
desc_overview

rates<-rate_table(base$time, base$ev_code, cause_codes, cause_labels, time_per_year)
rates
fwrite(rates, safe_path(paste0("Descriptive_rates and follow-up_", Sys.Date(), ".csv")), sep=";")

prev<-prevalence_table(as.integer(base_raw$ascites==1), "Ascites at inclusion") # replace by your variable
prev
fwrite(prev, safe_path(paste0("Descriptive_prevalence_", Sys.Date(), ".csv")), sep=";")

## Measurement process: number of measurements per individual and time between two measurements (the model assumes that the timing of the visits does not depend on health)
n_meas<-meas %>% group_by(id) %>% summarise(n=n(), .groups="drop")
gaps<-meas %>% arrange(id, day) %>% group_by(id) %>% mutate(gap=day - lag(day)) %>% ungroup() %>% filter(!is.na(gap))
desc_meas<-tibble(measurements_per_individual_median=median(n_meas$n), q1=quantile(n_meas$n, 0.25), q3=quantile(n_meas$n, 0.75),
                  gap_median=median(gaps$gap), gap_q1=quantile(gaps$gap, 0.25), gap_q3=quantile(gaps$gap, 0.75))
desc_meas

## Missing values and age of the value used, by covariate (share of the person-time)
pt<-long_main$tstop - long_main$tstart
desc_tv<-bind_rows(lapply(seq_len(nrow(tv_dict)), function(i){
  v<-tv_dict$model_name[i]
  s<-long_main[[paste0("since_", v)]]
  tibble(covariate=tv_dict$label[i], pct_person_time_missing=100 * sum(pt[is.na(long_main[[v]])])/sum(pt),
         median_age_of_value=median(s, na.rm=TRUE), q3_age_of_value=quantile(s, 0.75, na.rm=TRUE))
}))
desc_tv

fwrite(desc_overview, safe_path(paste0("Descriptive_overview_", Sys.Date(), ".csv")), sep=";")
fwrite(desc_meas, safe_path(paste0("Descriptive_measurement process_", Sys.Date(), ".csv")), sep=";")
fwrite(desc_tv, safe_path(paste0("Descriptive_time-varying covariates_", Sys.Date(), ".csv")), sep=";")

#----------------------------------------------------------------
#### Cox models with time-varying covariates ####

## Datasets of the sensitivity analyses
long_base<-long_main %>% group_by(id) %>% mutate(across(all_of(tv_names), ~ first(.x))) %>% ungroup() # value of the first interval, kept fixed (no future information)

long_dec<-long_main %>% group_by(id) %>% mutate(across(all_of(tv_names), ~ first(.x), .names="{.col}_base")) %>% ungroup()
for(v in tv_names) long_dec[[paste0(v, "_change")]]<-long_dec[[v]] - long_dec[[paste0(v, "_base")]] # change since the first interval
dec_terms<-c(fixed_terms, paste0(tv_names, "_base"), paste0(tv_names, "_change"))

long_lag<-build_long(base, meas, tv_names, lag_sens, max_carry_main) %>% add_time_scale()
long_old<-build_long(base, meas, tv_names, lag_main, max_carry_sens) %>% add_time_scale()

## Models for the cause of interest
models<-list("Current value (main)"=fit_cox(long_main, main_terms, cause_main),
             "Baseline value only"=fit_cox(long_base, main_terms, cause_main),
             "Baseline value + change"=fit_cox(long_dec, dec_terms, cause_main),
             "Lagged value"=fit_cox(long_lag, main_terms, cause_main),
             "Old values excluded"=fit_cox(long_old, main_terms, cause_main))
cox_main<-models[[1]]
summary(cox_main)

## Other causes (cause-specific hazards of the competing events), current value
for(k in cause_codes[-1]) models[[paste0("Current value (", cause_labels[as.character(k)], ")")]]<-fit_cox(long_main, main_terms, k)

## Number of individuals and events actually used by each model (complete cases)
model_n<-tibble(model=names(models), n_rows=sapply(models, function(m) m$n), n_events=sapply(models, function(m) m$nevent))
model_n

## Tables
res_all<-bind_rows(lapply(names(models), function(nm) tidy_cox(models[[nm]], nm)))
fwrite(res_all, safe_path(paste0("Cox models_hazard ratios_", Sys.Date(), ".csv")), sep=";")
fwrite(model_n, safe_path(paste0("Cox models_sample sizes_", Sys.Date(), ".csv")), sep=";")
qs_save(models, safe_path(paste0("Cox models_", Sys.Date(), ".qs2")))

#----------------------------------------------------------------
#### Assumption testing (main model) ####

# The tests below use the main model fitted without cluster(id): the same coefficients, the tests are not defined with robust variances
cox_nc<-fit_cox(long_main, main_terms, cause_main, cluster=FALSE)

#- - - - - - - - - - - - - - - -
## Independence assumption: an individual contributes several rows, which are handled by the robust variance (cluster(id)); compare the standard errors
se_compare<-tibble(term=names(coef(cox_main)), se_robust=sqrt(diag(vcov(cox_main))), se_naive=sqrt(diag(vcov(cox_nc)))) %>% mutate(ratio=se_robust/se_naive)
se_compare # a ratio far from 1 means that the within-individual correlation matters

#- - - - - - - - - - - - - - - -
## Proportional hazards (scaled Schoenfeld residuals), global and per covariate; p < alpha_ph = violation

PH_test<-cox.zph(cox_nc)
ph_summary<-tidy_ph(PH_test)
ph_summary
ph_violated_vars<-ph_summary %>% filter(PH_violated, term!="GLOBAL") %>% pull(term)

png(file=safe_path(paste0("PH test_scaled Schoenfeld residuals_", Sys.Date(), ".png")), res=300, width=3000, height=2400)
print(ggcoxzph(PH_test, font.main=11))
dev.off()
fwrite(ph_summary, safe_path(paste0("PH test_", Sys.Date(), ".csv")), sep=";")

#- - - - - - - - - - - - - - - -
## Linearity of the continuous covariates: likelihood-ratio test of a natural spline (spline_df degrees of freedom) against a linear term
#
# - plot: log hazard ratio against the covariate value (spline fit, 95% CI)

cont_terms<-main_terms[sapply(main_terms, function(v) is.numeric(long_main[[v]]) && length(unique(na.omit(long_main[[v]])))>5)]
lin_res<-list()
lin_plots<-list()
for(v in cont_terms){
  rhs_s<-ifelse(main_terms==v, paste0("splines::ns(", v, ", df=", spline_df, ")"), main_terms)
  d<-long_main %>% drop_na(all_of(main_terms))
  m_lin<-fit_cox(d, main_terms, cause_main, cluster=FALSE)
  m_ns<-fit_cox(d, rhs_s, cause_main, cluster=FALSE)
  lrt<-anova(m_lin, m_ns)
  lin_res[[v]]<-tibble(term=v, chisq=lrt$Chisq[2], df=lrt$Df[2], p_value=lrt$`Pr(>|Chi|)`[2], non_linear=lrt$`Pr(>|Chi|)`[2]<alpha_ph)
  pr<-predict(m_ns, type="terms", se.fit=TRUE)
  j<-grep(paste0("ns\\(", v), colnames(pr$fit))[1]
  dd<-tibble(x=d[[v]], fit=pr$fit[, j], se=pr$se.fit[, j]) %>% arrange(x)
  lin_plots[[v]]<-ggplot(dd, aes(x=x, y=fit)) +
    geom_ribbon(aes(ymin=fit - 1.96 * se, ymax=fit + 1.96 * se), fill="#CAE1FF") +
    geom_line(color="#C35C33", linewidth=1) +
    geom_hline(yintercept=0, linetype="dashed") +
    labs(x=v, y="Log HR (relative to the mean)") +
    theme_Gaia()
}
lin_res<-bind_rows(lin_res)
lin_res
fwrite(lin_res, safe_path(paste0("Linearity test_", Sys.Date(), ".csv")), sep=";")
if(length(lin_plots)>0) save_plot(lin_plots, "Linearity plots", ncol=2, ext="png")

#- - - - - - - - - - - - - - - -
## Multicollinearity: VIF overall and by quarter of follow-up (the correlation of time-varying covariates can change over time)

vif_res<-vif_by_period(long_main, main_terms)
vif_res
fwrite(vif_res, safe_path(paste0("Multicollinearity_", Sys.Date(), ".csv")), sep=";")

#- - - - - - - - - - - - - - - -
## Effect modification
#
# - likelihood-ratio test of each covariate x interaction_var term (one at a time)
# - global test (all interactions together),
# - hazard ratios in each stratum of interaction_var

if(!is.null(interaction_var) && interaction_var %in% main_terms){
  d<-long_main %>% drop_na(all_of(main_terms))
  others<-setdiff(main_terms, interaction_var)
  m0<-fit_cox(d, main_terms, cause_main, cluster=FALSE)
  
  # one interaction at a time
  interaction_res<-bind_rows(lapply(others, function(cv){
    m1<-fit_cox(d, c(main_terms, paste0(cv, ":", interaction_var)), cause_main, cluster=FALSE)
    lrt<-anova(m0, m1)
    tibble(covariate=cv, interaction=paste0(cv, " x ", interaction_var), chisq=lrt$Chisq[2], df=lrt$Df[2], p_value=lrt$`Pr(>|Chi|)`[2], modification=lrt$`Pr(>|Chi|)`[2]<0.05)
  }))
  interaction_res
  
  # global test: all the interactions together
  m_all<-fit_cox(d, c(main_terms, paste0(others, ":", interaction_var)), cause_main, cluster=FALSE)
  g<-anova(m0, m_all)
  interaction_global<-tibble(test="All interactions together", chisq=g$Chisq[2], df=g$Df[2], p_value=g$`Pr(>|Chi|)`[2])
  interaction_global
  
  # hazard ratios in each stratum
  lev<-levels(factor(d[[interaction_var]]))
  strat_res<-bind_rows(lapply(lev, function(l){
    ds<-d[d[[interaction_var]]==l, ]
    tidy_cox(fit_cox(ds, others, cause_main), paste0(interaction_var, " = ", l))
  }))
  strat_res
  
  fwrite(interaction_res, safe_path(paste0("Interaction test_", interaction_var, "_", Sys.Date(), ".csv")), sep=";")
  fwrite(interaction_global, safe_path(paste0("Interaction test_global_", interaction_var, "_", Sys.Date(), ".csv")), sep=";")
  fwrite(strat_res, safe_path(paste0("Cox models_stratified by_", interaction_var, "_", Sys.Date(), ".csv")), sep=";")
}

#----------------------------------------------------------------
#### In case of violation of the proportional hazards assumption ####

# The effect of each violating covariate is allowed to change with log(time): covariate x log(time) interaction (tt() term)

if(length(ph_violated_vars)>0){
  tt_terms<-ifelse(main_terms %in% ph_violated_vars, paste0("tt(", main_terms, ")"), main_terms)
  d<-long_main
  d$event_cs<-as.integer(d$endpt==cause_main)
  cox_tt<-coxph(as.formula(paste0("Surv(start, stop, event_cs) ~ ", paste(tt_terms, collapse=" + "), " + cluster(id)")), data=d,
                tt=function(x, t, ...) x * log(t), ties="efron")
  summary(cox_tt)
  # tt() terms give the log hazard ratio at log(time) = 0 and the change per unit of log(time): the HR at a given time t is exp(beta + gamma x log(t))
  res_tt<-tidy_cox(cox_tt, "Time-varying effect (covariate x log time)")
  fwrite(res_tt, safe_path(paste0("Cox model_time-varying effects_", Sys.Date(), ".csv")), sep=";")
  models[["Time-varying effect (covariate x log time)"]]<-cox_tt
  data.frame(model=c("Main","With tt()"), AIC=c(AIC(cox_nc), AIC(cox_tt))) # AIC of a model with cluster() is not strictly comparable: indicative only
}

#----------------------------------------------------------------
#### Competing risks (only if there are several causes) ####

if(length(cause_codes)>1){
  
  # data: one row per individual, covariates at inclusion (value of the first interval)
  bl<-long_base %>% group_by(id) %>% slice(1) %>% ungroup() %>% select(id, all_of(main_terms))
  d_cr<-base %>% select(id, time, ev_code) %>% inner_join(bl, by="id") %>% drop_na()
  d_cr$status_f<-factor(d_cr$ev_code, levels=c(0, cause_codes), labels=c("censored", paste0("c", cause_codes)))
  
  # cumulative incidence function of each cause (Aalen-Johansen), overall (time in years)
  d_cr$time_y<-d_cr$time/time_per_year
  cif_all<-tidycmprsk::cuminc(Surv(time_y, status_f) ~ 1, data=d_cr)
  p_cif<-ggsurvfit::ggcuminc(cif_all, outcome=paste0("c", cause_codes)) +
    ggsurvfit::add_confidence_interval() +
    labs(x="Time (years)", y="Cumulative incidence") +
    theme_Gaia()
  ggplot2::ggsave(paste0("Cumulative incidence function_", Sys.Date(), ".png"), p_cif, width=30, height=18, units="cm", dpi=300)
  
  # CIF by interaction_var (fixed covariate)
  if(!is.null(interaction_var) && interaction_var %in% names(d_cr)){
    cif_s<-tidycmprsk::cuminc(as.formula(paste0("Surv(time_y, status_f) ~ ", interaction_var)), data=d_cr)
    p_cif_s<-ggsurvfit::ggcuminc(cif_s, outcome=paste0("c", cause_main)) +
      ggsurvfit::add_confidence_interval() +
      labs(x="Time (years)", y=paste0("Cumulative incidence of ", cause_labels[[as.character(cause_main)]]), color=interaction_var, fill=interaction_var) +
      theme_Gaia()
    ggplot2::ggsave(paste0("Cumulative incidence function_by ", interaction_var, "_", Sys.Date(), ".png"), p_cif_s, width=30, height=18, units="cm", dpi=300)
  }
  
  # Fine-Gray model on the covariates at inclusion (sensitivity analysis): subdistribution HRs
  fg<-survival::finegray(Surv(time, status_f) ~ ., data=d_cr[, c("id", "time", "status_f", main_terms)], etype=paste0("c", cause_main))
  cox_fg<-coxph(as.formula(paste0("Surv(fgstart, fgstop, fgstatus) ~ ", paste(main_terms, collapse=" + "), " + cluster(id)")), data=fg, weights=fgwt)
  res_fg<-tidy_cox(cox_fg, "Fine-Gray (covariates at inclusion)")
  res_fg
  fwrite(res_fg, safe_path(paste0("Fine-Gray model_covariates at inclusion_", Sys.Date(), ".csv")), sep=";")
  models[["Fine-Gray (covariates at inclusion)"]]<-cox_fg
}

#----------------------------------------------------------------
#### Visualization ####

#- - - - - - - - - - - - - - - -
## Forest plot of each model (hazard ratio per 1 unit of the covariate, log scale; the label of each bar comes from the model itself)

for(nm in names(models)){
  d<-tidy_cox(models[[nm]], nm)
  d$term_lab<-factor(paste0(d$term, "\n", d$HR), levels=rev(paste0(d$term, "\n", d$HR))) # HR written under the name of the covariate
  d$direction<-factor(d$direction, levels=c("positive association","negative association","uncertain association"))
  p<-ggplot(d, aes(x=estimate, y=term_lab, color=direction)) +
    geom_vline(xintercept=1, linetype="dashed", color="grey40") +
    geom_errorbar(aes(xmin=conf.low, xmax=conf.high), orientation="y", width=0.25) +
    geom_point(size=2.5) +
    scale_x_log10() +
    scale_color_manual("Association:", values=c("positive association"="red", "negative association"="#40B696", "uncertain association"="black"), drop=FALSE) +
    labs(x="Hazard ratio per 1 unit (log scale)", y="") +
    theme_Gaia()
  save_plot(p, paste0("Forest plot_", to_filename(nm)), ext="png", legend_position="top")
}

#- - - - - - - - - - - - - - - -
## Comparison of the models for the cause of interest: terms present in all the models of the sensitivity analyses

cmp_models<-names(models)[1:5]
cmp<-bind_rows(lapply(cmp_models, function(nm) tidy_cox(models[[nm]], nm))) %>% filter(term %in% main_terms) %>%
  mutate(model=factor(model, levels=cmp_models), term=factor(term, levels=rev(main_terms)))
p_cmp<-ggplot(cmp, aes(x=estimate, y=term, color=model)) +
  geom_vline(xintercept=1, linetype="dashed", color="grey40") +
  geom_errorbar(aes(xmin=conf.low, xmax=conf.high), orientation="y", width=0.3, position=position_dodge(width=0.7)) +
  geom_point(size=2, position=position_dodge(width=0.7)) +
  scale_x_log10() +
  labs(x="Hazard ratio per 1 unit (log scale)", y="", color="Model:") +
  theme_Gaia()
save_plot(p_cmp, "Forest plot_comparison of the models", ext="png", legend_position="top")

#- - - - - - - - - - - - - - - -
## Trajectories of the time-varying covariates according to the final status of the individual (smooth of the repeated measurements over time)

traj<-meas %>%
  left_join(base %>% select(id, ev_code), by="id") %>%
  mutate(final_status=factor(ev_code, levels=c(0, cause_codes), labels=c("Censored", unname(cause_labels)))) %>%
  pivot_longer(all_of(tv_names), names_to="covariate", values_to="value") %>%
  left_join(tv_dict %>% select(model_name, label), by=c("covariate"="model_name")) %>%
  drop_na(value)
p_traj<-ggplot(traj, aes(x=day, y=value, color=final_status)) +
  geom_smooth(method="loess", formula=y ~ x, se=TRUE, linewidth=1) +
  facet_wrap(~label, scales="free_y") +
  labs(x=paste0("Time (", meas_time_var, ")"), y="Value of the covariate (transformed scale)", color="Final status:") +
  theme_Gaia()
save_plot(p_traj, "Trajectories of the time-varying covariates", ext="png", legend_position="top")

#- - - - - - - - - - - - - - - -
## Survival, cumulative events and cumulative hazard (Kaplan-Meier on the baseline table: one row per individual, no covariate)

save_surv<-function(p, name){
  tryCatch(save_plot(p, name, ext="png", width=30, height=18),
           error=function(e) ggplot2::ggsave(paste0(name, "_", Sys.Date(), ".png"), p, width=30, height=18, units="cm", dpi=300))
}

km_data<-base %>% mutate(time_y=time/time_per_year, ev_main=as.integer(ev_code==cause_main))
curve_types<-list(list(type="survival", y="Survival probability", nm="Survival"),
                  list(type="risk", y="Cumulative event probability", nm="Cumulative events"),
                  list(type="cumhaz", y="Cumulative hazard", nm="Cumulative hazard"))

for(cfg in curve_types){
  # overall
  p<-ggsurvfit::survfit2(Surv(time_y, ev_main) ~ 1, data=km_data) %>%
    ggsurvfit::ggsurvfit(type=cfg$type, color="#40B696") +
    ggsurvfit::add_confidence_interval(fill="#40B696") +
    ggsurvfit::add_risktable() +
    labs(x="Time (years)", y=cfg$y) +
    theme_Gaia()
  save_surv(p, paste0(cfg$nm, "_overall"))
  
  # by interaction_var (fixed covariate), with log-rank test
  if(!is.null(interaction_var) && interaction_var %in% names(km_data)){
    ps<-ggsurvfit::survfit2(as.formula(paste0("Surv(time_y, ev_main) ~ ", interaction_var)), data=km_data) %>%
      ggsurvfit::ggsurvfit(type=cfg$type) +
      ggsurvfit::add_confidence_interval() +
      ggsurvfit::add_risktable() +
      ggsurvfit::add_pvalue(caption="Log-rank {p.value}") +
      labs(x="Time (years)", y=cfg$y, color=interaction_var, fill=interaction_var) +
      theme_Gaia()
    save_surv(ps, paste0(cfg$nm, "_by ", interaction_var))
  }
}

#----------------------------------------------------------------
#### Session information (reproducibility) ####

writeLines(capture.output(sessionInfo()), safe_path(paste0("Session info_", Sys.Date(), ".txt")))
