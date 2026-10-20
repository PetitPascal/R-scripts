# Creating a function to compute the geometric mean
gmean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}
gmean<<-gmean

# Example of use
gmean(c(1.5,9.46,20.3,10.2))
