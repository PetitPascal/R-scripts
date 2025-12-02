# Creating a function to compute the geometric standard deviation
gsd <- function(data) { 
  log_data <- log(data)
  gsd <- exp(sd(log_data[is.finite(log_data)]))
  return(gsd)
}
gsd<<-gsd

# Example of use
gsd(c(1.5,9.46,20.3,10.2))
