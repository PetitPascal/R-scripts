# Requires Rcpp for efficient computation of pairwise counts across bins
if (!requireNamespace("Rcpp", quietly = TRUE)) install.packages("Rcpp")
library(Rcpp)

# C++ function to efficiently count pairwise SHAP relationships between bins
# This avoids large O(n^2) memory usage in R by using double precision for counts
cppFunction('
Rcpp::List bin_pairwise_counts(NumericVector bx, NumericVector bs, NumericVector bc) {
  int B = bx.size();
  double pos = 0.0;  // count of pairs where higher feature -> higher SHAP
  double neg = 0.0;  // count of pairs where higher feature -> lower SHAP
  double total = 0.0; // total number of valid pairs

  for (int i = 0; i < B; ++i) {
    for (int j = i + 1; j < B; ++j) {
      double ci = bc[i];
      double cj = bc[j];
      if (ci <= 0.0 || cj <= 0.0) continue; // skip empty bins
      double pairs = ci * cj; // pair count as double to avoid overflow

      // Compare feature bin centers
      if (bx[i] > bx[j]) {
        total += pairs;
        if (bs[i] > bs[j]) pos += pairs;
        else if (bs[i] < bs[j]) neg += pairs;
      } else if (bx[j] > bx[i]) {
        total += pairs;
        if (bs[j] > bs[i]) pos += pairs;
        else if (bs[j] < bs[i]) neg += pairs;
      }
      // if equal bin centers, skip as no ordering information
    }
  }

  return Rcpp::List::create(Rcpp::Named("pos") = pos,
                            Rcpp::Named("neg") = neg,
                            Rcpp::Named("total") = total);
}
', depends = "Rcpp")

# Main function to determine the SHAP direction for a given feature
determine_SHAP_direction <- function(df, feature_col, shap_col,
                                     majority_threshold = 0.55,
                                     n_quantile_bins = 200) {
  # Extracting vectors
  x <- df[[feature_col]]
  s <- df[[shap_col]]
  
  # Removing missing values
  valid <- complete.cases(x, s)
  x <- x[valid]; s <- s[valid]
  
  # Returning "undefined" if feature or SHAP has no variation
  if (length(unique(x)) <= 1 || length(unique(s)) <= 1) return("undefined")
  
  # ---- Binary or sparse features ----
  if (length(unique(x)) <= 2 || quantile(x, 0.75, na.rm = TRUE) == 0) {
    low_val <- min(x, na.rm = TRUE)
    high_val <- max(x, na.rm = TRUE)
    g0 <- s[x == low_val]
    g1 <- s[x == high_val]
    n0 <- length(g0); n1 <- length(g1)
    if (n0 == 0 || n1 == 0) return("neutral")
    
    # Pairwise counting using findInterval
    g0_sorted <- sort(g0)
    pos_count <- sum(findInterval(g1, g0_sorted, left.open = TRUE))
    total_pairs <- as.double(n0) * as.double(n1)
    prop_pos <- pos_count / total_pairs
    
    if (prop_pos >= majority_threshold) return("promoting")
    if (prop_pos <= (1 - majority_threshold)) return("mitigating")
    return("neutral")
  }
  
  # ---- Continuous or ordered features ----
  # Quantile-based binning to reduce computation and memory usage
  probs <- seq(0, 1, length.out = n_quantile_bins + 1)
  breaks <- unique(quantile(x, probs = probs, na.rm = TRUE, type = 7))
  if (length(breaks) <= 2) breaks <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length.out = 3)
  bins <- cut(x, breaks = breaks, include.lowest = TRUE)
  
  # Computing bin-level statistics: mean feature, mean SHAP, and count
  bx <- tapply(x, bins, mean, na.rm = TRUE)
  bs <- tapply(s, bins, mean, na.rm = TRUE)
  bc <- tapply(s, bins, length)
  keep <- !is.na(bx) & !is.na(bs) & !is.na(bc)
  bx <- as.numeric(bx[keep])
  bs <- as.numeric(bs[keep])
  bc <- as.numeric(bc[keep])
  
  if (length(bx) < 2) return("neutral") # Not enough bins to compare
  
  # Using Rcpp function to count pairwise relationships across bins
  cnts <- bin_pairwise_counts(bx, bs, bc)
  pos_pairs <- as.numeric(cnts$pos)
  neg_pairs <- as.numeric(cnts$neg)
  total_pairs <- as.numeric(cnts$total)
  if (total_pairs <= 0) return("neutral")
  
  # Computing proportions and assigning direction
  prop_pos <- pos_pairs / total_pairs
  prop_neg <- neg_pairs / total_pairs
  
  if (prop_pos >= majority_threshold) return("promoting")
  if (prop_neg >= majority_threshold) return("mitigating")
  return("neutral")
}
