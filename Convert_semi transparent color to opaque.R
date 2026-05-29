# ------------------------------------------------------------
# Function: de_alpha
# Purpose:  Convert a semi-transparent color to the equivalent
#           opaque RGB color that would appear visually similar
#           on a white background.
#
# Args:
#   col : A color string (named color, hex, or rgba) that may 
#         contain transparency (alpha).
#
# Returns:
#   A hex color string representing the fully opaque version.
# ------------------------------------------------------------

## Installing and loading packages
pack_needed<-c("grDevices","scales","colorspace")
is_installed<-pack_needed %in% rownames(installed.packages(all.available=TRUE))
if(any(is_installed == FALSE)){
  install.packages(pack_needed[!is_installed],repos = "http://cran.us.r-project.org")
}
invisible(lapply(pack_needed, library, character.only = TRUE))

## Creating the function to remove transparency
de_alpha<-function(col) {
  
  # Converting to RGB with alpha (0–255)
  rgba <- grDevices::col2rgb(col, alpha = TRUE)
  
  # Separating alpha and RGB
  alpha <- rgba["alpha", ] / 255
  rgb   <- rgba[1:3, , drop = FALSE]
  
  # Computing the effect of transparency on a white background
  # visual_rgb = alpha * rgb + (1 - alpha) * 255
  visual_rgb <- rgb * alpha + (1 - alpha) * 255
  
  # Returning as hex color
  grDevices::rgb(red=visual_rgb[1, ] / 255,
    green=visual_rgb[2, ] / 255,
    blue=visual_rgb[3, ] / 255)
}

# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------

# Defining a semi-transparent color
transparent_col <- scales::alpha("#009E73", 0.35)

# Converting to its fully opaque equivalent
nontransparent_col<-de_alpha(transparent_col)
nontransparent_col

# Darkening a color using colorspace
colorspace::darken("#00a064",0.4)
