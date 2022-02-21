# Accessibility Mapping in R
# 
# Dan Weiss
# Telethon Kids Institute, Perth 
# Malaria Atlas Project, University of Oxford
# 2020-08-06
#
# This script requires the gdistance package (van Etten, J. R Package gdistance: Distances and Routes on Geographical Grids. Journal of Statistical Software 76, 1-21)
# https://cran.r-project.org/web/packages/gdistance/index.html
#
# This script requires the two user supplied datasets:
# (a) A friction surface, two of which are available here: https://malariaatlas.org/research-project/accessibility_to_healthcare/
# (b) A user-supplied .csv of points (i.e., known geographic coordinates) 
#
# Notes:
# (a) All file paths and names should be changed as needed.
# (b) Important runtime details can be found in the comments.
# (c) This script is suitable only for analyses of moderately sized areas (e.g., up to 10 million km^2 in lower latitude settings - GLOBAL RUNS WILL NOT WORK).
#     We recommend using Google Earth Engine for larger areas, with the exception of high-latitude areas where custom approaches are typically required.
#
# Citations: 
#
# D. J. Weiss, A. Nelson, C. A. Vargas-Ruiz, K. Gligoric, S. Bavadekar, E. Gabrilovich, A. Bertozzi-Villa, J. Rozier, H. S. Gibson, T. Shekel, C. Kamath, A. Lieber, K. Schulman,
# Y. Shao, V. Qarkaxhija, A. K. Nandi, S. H. Keddie, S. Rumisha, P. Amratia, R. Arambepola, E. G. Chestnutt, J. J. Millar, T. L. Symons, E. Cameron, K. E. Battle, S. Bhatt, 
# and P. W. Gething. Global maps of travel time to healthcare facilities. (2020) Nature Medicine.
#
# A. Nelson, D. J. Weiss, J. van Etten, A. Cattaneo, T. S. McMenomy,and J. Koo. A suite of global accessibility indicators. (2019). Nature Scientific Data. doi.org/10.1038/s41597-019-0265-5
#
# D. J. Weiss, A. Nelson, H.S. Gibson, W. Temperley, S. Peedell, A. Lieber, M. Hancher, E. Poyart, S. Belchior, N. Fullman, B. Mappin, U. Dalrymple, J. Rozier, 
# T.C.D. Lucas, R.E. Howes, L.S. Tusting, S.Y. Kang, E. Cameron, D. Bisanzio, K.E. Battle, S. Bhatt, and P.W. Gething. A global map of travel time to cities to assess 
# inequalities in accessibility in 2015. (2018). Nature. doi:10.1038/nature25181.
# 

## Required Packages
require(gdistance)

# User Defined Variables - used if clipping from the global layer, if no clipping is needed, see lines 66-67 (currently commented out).
# This could also be accomplished by importing a shapefile (for example) 
# Geographic Coordinates (WGS84)
left   <- -2.0
right  <- 0.0
bottom <- 50.0
top    <- 52.0
transition.matrix.exists.flag <- 0 # if the geo-corrected graph has already been made, this can save time.  Uses the same T.GC.filename as specified using the T.GC.filename variable.

# Input Files
# Note: the alternate, 'walking_only' friction surface is named friction_surface_2019_v51_walking_only.tif
friction.surface.filename <- 'E:\\accessibility\\friction_surface_2019_v51.tif.tif'
point.filename <- 'E:\\accessibility\\points.csv' # Just 2 columns.  Structured as [X_COORD, Y_COORD] aka [LONG, LAT].  Use a header.

# Output Files
T.filename <- 'E:\\accessibility\\study.area.T.rds'
T.GC.filename <- 'E:\\accessibility\\study.area.T.GC.rds'
output.filename <- 'E:\\accessibility\\study.area.accessibility.tif'

# Read in the points table
points <- read.csv(file = point.filename)

# Fetch the number of points
temp <- dim(points)
n.points <- temp[1]

#  Define the spatial template
friction <- raster(friction.surface.filename)
fs1 <- crop(friction, extent(left, right, bottom, top))
# Use the following line instead of the preceding 2 if clipping is not needed (i.e., to run globally), but be warned that trying this will far exceed the computational capacity available to most users.
# fs1 <- raster(friction.surface.filename) 

# Make the graph and the geocorrected version of the graph (or read in the latter).
if (transition.matrix.exists.flag == 1) {
  # Read in the transition matrix object if it has been pre-computed
  T.GC <- readRDS(T.GC.filename)
} else {
  # Make and geocorrect the transition matrix (i.e., the graph)
  T <- transition(fs1, function(x) 1/mean(x), 8) # RAM intensive, can be very slow for large areas
  saveRDS(T, T.filename)
  T.GC <- geoCorrection(T)                    
  saveRDS(T.GC, T.GC.filename)
}

# Convert the points into a matrix
xy.data.frame <- data.frame()
xy.data.frame[1:n.points,1] <- points[,1]
xy.data.frame[1:n.points,2] <- points[,2]
xy.matrix <- as.matrix(xy.data.frame)

# Run the accumulated cost algorithm to make the final output map. This can be quite slow (potentially hours).
temp.raster <- accCost(T.GC, xy.matrix)

# Write the resulting raster
writeRaster(temp.raster, output.filename)