## Libraries
library(ggplot2)
library(rgdal)
library(spdep)

## plotAreaCol
# This functions plots the values and saves the figure to the specified file.
# It is advisable to not plot directly in R since it will take a long time.
# The arguments are:
#   fNamme: file name for saving figure
#   width: width of figure in inches
#   height: height of figure in inches
#   estVal: the k values to be plotted on geoMap
#   geoMap: the map containing k regions
#   leg: name to use on top of legend
#   colLim: control lower limit and upper limit of color scale (colLim = c(lowVal, highVal))
plotAreaCol = function(fName, width, height, estVal, geoMap, leg, colLim = NULL){
  if(is.null(colLim)){
    colLim = range(estVal)
  }
  
  # Set up data object for plotting
  nigeriaMapTmp = geoMap
  nigeriaMapTmp$MCV1 = estVal
  nigeria.df = merge(fortify(nigeriaMapTmp), as.data.frame(nigeriaMapTmp), by.x = "id", by.y = 0)
  nigeria.df$Longitude = nigeria.df$long
  nigeria.df$Latitude  = nigeria.df$lat
  
  # Plot
  map = ggplot() +
    geom_polygon(data = nigeria.df,
                 aes(x = Longitude, y = Latitude, group = group, fill = MCV1),
                 color = 'gray', size = .2)+
    scale_fill_viridis_c(direction = 1,
                         begin = 1,
                         end = 0,
                         limit = colLim,
                         name = leg) + 
    coord_fixed() + 
    theme(text = element_text(size=40),
          legend.key.height = unit(4, 'cm'),
          legend.key.width  = unit(1.75, 'cm'))
  ggsave(filename = fName,
         plot = map,
         width = width, 
         height = height)
}

