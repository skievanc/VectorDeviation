#
#Functions for X components
#
#Sum of X north/south magnitudes
require(ggplot2)
require(gridExtra)
require(raster)
require(rgdal)

Xnsum = function(vdata){
  i=1
  Xnmsum = 0
  while(i<=length(vdata)){
    Xnmsum = Xnmsum+(vdata[i,1]*cos((vdata[i,2])/(180/pi)))
    i=i+1
  }
  return(as.double(Xnmsum))
}
#
#Sum of X east/west magnitudes
Xesum = function(vdata){
  i=1
  Xemsum = 0
  while(i<=length(vdata)){
    Xemsum = Xemsum+(vdata[i,1]*sin((vdata[i,2])/(180/pi)))
    i=i+1
  }
  return(as.double(Xemsum))
}
#
#Average of X magnitudes
#Calls on Xnsum and Xesum functions
XMbar = function(vdata){
  xbar=sqrt(((Xnsum(vdata))^2)+((Xesum(vdata))^2))/length(vdata)
  return(xbar)
}
#
#Average X direction
#Calls on Xnsum and Xesum functions
#If Xnsum and Xesum both return zero, will return NA
XDirAVG = function(vdata){
  XN=Xnsum(vdata)
  XE=Xesum(vdata)
  if(XN>0){
    return(atan(XE/XN)*(180/pi))
  }else if(XN<=0 & XE>0){
    return(90-atan(XN/XE)*(180/pi))
  }else if(XN<=0 & XE<0){
    return(-90-atan(XN/XE)*(180/pi))
  }else if(XN<0 & XE==0){
    return(180)
  }else{
    return(NA)
  }
}
#
##
#
#Functions for Y components
#
#Sum of Y north/south magnitudes
Ynsum = function(vdata){
  i=1
  Ynmsum = 0
  while(i<=length(vdata)){
    Ynmsum = Ynmsum+(vdata[i,3]*cos((vdata[i,4])/(180/pi)))
    i=i+1
  }
  return(as.double(Ynmsum))
}
#
#Sum of Y east/west magnitudes
Yesum = function(vdata){
  i=1
  Yemsum = 0
  while(i<=length(vdata)){
    Yemsum = Yemsum+(vdata[i,3]*sin((vdata[i,4])/(180/pi)))
    i=i+1
  }
  return(as.double(Yemsum))
}
#
#Average of Y magnitudes
#Calls on Ynsum and Yesum functions
YMbar = function(vdata){
  Ybar=sqrt(((Ynsum(vdata))^2)+((Yesum(vdata))^2))/length(vdata)
  return(Ybar)
}
#
#Average Y direction
#Calls on Ynsum and Yesum functions
#If Ynsum and Yesum both return zero, will return NA
YDirAVG = function(vdata){
  YN=Ynsum(vdata)
  YE=Yesum(vdata)
  if(YN>0){
    return(atan(YE/YN)*(180/pi))
  }else if(YN<=0 & YE>0){
    return(90-atan(YN/YE)*(180/pi))
  }else if(YN<=0 & YE<0){
    return(-90-atan(YN/YE)*(180/pi))
  }else if(YN<0 & YE==0){
    return(180)
  }else{
    return(NA)
  }
}
#
##
#
#Functions for final metrics
#
#Magnitude Mean Deviation
#Returns y magnitude average minus x magnitude average
#Calls on YMbar and XMbar functions

#'Magnitude Mean Deviation
#'
#'@description
#'Computes the mean deviation of magnitude between two vector data sets.
#'Uses the equation mean vector of \code{y} - mean vector of \code{x}.
#'
#'@usage
#'magnitude_md(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in \emph{Metrics that make a difference.}"
#'@export
magnitude_md = function(vdata){
  return(YMbar(vdata)-XMbar(vdata))
}
#
#Magnitude Quantity Component
#Returns absolute value of magnitude_md function
#'Magnitude Quantity Component
#'
#'@description
#'Computes the quantity component of magnitude deviation between two vector data sets.
#'Quantity component is calculated as the absolute value of magnitude mean deviation.
#'
#'@usage
#'magnitude_qc(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in \emph{Metrics that make a difference.}"
#'@export
magnitude_qc = function(vdata){
  return(abs(magnitude_md(vdata)))
}
#
#Returns the sum of absolute value difference averages
#Used in calculating allocation component and MAD
MagnitudeDiffAVG = function(vdata){
  i=1
  MagAVG=0
  while(i<=length(vdata)){
    MagAVG = MagAVG+(abs(vdata[i,3]-vdata[i,1])/length(vdata))
    i=i+1
  }
  return(as.double(MagAVG))
}
#
#Magnitude Allocation Component
#Returns greater of zero and MagnitudeDiffAVG function minus magnitude_qc function
#Calls MagnitudeDiffAVG and magnitude_qc functions
#'Magnitude Allocation Component
#'
#'@description
#'Computes the allocation component of magnitude deviation between two vector data sets.
#'Allocation component is calculated as the Mean Absolute Deviation minus the Quantity Component.
#'
#'@usage
#'magnitude_ac(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in \emph{Metrics that make a difference.}"
#'@export
magnitude_ac = function(vdata){
  AC=0
  AVGQC=MagnitudeDiffAVG(vdata)-magnitude_qc(vdata)
  if(AVGQC>0){
    AC=AVGQC
  }
  return(AC)
}
#
#Magnitude Mean Absolute Deviation
#Returns greater of MagnitudeDiffAVG and magnitude_qc functions
#'Mean Absolute Deviation of Magnitude
#'
#'@description
#'Computes the mean absolute deviation (MAD) of magnitude between two vector data sets.
#'Mean Absolute Deviation is calculated as the sum of Allocation and Quantity Components.
#'
#'@usage
#'magnitude_mad(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in \emph{Metrics that make a difference.}"
#'@export
magnitude_mad = function(vdata){
  if(magnitude_qc(vdata)>MagnitudeDiffAVG(vdata)){
    return(magnitude_qc(vdata))
  }else{
    return(MagnitudeDiffAVG(vdata))
  }
}
#
#Direction absolute value difference average
#Used in calculating direction allocation component and MAD

DirectionDiffAVG = function(vdata){
  i=1
  DirAVG=0
  if(!(vdata[i,1] == 0 | vdata[i,3] ==0)){
    while(i<=length(vdata)){
      DirAVG = DirAVG+(abs(vdata[i,2]-vdata[i,4])/length(vdata))
      i=i+1
    }
  }
  return(as.double(DirAVG))
}


#Direction Mean Deviation
#If X or Y avg direction is NA, returns zero
#Calls on XDirAVG and YDirAVG functions

#'Direction Mean Deviation
#'
#'Computes the mean deviation of direction between two vector data sets.
#'If X or Y avg direction is NA, returns zero.
#'
#'
#'@usage
#'direction_md(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in emph{Metrics that make a difference.}"
#'@export
direction_md = function(vdata){
  XD=XDirAVG(vdata)
  YD=YDirAVG(vdata)
  if(is.na(XD) | is.na(YD)){
    return(0)
  }else if(abs(YD-XD)<180){
    return(YD-XD)
  }else if((YD-XD)>180){
    return(YD-XD-360)
  }else if((YD-XD)<(-180)){
    return(YD-XD+360)
  }else if(abs(YD-XD)==180){
    return(180)
  }
}
#
#Direction Quantity Component
#Returns absolute value of direction_md function
#'Direction Quantity Component
#'
#'@description
#'Computes the quantity component of direction deviation between two vector data sets.
#'Quantity component is calculated as the absolute value of direction mean deviation.
#'
#'@usage
#'direction_qc(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in \emph{Metrics that make a difference.}"
#'@export
direction_qc = function(vdata){
  return(abs(direction_md(vdata)))
}
#
#Direction Allocation Component
#Returns greater of zero and DirectionDiffAVG function minus direction_qc function
#Calls DirectionDiffAVG and direction_qc functions
#'Direction Allocation Component
#'
#'@description
#'Computes the allocation component of direction deviation between two vector data sets.
#'Allocation component is calculated as the Mean Absolute Deviation minus the Quantity Component.
#'
#'@usage
#'direction_ac(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in \emph{Metrics that make a difference.}"
#'@export
direction_ac = function(vdata){
  AC=0
  AVGQC=DirectionDiffAVG(vdata)-direction_qc(vdata)
  if(AVGQC>0){
    AC=AVGQC
  }
  return(AC)
}
#
#Direction Mean Absolute Deviation
#Returns greater of DirectionDiffAVG and direction_qc functions
#'Mean Absolute Deviation of Direction
#'
#'@description
#'Computes the mean absolute deviation (MAD) of direction between two vector data sets.
#'Mean Absolute Deviation is calculated as the sum of Allocation and Quantity Components.
#'
#'@usage
#'direction_mad(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in \emph{Metrics that make a difference.}"
#'@export
direction_mad = function(vdata){
  if(direction_qc(vdata)>DirectionDiffAVG(vdata)){
    return(direction_qc(vdata))
  }else{
    return(DirectionDiffAVG(vdata))
  }
}

#
#




#
##
#Graphs
#
#
#'Mean Absolute Deviation Plot
#'
#'@description
#'Creates a stacked bar graph of quantity and allocation components of MAD
#'for both magnitude and direction of two vector variable sets.
#'@usage
#'mad_plot(vdata)
#'@param vdata a numeric data frame of four columns; in the format: x magnitude, x direction, y magnitude, y direction.
#'@details
#'Equations and method for this function are based on Robert Pontius' vector deviation equations in \emph{Metrics that make a difference.}"
#'@export
mad_plot = function(vdata){
  magDF = data.frame(
    "Type"=c("Magnitude","Magnitude"),
    "Component"=c("Quantity Component","Allocation Component"),
    "Value"=c(magnitude_qc(vdata),magnitude_ac(vdata)))

  magplot=ggplot(magDF, aes(fill =  Component, x = Type, y = Value)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip()+
    labs(y = "Mean Absolute Deviation" , x = "")

  dirDF = data.frame(
    "Type"=c("Direction","Direction"),
    "Component"=c("Quantity Component","Allocation Component"),
    "Value"=c(direction_qc(vdata),direction_ac(vdata)))

  dirplot=ggplot(dirDF, aes(fill =  Component, x = Type, y = Value)) +
    geom_bar(stat = "identity", position = "stack") +
    coord_flip()+
    labs(y = "Mean Absolute Deviation" , x = "")

  grid.arrange(magplot,dirplot,nrow=2)
}

#
##
#Section for using equations on maps and u/v components in general

#Turn a vector map into a table
#Input order mag1,dir1,mag2,dir2; output table with variables in same order
#'Raster to Table Input (vector)
#'
#'@description
#'Converts 4 raster files of 2 vector variable sets into table for use in deviation
#'calculations. Expects a magnitude and direction component raster for each vector set.
#'@usage raster_input_vector(m1,d1,m2,d2)
#'@param m1 raster of magnitude component of vector variable 1
#'@param d1 raster of direction component of vector variable 1
#'@param m2 raster of magnitude component of vector variable 2
#'@param d2 raster of direction component of vector variable 2
#'@details
#'Accepts any raster map file type supported by the raster package. Rasters must be
#'input in the order specified. Prepares data to be used in vector deviation functions.
raster_input_vector = function(m1,d1,m2,d2){
  mstack=stack(list(m1,d1,m2,d2))

}
#Turns 4 maps into table; still has u/v values; input order u1,v1,u2,v2
#For troubleshooting purposes
mapInputuv=function(m1u,m1v,m2u,m2v){
  mstack=stack(list(m1u,m1v,m2u,m2v))
  rdfMap=na.omit(as.data.frame(mstack))
  colnames(rdfMap) = c("u1","v1","u2","v2")
  return(rdfMap)
}
#Turns 4 maps into table with vector variables of mag1,dir1,mag2,dir2; input order u1,v1,u2,v2
#
#'Raster to Table Input (u/v component)
#'
#'@description
#'Converts 4 input maps of u and v components of wind into a table of vector components (direction and magnitude).
#'Prepares meteorological wind data from raster files for use in vector deviation calculations.
#'
#'@usage
#'raster_input_uv(m1u,m1v,m2u,m2v)
#'@param m1u raster of u component of variable 1
#'@param m1v raster of v component of variable 1
#'@param m2u raster of u component of variable 2
#'@param m2v raster of v component of variable 2
#'@details
#'Accepts any raster map file type supported by the raster package.
#'Conversion of u and v comonpent wind into vectors maintains and combines directions into
#'the direction attribute, and speed into the magnitude attribute.
raster_input_uv=function(m1u,m1v,m2u,m2v){
  mstack=stack(list(m1u,m1v,m2u,m2v))
  rdfMap=na.omit(as.data.frame(mstack))
  vMap=data.frame()
  x=1
  while(x<=nrow(rdfMap)){
    vMap[x,1]=sqrt((rdfMap[x,1]^2)+(rdfMap[x,2]^2))
    vMap[x,2]=180+((180/pi)*atan2(rdfMap[x,1],rdfMap[x,2]))
    vMap[x,3]=sqrt((rdfMap[x,3]^2)+(rdfMap[x,4]^2))
    vMap[x,4]=180+((180/pi)*atan2(rdfMap[x,3],rdfMap[x,4]))
    x=x+1
  }
  colnames(vMap) = c("xMag","xDir","yMag","yDir")
  return(vMap)
}
#Converts u/v table into vector variable table; from u1,v1,u2,v2 into mag1,dir1,mag2,dir2
#'Convert Table of u and v Components into Vector Components
#'
#'@description
#'Converts a table of u and v component winds of 2 variables into vector
#'components. Input table must be in form variable 1 u, variable 1 v,
#'variable 2 u, variable 2 v.
#'Prepares meteorological wind data from raster files for use in vector deviation calculations.
#'@usage
#'input_uv(data)
#'@param data table of u and v component of winds from 2 variables; must be in form variable 1 u, variable 1 v,
#'variable 2 u, variable 2 v.
#'@details
#'Conversion of u and v component wind into vectors maintains and combines directions into
#'the direction attribute, and speed into the magnitude attribute.
input_uv=function(data){
  vTab=data.frame()
  x=1
  while(x<=nrow(data)){
    vTab[x,1]=sqrt((data[x,1]^2)+(data[x,2]^2))
    vTab[x,2]=180+((180/pi)*atan2(data[x,1],data[x,2]))
    vTab[x,3]=sqrt((data[x,3]^2)+(data[x,4]^2))
    vTab[x,4]=180+((180/pi)*atan2(data[x,3],data[x,4]))
    x=x+1
  }
  colnames(vTab) = c("xMag","xDir","yMag","yDir")
  return(vTab)
}



