# Change Detection for Land Cover Mapping using Random Forest applied to Landsat 5 and Landsat 7 imagery.

This notebook provides a workflow to detect and classify the land cover changes in a preprocessed multi date image at 30m resolution. The model will detect the following Land Cover Changes: Water to Water, Forest to Non-Forest, Non-Forest to Forest and Forest to Forest.

## Setting up the environment 

###### Importing packages and checking versions
```
#Loading the required libraries
require(maptools) #Tools for handling spatial data
require(sp) #Methods for dealing with spatial data
require(randomForest) #implements Breiman's random forest algorithm for classification
require(raster) #Methods to create a RasterLayer object
require (rgdal) #permit spatial data to be associated with coordinate reference systems
```
###### Setting working directory
```
# Set working directory
setwd("~/Documents/Data/Notebook/gis_natgeo-master/test5")
# Name and path for the Shapefile (don't need the .shp extension)
shapefile <- "~/Documents/Data/Notebook/gis_natgeo-master/test5/training_.shp"
# Class numbers that you want to select training sample from
classNums <- c(11,12,21,22)
# For each land cover class the approximate number of training samples to be randomly selected 
# If a value is "0" then all pixels in all of the polygons for that class will be used 
classSampNums <- c(100,100,100,100)
# Name of the attribute that holds the integer land cover type identifer
attName <- 'id'
# No-data value for the input image
nd <- -9999
# Name and path for the input satellite image
inImageName <-"~/Documents/Data/Notebook/gis_natgeo-master/test5/mfnp_2014_2019_stacked_upper.tif"
```
###### Setting the probability threshhold
```
# Note that if this file exists the write will fail with the message "Creation of output file failed"  
outMarginFile <- 'margin.shp'
# Output classification image (enter TRUE or FALSE)
classImage <- TRUE
# Output probability image layer (enter TRUE or FALSE)
probImage <- FALSE
# Output classification layer and set pixels with probability less than "probThreshold" to 0 (enter TRUE or FALSE)
threshImage <- FALSE
# Enter threshold probability in percent (values must be between 0 and 100) only used if threshImage=TRUE
probThreshold <- 75
# Layer number (band number) for the X and Y axis of the feature space plot. 
# If you do not want to calculate a feature plot enter 0 as the layer number
xBand <- 2
yBand <- 3
```
## Start Processing

###### This step also loads the shapefile with the training polygons

```
startTime <- Sys.time()
cat("Start time", format(startTime))

# Read the Shapefile
#vec <- readShapePoly(shapefile): This was the old code, and was updated for the new versions of the R pacakges 
vec<-readOGR(dsn = ".", layer = "training_")
```
```
# Load the image then flag all no-data values(nd) so they are not processed
satImage <- brick(inImageName)
NAvalue(satImage) <- nd
for (b in 1:nlayers(satImage)) { NAvalue(satImage[[b]]) <- nd }
```
## Create a vector of unique land cover types

```
# Create vector of unique land cover attribute values
allAtt <- vec@data
tabAtt <-table(allAtt[[attName]])
uniqueAtt <-as.numeric(names(tabAtt))
```
```
# Check if lenght of classNums and classSampNums is equal
if (length(classNums) != length(classSampNums)) {
  cat("length of classNums and classSampNums no equal")
  stop("Check the classNums and classSampNums variable", call.=FALSE)
}

# Check if all classNums exist in uniqueAtt
    #### CHECK THIS FUNCTION TO SEE IF classNums ARE IN uniqueAtt  ################
if (sum(classNums %in% uniqueAtt) != length(uniqueAtt)) {
  cat("not all classes in classNums are defined in the vecotr file")
  stop("Check classNums and vector attribute table", call.=FALSE)
}
    
# Create input data from a Shapefile using all training data 
cat("Create training data using all pixels in training polygons")
predictors <- data.frame()
response <- numeric()
xyCoords <- data.frame()
```
## Create training data to train model
```
cat("Create training data to train model")
# If all pixels in a polygon are to be used process this block
for (n in 1:length(classNums)) {
  if (classSampNums[n] == 0) {
    # Get the metadata for all polygons for this particular class
    class_data<- vec[vec[[attName]]==classNums[n],]
    # Extract and combine predictor and response variables for each polygon within a class
    for (i in 1:dim(class_data)[1]) {
      satValues <- extract(satImage, class_data[i,], cellnumbers=TRUE, df=TRUE)
     ## satValues <- as.data.frame(do.call(rbind,satValues))
      attributeVector <- rep.int(classNums[n],nrow(satValues))
      xyCoords <- rbind(xyCoords, xyFromCell(satImage, satValues[,2]))
      predictors <- rbind(predictors, satValues[,-1:-2])
      response <- c(response, attributeVector)
      
    }
  } else {
    # Create input data from a Shapefile by sampling training data polygons
    # Get the metadata for all polygons for a particular class (based on the uniqueAtt variable)
    class_data<- vec[vec[[attName]]==classNums[n],]
    # Get the area of each polygon for a particular class
    areas <- sapply(slot(class_data, "polygons"), slot, "area")
    # Calculate the number of samples for each polygon based on the area in proportion to total area for a class
    nsamps <- ceiling(classSampNums[n]*(areas/sum(areas)))
    # Use random sampling to select training points (proportial based on area) from each polygon for a given class 
    for (i in 1:dim(class_data)[1]) {
      xy_class <- spsample(class_data[i,], type="random", n=nsamps[i])
      # Add coordinates to create a list of random points for all polygons
      if (i == 1) cpts <- xy_class
      else cpts <- rbind(cpts, xy_class)
    }
    # The number of points might not match numsamps exactly.
    xy_ForClass <- cpts
    xyCoords <- rbind(xyCoords, xy_ForClass@coords)

  # Get class number for each sample point for responce variable
  response <- c(response, over(xy_ForClass, vec)[[attName]])
  # Get pixel DNs from the image for each sample point
  predictors <- rbind(predictors, extract(satImage, xy_ForClass))
  }
}
    
trainvals <- cbind(response, predictors)    
```
## Test if feature space plot is needed, if yes; a feature space plot will be plotted

```
# Test if feature space plot is needed
if (xBand != 0 & yBand != 0) {
  #Plot feature space and samples
  continue <- "c"
  while (continue == "c") {
    plotImage <- stack(satImage[[xBand]], satImage[[yBand]])
    # Get pixel values from the image under each sample point and create a table with 
    # observed and predicted values
    cat("Getting pixel values to create feature space plot")
    featurePlotPoints <- sampleRegular(plotImage,100000 )
  
    # Remove NA values from trainvals table created above
    featurePlotPoints <- na.omit(featurePlotPoints)
  
    minBand1 <- min(featurePlotPoints[,1])
    maxBand1 <- max(featurePlotPoints[,1])
    minBand2 <- min(featurePlotPoints[,2])
    maxBand2 <- max(featurePlotPoints[,2])
    rangeBand1 <- maxBand1 - minBand1 + 1
    rangeBand2 <- maxBand2 - minBand2 + 1
  
    xAxisLabel <- paste("Layer", xBand, sep=" ")
    yAxisLabel <- paste("Layer", yBand, sep=" ")
  
    plot(featurePlotPoints[,1], featurePlotPoints[,2], col="lightgrey", xlab=xAxisLabel, ylab=yAxisLabel)
  
    uniqueValues <- unique(trainvals[,1])
    for (v in 1:length(uniqueValues)) {
      points(trainvals[which(trainvals[,1]==uniqueValues[v]), xBand+1], trainvals[which(trainvals[,1]==uniqueValues[v]), yBand+1], col=v, pch=20)
    }
  
    legend(minBand1, maxBand2, col=1:v, pch=20, title="Classes", legend=as.character(uniqueValues))
  
    continue <- readline(prompt="Type n to stop, c to change feature space bands, s to define a rectangle to locate gaps in feature space, or any other key to continue with randome forests model creation and prediciton, type the key in the console and press enter: ")
  
    if (substr(continue, 1,1) == "n") {
      stop("Processing stopped at users request", call.=FALSE)
    }
    if (substr(continue, 1,1) == "s") {
      cat("Click two points to define the area on the feature space plot that you want to highlight")
      coords <- locator(n=2)
      coords <- unlist(coords)
      xvals <- coords[1:2]
      yvals <- coords[3:4]
      
      # Print out the corner coordinates for the rectangle
      cat("min X =", min(xvals), "\n")
      cat("max X =", max(xvals), "\n")
      cat("min y =", min(yvals), "\n")
      cat("max y =", max(yvals), "\n")
      
      # Draw the rectangle on the feature space plot
      rectangle <- matrix(nrow=5, ncol=2)
      rectangle[1,] <- c(min(xvals), max(yvals))
      rectangle[2,] <- c(max(xvals), max(yvals))
      rectangle[3,] <- c(max(xvals), min(yvals))
      rectangle[4,] <- c(min(xvals), min(yvals))
      rectangle[5,] <- c(min(xvals), max(yvals))
      lines(rectangle[,1], rectangle[,2])
      
      # Get the bands used to calculate the feature space plot
      b1 <- raster(plotImage, layer=1)
      b2 <- raster(plotImage, layer=2)
      
      # Threshold satImage so all values selected in the rectangle on the feature space plot are set to 255
      satImage[(b1 > min(xvals)) & (b1 < max(xvals)) & (b2 > min(yvals)) & (b2 < max(yvals))] <- 255
      
      
      cat("White pixels in the plotted image were selected in the rectangle drawn on the feature space plot")
      stop("Add new training data and re-run the script", call.=FALSE)
    }
    if (substr(continue, 1,1) == "c") {
      xBand <- as.numeric(readline(prompt="Enter the band number for the x axis: "))
      yBand <- as.numeric(readline(prompt="Enter the band number for the y axis: "))
    }
  }
}

# Remove NA values 
trainvals <- na.omit(trainvals)
```
## Verify that shapefile with training polygon and input image are in the same projection

```
## Check to make sure Shapefile and input image are in the same projection
###### Then run the random forest model and start predictions. The output rasters will be created in this block of code

if (nrow(trainvals) == 0) {
  cat("No training data found")
  stop("It is possible the projection of the Shapefile with training data and input image are different </p>Check projections and run again", call.=FALSE)
}

# Run Random Forest
cat("Calculating random forest object")
randfor <- randomForest(as.factor(response) ~., data=trainvals, importance=TRUE, na.action=na.omit)

# Start predictions
cat("Starting predictions")
# Calculate the image block size for processing
bs <- blockSize(satImage)

extensionName <- unlist(strsplit(inImageName, "\\."))[length(unlist(strsplit(inImageName, "\\.")))]
outFileBaseName <- unlist(strsplit(inImageName, paste("\\.", extensionName, sep="")))[1]

# Create the output rasters
if (classImage) {
  outClassImage <- raster(satImage)
  outClassImage <- writeStart(outClassImage, filename=paste(outFileBaseName, "_Class.tif", sep=""), navalue=0, progress='text', format='GTiff', datatype='INT1U', overwrite=TRUE)
}
if (probImage) {
  outProbImage <- raster(satImage)
  outProbImage <- writeStart(outProbImage, filename=paste(outFileBaseName, "_Prob.tif", sep=""), navalue=0, progress='text', format='GTiff', datatype='INT1U', overwrite=TRUE)
}
if (threshImage) {
  outThreshImage <- raster(satImage)
  outThreshImage <- writeStart(outThreshImage, filename=paste(outFileBaseName, "_Thresh.tif", sep=""), navalue=0, progress='text', format='GTiff', datatype='INT1U', overwrite=TRUE)
}

# Loop though each of the image blocks to calculate the output layers selected in the variables section
for (i in 1:bs$n) {
  cat("processing block", i, "of", bs$n, "\r")
  imageBlock <-  getValuesBlock(satImage, row=bs$row[i], nrows=bs$nrows[i])
  predValues <- predict(randfor, imageBlock, type='response')
  classValues <- as.numeric(levels(predValues))[predValues]
  
  if (classImage) {
    #outClassMatrix <- matrix(classValues, nrow=nrow(imageBlock), ncol=1)
    outClassImage <- writeValues(outClassImage, classValues, bs$row[i])
  }
  if (probImage || threshImage) { 
    predProbs <- as.data.frame(predict(randfor, imageBlock, type='prob'))
    maxProb <- round(apply(predProbs, 1, max) * 100)
    if (probImage) { 
      #outProbMatrix <- matrix(maxProb, nrow=nrow(imageBlock), ncol=1)
      outProbImage <- writeValues(outProbImage, maxProb, bs$row[i])
    }
    if (threshImage) {
      threshValues <- classValues
      threshValues[which(maxProb <= probThreshold)] <- 0
      #outThreshMatrix <- matrix(threshValues, nrow=nrow(imageBlock), ncol=1)
      outThreshImage <- writeValues(outThreshImage, threshValues, bs$row[i])
    }
  }
}

# Stop writing and close the file
if (classImage) {
  outClassImage <- writeStop(outClassImage)
}
if (probImage) {
  outProbImage <- writeStop(outProbImage)
}
if (threshImage) {
  outThreshImage <- writeStop(outThreshImage)
}

```

## Plot output images

###### One of three possible images will be ploted if the if/else satement is true. The three images are probability image, threshold image or class image.

```
#Plot the image output if probImage is True then the outProbImage is plotted, 
#else if ThreshImage is True then the outThreshImage is plotted 
#and if ClassImage is True then the outclassImage is plotted

if (probImage==TRUE) {
  plot(outProbImage)
  print("This is a ProbImage")
} else if (threshImage==TRUE) {
  plot(outThreshImage)
  print("This is a ThreshImage")
}  else if(classImage==TRUE) {
  plot(outClassImage)
  print("This is a classImage")
}
```
## Print error rate and confusion matrix for the classification

```
# Print error rate and confusion matrix for this classification
confMatrix <- randfor$confusion
cat("OOB error rate estimate\n", 1 - (sum(diag(confMatrix)) / sum(confMatrix[,1:ncol(confMatrix)-1])), "%\n\n", sep="")
print(randfor$confusion)
```
```{r Creating the marginfile}
if (outMarginFile != "") {
  # Calculate margin (proportion of votes for correct class minus maximum proportion of votes for other classes)
  marginData <- margin(randfor)
  trainingAccuracy <- cbind(marginData[order(marginData)], trainvals[order(marginData),1])
  
  
  # Add column names to attributes table
  colnames(trainingAccuracy) <- c("margin", "classNum")  
  # Order X and Y coordinates 
  xyCoords <- xyCoords[order(marginData),]
  
  # Create and write point Shapefile with margin information to help improve training data
  pointVector <- SpatialPointsDataFrame(xyCoords, as.data.frame(trainingAccuracy), coords.nrs = numeric(0), proj4string = satImage@crs)
  writeOGR(pointVector, outMarginFile, "layer", driver="ESRI Shapefile", check_exists=TRUE)
}
```
```{r Ploting Variable of Importance}
# Plotting variable importance plot
varImpPlot(randfor)
```
## Print computation time

```
# Calculate processing time
timeDiff <- Sys.time() - startTime
cat("\nProcessing time", format(timeDiff), "\n")
```







