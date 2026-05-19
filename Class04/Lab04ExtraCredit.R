# To get and download a table or data online open it and read it
library(ggplot2)
source("http://thegrantlab.org/misc/cdc.R")
View(cdc)

# Head print out only 6 first entries of the table with default
head(cdc$height)

# Tail print out the last 6 entries by default but here we choose last 20 values of wieght column
tail(cdc$weight, n=20)

# Make a scatterplot of height vs weight
# We can use just plot or ggplot
# ggplot is easier and good for large tables like this


ggplot(data=cdc, mapping = aes(x = height, y = weight)) + geom_point()

# To find the correlation of the height vs weight we use the cor()
 cor(cdc$height, cdc$weight)

# to convert height in to m we simply in a new vecotr tome it to 0.0254
 height_m <- cdc$height * 0.0254
# To convert weight to kg we simply times 0.454
 weight_kg <- cdc$weight * 0.454

# BMI is calculated as weight in kg devided by height in meters squared 
BMI <- (weight_kg)/(height_m^2)
 plot(cdc$height, BMI)

 # correlation between the Height and BMI
 cor(cdc$height, BMI) 
 
 #Q8
 # Now to find patient with BMI of 30 or above
 head(BMI >=30, n=100)

 sum(BMI >=30)

 # To find the proportion of obese individuals we can use
 sum(BMI>=30)/length(BMI)
# or the precent value
 sum(BMI>=30)/length(BMI) *100
# round it to 1 decimal place
 round(sum(BMI>=30)/length(BMI) *100, 1)
 
 #Use the bracket notion to make a scatterplot of height and weight for the first 100 respondents
 plot(cdc[1:100, "height"], cdc[1:100, "weight"], xlab = "Height", ylab = "Weight", main = "Height vs Weight (First 100 Respondents)")

#Q10
 #How many obese individuals are male in the full dataset
 
g




