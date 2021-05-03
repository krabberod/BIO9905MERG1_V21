# This is a simple R script.
# Consider using R Notebook if you want to do markdown style commenting.
# Use hashtags to make comments. R will ignore everything after hashtag.
#' Hastags with a single quote will be written as "normal text" when knitted (more on that later)
#
# Make comments! Lots of comments to your script. It's important to document code
# Both so you can remember what you are doing, and so others can understand you what you have done
# https://medium.com/@andrewgoldis/how-to-document-source-code-responsibly-2b2f303aa525
# Execute lines by pressing cmd+enter (mac), or ctrl+enter (windows)

# To check the working directory:
getwd()

# Changing the working directory
setwd(dir = "./../Setup/")

# Save workspace

setwd(dir = "./../Intro_to_Rstudio//")
#### Several hashtags makes a section  ####
library(ggplot2)

#save.image("Some_data.Rdata")


#### Start a new section ####
#' do some stats:   
#' Example from http://r-statistics.co/Top50-Ggplot2-Visualizations-MasterList-R-Code.html  
#' Notice the different comment-mark?   
#' It will make a difference when writing a report (aka. knit the document)
#' A package called knitr exists, but will not be covered in this course. 
#'   
#'      
theme_set(theme_bw())  

# Data Prep
data("mtcars")  # load data
mtcars$`car name` <- rownames(mtcars)  # create new column for car names
mtcars$mpg_z <- round((mtcars$mpg - mean(mtcars$mpg))/sd(mtcars$mpg), 2)  # compute normalized mpg
mtcars$mpg_type <- ifelse(mtcars$mpg_z < 0, "below", "above")  # above / below avg flag
mtcars <- mtcars[order(mtcars$mpg_z), ]  # sort
mtcars$`car name` <- factor(mtcars$`car name`, levels = mtcars$`car name`)  # convert to factor to retain sorted order in plot.

# Diverging Barcharts
ggplot(mtcars, aes(x=`car name`, y=mpg_z, label=mpg_z)) + 
  geom_bar(stat='identity', aes(fill=mpg_type), width=.5)  +
  scale_fill_manual(name="Mileage", 
                    labels = c("Above Average", "Below Average"), 
                    values = c("above"="#00ba38", "below"="#f8766d")) + 
  labs(subtitle="Normalised mileage from 'mtcars'", 
       title= "Diverging Bars") + 
  coord_flip()

#' ### NOTE
#' the section header  in an R-script does not work will with making a simple report
#' If you want to make "better" reports make a Rnotebook file instead of a script
#' It depends on what functionallity you want...
#' 

