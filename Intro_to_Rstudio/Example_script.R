# This is a simple R scritp.
# Consider using R Notebook if you want to do markdown style commenting.
# Use hastags to make comments. R will ingore everything after hashtag
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
save.image("Some_data.Rdata")
