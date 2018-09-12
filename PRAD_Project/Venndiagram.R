#R-script for making both weighted as well as unweighted venn diagrams form upto 8 sets of gene lists.
#The gene list should be in .xlsx format with each condition occupying only one column)
#Duplicate entries in each sets will be ignored.
#Requires Perl and packages gdata, vennerable and gplots.
#To install vennerable run the following commands:
#   source("https://bioconductor.org/biocLite.R"); biocLite(c("RBGL","graph"))
#   install.packages("devtools"); library(devtools);
#   install_github("js229/Vennerable")
#Other packages can be installed from cran.

library (gdata)
library (Vennerable)
library (gplots)


#Reading the file-name and subsequently reading the gene list
#The file should be in the working directory
x <- readline("Enter the file name:") 
d <- paste(x,".xlsx",sep="")

geneLists <- read.xls(d, sheet = 1, stringsAsFactors = FALSE, header = FALSE)
head(geneLists)

#Reading the number of conditions i.e the number of sets to be plotted.

y <- readline("Enter the number of conditions:") #This module can plot upto 6 conditions

# There will be empty strings in genelist from empty cells in the excel file
tail(geneLists)
# Function to remove empty strings from the gene list 

removeEMPTYstrings <- function(x) {

    newVectorWOstrings <- x[x != ""]
    return(newVectorWOstrings)

    }


#Intialinzing the gene lists
#Note: 8 lists will be initialized by default but only those lists pointed out by the user will be used for plotting

genelist <- lapply(as.list(geneLists), removeEMPTYstrings)
genelist2 <- lapply(as.list(geneLists), removeEMPTYstrings)
genelist3 <- lapply(as.list(geneLists), removeEMPTYstrings)
genelist4 <- lapply(as.list(geneLists), removeEMPTYstrings)
genelist5 <- lapply(as.list(geneLists), removeEMPTYstrings)
genelist6 <- lapply(as.list(geneLists), removeEMPTYstrings)
genelist7 <- lapply(as.list(geneLists), removeEMPTYstrings)
genelist8 <- lapply(as.list(geneLists), removeEMPTYstrings)


#Array to assign the name of conditions to each list
con <- array(, dim = c(y, 1, 1)) #This array contains the names of the conditions
print ("Enter the name of conditions:")

z = 1
while (z <= y) {
  
    a <- "Condition"
    a <- paste(a,z,":",sep=" ")
    con[z, 1, 1] <- readline(a)
    z = z + 1

    }


# We rename pre-intialized data with the desired names here
names(genelist) <- c(con)


#Plotting the venn diagram and accesing the elements in each section of the diagram

v <- readline("Do you want weighted venn diagram?(yes/no): ")
if (v=="yes") {
  v = TRUE
}
if (v=="no") {
  v=FALSE
}
elementlist <- venn(genelist[1:y], show.plot=FALSE)

plot(Venn(Sets = genelist[1:y]), doWeight = v)


#Printing the elements list into a txt file
#The txt file will be saved in the current working directory

t <-  paste(x,".txt",sep="")

sink(file = t, append = TRUE, type = c("output", "message"), split = FALSE)
    print (elementlist)
sink()

