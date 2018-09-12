
# remove any previously loaded functions
rm(list=ls())

# suppress warnings
oldw <- getOption("warn")
options(warn = -1)

# Load the required packages
list.of.packages <- c("optparse", "Hmisc")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(optparse)
library(Hmisc)


# set arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.05, 
              help="input pvalue between 0 and 1", metavar="character"),
  make_option(c("-o", "--out"), type="character", 
              default="output_coexpressed_R.csv", 
              help="output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# To run interactively, uncomment the following two lines
setwd("C:\\Users\\alokum\\Desktop\\Menuscript\\PNAS\\Menuscript\\ADEX")
infile1 <- read.table(file = "all_aml_train.preprocessed.gct", header = T, sep = "\t", strip.white = T, quote = "\"", stringsAsFactors = F, row.names = 1)

# ~~~~~~~~~~~~~
# ~~ Find co-expressed genes by Pearson correlation, output table
# ~~~~~~~~~~~~~

# function to calculate the correlations, p-values and make bonferroni corrections
PearsonCorrelationTable <- function(infile, pvalue, outfile) {

  # read in data matrix
  infile <- read.csv(file = opt$input, header = T, sep = "\t", 
                     strip.white = T, quote = "\"", stringsAsFactors = F, row.names = 1)
  
  # pearson correlation excluding name columns using Hmisc
  r <- rcorr(t(infile[,2:ncol(infile)]), type = "pearson")
  
  # reformat matrix to two columns of interactions with values
  df <- flattenCorrMatrix(r$r, r$P)
  
  # correct for multiple hypothesis testing
  df$bonferroni <- p.adjust(df$p, method = "bonferroni")
  
  # calculate the t-statistic to include in the final output
  df$t.statistic <- df$cor / (sqrt( (1-(df$cor^2)/(ncol(r$r)-2) )))
  
  # subset the result to only significant correlations
  significants <- subset(df, bonferroni <= pvalue)
  
  # reformat output
  output <- data.frame(significants$row, significants$column, significants$cor, 
                      significants$p, significants$bonferroni)
  colnames(output) <- c("InteractorA", "InteractorB", "Correlation",
                        "p-value", "Bonferroni-adjusted p-value")
  output <- output[order(output$`p-value`),]
  
  # write to file
  write.csv(output, outfile, 
            row.names = F, quote = F)
}

# function to reformat the rcorr output matrix to two columns of interactions with values

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = t(cormat)[ut],
    p = pmat[ut]
  )
}


# Run the function given the user inputs 

PearsonCorrelationTable(opt$input, opt$pvalue, opt$out)

#Excluded from final version, include if using unpreprocessed data
# ~~~~~~~~~~~~~ 
# ~~ Apply low and high thresholds and remove values that never change
# ~~~~~~~~~~~~~
# This requires first cleaning the data. As done using GenePattern earlier in this assignment, values <20 will be set to 20, values > 20000 will be set to 20000, and values that change by less than 3 fold between any sample will be removed.
# 
# # set values lower than 20 to 20
# LowThreshold <- function(x) { x[x<20] <- 20; x }
# df.low <- as.data.frame(LowThreshold(input))
# 
# #set values higher than 20000 to 20000
# HighThreshold <- function(x) { x[x>20000] <- 20000; x }
# df.low.high <- as.data.frame(HighThreshold(df.low[3:40]))
# 
# # Remove rows that vary by less than 3 fold in any sample
# df.low.high$folddiff <- apply(df.low.high, 1, FUN = function(x) {max(x)/min(x)})
# df.desc <- cbind(df.low$Name, df.low$Description, df.low.high)
# df.clean <- df.desc[df.low.high$folddiff > 3,]
# df.clean$folddiff <- NULL

# return warning setting
options(warn = oldw)
