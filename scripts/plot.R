
setwd("/home/jmht/Documents/test/CC/run3")
data <- read.table(file="results_obj.csv",sep=',', header=T)
data$success <- as.numeric( data$shelxeCC >= 25 & data$shelxeAvgChainLength >= 10 )
data$success <- replace( data$success, is.na(data$success), 0 )

# Count # jobs
njobs <- aggregate( data$success, by=list( data$pdbCode ), FUN=length)[2]
# Count succcess
nsuccess <- aggregate( data$success, by=list( data$pdbCode ), FUN=sum)[2]
# Failed models
nfail <- njobs-nsuccess

# Get a data frame just with the successk
# Order gives us the indices of the frame rows ordered by the values we've given
# duplciated gives us an array of bools for which items are duplicated


#http://stackoverflow.com/questions/6289538/aggregate-a-dataframe-on-a-given-column-and-display-another-column
#http://stackoverflow.com/questions/2822156/select-rows-with-largest-value-of-variable-within-a-group-in-r

# For each case need to get the best success - i.e. success with max CC
# select success
x <- data[ data$success==1, ]
# Order by pdbCode and CC
x <- x[ order( x$pdbCode, x$shelxeCC, decreasing=TRUE ),   ]
# Select top by selecting not duplcates on pdbCode
x <- x[ !duplicated(x$pdbCode), ]


q()
png(filename='plot.png')

data <- read.table(file="results_obj.csv",sep=',', header=T)

# Fix missing values
#completeEnsemble <- complete.cases(data$ensembleNumResidues)
#completeHelix <- complete.cases(data$lenHelix)

#pdata <- subset( data, Shelxe.CC >= 25 & Shelxe.avg..chain.length >= 10 )
#pdata <- subset( data, shelxeCC >= 25 & shelxeAvgChainLength >= 10 )
pdata <-  data[ data$shelxeCC >= 25 & data$shelxeAvgChainLength >= 10, ]

dim( data )
dim( pdata )

pdbAll <- unique( data$pdbCode )
sprintf("Len all PDBs %s", length(pdbAll) )
pdbSuccess <- unique( pdata$pdbCode )
sprintf("Len success PDBs %s", length(pdbSuccess) )

q()


par(bg="white")
amax=max( na.omit(pdata$Ensemble.num.residues), na.omit( pdata$Helix.length ) )
plot( pdata$Ensemble.num.residues, pdata$Helix.length,
xlim=c(0, amax ), ylim=c(0, amax ),
xlab="Ensemble Num Residues", ylab="Helix length"  )

#guideline
abline( 0, 1, col="red" )

#main="Num residues vs helix length",
title("Num residues vs helix length")

dev.off()
