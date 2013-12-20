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
