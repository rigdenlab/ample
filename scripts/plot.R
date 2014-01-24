

library(ggplot2)
#setwd("/home/jmht/Documents/test/CC/run3")
setwd("/Users/jmht/Documents/AMPLE/data/coiled-coils/ensemble")
data <- read.table(file="results_obj.csv",sep=',', header=T)
data$success <- as.numeric( data$shelxeCC >= 25 & data$shelxeAvgChainLength >= 10 )
data$success <- replace( data$success, is.na(data$success), 0 )

# Object to hold successes
#sdata <- data[ data$success == 1, ]

# Data to calculate
# NEED TO ADD TIMING DATA OF OVERALL AMPLE RUN
# * mean resolution of failing/succesful cases
mean(  data[ data$success == 1, ]$resolution ) # 1.949077
mean(  data[ data$success == 0, ]$resolution ) # 2.025158

# * mean of solvent content of failing/succesful cases
unique( data[ is.na( data$solventContent ),]$pdbCode ) # 1G1J 1KYC 1P9I 3CVF
mean( data[ data$success == 1, ]$solventContent )
mean( data[ data$success == 0, ]$solventContent )

# most common number of residues/proportion of target in search model
median( data$ensembleNumResidues ) # 30
median( data$ensemblePercentModel ) # 53

# max/min of number of models in ensembles
max( data$ensembleNumModels ) # 30
min( data$ensembleNumModels ) # 2

# Plots

# length distribution of successful search models both as number of residues 
# and as fraction of target chain length remaining in the search model.
#p <- ggplot(data=data[ data$success == 1, ], aes(ensemblePercentModel) )
#p + layer(geom="histogram") +

## Below for side-by-side but grid lines are off
#p + geom_bar( position = "dodge", binwidth = 5 ) +
#scale_fill_manual( values=c("#FF0000", "#00FF00"),
#		name="Success/Failure",
#		labels=c("Failure", "Success")
#		) +

scolour = "#FF0000"
fcolour = "#00FF00"

# Try showing success and failure together
p <-ggplot(data=data, aes(x=ensembleNumResidues, fill=factor(success) ) )
p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
		scale_fill_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
	xlab("Number of residues in ensemble") +
	ylab("Number of ensembles") +
	ggtitle("Number of residues in ensemble for successful and failing cases")
ggsave("ResiduesVsEnsemble2.png")

p <-ggplot(data=data, aes(x=ensemblePercentModel, fill=factor(success) ) )
p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
		scale_fill_manual( values=c("#FF0000", "#00FF00"),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Percent of residues in ensemble") +
		ylab("Number of ensembles") +
		ggtitle("Percentage of residues in ensemble for successful and failing cases")
ggsave("PercentVsEnsemble2.png")

#ÊComparison of the RIO figures for successful search models with the number of in-register 
# residues similarly defined.
# goodContacts, inregisterContacts
p <- ggplot(data=data[ data$symmatchOriginOk == "True", ], 
			aes(reforiginRMSD,
				inregisterContacts,
				color=factor(success)
				)
			)

			
odata = data[  data$floatingOrigin == 'False' | ( data$floatingOrigin == 'True' & data$csymmatchOriginOk=='True' ) , ]

x <-  odata[ odata$success == 1 & odata$numContacts == 0 , ,c("pdbCode","ensembleName") ]

#p <- ggplot(data=odata, aes(reforiginRMSD, goodContacts/fastaLength, color=factor(success)) )
p <- ggplot(data=odata, aes(reforiginRMSD, goodContacts, color=factor(success)) )
p + layer(geom="point")

# histrogram
#p <- ggplot(data=odata, aes(x=goodContacts, fill=factor(success)) )
p <- ggplot(data=odata, aes(x=nocatContacts, fill=factor(success)) )
p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 )

p <- ggplot(data=odata, aes(reforiginRMSD, inregisterContacts, color=factor(success)) )
p + layer(geom="point")


# As Fig ? shows, many successes are achieved with LLG <0

qplot(shelxeCC, resolution, data=data, colour=success)

p <- ggplot( data, aes(fastaLength, resolution, color=factor(success) ) )

#p + geom_point()
p + layer(geom="point")

p + layer(geom="point") + facet_wrap(~success)

# Can add colour to each layer separatly or to the aes added to the ggplot object
# in which case all layers inherit it
# geom_point(aes(color = factor(cyl)))

#p + xlab("Body Weight") + ylab("Total Hours Sleep") + ggtitle("Some Sleep Data")

q()

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
# Select top by selecting not duplicates on pdbCode
x <- x[ !duplicated(x$pdbCode), ]
# NB LOOK AT REORDER




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
