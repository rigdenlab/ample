

library(ggplot2)
scolour = "#FF0000"
fcolour = "#00FF00"
#setwd("/home/jmht/Documents/test/CC/run3")
setwd("/Users/jmht/Documents/AMPLE/data/coiled-coils/ensemble")
data <- read.table(file="results_bucc.csv",sep=',', header=T)
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
mean( data[ data$success == 1, ]$solventContent, na.rm=TRUE ) # 46.52952
mean( data[ data$success == 0, ]$solventContent, na.rm=TRUE ) # 51.3326

# most common number of residues/proportion of target in search model
median( data$ensembleNumResidues ) # 30
median( data$ensemblePercentModel ) # 53

# max/min of number of models in ensembles
max( data$ensembleNumModels ) # 30
min( data$ensembleNumModels ) # 2

# Most frequent number of models in a succesful ensemble
median( data[ data$success == 1, ]$ensembleNumResidues ) # 27
median( data[ data$success == 1, ]$ensemblePercentModel ) # 53

# Overall quality of the ab inito models
ensembleNativeTM
ensembleNativeRMSD


# Plots

# length distribution of successful search models both as number of residues 
# and as fraction of target chain length remaining in the search model.

# Try showing success and failure together
p <-ggplot(data=data, aes(x=ensembleNumResidues, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
	xlab("Number of residues in ensemble") +
	ylab("Number of ensembles") +
	ggtitle("Number of residues in ensemble for successful and failing cases")
ggsave("ResiduesVsEnsemble.png")

p <-ggplot(data=data, aes(x=ensemblePercentModel, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c("#FF0000", "#00FF00"),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Percent of residues in ensemble") +
		ylab("Number of ensembles") +
		ggtitle("Percentage of residues in ensemble for successful and failing cases")
ggsave("PercentVsEnsemble.png")

# CC vs final rFree coloured by success
# Need to remove nolog values from buccFinalRfree

# Horrible hack to get Rfree back to a number...
data$buccFinalRfree[ data$buccFinalRfree == "nolog" ] <- NA
data$buccFinalRfree <- as.numeric( as.character(data$buccFinalRfree) )

#x <-  data[ data$buccFinalRfree != 'nolog' & ! is.na( data$shelxeCC  ),
#		names(data) %in% c("shelxeCC","buccFinalRfree", "success") ]
# old skool...
#plot(x$buccFinalRfree,
#		x$shelxeCC,
#		pch=20,
#		cex=1,
#		col=c('red','green')[x$success+1],
#		xlab="Buccaneer Rfree",
#		ylab="ShelxeCC"  )
#write.csv(x, "x.csv", row.names=FALSE)

# another way...
#qplot( shelxeCC, buccFinalRfree, data=x, colour=factor(success) )

p <-ggplot(data=data, aes(shelxeCC, buccFinalRfree, colour=factor(success) ) )
p + geom_point() +
		scale_colour_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Shelxe CC score") +
		ylab("Final Buccaneer Rfree") +
		ggtitle("Shelxe vs Buccanner RFree score")
ggsave("CCVsRfree.png")

# Comparison of model quality
p <-ggplot(data=data, aes(x=ensembleNativeTM, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 0.01 ) +
		scale_fill_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("TM score of model to native") +
		ggtitle("Histogram of TM score for models for successful and failing cases")
ggsave("ModelTM.png")

p <-ggplot(data=data, aes(x=ensembleNativeRMSD, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("RMSD score of model to native") +
		ggtitle("Histogram of RMSD score for models for successful and failing cases")
ggsave("ModelRMSD.png")


# Plot of Reforigin RMSD for successes and failures
p <-ggplot(data=data, aes(x=reforiginRMSD, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("reforigin RMSD of placed model to native") +
		ggtitle("Histogram of reforigin RMSD scores by success")
ggsave("reforiginRMSD.png")

#ÊComparison of the RIO figures for successful search models with the number of in-register 
# residues similarly defined.
rdata = data[ data$floatingOrigin == 'False', ]

p <-ggplot(data=rdata, aes(x=nrGoodContacts, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
scale_fill_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Number of good contacts.") +
		ggtitle("Histogram of \"good\" contacts for non-floating origins.")
ggsave("goodContacts.png")

p <-ggplot(data=rdata, aes(x=nrNumContacts - nrGoodContacts, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("All contacts - good contacts ") +
		ggtitle("Histogram of \"uncategorised\" contacts for non-floating origins.")
ggsave("badContacts.png")

p <-ggplot(data=rdata, aes(x=nrGoodContacts / (nrNumContacts - nrGoodContacts), fill=factor(success) ) )
p + geom_histogram( position = 'dodge' ) +
		scale_fill_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Good contacts / \"uncateogorised\" contacts ") +
		ggtitle("Histogram of ratio of good/\"uncategorised\" contacts for non-floating origins.")
ggsave("contactsRatio.png")

p <-ggplot(data=rdata, aes(x=nrGoodContacts, y=nrInRegisterContacts, colour=factor(success) ) )
p + geom_point() +
		scale_colour_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Good Contacts") +
		ylab("In-register Contacts") +
		ggtitle("Good contacts vs in-register contacts for non-floating origins")
ggsave("goodVsInRegister.png")


# plot of nrInRegisterContacts vs nrGoodContacts
p <-ggplot(data=rdata, aes(x=shelxeCC, y=nrGoodContacts, colour=factor(success) ) )
p + geom_point() + scale_colour_manual( values=c(scolour, fcolour),
		name="Success/Failure",
		labels=c("Failure", "Success")
) +
xlab("Shelxe CC") +
ylab("Good Contacts") +
ggtitle("Good contacts vs shelxe CC for non-floating origins")
ggsave("CCVsInRegister.png")


# TFZ/LLG
p <-ggplot(data=data, aes(x=phaserTFZ, y=phaserLLG, colour=factor(success) ) )
p + geom_point() + scale_colour_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Phaser TFZ") +
		ylab("Phaser LLG") +
		ggtitle("Phaser LLG vs Phaser TFZ")
ggsave("LLGvsTFZ.png")

p <-ggplot(data=data, aes(x=phaserTFZ, y=phaserLLG, colour=factor(success) ) )
p + geom_point() +
		scale_colour_manual( values=c(scolour, fcolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		scale_y_continuous( limits=c(-3000, max(data$phaserLLG, na.rm=TRUE) )  )+
		xlab("Phaser TFZ") +
		ylab("Phaser LLG") +
		ggtitle("Phaser LLG vs Phaser TFZ")
ggsave("LLGvsTFZtrunc.png")


write.csv(odata, "rdata.csv", row.names=FALSE)

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
