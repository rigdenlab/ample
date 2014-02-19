

library(ggplot2)
scolour="#3333FF"
fcolour="#FF0000"
setwd("/Users/jmht/Documents/AMPLE/data/coiled-coils/ensemble")
#setwd("/home/jmht/Documents/test/CC/contacts")
#data <- read.table(file="results_bucc.csv",sep=',', header=T)
data <- read.table(file="results_bucc.csv",sep=',', header=T)

# Categorise successes
data$success <- as.numeric( data$shelxeCC >= 25 & data$shelxeAvgChainLength >= 10 )
data$success <- replace( data$success, is.na(data$success), 0 )

# Need to remove nolog values from buccFinalRfree
# Horrible hack to get Rfree back to a number...
data$buccFinalRfree[ data$buccFinalRfree == "nolog" ] <- NA
data$buccFinalRfree <- as.numeric( as.character(data$buccFinalRfree) )

# Calculate number of copies of decoy that were placed
data$numPlacedCopies <- data$numPlacedAtoms / data$ensembleNumAtoms

# Need to select the "best" origin and collate the properties

# CAN ONLY COMPARE NON_FLOATING ORIGINS and where a phaser model was produced and where we could find an origin
odata = data[ data$floatingOrigin == 'False' & !is.na( data$phaserLLG), ]
# HACK - set NA aoNumContacts to zero
#odata$aoNumContacts[ is.na( odata$aoNumContacts ) ] <- 0

# Find where origins not same
x <- odata[ ! ( is.na( odata$aoOrigin ) & is.na( odata$roOrigin ) ) & 
			! ( odata$aoOrigin == "" | odata$roOrigin == "" ) & 
				odata$aoOrigin != odata$roOrigin,
		c("aoOrigin","roOrigin", "pdbCode", "ensembleName") ]

aoOrigin
aoNumContacts
aoNumRio
aoRioInregister
aoRioOoRegister
aoRioBackwards
aoRioGood
aoRioNocat
roOrigin
roNumContacts
roNumRi
roRioInregister
roRioOoRegister
roRioBackwards
roRioGood
roRioNocat




# Object to hold successes
#sdata <- data[ data$success == 1, ]

# * mean resolution of failing/succesful cases
x <- mean(  data[ data$success == 1, ]$resolution ) # 1.949077
cat( "Mean resolution for successes is ",x,"\n")
x <- mean(  data[ data$success == 0, ]$resolution ) # 2.025158
cat( "Mean resolution for failures is ",x,"\n")

# * mean of solvent content of failing/succesful cases
#unique( data[ is.na( data$solventContent ),]$pdbCode ) # 1G1J 1KYC 1P9I 3CVF
x <-mean( data[ data$success == 1, ]$solventContent, na.rm=TRUE ) # 46.52952
cat( "Mean solvent content for successes is ",x,"\n")
x <-mean( data[ data$success == 0, ]$solventContent, na.rm=TRUE ) # 51.3326
cat( "Mean solvent content for failures is ",x,"\n")

# most common number of residues/proportion of target in search model
median( data$ensembleNumResidues ) # 30
median( data$ensemblePercentModel ) # 53

# max/min of number of models in ensembles
max( data$ensembleNumModels ) # 30
min( data$ensembleNumModels ) # 2

# Most frequent number of models in a succesful ensemble
median( data[ data$success == 1, ]$ensembleNumResidues ) # 27
median( data[ data$success == 1, ]$ensemblePercentModel ) # 53


# Plots

# Timings
# proportion for all runs where there was a shelxe build

# For some phaser runs the logs are not complete so we have no timing data
#sdata = data[ data$shelxeTime > 0 & ! is.na( data$phaserTime ), ]
#mean( sdata$fragmentTime, na.rm=TRUE )
#tlabels = c("fragmentTime", "modelTime", "ensembleTime", "phaserTime", "shelxeTime")
#aggregate( sdata[ , tlabels ], by=list( pdbCode = sdata$pdbCode ), FUN=mean)

#x <- data.frame( mean(sdata$fragmentTime), mean(sdata$modelTime), mean(sdata$ensembleTime), mean(sdata$phaserTime), mean(sdata$shelxeTime) )
# horrible...
#png("allSuccessTimingsPie.png")
#pie( as.numeric( x ), labels=tlabels, main="Timings for all runs inc. shelxe" )
#dev.off()

# length distribution of successful search models both as number of residues 
# and as fraction of target chain length remaining in the search model.

# Distribution of all ensembles
p <-ggplot(data=data, aes(x=ensembleNumResidues, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Number of residues in ensemble") +
		ylab("Number of ensembles") +
		ggtitle("Number of residues in ensemble for successful and failing cases")
ggsave("ResiduesVsEnsemble.png")



# Try showing success and failure together
p <-ggplot(data=data, aes(x=ensembleNumResidues, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_colour_manual( values=c(fcolour,scolour),
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
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Percent of residues in ensemble") +
		ylab("Number of ensembles") +
		ggtitle("Percentage of residues in ensemble for successful and failing cases")
ggsave("PercentVsEnsemble.png")

# CC vs final rFree coloured by success
p <-ggplot(data=data, aes(shelxeCC, buccFinalRfree, colour=factor(success) ) )
p + geom_point() +
		scale_colour_manual( values=c(fcolour,scolour),
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
		scale_colour_manual( values=c(fcolour,scolour),
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
		scale_colour_manual( values=c(fcolour,scolour),
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
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("reforigin RMSD of placed model to native") +
		ggtitle("Histogram of reforigin RMSD scores by success")
ggsave("reforiginRMSD.png")

#Comparison of the RIO figures for successful search models with the number of in-register 
# residues similarly defined.

p <-ggplot(data=odata, aes(x=nrGoodContacts, y=nrInRegisterContacts, colour=factor(success) ) )
p + geom_point() +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Good Contacts") +
		ylab("In-register Contacts") +
		ggtitle("Good contacts vs in-register contacts for non-floating origins")
ggsave("goodVsInRegister.png")


p <-ggplot(data=odata, aes(x=nrInRegisterContacts, y=nrOoRegisterContacts, colour=factor(success) ) )
p + geom_point( size=1.5 ) + stat_sum() +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("In-register Contacts") +
		ylab("Out-of-register Contacts") +
		ggtitle("In- vs out-of-register contacts for non-floating origins")
ggsave("InVsOutContacts.png")

p <-ggplot(data=odata, aes(x=nrGoodContacts, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Number of good contacts.") +
		ggtitle("Histogram of good contacts (non-floating origins)")
ggsave("goodContacts.png")

p <-ggplot(data=odata, aes(x=nrNumContacts - nrGoodContacts, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("All contacts - good contacts ") +
		ggtitle("Histogram of uncategorised contacts (non-floating origins)")
ggsave("badContacts.png")

p <-ggplot(data=odata, aes(x=nrGoodContacts / (nrNumContacts - nrGoodContacts), fill=factor(success) ) )
p + geom_histogram( position = 'dodge' ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Good contacts / uncateogorised contacts ") +
		ggtitle("Histogram of ratio of good/uncategorised contacts (non-floating origins)")
ggsave("contactsRatio.png")

p <-ggplot(data=odata, aes(x=nrGoodContacts - (nrNumContacts - nrGoodContacts), fill=factor(success) ) )
p + geom_histogram( position = 'dodge' ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Good contacts - uncateogorised contacts ") +
		ggtitle("Histogram of Good minus Uncategorised contacts (non-floating origins)")
ggsave("contactsDifference.png")

p <-ggplot(data=odata, aes(x=nrGoodContacts, y=nrInRegisterContacts, colour=factor(success) ) )
p + geom_point() +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Good Contacts") +
		ylab("In-register Contacts") +
		ggtitle("Good contacts vs in-register contacts for non-floating origins")
ggsave("goodVsInRegister.png")

# plot of nrInRegisterContacts vs nrGoodContacts
p <-ggplot(data=odata, aes(x=shelxeCC, y=nrGoodContacts, colour=factor(success) ) )
p + geom_point() + scale_colour_manual( values=c(fcolour,scolour),
		name="Success/Failure",
		labels=c("Failure", "Success")
) +
xlab("Shelxe CC") +
ylab("Good Contacts") +
ggtitle("Good contacts vs shelxe CC for non-floating origins")
ggsave("CCVsInRegister.png")


p <-ggplot(data=odata, aes(x=numPlacedAtoms-aoNumContacts, y=aoNumContacts, colour=factor(success) ) )
p + geom_point( size=1 ) +
		scale_size() +
		facet_grid( ensembleSideChainTreatment ~ .) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("Num. unmatched atoms") +
		ylab("All contacts (< 0.5A)") +
		ggtitle("Num coincident vs unmatched atoms non-floating origins")
ggsave("CoindicentVsUnmatched.png")

p <-ggplot(data=odata[ odata$aoOrigin == odata$roOrigin , ], aes(x=roRioGood, y=aoNumContacts, colour=factor(success) ) )
p + geom_point( size=1 ) +
		scale_size() +
		facet_grid( ensembleSideChainTreatment ~ .) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("RIO score (in- + out-of-register)") +
		ylab("All contacts (< 0.5A)") +
		ggtitle("Num coincident vs RIO C-alpha non-floating & matching origins")
ggsave("CoindicentVsRIO.png")


# plot of allContacts against numAtoms
#	scale_y_continuous( limits=c(0, max( odata$numPlacedAtoms - odata$aoNumContacts , na.rm=TRUE) ) ) +
#	stat_sum( aes(size = ..n..) ) +
p <-ggplot(data=odata, aes(x=numPlacedAtoms-aoNumContacts, y=aoNumContacts, colour=factor(success) ) )
p + geom_point() +
	scale_size() +
	facet_grid( ensembleSideChainTreatment ~ .) +
	scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
	xlab("Num. unmatched atoms") +
	ylab("All contacts (< 0.5A)") +
	ggtitle("Num coincident vs unmatched atoms non-floating origins")
ggsave("CoindicentVsUnmatched.png")

#odata[ odata$numPlacedAtoms - odata$aoNumContacts > 3000, ]
#		scale_y_continuous( limits=c(0, max( odata$numPlacedAtoms, na.rm=TRUE) )  ) +
#	scale_x_continuous( limits=c(-10,100) ) +
#	scale_y_continuous( limits=c(-10,100) ) +

p <-ggplot(data=odata[ odata$success == 0,], aes(x=numAllAtoms-numAllContacts, y=numAllContacts ) )
p + geom_point() +
		scale_y_continuous( limits=c(0, 1000)  ) +
		xlab("Num. unmatched atoms") +
		ylab("All contacts (< 0.5A)") +
		ggtitle("Num contacts vs num unmatched atoms non-floating origins")


# TFZ/LLG
p <-ggplot(data=data, aes(x=phaserTFZ, y=phaserLLG, colour=factor(success) ) )
p + geom_point() + scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Phaser TFZ") +
		ylab("Phaser LLG") +
		ggtitle("Phaser LLG vs Phaser TFZ")
ggsave("LLGvsTFZ.png")

p <-ggplot(data=data, aes(x=phaserTFZ, y=phaserLLG, colour=factor(success) ) )
p + geom_point() +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		scale_y_continuous( limits=c(-3000, max(data$phaserLLG, na.rm=TRUE) )  )+
		xlab("Phaser TFZ") +
		ylab("Phaser LLG") +
		ggtitle("Phaser LLG vs Phaser TFZ")
ggsave("LLGvsTFZtrunc.png")


# CC distribution
p <-ggplot(data=data, aes(x=shelxeCC, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("CC score") +
		ggtitle("Histogram of CC scores by success")
ggsave("CCscores.png")


# Resolution vs length
p <-ggplot(data=data, aes(x=resolution, y=fastaLength, colour=factor(success) ) )
p + geom_point() +
		stat_sum() +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Resolution (A)") +
		ylab("Protein length (residues)") +
		ggtitle("Resolution vs protein length")
ggsave("resolutionVsLength.png")



q()
#write.csv(odata, "odata.csv", row.names=FALSE)


#p + layer(geom="point") + facet_wrap(~success)

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
x <- x[ order( x$pdbCode, x$shelxeCC, decreasing=TRUE ), ]
# Select top by selecting not duplicates on pdbCode
x <- x[ !duplicated(x$pdbCode), ]
# NB LOOK AT REORDER

# Fix missing values
#completeEnsemble <- complete.cases(data$ensembleNumResidues)
#completeHelix <- complete.cases(data$lenHelix)

pdbAll <- unique( data$pdbCode )
sprintf("Len all PDBs %s", length(pdbAll) )
pdbSuccess <- unique( pdata$pdbCode )
sprintf("Len success PDBs %s", length(pdbSuccess) )

q()
