

library(ggplot2)
scolour="#3333FF"
fcolour="#FF0000"
setwd("/Users/jmht/Documents/AMPLE/data/coiled-coils/ensemble")
#setwd("/home/jmht/Documents/test/CC/contacts")
#data <- read.table(file="results_bucc.csv",sep=',', header=T)
data <- read.table(file="results.csv",sep=',', header=T)

# Categorise successes
data$success <- as.numeric( data$shelxeCC >= 25 & data$shelxeAvgChainLength >= 10 )
data$success <- replace( data$success, is.na(data$success), 0 )

# Need to remove nolog values from buccFinalRfree
# Horrible hack to get Rfree back to a number...
data$buccFinalRfree[ data$buccFinalRfree == "nolog" ] <- NA
data$buccFinalRfree <- as.numeric( as.character(data$buccFinalRfree) )

# Calculate number of copies of decoy that were placed
data$numPlacedChains <- data$numPlacedAtoms / data$ensembleNumAtoms

#
# Summary information
#

pdbAll <- unique( data$pdbCode )
sprintf("Len all PDBs %s", length(pdbAll) )
pdbSuccess <- unique( data[ data$success==1, ]$pdbCode )
sprintf("Len success PDBs %s", length(pdbSuccess) )

# Failures
#setdiff( pdbAll, pdbSuccess )


# * mean resolution of failing/succesful cases
x <- mean(  data[ data$success == 1, ]$resolution ) # 1.969702
cat( "Mean resolution for successes is ",x,"\n")
x <- mean(  data[ data$success == 0, ]$resolution ) # 2.007124
cat( "Mean resolution for failures is ",x,"\n")

# * mean of solvent content of failing/succesful cases
#unique( data[ is.na( data$solventContent ),]$pdbCode ) # 1G1J 1KYC 1P9I 3CVF
x <-mean( data[ data$success == 1, ]$solventContent, na.rm=TRUE ) # 47.20036
cat( "Mean solvent content for successes is ",x,"\n")
x <-mean( data[ data$success == 0, ]$solventContent, na.rm=TRUE ) # 50.47158
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



#
# Summary table
#
#http://stackoverflow.com/questions/6289538/aggregate-a-dataframe-on-a-given-column-and-display-another-column
#http://stackoverflow.com/questions/2822156/select-rows-with-largest-value-of-variable-within-a-group-in-r

# For each case need to get the best success - i.e. success with max CC
# Order by pdbCode and CC
x <- data[ order( data$pdbCode, data$success, data$shelxeCC, decreasing=TRUE ), ]
# Select top by selecting not duplicates on pdbCode
x <- x[ !duplicated(x$pdbCode), c("pdbCode", "success","fastaLength","resolution","numChains","numPlacedChains", "shelxeCC", "shelxeAvgChainLength")  ]
# Now put in alphabetical order
x <- x[ order( x$pdbCode ), ]

# Need to get numbers of success
# This gets the stats  - we can join because the by function is the pdbCode which is in similar alphabetic order
x["successfulModels"] <- aggregate( data$success, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["failedModels"] <- aggregate( data$success, by=list(data$pdbCode), FUN=function(x){ sum( x == 0 ) } )[2]

x["successPolyAla"]  <- aggregate( 
		data[ data$ensembleSideChain == "poly_ala", ]$success,
		by=list( data[ data$ensembleSideChain == "poly_ala", ]$pdbCode ),
		FUN=function(x){ sum( x == 1 ) } )[2]
x["successScwrl"]  <- aggregate(
		data[ data$ensembleSideChain == "SCWRL_reliable_sidechains", ]$success,
		by=list( data[ data$ensembleSideChain == "SCWRL_reliable_sidechains", ]$pdbCode ),
		FUN=function(x){ sum( x == 1 ) } )[2]
x["successAllAtom"]  <- aggregate(
		data[ data$ensembleSideChain == "All_atom", ]$success,
		by=list( data[ data$ensembleSideChain == "All_atom", ]$pdbCode ),
		FUN=function(x){ sum( x == 1 ) } )[2]
# NB LOOK AT REORDER
summary <- x[ order( x$success, decreasing=TRUE ), ]
write.csv(summary, "summary.csv", row.names=FALSE)


# highest percentage of model that solved - selecting the largest
# Could use which.max but want to select largest
#sdata = data[ data$success == 1, ]
#x <- sdata[ sdata$ensemblePercentModel == max( sdata$ensemblePercentModel ), ]
#x[ order( x$fastaLength, decreasing=TRUE ), ][1,]
# 3K29 SCWRL_reliable_sidechains_trunc_168.103952_rad_3 100% 169 residues

# Min
#x <- sdata[ sdata$ensemblePercentModel == min( sdata$ensemblePercentModel ), ]
#x[ order( x$fastaLength, decreasing=TRUE ), ][1,]
# 3OKQ poly_ala_trunc_0.202334_rad_1 4% 7 residues

# Longest model that solved
#x <- sdata[ sdata$ensembleNumResidues == max( sdata$ensembleNumResidues ), ]
#x[ order( x$fastaLength, decreasing=TRUE ), ][1,]
# 3K29 SCWRL_reliable_sidechains_trunc_168.103952_rad_3 169 residues

# Shortest model that solved
#x <- sdata[ sdata$ensembleNumResidues == min( sdata$ensembleNumResidues ), ]
#x[ order( x$fastaLength, decreasing=TRUE ), ][1,]
#  2OVC SCWRL_reliable_sidechains_trunc_0.012153_rad_1  5 residues fasta: 33



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



###############################################################################################################################
#
# Overall/Ensemble info
#
###############################################################################################################################

# length distribution of successful search models both as number of residues 
# and as fraction of target chain length remaining in the search model.

# Distribution of all ensembles
p <-ggplot(data=data, aes(x=ensembleNumResidues, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(fcolour,scolour),
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
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Percent of residues in ensemble") +
		ylab("Number of ensembles") +
		ggtitle("Percentage of residues in ensemble for successful and failing cases")
ggsave("PercentVsEnsemble.png")

# Number that solved under each side-chain treatment
p <-ggplot( data=data, aes( factor(ensembleSideChainTreatment), fill=factor(success) ) ) 
p + geom_histogram( position = 'dodge' ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Side-chain Treatment") +
		ggtitle("Histogram of side-chain treatment for successful and failing cases")
ggsave("sideChain.png")


# As a table
#x <- aggregate( data$success==0, by=list( SideChainTreatment = data$ensembleSideChainTreatment ), FUN=sum)
#names(x) <- sub("^x$", "failures", names(x))
#x$success <- aggregate( data$success==1, by=list( SideChainTreatment = data$ensembleSideChainTreatment ), FUN=sum)$x

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

# CC distribution
p <-ggplot(data=data, aes(x=shelxeCC, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("CC score") +
		ggtitle("Histogram of CC scores by success")
ggsave("CCscores.png")


# Resolution vs length
#		stat_sum( group=c(1,2) ) +
#		stat_sum( aes(size=..n..) ) +
p <-ggplot(data=data, aes(x=resolution, y=fastaLength, colour=factor(success) ) )
p + geom_point() +
		stat_sum( aes(size=..n..) ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Resolution (A)") +
		ylab("Protein length (residues)") +
		ggtitle("Resolution vs protein length")
ggsave("resolutionVsLength.png")

# Binned (20) bar graph with length along the bottom and the bars dividied/coloured by success
# NB: Uuses data from summary
p <-ggplot( data=summary, aes( fastaLength, fill=factor(success) ) )
p + geom_histogram( binwidth=20 ) +
	scale_x_continuous( breaks=seq(0,260,20) ) +
	scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
	) +
	xlab("Length in residues") +
	ylab("Number of targets") +
	ggtitle("Number of successes/failures by protein length")
ggsave("targetByLength.png")

#Distribution of chain lengths of all ensembles

# percentage traced vs shelxeCC
# FIX _ ASSUMES ALL CHAINS ARE THE SAME!
#p <-ggplot(data=data, aes(x= ( ( shelxeAvgChainLength * shelxeNumChains ) / ( fastaLength * numChains ) ) * 100,
#				y=buccFinalRfree, colour=factor(success) ) )
#p + geom_point( size= 1) +
#		stat_sum( aes(size=..n..) ) +
#		scale_colour_manual( values=c(fcolour,scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success")
#		) +
#		xlab("% Shelxe Trace") +
#		ylab("Shelxe CC") +
#		ggtitle("Percentage traced by Shelxe vs Shelxe CC score")
#ggsave("shelxePropVsCC.png")


###############################################################################################################################
#
# Comparison of model quality
#
###############################################################################################################################
p <-ggplot(data=data, aes(x=ensembleNativeTM, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 0.01 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("TM score of model to native") +
		ggtitle("Histogram of TM score for models for successful and failing cases")
ggsave("modelTM.png")

p <-ggplot(data=data, aes(x=ensembleNativeRMSD, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("RMSD score of model to native") +
		ggtitle("Histogram of RMSD score for models for successful and failing cases")
ggsave("modelRMSD.png")

# Plot of Reforigin RMSD for successes and failures
p <-ggplot(data=data, aes(x=reforiginRMSD, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("reforigin RMSD of placed model to native") +
		ggtitle("Histogram of reforigin RMSD scores by success")
ggsave("reforiginRMSD.png")

###############################################################################################################################
#
# RIO data
#
###############################################################################################################################

#
# Can only compare where a phaser model was produced and where we could find an origin
#
odata = data[ ! is.na(data$ccmtzOrigin) & ! is.na( data$phaserLLG), ]
## HACK - set NA aoNumContacts to zero - we need to do this as we set NA when they were 0 when looking for origins
#missing <- is.na( odata$aoNumContacts )
#odata$aoNumContacts[ missing ] <- 0
#odata$aoNumRio[ missing ] <- 0
#odata$aoRioInregister[ missing ] <- 0
#odata$aoRioOoRegister[ missing ] <- 0
#odata$aoRioBackwards[ missing ] <- 0
#odata$aoRioGood[ missing ] <- 0
#odata$aoRioNocat[ missing ] <- 0

## Find where origins are the same
#x <- odata[ ! ( is.na( odata$aoOrigin ) & is.na( odata$roOrigin ) ) & 
#				! ( odata$aoOrigin == "" & odata$roOrigin == "" ) & 
#				odata$aoOrigin == odata$roOrigin, ]
# Find where they have been found and data is different - only 90
#x <- odata[ ! is.na( odata$aoOrigin ) & ! is.na( odata$roOrigin )  & 
#				! odata$aoOrigin == "" & ! odata$roOrigin == ""  & 
#				odata$aoOrigin != odata$roOrigin, ]
#x <- odata[ ! is.na( odata$aoOrigin ) & ! is.na( odata$roOrigin )  & 
#				! odata$aoOrigin == "" & ! odata$roOrigin == ""  & 
#				odata$aoOrigin != odata$roOrigin, ]

#Comparison of the RIO figures for successful search models with the number of in-register 
# residues similarly defined.
p <-ggplot(data=odata, aes(x=ccmtzRioOoRegister, y=ccmtzRioInregister, colour=factor(success) ) ) 
p + geom_point() +
		stat_sum( aes(size=..n..) ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Out-of-register Contacts") +
		ylab("In-register Contacts") +
		ggtitle("In- vs, Out-of-register")
ggsave("ooRegisterVsInRegister.png")

p <-ggplot(data=odata, aes(x=ccmtzRioGood, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Number of good contacts.") +
		ggtitle("Histogram of in- plus out-of-register contacts.")
ggsave("goodContacts.png")

p <-ggplot(data=odata, aes(x=ccmtzRioNoCat, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("All contacts - good contacts ") +
		ggtitle("Histogram of uncategorised contacts")
ggsave("badContacts.png")

# plot of nrInRegisterContacts vs nrGoodContacts
p <-ggplot(data=odata, aes(x=shelxeCC, y=ccmtzRioGood, colour=factor(success) ) )
p + geom_point() +
	scale_colour_manual( values=c(fcolour,scolour),
	name="Success/Failure",
	labels=c("Failure", "Success") ) +
	xlab("Shelxe CC") +
	ylab("Good Contacts") +
	ggtitle("Good contacts vs shelxe CC for non-floating origins")
ggsave("CCVsInRegister.png")

#
# coincident vs unmatched
# ronan's idea to divide by number of unmatched
#
##p <-ggplot(data=odata, aes(x=numPlacedAtoms-aoNumContacts, y=aoNumContacts, colour=factor(success) ) )
#p <-ggplot(data=odata, aes(x=(numPlacedAtoms-aoNumContacts)/numPlacedAtoms, y=aoNumContacts/numPlacedAtoms, colour=factor(success) ) )
#		#scale_size() +
#p + geom_point( size=2 ) +
#		stat_sum( aes(size=..n..) ) +
#		facet_grid( ensembleSideChainTreatment ~ .) +
#		scale_colour_manual( values=c(fcolour, scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success") ) +
#		xlab("% unmatched atoms") +
#		ylab("% contacts (< 0.5A)") +
#		ggtitle("Num coincident vs unmatched atoms non-floating origins")
#ggsave("CoincidentVsUnmatched.png")

# facet_grid( ensembleSideChainTreatment ~ .) +
# label_both
l = function( variable, value ) {
	if ( variable == "success" ) {
		return( c("Failure","Success"))
	} else if ( variable == "ensembleSideChainTreatment" ) {
		return( c("All atom","Poly-alanine","SCWRL reliable") )
	} else {
		return ("foo")
	}
}

p <-ggplot(data=odata,
		aes(x=ccmtzRioGood, y=numPlacedAtoms-ccmtzAaNumContacts,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..)) +
		facet_grid( ensembleSideChainTreatment ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("RIO score (in- + out-of-register)") +
		ylab("Num non-coincident atoms (< 0.5A)") +
		ggtitle("Rio score vs unmatched atoms")
ggsave("RIOVsUnmatched.png")


p <-ggplot(data=odata,
		aes(x=numPlacedAtoms-ccmtzAaNumContacts, y=ccmtzAaNumContacts,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..)) +
		facet_grid( ensembleSideChainTreatment ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("Num unmatched atoms") +
		ylab("Num contacts (< 0.5A)") +
		ggtitle("Number of coincident vs unmatched atoms")
ggsave("coincidentVsUnmatched.png")

p <-ggplot(data=odata, aes(x=ccmtzRioGood/numPlacedCA, y=ccmtzAaNumContacts/numPlacedAtoms, colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( ensembleSideChainTreatment ~ success, labeller=l ) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("% RIO score (in- + out-of-register)") +
		ylab("% contacts (< 0.5A)") +
		ggtitle("% coincident vs RIO C-alpha")
ggsave("percentCoincidentVsRIO.png")


p <-ggplot(data=odata,
		aes(x=(numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms, y=ccmtzAaNumContacts/numPlacedAtoms,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..)) +
		facet_grid( ensembleSideChainTreatment ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("% unmatched atoms") +
		ylab("% contacts (< 0.5A)") +
		ggtitle("Percentage coincident vs unmatched atoms")
ggsave("percentCoincidentVsUnmatched.png")

p <-ggplot(data=odata, aes(x=ccmtzRioGood/numPlacedCA, y=ccmtzAaNumContacts/numPlacedAtoms, colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( ensembleSideChainTreatment ~ success, labeller=l ) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("% RIO score (in- + out-of-register)") +
		ylab("% contacts (< 0.5A)") +
		ggtitle("% coincident vs RIO C-alpha")
ggsave("percentCoincidentVsRIO.png")

# plot of allContacts against numAtoms
#	scale_y_continuous( limits=c(0, max( odata$numPlacedAtoms - odata$aoNumContacts , na.rm=TRUE) ) ) +
#	stat_sum( aes(size = ..n..) ) +
p <-ggplot(data=odata, aes(x=numPlacedAtoms-ccmtzAaNumContacts, y=ccmtzAaNumContacts, colour=factor(success) ) )
p + geom_point() +
	scale_size() +
	facet_grid( ensembleSideChainTreatment ~ .) +
	scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
	xlab("Num. unmatched atoms") +
	ylab("All contacts (< 0.5A)") +
	ggtitle("Num coincident vs unmatched atoms")
ggsave("coindicentVsUnmatched.png")

###############################################################################################################################
#
# TFZ/LLG
#
###############################################################################################################################
#
p <-ggplot(data=data[ ! is.na( data$phaserLLG) & ! is.na( data$phaserTFZ),], 
		aes(x=phaserTFZ, y=phaserLLG, colour=factor(success) ) )
p + geom_point() + 
		scale_colour_manual( values=c(fcolour,scolour),
							name="Success/Failure",
							labels=c("Failure", "Success") ) +
		xlab("Phaser TFZ") +
		ylab("Phaser LLG") +
		ggtitle("Phaser LLG vs Phaser TFZ")
ggsave("LLGvsTFZ.png")
# Below for truncating the y-axis
# scale_y_continuous( limits=c(-3000, max(data$phaserLLG, na.rm=TRUE) )  )+
