library(ggplot2)
scolour="#3333FF"
fcolour="#FF0000"
#data <- read.table(file="results_all_arpwarp.csv",sep=',', header=T)
#data <- read.table(file="results_all_arpwarp.csv",sep=',', header=T)
data <- read.table(file="results_all_arpwarpH.csv",sep=',', header=T)

# Categorise successes
data$successShelxe <- as.numeric( data$shelxeCC >= 25 & data$shelxeAvgChainLength >= 10 )
data$successShelxe <- replace( data$successShelxe, is.na(data$successShelxe), 0 )

# Definition of success where at least one rebuild worked 
data$successRebuild <- as.numeric( data$arpWarpFinalRfree <= 0.45 | data$buccFinalRfree <= 0.45 )
#data$successRebuild <- as.numeric(  data$buccFinalRfree <= 0.45 )
data$successRebuild <- replace( data$successRebuild, is.na(data$successRebuild), 0 )

# Data where MR seemingly worked according to refmac rfree
data$successRefmac <- as.numeric( data$rfree <= 0.45 )
data$successRefmac <- replace( data$successRefmac, is.na(data$successRefmac), 0 )

# Data where MR seemingly worked according to phaser
data$successPhaser <- as.numeric( data$phaserTFZ > 8.0 | data$phaserLLG > 120 )
data$successPhaser <- replace( data$successPhaser, is.na(data$successPhaser), 0 )

# Gold standard
data$success <- as.numeric( data$successShelxe == 1 & data$successRebuild == 1 )
data$success <- replace( data$success, is.na(data$success), 0 )

# single-model success
data$successSingleModel <- as.numeric( data$single_model_shelxeCC >= 25 & 
				data$single_model_shelxeAvgChainLength >= 10 &
				( data$single_model_buccFinalRfree <= 0.45 | data$single_model_arpWarpFinalRfree  <= 0.45 ) )
data$successSingleModel <- replace( data$successSingleModel, is.na(data$successSingleModel), 0 )

# helix success
#ok <- (data$helix_All_atom_shelxeCC >= 25 & data$helix_All_atom_shelxeAvgChainLength >= 10 & (data$helix_All_atom_buccFinalRfree <= 0.45 | data$helix_All_atom_arpWarpFinalRfree <= 0.45))| 
#(data$helix_SCWRL_reliable_sidechains_shelxeCC >= 25 & data$helix_SCWRL_reliable_sidechains_shelxeAvgChainLength >= 10 & (data$helix_SCWRL_reliable_sidechains_buccFinalRfree <= 0.45 | data$helix_SCWRL_reliable_sidechains_arpWarpFinalRfree <= 0.45))| 
#(data$helix_poly_ala_shelxeCC >= 25 & data$helix_poly_ala_shelxeAvgChainLength >= 10 & (data$helix_poly_ala_buccFinalRfree <= 0.45 | data$helix_poly_ala_arpWarpFinalRfree <= 0.45) )
ok <-(data$helix_poly_ala_shelxeCC >= 25 & data$helix_poly_ala_shelxeAvgChainLength >= 10 & (data$helix_poly_ala_buccFinalRfree <= 0.45 | data$helix_poly_ala_arpWarpFinalRfree <= 0.45) )
data$successHelix <- as.numeric( ok )
data$successHelix <- replace( data$successHelix, is.na(data$successHelix), 0 )


# Calculate number of copies of decoy that were placed
data$numPlacedChains <- data$numPlacedAtoms / data$ensembleNumAtoms

#write(data,"foo.csv",row.NAMES=FALSE)

# Data for which there are placed models and we have an origin
odata = data[ ! is.na(data$ccmtzOrigin) & ! is.na( data$phaserLLG), ]

# Categories for resolution
data$resCat <- 0
data$resCat[ data$resolution < 2.0 ] <- 1
data$resCat[  data$resolution >= 2.0 & data$resolution < 2.2 ] <- 2
data$resCat[  data$resolution >= 2.2 & data$resolution < 2.4 ] <- 3
data$resCat[  data$resolution >= 2.4 & data$resolution < 2.6 ] <- 4
data$resCat[data$resolution >= 2.6 ] <- 5

# Horrible hack - use label as factor - needed as no labeller attribute to facet_wrap
data$resCatw <- 0
data$resCatw[ data$resolution < 2.0 ] <- "1: res. < 2.0"
data$resCatw[  data$resolution >= 2.0 & data$resolution < 2.2 ] <- "2: 2.0 < res. < 2.2"
data$resCatw[  data$resolution >= 2.2 & data$resolution < 2.4 ] <- "3: 2.2 < res. < 2.4"
data$resCatw[  data$resolution >= 2.4 & data$resolution < 2.6 ] <- "4: 2.4 < res. < 2.6"
data$resCatw[data$resolution >= 2.6 ] <- "5: res. > 2.6"

l = function( variable, value ) {
	if ( variable == "success" ) {
		return( c("Failure","Success"))
	} else if ( variable == "ensembleSideChainTreatment" ) {
		return( c("All atom","Poly-alanine","SCWRL reliable") )
	} else if ( variable == "resCat" ) {
		return( c("res. < 2.0","2.0 < res. < 2.2","2.2 < res. < 2.4","2.4 < res. < 2.6","res. > 2.6") )
	} else {
		return ("foo")
	}
}

#http://stackoverflow.com/questions/19357668/r-ggplot2-facetting-keep-ratio-but-override-define-output-plot-size

# Comparison of Ensemble with Single model
smdata <- data[ data$single_model_ran == "True", ]
x <- smdata[ order( smdata$pdbCode, smdata$success, smdata$shelxeCC, decreasing=TRUE ), ]

# Select top by selecting not duplicates on pdbCode
x <- x[ !duplicated(x$pdbCode), c("pdbCode", "fastaLength","resolution","numChains","numPlacedChains", "shelxeCC", "shelxeAvgChainLength")  ]

# Now put in alphabetical order
x <- x[ order( x$pdbCode ), ]

# Need to get numbers of success
# This gets the stats  - we can join because the by function is the pdbCode which is in similar alphabetic order
x["numModels"] <- aggregate( smdata$pdbCode, by=list(smdata$pdbCode), FUN=length )[2]
x["worked"] <- aggregate( smdata$success, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x$worked <- replace( x$worked, x$worked > 0, 1 )

# Different measures of success
x["success"] <- aggregate( smdata$success, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successShelxe"] <- aggregate( smdata$successShelxe, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successRefmac"] <- aggregate( smdata$successRefmac, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successPhaser"] <- aggregate( smdata$successPhaser, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successRebuild"] <- aggregate( smdata$successRebuild, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]

x["successSingleModel"] <- aggregate( smdata$successSingleModel, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successHelix"] <- aggregate( smdata$successHelix, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["helixWorked"] <- x$successHelix
x$helixWorked <- replace( x$helixWorked, x$helixWorked > 0, 1 )


# NB LOOK AT REORDER
# Now put in order by success, resolution
summaryData <- x[ order( -x$worked, x$resolution ), ]
#x <- x[ order( x$success, decreasing=TRUE ), ]
write.csv(summaryData, "summary1.csv", row.names=FALSE)

#
# Summary information
#

pdbAll <- unique( data$pdbCode )
sprintf("Len all PDBs %s", length(pdbAll) )
pdbSuccess <- unique( data[ data$success==1, ]$pdbCode )
sprintf("Len success PDBs %s", length(pdbSuccess) )

smSuccess <- unique( data[ data$successSingleModel==1, ]$pdbCode )
sprintf("Len single-model success PDBs %s", length(smSuccess) )

hSuccess <- unique( data[ data$helix_success==1, ]$pdbCode )
sprintf("Len helix success PDBs %s", length(hSuccess) )

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
x <- x[ !duplicated(x$pdbCode), c("pdbCode", "fastaLength","resolution","numChains","numPlacedChains", "shelxeCC", "shelxeAvgChainLength")  ]

# Now put in alphabetical order
x <- x[ order( x$pdbCode ), ]

# Need to get numbers of success
# This gets the stats  - we can join because the by function is the pdbCode which is in similar alphabetic order
x["numModels"] <- aggregate( data$pdbCode, by=list(data$pdbCode), FUN=length )[2]
x["worked"] <- aggregate( data$success, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x$worked <- replace( x$worked, x$worked > 0, 1 )


# Different measures of success
x["success"] <- aggregate( data$success, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successShelxe"] <- aggregate( data$successShelxe, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successRefmac"] <- aggregate( data$successRefmac, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successPhaser"] <- aggregate( data$successPhaser, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successRebuild"] <- aggregate( data$successRebuild, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]

x["successSingleModel"] <- aggregate( data$successSingleModel, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
x["successHelix"] <- aggregate( data$successHelix, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]

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


# ADD SINGLE MODEL AND HELIX
#x["successfulSingleModel"] <- aggregate( data$successSingleModel, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
#x$successfulSingleModel[ x$successfulSingleModel >= 1 ] <- 1
#x["successfulHelix"] <- aggregate( data$helix_success, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
#x$successfulHelix[ x$successfulHelix >= 1 ] <- 1

# NB LOOK AT REORDER
# Now put in order by success, resolution
summaryData <- x[ order( -x$worked, x$resolution ), ]
#x <- x[ order( x$success, decreasing=TRUE ), ]
write.csv(summaryData, "summary.csv", row.names=FALSE)

sdata = data[ data$success == 1, ]
# Most accurate models that solved 
# Highest number RIO
sdata[ sdata$ccmtzRioGood == max(sdata$ccmtzRioGood,na.rm=TRUE), c("pdbCode","ensembleName") ]
# 2IC9 poly_ala_trunc_6.098002_rad_1

# Best Reforigin RMSD reforiginRMSD
sdata[ sdata$reforiginRMSD == min(sdata$reforiginRMSD,na.rm=TRUE), c("pdbCode","ensembleName") ]
# 3LJM All_atom_trunc_0.000847_rad_1

# Longest that solved with zero RIO score
x <- sdata[ sdata$ccmtzRioGood == 0, ]
x[ x$fastaLength == max(x$fastaLength), c("pdbCode","ensembleName") ]
# 2QIH poly_ala_trunc_24.292064_rad_2

# Longest that solved with zero AIO score
x <- sdata[ sdata$ccmtzAaNumContacts== 0, ]
x[ x$fastaLength == max(x$fastaLength), c("pdbCode","ensembleName") ]
# 3M91 SCWRL_reliable_sidechains_trunc_0.066638_rad_1

# highest percentage of model that solved - selecting the largest
# Could use which.max but want to select largest
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


#
# Plots
#

# Ensemble vs single model
p <-ggplot(data=smdata,
		aes(shelxeCC,single_model_shelxeCC,colour=factor(success),shape=factor(successSingleModel)))
p + geom_point() +
scale_colour_manual( values=c(fcolour,scolour),
		name="Ensemble Success",
		labels=c("Failure", "Success")
) +
scale_shape_manual( values=c(0,2),
		name="Single Model Success",
		labels=c("Failure", "Success")
) +
xlab("Ensemble Shelxe CC score") +
ylab("Single Model Shelxe CC score") +
ggtitle("Comparison of Ensemble and Single Model CC scores")
ggsave("SMEShelxe.png")

# Helix vs Ensemble
hdata = data[ data$helix_ran == "True", ]
p <-ggplot(data=hdata,
		aes(shelxeCC,helix_poly_ala_shelxeCC,colour=factor(success),shape=factor(successHelix)))
p + geom_point() +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Ensemble Success",
				labels=c("Failure", "Success")
		) +
		scale_shape_manual( values=c(0,2),
				name="Idealised helix Success",
				labels=c("Failure", "Success")
		) +
		scale_x_continuous(limits=c(0,80)) +
		scale_y_continuous(limits=c(0,80)) +
		xlab("Ensemble Shelxe CC score") +
		ylab("Idealised Helix Shelxe CC score") +
		ggtitle("Comparison of Idealised Helices and Single Model CC scores")
ggsave("SMHShelxe.png")


# Ensemble vs helix
smdata$helixCC <- max(smdata$helix_All_atom_shelxeCC,smdata$helix_SCWRL_reliable_sidechains_shelxeCC,smdata$helix_poly_ala_shelxeCC)
p <-ggplot(data=smdata[ smdata$helix_ran == "True",],
		aes(shelxeCC,helixCC))
p + geom_point() +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Ensemble Success",
				labels=c("Failure", "Success")
		) +
		scale_shape_manual( values=c(0,2),
				name="Single Model Success",
				labels=c("Failure", "Success")
		) +
		xlab("Ensemble Shelxe CC score") +
		ylab("Single Model Shelxe CC score") +
		ggtitle("Comparison of Ensemble and Single Model CC scores")
ggsave("SMEShelxe.png")


# How the different treatments fared - facet by resolution
# Number of success

#p <-ggplot(data=summaryData)
#p + geom_point(mapping=aes(x=pdbCode,y=success), colour="red") +
#geom_point(mapping=aes(x=pdbCode,y=successSingleModel), colour="blue") +
#geom_point(mapping=aes(x=pdbCode,y=successHelix), colour="green")
#
#
#p <-ggplot(data=summaryData)
#p + geom_point(mapping=aes(x=resolution,y=success), colour="red") +
#geom_point(mapping=aes(x=resolution,y=successSingleModel ), colour="blue") +
#geom_point(mapping=aes(x=resolution,y=successHelix), colour="green")



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
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
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

p <-ggplot(data=data[ data$success==1, ], aes(x=buccFinalRfree, y=arpWarpFinalRfree) )
p + geom_point() +
		#facet_grid( .~resCat, space="free", scales="free", labeller=l) +
		facet_wrap( ~resCatw) +
		#stat_sum( aes(size=..n..) ) +
		xlab("buccFinalRfree") +
		ylab("arpWarpFinalRfree") +
		ggtitle("Bucc vs Arpwarp")
ggsave("buccVsArp.png")

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
p <-ggplot( data=data, aes( fastaLength, fill=factor(success) ) )
p + geom_histogram( binwidth=5, position = 'dodge' ) +
	scale_x_continuous( breaks=seq(0,260,20) ) +
	scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")) +
	xlab("Length in residues") +
	ylab("Number of targets") +
	ggtitle("Successes/failures by target protein length")
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

#
# Histograms
#
p <-ggplot(data=odata, aes(x=ccmtzRioGood, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("RIO score") +
		ggtitle("Histogram of RIO scores")
ggsave("rioHistogram.png")

p <-ggplot(data=odata, aes(x=ccmtzRioNoCat, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Number non-RIO C-alpha") +
		ggtitle("non-RIO C-alpha atoms")
ggsave("noRioHistogram.png")


# Proportion RIO of model
p <-ggplot(data=odata, aes(x=ccmtzRioGood/numPlacedCA, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.05 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Proportion RIO") +
		ggtitle("RIO as proportion of model")
ggsave("rioModelHistogramProp.png")

#Calculate proportion of in/out in rioGood
#prop.table(as.matrix(x[-1]),margin=1)
odata$propRioIn <- prop.table(as.matrix( odata[, c("ccmtzRioInregister","ccmtzRioOoRegister") ]),margin=1)[,1]
odata$propRioOut <- prop.table(as.matrix( odata[, c("ccmtzRioInregister","ccmtzRioOoRegister") ]),margin=1)[,2]
#odata$propRioIn <- replace( odata$propRioIn, is.na(odata$propRioIn), -1 )

#p <-ggplot(data=odata[ odata$ccmtzRioGood > 0, ], aes(x=propRioIn, fill=factor(success)) )
#p + geom_histogram( binwidth = 0.05, position='dodge' ) 
p <-ggplot(data=odata[ odata$ccmtzRioGood > 0, ],
		aes(x=propRioIn,
			y=numResidues,
			colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( .~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("In-register proportion of RIO") +
		ylab("Number of residues in ASU") +
		ggtitle("In-register CA as a prop. of RIO vs num. residues for RIO > 0")
ggsave("rioInRegisterPropVsLength.png")

# Proportion RIO of native
p <-ggplot(data=odata, aes(x=ccmtzRioGood/numResidues, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.05 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Proportion RIO") +
		ggtitle("RIO as proportion of native")
ggsave("rioNativeHistogramProp.png")

# Proportion in density of model
p <-ggplot(data=odata, aes(x=ccmtzAaNumContacts/numPlacedAtoms, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.01 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Proportion in density") +
		ggtitle("Atoms in density as a proportion of the model")
ggsave("inDensityModelHistogramProp.png")

# Proportion in density of native
p <-ggplot(data=odata, aes(x=ccmtzAaNumContacts/numAtoms, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.001 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Proportion in density") +
		ggtitle("Atoms in density as proportion of native")
ggsave("inDensityNativeHistogramProp.png")

# Number in density
p <-ggplot(data=odata, aes(x=ccmtzAaNumContacts, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Number in density") +
		ggtitle("Number of atoms in density")
ggsave("inDensityHistogramNum.png")


# Number out of density
p <-ggplot(data=odata, aes(x=numPlacedAtoms-ccmtzAaNumContacts, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Number out of density") +
		ggtitle("Number of atoms out of density")
ggsave("outOfDensityHistogramNum.png")

# Out of density as proportion of model
p <-ggplot(data=odata, aes(x=(numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.01 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Proportion out of density") +
		ggtitle("Atoms out of density as proportion of model")
ggsave("outOfDensityHistogramModel.png")

# Out of density as proportion of native
p <-ggplot(data=odata, aes(x=(numPlacedAtoms-ccmtzAaNumContacts)/numAtoms, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.01 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Proportion out of density") +
		ggtitle("Atoms out of density as proportion of native")
ggsave("outOfDensityHistogramNative.png")

#
# Graphs
#
#Comparison of the RIO figures for successful search models with the number of in-register 
# residues similarly defined.
p <-ggplot(data=odata, aes(x=ccmtzRioOoRegister, y=ccmtzRioInregister, colour=factor(success) ) ) 
p + geom_point() +
		stat_sum( aes(size=..n..) ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		xlab("Out-of-register C-alpha") +
		ylab("In-register C-alpha") +
		ggtitle("In- vs, Out-of-register C-alphas")
ggsave("ooRegisterVsInRegister.png")


# plot of nrInRegisterContacts vs nrGoodContacts
p <-ggplot(data=odata, aes(x=shelxeCC, y=ccmtzRioGood, colour=factor(success) ) )
p + geom_point() +
	scale_colour_manual( values=c(fcolour,scolour),
	name="Success/Failure",
	labels=c("Failure", "Success") ) +
	xlab("Shelxe CC") +
	ylab("RIO score") +
	ggtitle("RIO score vs shelxe CC")
ggsave("CCVsRIO.png")

#
# inDensity vs unmatched
# ronan's idea to divide by number of unmatched
#
# facet_grid( ensembleSideChainTreatment ~ .) +
# label_both


# Proportion RIO vs proportion in density
p <-ggplot(data=odata,
		aes(x=ccmtzAaNumContacts/numAtoms,
				y=ccmtzRioGood/numResidues,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("Proportion  atoms in density (< 0.5A)") +
		ylab("RIO C-alphas as prop. of native (< 1.5A)") +
		ggtitle("RIO vs in-density as proportion of native")
ggsave("inDensityVsRIObyRes.png")


# Measure of how many correct vs number/proportion misplaced
# Misplaced atoms of model measured by: numPlacedAtoms - ccmtzAaNumContacts - proportion is (numPlacedAtoms-ccmtzAaNumContacts)/numAtoms
# Proportion correct measured by: ccmtzAaNumContacts/numAtoms
#                             or: ccmtzRioGood/numResidues
# so...
p <-ggplot(data=odata,
		aes(x=ccmtzAaNumContacts/numAtoms,
			y=(numPlacedAtoms-ccmtzAaNumContacts)/numAtoms,
			colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("Atoms in density as prop. of native (< 0.5A)") +
		ylab("Atoms outside of density as prop. of native (> 0.5A)") +
		ggtitle("In- vs out-of-density as prop. of native")
ggsave("inDensityVsMisplacedPropNative.png")

p <-ggplot(data=odata,
		aes(x=ccmtzAaNumContacts/numPlacedAtoms,
				y=(numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("Atoms in density as prop. of model (< 0.5A)") +
		ylab("Atoms outside of density as prop. of model (> 0.5A)") +
		ggtitle("In- vs out-of-density as prop. of model")
ggsave("inDensityVsMisplacedPropModel.png")

# Number outside of density as proportion of model shows negative signal as proportion of the input signal
# Number in density shows what actually matched
p <-ggplot(data=odata,
		aes(x=ccmtzAaNumContacts/numAtoms,
				y=(numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("Atoms in density as prop. of native (< 0.5A)") +
		ylab("Atoms outside of density as prop. of model (> 0.5A)") +
		ggtitle("In- vs out-of-density as native/model")
ggsave("inDensityVsMisplacedPropBoth.png")

p <-ggplot(data=odata,
		aes(x=ccmtzAaNumContacts,
				y=numPlacedAtoms-ccmtzAaNumContacts,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("Num. atoms in density (< 0.5A)") +
		ylab("Num. out of density (> 0.5A)") +
		ggtitle("Number placed in- vs out-of-density")
ggsave("inDensityVsMisplacedNum.png")


p <-ggplot(data=odata,
		aes(x=ccmtzRioGood/numResidues,
				y=(numPlacedAtoms-ccmtzAaNumContacts)/numAtoms,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("RIO as proportion of native") +
		ylab("Num. outside density as prop. native (> 0.5A)") +
		ggtitle("Out-of-density as prop. of native vs RIO as prop. of native")
ggsave("correctlyVsMisplacedProp.png")

p <-ggplot(data=odata,
		aes(x=ccmtzRioGood/numResidues,
				y=(numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("RIO as proportion of native") +
		ylab("Num. outside density as prop. model (> 0.5A)") +
		ggtitle("Out-of-density as prop. of model vs RIO as prop. of native")
ggsave("correctlyVsMisplacedPropModel.png")

p <-ggplot(data=odata,
		aes(x=ccmtzRioGood/numResidues,
				y=numPlacedAtoms-ccmtzAaNumContacts,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success") ) +
		xlab("RIO as proportion of native") +
		ylab("Numer placed atoms - AIO score") +
		ggtitle("RIO as proportion of native vs num. of misplaced atoms.")
ggsave("correctlyVsMisplacedNum.png")
##
##
##
# Selecting interesting cases
# Successes with RIO == 0
# odata[ odata$ccmtzRioGood == 0 & odata$success == 1, c("pdbCode","ensembleName") ]
# Failures with RIO as prop of native > 0.5
# odata[ odata$ccmtzRioGood/odata$numResidues > 0.5 & odata$success == 0, c("pdbCode","ensembleName") ]

#p <-ggplot(data=odata,
#		aes(x=numPlacedAtoms-ccmtzAaNumContacts, y=ccmtzAaNumContacts,
#				colour=factor(success) ) )
#p + geom_point( size=1 ) +
#		stat_sum( aes(size=..n..)) +
#		#facet_grid( ensembleSideChainTreatment ~ success, labeller=l) +
#		facet_grid( resCat ~ success, labeller=l) +
#		scale_colour_manual( values=c(fcolour, scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success") ) +
#		xlab("Num unmatched atoms") +
#		ylab("Num contacts (< 0.5A)") +
#		ggtitle("Number of coincident vs unmatched atoms")
#ggsave("coincidentVsUnmatchedRes.png")
#
## plot of allContacts against numAtoms
##	scale_y_continuous( limits=c(0, max( odata$numPlacedAtoms - odata$aoNumContacts , na.rm=TRUE) ) ) +
##	stat_sum( aes(size = ..n..) ) +
#p + geom_point() +
#		scale_size() +
#		#facet_grid( ensembleSideChainTreatment ~ .) +
#		facet_grid( ensembleSideChainTreatment ~ success, labeller=l) +
#		scale_colour_manual( values=c(fcolour, scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success") ) +
#		xlab("Num. unmatched atoms") +
#		ylab("All contacts (< 0.5A)") +
#		ggtitle("Num coincident vs unmatched atoms")
#ggsave("coincidentVsUnmatchedResSideChain.png")
#
#p <-ggplot(data=odata, aes(x=ccmtzRioGood/numPlacedCA, y=ccmtzAaNumContacts/numPlacedAtoms, colour=factor(success) ) )
#p + geom_point( size=1 ) +
#		stat_sum( aes(size=..n..) ) +
#		facet_grid( ensembleSideChainTreatment ~ success, labeller=l ) +
#		scale_colour_manual( values=c(fcolour, scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success") ) +
#		xlab("% RIO score (in- + out-of-register)") +
#		ylab("% contacts (< 0.5A)") +
#		ggtitle("% coincident vs RIO C-alpha")
#ggsave("percentCoincidentVsRIO.png")
#
#
#p <-ggplot(data=odata,
#		aes(x=(numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms, y=ccmtzAaNumContacts/numPlacedAtoms,
#				colour=factor(success) ) )
#p + geom_point( size=1 ) +
#		stat_sum( aes(size=..n..)) +
#		facet_grid( ensembleSideChainTreatment ~ success, labeller=l) +
#		scale_colour_manual( values=c(fcolour, scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success") ) +
#		xlab("% unmatched atoms") +
#		ylab("% contacts (< 0.5A)") +
#		ggtitle("Percentage coincident vs unmatched atoms")
#ggsave("percentCoincidentVsUnmatched.png")
#
#p <-ggplot(data=odata, aes(x=ccmtzRioGood/numPlacedCA, y=ccmtzAaNumContacts/numPlacedAtoms, colour=factor(success) ) )
#p + geom_point( size=1 ) +
#		stat_sum( aes(size=..n..) ) +
#		facet_grid( ensembleSideChainTreatment ~ success, labeller=l ) +
#		scale_colour_manual( values=c(fcolour, scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success") ) +
#		xlab("% RIO score (in- + out-of-register)") +
#		ylab("% contacts (< 0.5A)") +
#		ggtitle("% coincident vs RIO C-alpha")
#ggsave("percentCoincidentVsRIO.png")


# Length distribution of helcies
p <-ggplot(data=odata, aes(x=rioLenHelix, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1  ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success")
		) +
		ylab("Number of cases") +
		xlab("Length of RIO helix") +
		ggtitle("Distribution of RIO helices")
ggsave("rioHelices.png")


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
