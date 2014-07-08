library(ggplot2)
scolour="#3333FF"
fcolour="#FF0000"
comparison=FALSE
#data <- read.table(file="results_all_arpwarpH.csv",sep=',', header=T)
#data <- read.table(file="comparison_results.csv",sep=',', header=T)
data <- read.table(file="final_results.csv",sep=',', header=T)

# Categorise successes
data$successShelxe <- as.numeric( data$shelxeCC >= 25 & data$shelxeAvgChainLength >= 10 )
data$successShelxe <- replace( data$successShelxe, is.na(data$successShelxe), 0 )

# Definition of success where at least one rebuild worked 
data$minRfree <- pmin( data$arpWarpFinalRfree, data$buccFinalRfree, na.rm=TRUE )
data$successRebuild <- as.numeric( data$minRfree <= 0.45 )
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

if (comparison){
	# single-structure success
	data$successSingleStructure <- as.numeric( data$single_model_shelxeCC >= 25 & 
					data$single_model_shelxeAvgChainLength >= 10 &
					( data$single_model_buccFinalRfree <= 0.45 | data$single_model_arpWarpFinalRfree  <= 0.45 ) )
	data$successSingleStructure <- replace( data$successSingleStructure, is.na(data$successSingleStructure), 0 )
	
	# helix success
	#ok <- (data$helix_All_atom_shelxeCC >= 25 & data$helix_All_atom_shelxeAvgChainLength >= 10 & (data$helix_All_atom_buccFinalRfree <= 0.45 | data$helix_All_atom_arpWarpFinalRfree <= 0.45))| 
	#(data$helix_SCWRL_reliable_sidechains_shelxeCC >= 25 & data$helix_SCWRL_reliable_sidechains_shelxeAvgChainLength >= 10 & (data$helix_SCWRL_reliable_sidechains_buccFinalRfree <= 0.45 | data$helix_SCWRL_reliable_sidechains_arpWarpFinalRfree <= 0.45))| 
	#(data$helix_poly_ala_shelxeCC >= 25 & data$helix_poly_ala_shelxeAvgChainLength >= 10 & (data$helix_poly_ala_buccFinalRfree <= 0.45 | data$helix_poly_ala_arpWarpFinalRfree <= 0.45) )
	#ok <-(data$helix_poly_ala_shelxeCC >= 25 & data$helix_poly_ala_shelxeAvgChainLength >= 10 & (data$helix_poly_ala_buccFinalRfree <= 0.45 | data$helix_poly_ala_arpWarpFinalRfree <= 0.45) )
	#data$successHelix <- as.numeric( ok )
	#data$successHelix <- replace( data$successHelix, is.na(data$successHelix), 0 )
	ok <-(data$helix_poly_ala_shelxeCC >= 25 & data$helix_poly_ala_shelxeAvgChainLength >= 10 & (data$helix_poly_ala_buccFinalRfree <= 0.45 | data$helix_poly_ala_arpWarpFinalRfree <= 0.45) )
	ok <- as.numeric( ok )
	data$successHelix1 <- replace( ok, is.na(ok), 0 )
	ok <-(data$helix_SCWRL_reliable_sidechains_shelxeCC >= 25 & data$helix_SCWRL_reliable_sidechains_shelxeAvgChainLength >= 10 & (data$helix_SCWRL_reliable_sidechains_buccFinalRfree <= 0.45 | data$helix_SCWRL_reliable_sidechains_arpWarpFinalRfree <= 0.45) )
	ok <- as.numeric( ok )
	data$successHelix2 <- replace( ok, is.na(ok), 0 )
	ok <-(data$helix_All_atom_shelxeCC >= 25 & data$helix_All_atom_shelxeAvgChainLength >= 10 & (data$helix_All_atom_buccFinalRfree <= 0.45 | data$helix_All_atom_arpWarpFinalRfree <= 0.45) )
	ok <- as.numeric( ok )
	data$successHelix3 <- replace( ok, is.na(ok), 0 )
	
	# CHECK
	data$helixWorked <- as.numeric( data$successHelix1 | data$successHelix2 | data$successHelix3 )
} # End comparison

# Calculate number of copies of decoy that were placed
data$numPlacedChains <- data$numPlacedAtoms / data$ensembleNumAtoms

#write(data,"foo.csv",row.NAMES=FALSE)


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


# Data for which there are placed models and we have an origin
odata = data[ ! is.na(data$ccmtzOrigin) & ! is.na( data$phaserLLG), ]

# Data for those that worked
sdata = data[ data$success == 1, ]

if (comparison){
	# Data for which we can compare the single-model casse
	smdata <- data[ data$single_model_ran == "True", ]
	
	#http://stackoverflow.com/questions/19357668/r-ggplot2-facetting-keep-ratio-but-override-define-output-plot-size
	# Comparison of Ensemble with Single model
	x <- smdata[ order( smdata$pdbCode, smdata$success, smdata$shelxeCC, decreasing=TRUE ), ]
	
	# Select top by selecting not duplicates on pdbCode
	x <- x[ !duplicated(x$pdbCode), c("pdbCode", "fastaLength","resolution","numChains","numPlacedChains", "shelxeCC", "shelxeAvgChainLength")  ]
	
	# Now put in alphabetical order
	x <- x[ order( x$pdbCode ), ]
	
	# Need to get numbers of success
	# This gets the stats  - we can join because the by function is the pdbCode which is in similar alphabetic order
	x["numModels"] <- aggregate( smdata$pdbCode, by=list(smdata$pdbCode), FUN=length )[2]
	
	# Different measures of success
	x["success"] <- aggregate( smdata$success, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	#x["worked"] <- aggregate( smdata$success, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	#x$worked <- replace( x$worked, x$worked > 0, 1 )
	x["worked"] <- replace( x$success, x$success > 0, 1 )
	x["successShelxe"] <- aggregate( smdata$successShelxe, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	x["successRefmac"] <- aggregate( smdata$successRefmac, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	x["successPhaser"] <- aggregate( smdata$successPhaser, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	x["successRebuild"] <- aggregate( smdata$successRebuild, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	
	x["successSingleStructure"] <- aggregate( smdata$successSingleStructure, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	#x["successHelix"] <- aggregate( smdata$successHelix, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	#x["helixWorked"] <- x$successHelix
	#x$helixWorked <- replace( x$helixWorked, x$helixWorked > 0, 1 )
	x["successHelix"] <- aggregate( smdata$helixWorked, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	x["helixWorked"] <- replace( x$successHelix, x$successHelix > 0, 1 )
	x["successHelix1"] <- aggregate( smdata$successHelix1, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	x["successHelix2"] <- aggregate( smdata$successHelix2, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	x["successHelix3"] <- aggregate( smdata$successHelix3, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	
	
	# NB LOOK AT REORDER
	# Now put in order by success, resolution
	summaryData <- x[ order( -x$worked, x$resolution ), ]
	#x <- x[ order( x$success, decreasing=TRUE ), ]
	write.csv(summaryData, "summaryComparison.csv", row.names=FALSE)

} # End comparison

#
# Summary information
#

pdbAll <- unique( data$pdbCode )
sprintf("Len all PDBs %s", length(pdbAll) )
pdbSuccess <- unique( data[ data$success==1, ]$pdbCode )
sprintf("Len success PDBs %s", length(pdbSuccess) )

smSuccess <- unique( data[ data$successSingleStructure==1, ]$pdbCode )
sprintf("Len single-model success PDBs %s", length(smSuccess) )

hSuccess <- unique( data[ data$helix_success==1, ]$pdbCode )
sprintf("Len helix success PDBs %s", length(hSuccess) )

# Failures
#setdiff( pdbAll, pdbSuccess )


# * mean resolution of failing/succesful cases
x <- mean(  data[ data$success == 1, ]$resolution ) # 1.882012
cat( "Mean resolution for successes is ",x,"\n")
x <- mean(  data[ data$success == 0, ]$resolution ) # 2.034655
cat( "Mean resolution for failures is ",x,"\n")

# * mean of solvent content of failing/succesful cases
#unique( data[ is.na( data$solventContent ),]$pdbCode ) # 1G1J 1KYC 1P9I 3CVF
x <-mean( data[ data$success == 1, ]$solventContent, na.rm=TRUE ) # 46.58613
cat( "Mean solvent content for successes is ",x,"\n")
x <-mean( data[ data$success == 0, ]$solventContent, na.rm=TRUE ) # 50.39069
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
summaryData <- data[ order( data$pdbCode, data$success, data$shelxeCC, decreasing=TRUE ), ]

# Select top by selecting not duplicates on pdbCode
summaryData <- summaryData[ !duplicated(summaryData$pdbCode), c("pdbCode","fastaLength","resolution","numChains",
				"numResidues", "shelxeCC", "shelxeAvgChainLength","shelxeNumChains")  ]

# Now put in alphabetical order
summaryData <- summaryData[ order( summaryData$pdbCode ), ]

# Need to get numbers of success
# This gets the stats  - we can join because the by function is the pdbCode which is in similar alphabetic order
summaryData["numModels"] <- aggregate( data$pdbCode, by=list(data$pdbCode), FUN=length )[2]
summaryData["worked"] <- aggregate( data$success, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData$worked <- replace( summaryData$worked, summaryData$worked > 0, 1 )


# Different measures of success
summaryData["success"] <- aggregate( data$success, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["successShelxe"] <- aggregate( data$successShelxe, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["successRefmac"] <- aggregate( data$successRefmac, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["successPhaser"] <- aggregate( data$successPhaser, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["successRebuild"] <- aggregate( data$successRebuild, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]

if (comparison){
	summaryData["successSingleStructure"] <- aggregate( data$successSingleStructure, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	summaryData["successHelix"] <- aggregate( data$helixWorked, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
}

# side-chain treatment
summaryData["successPolyAla"]  <- aggregate( 
		data[ data$ensembleSideChain == "poly_ala", ]$success,
		by=list( data[ data$ensembleSideChain == "poly_ala", ]$pdbCode ),
		FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["successScwrl"]  <- aggregate(
		data[ data$ensembleSideChain == "SCWRL_reliable_sidechains", ]$success,
		by=list( data[ data$ensembleSideChain == "SCWRL_reliable_sidechains", ]$pdbCode ),
		FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["successAllAtom"]  <- aggregate(
		data[ data$ensembleSideChain == "All_atom", ]$success,
		by=list( data[ data$ensembleSideChain == "All_atom", ]$pdbCode ),
		FUN=function(x){ sum( x == 1 ) } )[2]

# Subclustering radii
pdbCodes <- unique( data$pdbCode )

x <- aggregate( success~pdbCode,
		data=data,
		FUN=function(x){ sum( x == 1 ) },
		subset=data$ensembleRadiusThreshold==1 )
missing <- setdiff( pdbCodes, x$pdbCode)
y<-rep(0,length(missing)) # Make vector of zeros as long as the missing codes
y <- data.frame( missing, y) # Create a data frame and name it
names(y) <- c("pdbCode","success")
y<-rbind(x,y)
y <- y[ order(y$pdbCode), ]
summaryData$subclustering1A <- y$success

x <- aggregate( success~pdbCode,
		data=data,
		FUN=function(x){ sum( x == 1 ) },
		subset=data$ensembleRadiusThreshold==2 )
missing <- setdiff( pdbCodes, x$pdbCode)
y<-rep(0,length(missing)) # Make vector of zeros as long as the missing codes
y <- data.frame( missing, y) # Create a data frame and name it
names(y) <- c("pdbCode","success")
y<-rbind(x,y)
y <- y[ order(y$pdbCode), ]
summaryData$subclustering2A <- y$success

x <- aggregate( success~pdbCode,
		data=data,
		FUN=function(x){ sum( x == 1 ) },
		subset=data$ensembleRadiusThreshold==3 )
missing <- setdiff( pdbCodes, x$pdbCode)
y<-rep(0,length(missing)) # Make vector of zeros as long as the missing codes
y <- data.frame( missing, y) # Create a data frame and name it
names(y) <- c("pdbCode","success")
y<-rbind(x,y)
y <- y[ order(y$pdbCode), ]
summaryData$subclustering3A <- y$success

# ADD SINGLE MODEL AND HELIX
#x["successfulSingleModel"] <- aggregate( data$successSingleStructure, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
#x$successfulSingleModel[ x$successfulSingleModel >= 1 ] <- 1
#x["successfulHelix"] <- aggregate( data$helix_success, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
#x$successfulHelix[ x$successfulHelix >= 1 ] <- 1

# NB LOOK AT REORDER
# Now put in order by success, resolution
summaryData <- summaryData[ order( -summaryData$worked, summaryData$resolution ), ]
#x <- x[ order( x$success, decreasing=TRUE ), ]
write.csv(summaryData, "summary.csv", row.names=FALSE)

# Plot of targets by chain length
p <-ggplot( data=summaryData, aes( fastaLength, fill=factor(worked) ) ) 
p + geom_histogram( position = 'stack', binwidth = 10 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of targets") +
		xlab("Target chain length") +
		ggtitle("Histogram of target chain length for successful and failing cases")
ggsave("targetsByLength.png")

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
scols <- c("0" = "red","1" = "yellow","2" = "orange", "3" = "blue")

if (comparison){
	
	# Set all NA CC to 0
	smdata$shelxeCC <- replace( smdata$shelxeCC,is.na(smdata$shelxeCC),0)
	smdata$single_model_shelxeCC <- replace( smdata$single_model_shelxeCC,is.na(smdata$single_model_shelxeCC),0)
	smdata[,"ESresult"] <- NA
	smdata$ESresult[ smdata$success == 0 & smdata$successSingleStructure == 0 ] <- 0
	smdata$ESresult[ smdata$success == 1 & smdata$successSingleStructure == 0 ] <- 1
	smdata$ESresult[ smdata$success == 0 & smdata$successSingleStructure == 1 ] <- 2
	smdata$ESresult[ smdata$success == 1 & smdata$successSingleStructure == 1 ] <- 3
	
	# 1006 both succeed
	# 488 ensemble success
	# 165 single model success
	# 3627 both fail
	
	# Ensemble vs single model
	p <-ggplot(data=smdata,
			aes(shelxeCC,single_model_shelxeCC,colour=factor(ESresult)))
	p + geom_point(size=1.5) +
	stat_sum( aes(size=..n..) ) +
	scale_colour_manual( name="Result",
			values=scols,
			labels=c("Both Failed","Ensemble success", "Single-structure success", "Both succeeed")
	) +
	xlab("Ensemble Shelxe CC score") +
	ylab("Single centroid structure CC score") +
	ggtitle("Ensemble and Single-structure CC scores")
	ggsave("SSEShelxe.png")
	
	# Helix vs Single model
	# Need to pull out the different side chain treamtments
	# hacky way to get names to match - must be better way of doing this
	hdata <- data[ data$helix_ran == "True" & data$ensembleSideChainTreatment == "poly_ala",
			c("single_model_shelxeCC","successSingleStructure","helix_poly_ala_shelxeCC","successHelix1") ]
	names(hdata) <- c("single_model_shelxeCC","successSingleStructure","helix_shelxeCC","successHelix")
	
	x <- data[ data$helix_ran == "True" & data$ensembleSideChainTreatment == "SCWRL_reliable_sidechains",
			c("single_model_shelxeCC","successSingleStructure","helix_SCWRL_reliable_sidechains_shelxeCC","successHelix2") ]
	names(x) <- c("single_model_shelxeCC","successSingleStructure","helix_shelxeCC","successHelix")
	hdata <- rbind(hdata,x)
	
	x <- data[ data$helix_ran == "True" & data$ensembleSideChainTreatment == "All_atom",
			c("single_model_shelxeCC","successSingleStructure","helix_All_atom_shelxeCC","successHelix3") ]
	names(x) <- c("single_model_shelxeCC","successSingleStructure","helix_shelxeCC","successHelix")
	hdata <- rbind(hdata,x)
	
	
	hdata$single_model_shelxeCC <- replace( hdata$single_model_shelxeCC,is.na(hdata$single_model_shelxeCC),0)
	hdata$helix_shelxeCC <- replace( hdata$helix_shelxeCC,is.na(hdata$helix_shelxeCC),0)
	hdata[,"result"] <- NA
	hdata$result[ hdata$successSingleStructure == 0 & hdata$successHelix == 0 ] <- 0
	hdata$result[ hdata$successSingleStructure == 1 & hdata$successHelix == 0 ] <- 1
	hdata$result[ hdata$successSingleStructure == 0 & hdata$successHelix == 1 ] <- 2
	hdata$result[ hdata$successSingleStructure == 1 & hdata$successHelix == 1 ] <- 3
	
	
	# 1637 in total - dim(hdata)
	# 398 single model success
	# 127 helix worked
	# 617 both worked dim( hdata[ hdata$successSingleStructure==1 & hdata$successHelix==1,])
	
	p <-ggplot(data=hdata,
			aes(single_model_shelxeCC,helix_shelxeCC,colour=factor(result)))
	p + geom_point( size=1.5 ) +
			stat_sum( aes(size=..n..) ) +
			scale_colour_manual( name="Result",
					values=scols,
					labels=c("Both Fail","Single-structure success","Helix success","Both succeeed")
			) +
			scale_x_continuous(limits=c(0,80)) +
			scale_y_continuous(limits=c(0,80)) +
			xlab("Single-structure Shelxe CC score") +
			ylab("Idealised Helix Shelxe CC score") +
			ggtitle("Idealised Helices and Single-structure CC scores")
	
	ggsave("SSHShelxe.png")

} # END comparison

# How the different treatments fared - facet by resolution
# Number of success

#p <-ggplot(data=summaryData)
#p + geom_point(mapping=aes(x=pdbCode,y=success), colour="red") +
#geom_point(mapping=aes(x=pdbCode,y=successSingleStructure), colour="blue") +
#geom_point(mapping=aes(x=pdbCode,y=successHelix), colour="green")
#
#
#p <-ggplot(data=summaryData)
#p + geom_point(mapping=aes(x=resolution,y=success), colour="red") +
#geom_point(mapping=aes(x=resolution,y=successSingleStructure ), colour="blue") +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("Number of residues in ensemble") +
		ylab("Number of ensembles") +
		ggtitle("Number of residues in ensemble for successful and failing cases")
ggsave("ResiduesVsEnsemble.png")

p <-ggplot(data=data, aes(x=ensemblePercentModel, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("Percent of residues in ensemble") +
		ylab("Number of ensembles") +
		ggtitle("Percentage of residues in ensemble for successful and failing cases")
ggsave("PercentVsEnsemble.png")

# Number that solved under each side-chain treatment
p <-ggplot( data=data, aes( factor(ensembleSideChainTreatment), fill=factor(success) ) ) 
p + geom_histogram( position = 'dodge' ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Side-chain Treatment") +
		ggtitle("Histogram of side-chain treatment for successful and failing cases")
ggsave("sideChain.png")


# Number that solved under each sub-clustering radius
p <-ggplot( data=data, aes( factor(ensembleRadiusThreshold), fill=factor(success) ) ) 
p + geom_histogram( position = 'dodge' ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Sub-clustering radius") +
		ggtitle("Histogram of sub-clustering for successful and failing cases")
ggsave("subClusteringRadius.png")


# As a table
#x <- aggregate( data$success==0, by=list( SideChainTreatment = data$ensembleSideChainTreatment ), FUN=sum)
#names(x) <- sub("^x$", "failures", names(x))
#x$success <- aggregate( data$success==1, by=list( SideChainTreatment = data$ensembleSideChainTreatment ), FUN=sum)$x

# CC vs final rFree coloured by success
p <-ggplot(data=data, aes(shelxeCC, buccFinalRfree, colour=factor(successShelxe) ) )
#p <-ggplot(data=data, aes(shelxeCC, minRfree, colour=factor(successShelxe) ) )
p + geom_point() +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Shelxe Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		facet_wrap( ~resCatw) +
		xlab("Shelxe CC score") +
		ylab("Buccaneer Rfree") +
		ggtitle("Shelxe Buccanner RFree score")
ggsave("CCVsRfree.png")

# CC distribution
p <-ggplot(data=data, aes(x=shelxeCC, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("CC score") +
		ggtitle("Histogram of CC scores by success")
ggsave("CCscores.png")

		#coord_fixed() +
p <-ggplot(data=data[ data$success==1, ], aes(x=buccFinalRfree, y=arpWarpFinalRfree) )
p + geom_point() +
		#scale_x_continuous(limits=c(0,1)) +
		#scale_y_continuous(limits=c(0,1)) +
		#facet_grid( .~resCat, space="fixed", scales="fixed", labeller=l) +
		facet_wrap( ~resCatw, ncol=1) +
		coord_equal() +
		#theme(aspect.ratio = 1) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("Resolution (\uc5)") +
		ylab("Chain length (residues)") +
		ggtitle("Resolution vs chain length")
ggsave("resolutionVsLength.png")

# Binned (20) bar graph with length along the bottom and the bars dividied/coloured by success
p <-ggplot( data=data, aes( fastaLength, fill=factor(success) ) )
p + geom_histogram( binwidth=5, position = 'dodge' ) +
	scale_x_continuous( breaks=seq(0,260,20) ) +
	scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
	xlab("Length in residues") +
	ylab("Number of targets") +
	ggtitle("Successes/failures by target chain length")
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("TM-score of model to native") +
		ggtitle("Histogram of TM-score for models for successful and failing cases")
ggsave("modelTM.png")

p <-ggplot(data=data, aes(x=ensembleNativeRMSD, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("RIO score") +
		ggtitle("Histogram of RIO scores")
ggsave("rioHistogram.png")

p <-ggplot(data=odata, aes(x=ccmtzRioNoCat, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Number non-RIO C-alpha") +
		ggtitle("non-RIO contacts")
ggsave("noRioHistogram.png")


# Proportion RIO of model
p <-ggplot(data=odata, aes(x=ccmtzRioGood/numPlacedCA, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.05 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("In-register proportion of RIO") +
		ylab("Number of residues in ASU") +
		ggtitle("In-register CA as a prop. of RIO vs num. residues for RIO > 0")
ggsave("rioInRegisterPropVsLength.png")

# Proportion RIO of native
p <-ggplot(data=odata, aes(x=ccmtzRioGood/numResidues, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.05 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Proportion RIO") +
		ggtitle("RIO as proportion of native")
ggsave("rioNativeHistogramProp.png")

# Proportion in density of model
p <-ggplot(data=odata, aes(x=ccmtzAaNumContacts/numPlacedAtoms, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.01 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Proportion in density") +
		ggtitle("Atoms in density as a proportion of the model")
ggsave("inDensityModelHistogramProp.png")

# Proportion in density of native
p <-ggplot(data=odata, aes(x=ccmtzAaNumContacts/numAtoms, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.001 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Proportion in density") +
		ggtitle("Atoms in density as proportion of native")
ggsave("inDensityNativeHistogramProp.png")

# Number in density
p <-ggplot(data=odata, aes(x=ccmtzAaNumContacts, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Number in density") +
		ggtitle("Number of atoms in density")
ggsave("inDensityHistogramNum.png")


# Number out of density
p <-ggplot(data=odata, aes(x=numPlacedAtoms-ccmtzAaNumContacts, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Number out of density") +
		ggtitle("Number of atoms out of density")
ggsave("outOfDensityHistogramNum.png")

# Out of density as proportion of model
p <-ggplot(data=odata, aes(x=(numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.01 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Proportion out of density") +
		ggtitle("Atoms out of density as proportion of model")
ggsave("outOfDensityHistogramModel.png")

# Out of density as proportion of native
p <-ggplot(data=odata, aes(x=(numPlacedAtoms-ccmtzAaNumContacts)/numAtoms, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 0.01 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Proportion out of density") +
		ggtitle("Atoms out of density as proportion of native")
ggsave("outOfDensityHistogramNative.png")

# Proportion of In vs out of register contacts for success
#3767 success
#3766 with RIO score > 0
#3209 in ==0 out-of-register >0 - 85.19%
#		130 out == 0, in > 0 - 3.45%
#		400 - some in & out - 10.62%
#		28 solved but no scores. - 0.74%
#		1364 - OO == RIO



#dim(odata[ odata$success==1 & odata$ccmtzRioOoRegister==odata$ccmtzRioNumContacts, ])
#dim(odata[ odata$success==1 & odata$ccmtzRioOoRegister>0 & odata$ccmtzRioInregister==0, ]) - 3209
#dim(odata[ odata$success==1 & odata$ccmtzRioOoRegister == odata$ccmtzRioNumContacts, ]) - 1364
#p <-ggplot(data=odata, aes(x=ccmtzRioOoRegister/ccmtzRioNumContacts, fill=factor(success) ) )
# odata[ odata$success==1,]$ccmtzRioOoRegister / odata[ odata$success==1,]$ccmtzRioNumContacts 
p <-ggplot(data=odata[ odata$success==1, ], aes(x=ccmtzRioOoRegister/ccmtzRioGood) )
p + geom_histogram( binwidth = 0.01, fill="#3333FF" ) +
		ylab("Number of cases") +
		xlab("Proportion out-of-register") +
		ggtitle("Prop. of Out-of-register contacts for successes")
ggsave("propOutRegister.png")

# Length distribution of RIO scores
p <-ggplot(data=odata, aes(x=ccmtzRioGood, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1  ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("RIO score") +
		ggtitle("Histogram of RIO scores")
ggsave("rioDistribution.png")

p <-ggplot(data=odata[ odata$ccmtzRioGood >0, ], aes(x=ccmtzRioGood, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1  ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("RIO score") +
		ggtitle("Histogram of RIO scores (RIO >0)")
ggsave("rioDistributionNoZero.png")


# Length distribution of ideal helcies
p <-ggplot(data=odata, aes(x=rioLenHelix, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1  ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Length of RIO helix") +
		ggtitle("Distribution of RIO helices")
ggsave("lengthIdealHelices.png")


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
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("Out-of-register C-alpha") +
		ylab("In-register C-alpha") +
		ggtitle("In- vs, Out-of-register C-alphas")
ggsave("ooRegisterVsInRegister.png")

p <-ggplot(data=odata[odata$ccmtzRioOoRegister + odata$ccmtzRioInregister > 0, ], aes(x=ccmtzRioOoRegister, y=ccmtzRioInregister, colour=factor(success) ) ) 
p + geom_point() +
		stat_sum( aes(size=..n..) ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("Out-of-register C-alpha") +
		ylab("In-register C-alpha") +
		ggtitle("In- vs, Out-of-register C-alphas (In+Out > 0)")
ggsave("ooRegisterVsInRegisterGtZero.png")

# plot of nrInRegisterContacts vs nrGoodContacts
p <-ggplot(data=odata, aes(x=shelxeCC, y=ccmtzRioGood, colour=factor(success) ) )
p + geom_point() +
	scale_colour_manual( values=c(fcolour,scolour),
	name="Success/Failure",
	labels=c("Failure", "Success"),
	guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("Proportion  atoms in density (< 0.5A)") +
		ylab("RIO C-alphas as prop. of native (< 1.5A)") +
		ggtitle("RIO vs in-density as proportion of native")
ggsave("inDensityVsRIObyRes.png")


# Measure of how many correct vs number/proportion misplaced
# Misplaced atoms of model measured by: numPlacedAtoms - ccmtzAaNumContacts - proportion is (numPlacedAtoms-ccmtzAaNumContacts)/numAtoms
# Proportion correct measured by: ccmtzAaNumContacts/numAtoms
#                             or: ccmtzRioGood/numResidues
# so...
#p <-ggplot(data=odata,
#		aes(x=ccmtzAaNumContacts/numAtoms,
#			y=(numPlacedAtoms-ccmtzAaNumContacts)/numAtoms,
#			colour=factor(success) ) )
#p + geom_point( size=1 ) +
#		stat_sum( aes(size=..n..) ) +
#		facet_grid( resCat ~ success, labeller=l) +
#		scale_colour_manual( values=c(fcolour, scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success"),
#				guide=FALSE) +
#		xlab("AIO as prop. of native") +
#		ylab("Num. native atoms - AIO") +
#		ggtitle("")
#ggsave("inDensityVsMisplacedPropNative.png")

p <-ggplot(data=odata,
		aes(x=ccmtzAaNumContacts/numPlacedAtoms,
				y=(numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=l) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
				labels=c("Failure", "Success"),
				guide=FALSE) +
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
							labels=c("Failure", "Success"),
							guide=FALSE) +
		xlab("Phaser TFZ") +
		ylab("Phaser LLG") +
		ggtitle("Phaser LLG vs Phaser TFZ")
ggsave("LLGvsTFZ.png")
# Below for truncating the y-axis
# scale_y_continuous( limits=c(-3000, max(data$phaserLLG, na.rm=TRUE) )  )+
