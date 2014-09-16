library(ggplot2)
# http://sape.inf.usi.ch/quick-reference/ggplot2/themes
#theme(panel.background = element_rect(colour = "#FFFFFF"))
#theme_set(theme_bw())
theme_set(theme_bw()+theme(legend.position="none",
				axis.line=element_line(size = 0.5),
				axis.title=element_text(size=15),
				axis.text=element_text(size=12) )
)



scolour="#3333FF"
fcolour="#FF0000"
pcolour='#FF00FF'

comparison=FALSE
#data <- read.table(file="results.csv",sep=',', header=T)

if (comparison){
	data <- read.table(file="/home/jmht/Documents/work/CC/all_results/results.csv",sep=',', header=T)
} else {
	data <- read.table(file="/media/data/shared/coiled-coils/ensemble/final_results/final_results.csv",sep=',', header=T)
}

# Calculate number of copies of decoy that were placed
data$numPlacedChains <- data$numPlacedAtoms / data$ensembleNumAtoms

# Kill redundant columns
todrop <- c(
		"helix_ran",
		"helix_All_atom_shelxeCC",
		"helix_All_atom_shelxeAvgChainLength",
		"helix_All_atom_buccFinalRfact",
		"helix_All_atom_buccFinalRfree",
		"helix_All_atom_arpWarpFinalRfact",
		"helix_All_atom_arpWarpFinalRfree",
		"helix_SCWRL_reliable_sidechains_shelxeCC",
		"helix_SCWRL_reliable_sidechains_shelxeAvgChainLength",
		"helix_SCWRL_reliable_sidechains_buccFinalRfact",
		"helix_SCWRL_reliable_sidechains_buccFinalRfree",
		"helix_SCWRL_reliable_sidechains_arpWarpFinalRfact",
		"helix_SCWRL_reliable_sidechains_arpWarpFinalRfree",
		"helix_poly_ala_shelxeCC",
		"helix_poly_ala_shelxeAvgChainLength",
		"helix_poly_ala_buccFinalRfact",
		"helix_poly_ala_buccFinalRfree",
		"helix_poly_ala_arpWarpFinalRfact",
		"helix_poly_ala_arpWarpFinalRfree")

data <- data[, !(colnames(data) %in% todrop)]

updateData <- function(data){
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
	
	return(data)
}

data <- updateData(data)

if (!comparison){
	# Hack to add in data for manual cases - see /media/data/shared/coiled-coils/ensemble/manual
	# 3BAS
	data[ data$pdbCode=="3BAS"& data$ensembleName=="poly_ala_trunc_51.929237_rad_2", ]$buccFinalRfact <- 0.2947
	data[ data$pdbCode=="3BAS"& data$ensembleName=="poly_ala_trunc_51.929237_rad_2", ]$buccFinalRfree <- 0.436
	data[ data$pdbCode=="3BAS"& data$ensembleName=="poly_ala_trunc_51.929237_rad_2", ]$successRebuild <- 1
	data[ data$pdbCode=="3BAS"& data$ensembleName=="poly_ala_trunc_51.929237_rad_2", ]$success <- 1
	#3CVF
	data[ data$pdbCode=="3CVF"& data$ensembleName=="poly_ala_trunc_3.742193_rad_1", ]$success <- 1
}
	

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


# Data for which there are placed models and we have an origin
odata = data[ ! is.na(data$ccmtzOrigin) & ! is.na( data$phaserLLG), ]

#Calculate proportion of in/out in rioGood
#prop.table(as.matrix(x[-1]),margin=1)
odata$propRioIn <- prop.table(as.matrix( odata[, c("ccmtzRioInregister","ccmtzRioOoRegister") ]),margin=1)[,1]
odata$propRioOut <- prop.table(as.matrix( odata[, c("ccmtzRioInregister","ccmtzRioOoRegister") ]),margin=1)[,2]

# Data for those that worked
sdata = data[ data$success == 1, ]

if (comparison){
	
	# Read in single-model data
	smdata <- read.table(file="/home/jmht/Documents/work/CC/single_model_results/single_model_results.csv",sep=',', header=T)
	smdata <- updateData(smdata)
	
	data$single_model_ran <- 0
	data$single_model_ran <- match(data$ensembleName,smdata$ensembleName, nomatch=0)
	data$single_model_ran <- replace( data$single_model_ran, data$single_model_ran > 0, 1 )
	
	# Use ensembleName as match as there are no duplicates
	data$single_model_shelxeCC <- -1
	data$single_model_shelxeCC <- smdata[ match(data$ensembleName,smdata$ensembleName),"shelxeCC"]
	data$single_model_shelxeAvgChainLength <- -1
	data$single_model_shelxeAvgChainLength <- smdata[ match(data$ensembleName,smdata$ensembleName),"shelxeAvgChainLength"]
	data$single_model_buccFinalRfree <- -1
	data$single_model_buccFinalRfree <- smdata[ match(data$ensembleName,smdata$ensembleName),"buccFinalRfree"]
	data$single_model_buccFinalRfact <- -1
	data$single_model_buccFinalRfact <- smdata[ match(data$ensembleName,smdata$ensembleName),"buccFinalRfact"]
	data$single_model_arpWarpFinalRfree <- -1
	data$single_model_arpWarpFinalRfree <- smdata[ match(data$ensembleName,smdata$ensembleName),"arpWarpFinalRfree"]
	data$single_model_arpWarpFinalRfact <- -1
	data$single_model_arpWarpFinalRfact <- smdata[ match(data$ensembleName,smdata$ensembleName),"arpWarpFinalRfact"]
	
	# single-structure success
	data$successSingleStructure <- as.numeric( data$single_model_shelxeCC >= 25 & 
					data$single_model_shelxeAvgChainLength >= 10 &
					( data$single_model_buccFinalRfree <= 0.45 | data$single_model_arpWarpFinalRfree  <= 0.45 ) )
	data$successSingleStructure <- replace( data$successSingleStructure, is.na(data$successSingleStructure), 0 )
	
	data$single_model_phaserTime <- NA
	data$single_model_phaserTime <- smdata[ match(data$ensembleName,smdata$ensembleName),"phaserTime"]
	
	# Suck in the helix data
	hdata <- read.table(file="/home/jmht/Documents/work/CC/polya_helices/results.csv",sep=',', header=T)
	hdata <- updateData(hdata)
	
	#
	
	#	# helix success
	#	ok <-(data$helix_poly_ala_shelxeCC >= 25 & data$helix_poly_ala_shelxeAvgChainLength >= 10 & (data$helix_poly_ala_buccFinalRfree <= 0.45 | data$helix_poly_ala_arpWarpFinalRfree <= 0.45) )
	#	ok <- as.numeric( ok )
	#	data$successHelix1 <- replace( ok, is.na(ok), 0 )
	#	ok <-(data$helix_SCWRL_reliable_sidechains_shelxeCC >= 25 & data$helix_SCWRL_reliable_sidechains_shelxeAvgChainLength >= 10 & (data$helix_SCWRL_reliable_sidechains_buccFinalRfree <= 0.45 | data$helix_SCWRL_reliable_sidechains_arpWarpFinalRfree <= 0.45) )
	#	ok <- as.numeric( ok )
	#	data$successHelix2 <- replace( ok, is.na(ok), 0 )
	#	ok <-(data$helix_All_atom_shelxeCC >= 25 & data$helix_All_atom_shelxeAvgChainLength >= 10 & (data$helix_All_atom_buccFinalRfree <= 0.45 | data$helix_All_atom_arpWarpFinalRfree <= 0.45) )
	#	ok <- as.numeric( ok )
	#	data$successHelix3 <- replace( ok, is.na(ok), 0 )
	#	# CHECK
	#	data$helixWorked <- as.numeric( data$successHelix1 | data$successHelix2 | data$successHelix3 )

	# Data for which we can compare the single-model casse
	smdata <- data[ data$single_model_ran == 1, ]
	
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
	
	x["workedSS"] <- replace( x$successSingleStructure, x$successSingleStructure > 0, 1 )
	
	x["successHelix"] <- aggregate( hdata$success, by=list(hdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	
	x["workedH"] <- replace( x$successHelix, x$successHelix > 0, 1 )
#	x["successHelix"] <- aggregate( smdata$helixWorked, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
#	x["helixWorked"] <- replace( x$successHelix, x$successHelix > 0, 1 )
#	x["successHelix1"] <- aggregate( smdata$successHelix1, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
#	x["successHelix2"] <- aggregate( smdata$successHelix2, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
#	x["successHelix3"] <- aggregate( smdata$successHelix3, by=list(smdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	
	
	# NB LOOK AT REORDER
	# Now put in order by success, resolution
	summaryData <- x[ order( -x$worked, x$resolution ), ]
	#x <- x[ order( x$success, decreasing=TRUE ), ]
	write.csv(summaryData, "summaryComparisonSM.csv", row.names=FALSE)

} # End comparison

#write.csv(data[data$success==1 & data$successSingleStructure==1,], "jens.csv", row.names=FALSE)

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
summaryData <- summaryData[ !duplicated(summaryData$pdbCode), c("pdbCode","fastaLength","resolution","spaceGroup","numChains",
				"estChainsASU","matthewsCoefficient","numResidues","numAtoms", "shelxeCC", "shelxeAvgChainLength","shelxeNumChains")  ]

# For helix
#summaryData <- summaryData[ !duplicated(summaryData$pdbCode), c("pdbCode","polyaLength","resolution","shelxeCC", "shelxeAvgChainLength","shelxeNumChains")  ]

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
	#summaryData["successHelix"] <- aggregate( data$helixWorked, by=list(data$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
	
	# Read in the helix summary
	#hsdata <- read.table(file="/home/jmht/Documents/work/CC/polya_helices/summary.csv",sep=',', header=T)
	#summaryData$successHelix <- 0
	#summaryData$successHelix <- hsdata[ match(summaryData$pdbCode,hsdata$pdbCode), "success"]
	summaryData["successHelix"] <- aggregate( hdata$success, by=list(hdata$pdbCode), FUN=function(x){ sum( x == 1 ) } )[2]
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

#
# Information on extreme models
#

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

# Longest that solved with zero AIO score (and RIO)
x <- sdata[ sdata$ccmtzAaNumContacts== 0, ]
x[ x$fastaLength == max(x$fastaLength), c("pdbCode","ensembleName") ]
# 3M91 SCWRL_reliable_sidechains_trunc_0.066638_rad_1

# highest percentage of model that solved - selecting the largest
# Could use which.max but want to select largest
#x <- sdata[ sdata$ensemblePercentModel == max( sdata$ensemblePercentModel ), ]
#x[ order( x$fastaLength, decreasing=TRUE ), ][1,]
# 3K29 SCWRL_reliable_sidechains_trunc_168.103952_rad_3 100% 169 residues

# Longest model that solved
#x <- sdata[ sdata$ensembleNumResidues == max( sdata$ensembleNumResidues ), ]
#x[ order( x$fastaLength, decreasing=TRUE ), ][1,]
# 3K29 SCWRL_reliable_sidechains_trunc_168.103952_rad_3 169 residues

# Min
#x <- sdata[ sdata$ensemblePercentModel == min( sdata$ensemblePercentModel ), ]
#x[ order( x$fastaLength, decreasing=TRUE ), ][1,]
# 3OKQ poly_ala_trunc_0.202334_rad_1 4% 7 residues

# Shortest model that solved
#x <- sdata[ sdata$ensembleNumResidues == min( sdata$ensembleNumResidues ), ]
#x[ order( x$fastaLength, decreasing=TRUE ), ][1,]
# 1YOD SCWRL_reliable_sidechains_trunc_0.000777_rad_1 

q()
#
# Plots
#
scols <- c("0" = "red","1" = "magenta","2" = "orange", "3" = "blue")

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
	
	# Ensemble vs single model - exclude those with zero scores
	#p <-ggplot(data=smdata[ smdata$shelxeCC >0 & smdata$single_model_shelxeCC >0,  ],
	#stat_sum( aes(size=..n..) ) + # skip as Olga didn't like
	p <-ggplot(data=smdata,
			aes(shelxeCC,single_model_shelxeCC,colour=factor(ESresult)))
	p + geom_point(size=1.5) +
	scale_colour_manual( name="Result",
			values=scols,
			labels=c("Both Failed","Ensemble success", "Single-structure success", "Both succeeed")
	) +
	xlab("Ensemble Shelxe CC score") +
	ylab("Single centroid structure CC score")
	#theme(legend.position="none")
	#ggtitle("Ensemble and Single-structure CC scores")
	ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_singleStructureVsEnsembleCC.eps",scale=1.5)
	ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_singleStructureVsEnsembleCC.png")
	
	# Length distribution of ideal helcies
	p <-ggplot(data=hdata[ hdata$success==1,], aes(x=polyaLength) )
	#p + geom_histogram( position = 'dodge', binwidth = 1, colour=scolour  ) +
	p + geom_histogram( binwidth=1, fill=scolour  ) +
#			scale_fill_manual( values=c(fcolour,scolour),
#					name="Success/Failure",
#					labels=c("Failure", "Success"),
#					guide=FALSE) +
			ylab("Number of cases") +
			xlab("Length of ideal helix") +
			#ggtitle("Distribution of polya helices")
	ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_lengthIdealHelices.eps",scale=1.5)
	ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_lengthIdealHelices.png")

#	# Helix vs Single model
#	# Need to pull out the different side chain treamtments
#	# hacky way to get names to match - must be better way of doing this
#	hdata <- data[ data$helix_ran == "True" & data$ensembleSideChainTreatment == "poly_ala",
#			c("single_model_shelxeCC","successSingleStructure","helix_poly_ala_shelxeCC","successHelix1") ]
#	names(hdata) <- c("single_model_shelxeCC","successSingleStructure","helix_shelxeCC","successHelix")
#	
#	x <- data[ data$helix_ran == "True" & data$ensembleSideChainTreatment == "SCWRL_reliable_sidechains",
#			c("single_model_shelxeCC","successSingleStructure","helix_SCWRL_reliable_sidechains_shelxeCC","successHelix2") ]
#	names(x) <- c("single_model_shelxeCC","successSingleStructure","helix_shelxeCC","successHelix")
#	hdata <- rbind(hdata,x)
#	
#	x <- data[ data$helix_ran == "True" & data$ensembleSideChainTreatment == "All_atom",
#			c("single_model_shelxeCC","successSingleStructure","helix_All_atom_shelxeCC","successHelix3") ]
#	names(x) <- c("single_model_shelxeCC","successSingleStructure","helix_shelxeCC","successHelix")
#	hdata <- rbind(hdata,x)
#	
#	
#	hdata$single_model_shelxeCC <- replace( hdata$single_model_shelxeCC,is.na(hdata$single_model_shelxeCC),0)
#	hdata$helix_shelxeCC <- replace( hdata$helix_shelxeCC,is.na(hdata$helix_shelxeCC),0)
#	hdata[,"result"] <- NA
#	hdata$result[ hdata$successSingleStructure == 0 & hdata$successHelix == 0 ] <- 0
#	hdata$result[ hdata$successSingleStructure == 1 & hdata$successHelix == 0 ] <- 1
#	hdata$result[ hdata$successSingleStructure == 0 & hdata$successHelix == 1 ] <- 2
#	hdata$result[ hdata$successSingleStructure == 1 & hdata$successHelix == 1 ] <- 3
#	
#	
#	# 1637 in total - dim(hdata)
#	# 398 single model success
#	# 127 helix worked
#	# 617 both worked dim( hdata[ hdata$successSingleStructure==1 & hdata$successHelix==1,])
#	
#	p <-ggplot(data=hdata,
#			aes(single_model_shelxeCC,helix_shelxeCC,colour=factor(result)))
#	p + geom_point( size=1.5 ) +
#			stat_sum( aes(size=..n..) ) +
#			scale_colour_manual( name="Result",
#					values=scols,
#					labels=c("Both Fail","Single-structure success","Helix success","Both succeeed")
#			) +
#			scale_x_continuous(limits=c(0,80)) +
#			scale_y_continuous(limits=c(0,80)) +
#			xlab("Single centroid structure Shelxe CC score") +
#			ylab("Idealised Helix Shelxe CC score") +
#			ggtitle("Idealised Helices and Single-structure CC scores")
#	
#	ggsave("SSHShelxe.eps",scale=1.5)

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

if (!comparison) {
	# skip when not comparison as didn't calcualte the timings
	# Timings
	# proportion for all runs where there was a shelxe build
	
	# For some phaser runs the logs are not complete so we have no timing data - we set these to 0
	data$phaserTime <- replace( data$phaserTime, is.na(data$phaserTime), 0 )
	data$shelxeTime <- replace( data$shelxeTime, is.na(data$shelxeTime), 0 )
	
	# modelling etc is same for all models so don't want to average across jobs, but targets
	tsum <- data[ !duplicated(data$pdbCode), c("pdbCode","fragmentTime","modelTime","ensembleTime")]
	mr <- aggregate( data[ , c("phaserTime", "shelxeTime") ], by=list( data$pdbCode ), FUN=mean)
	names(mr)[1] <- "pdbCode" # unnecessary - just for REM
	tsum$phaserTime <- mr$phaserTime
	tsum$shelxeTime <- mr$shelxeTime
	
	write.csv(tsum, "timingData.csv", row.names=FALSE)
	x <- c( mean(tsum$fragmentTime), mean(tsum$modelTime), mean(tsum$ensembleTime), mean(tsum$phaserTime), mean(tsum$shelxeTime) )
	
	lbls <- c("fragmentTime","modelTime","ensembleTime","phaserTime", "shelxeTime")
	lbls <- paste(lbls, round(x), sep="\n") # add times to labels 
	png("timingsPie.png")
	pie( as.numeric( x ),
			labels=lbls,
			col=rainbow(length(lbls))
	)
	#		main="Timings for all runs inc. shelxe" )
	dev.off()
}


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
		xlab("Number of residues per chain in search model") +
		ylab("Number of search models") +
		#ggtitle("Number of residues in ensemble for successful and failing cases")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_residuesPerEnsemble.png")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_residuesPerEnsemble.eps",scale=1.5)

p <-ggplot(data=data, aes(x=ensemblePercentModel, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 5 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("% of residues per chain in search model") +
		ylab("Number of search models") +
		#ggtitle("Percentage of residues in ensemble for successful and failing cases")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_residuesPerEnsemblePercent.png")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_residuesPerEnsemblePercent.eps",scale=1.5)

# Number that solved under each side-chain treatment
p <-ggplot( data=data, aes( factor(ensembleSideChainTreatment), fill=factor(success) ) ) 
p + geom_histogram( position = 'dodge', labeller=labeller ) +
		scale_x_discrete(limits=c("All_atom","SCWRL_reliable_sidechains","poly_ala"), # Change order
				labels=c("All_atom"="All atom", "poly_ala"="PolyA","SCWRL_reliable_sidechains"="Reliable")) + # Change labels
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Side-chain Treatment") +
		#ggtitle("Histogram of side-chain treatment for successful and failing cases")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_sideChainHistogram.png")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_sideChainHistogram.eps",scale=1.5)

# Number that solved under each sub-clustering radius
p <-ggplot( data=data, aes( factor(ensembleRadiusThreshold), fill=factor(success) ) ) 
p + geom_histogram( position = 'dodge' ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Sub-clustering radius") +
		#ggtitle("Histogram of sub-clustering for successful and failing cases")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_subClusteringRadiusHistogram.eps",scale=1.5)
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_subClusteringRadiusHistogram.png")


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
		#facet_wrap( ~resCatw) +
		xlab("Shelxe CC score") +
		ylab("Buccaneer Rfree") +
		ggtitle("Shelxe Buccanner Rfree score")
ggsave("buccaneerRfreeVsCC.eps",scale=1.5)

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
ggsave("CCHistogram.eps",scale=1.5)

#coord_fixed() +
p <-ggplot(data=data[ data$success==1, ], aes(x=buccFinalRfree, y=arpWarpFinalRfree) )
p + geom_point() +
		scale_x_continuous(limits=c(0,1)) +
		scale_y_continuous(limits=c(0,1)) +
		#facet_grid( .~resCat, space="fixed", scales="fixed", labeller=l) +
		#facet_wrap( ~resCatw, nrow=1) +
		#coord_equal() +
		#theme(aspect.ratio = 1) +
		#stat_sum( aes(size=..n..) ) +
		xlab("Buccaneer Rfree") +
		ylab("ARP/wARP Rfree") +
		ggtitle("Bucc vs Arpwarp")
ggsave("buccaneerVsArpwarpRfree.eps",scale=1.5)


## Resolution vs length1
#p <-ggplot(data=data, aes(x=resolution, y=fastaLength, colour=factor(success) ) )
#p + geom_point( size=0, alpha=0 ) +
#		stat_sum( aes(size=..n..), shape=1 ) +
#		scale_colour_manual( values=c(fcolour,scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success"),
#				guide=FALSE) +
#		xlab("Resolution (\uc5)") +
#		ylab("Chain length (residues)") +
#		ggtitle("Resolution vs chain length")
#ggsave("resolutionVsLength.eps",scale=1.5)

# Resolution vs length2

# Split into categories
# - all succcess
# - all failures
# above as filled blobs
# targets where success < failure - failures as big red circules - successs as filled blue circles
# targets where success > failure - succcess as big red circules - failures as filled red circles
summaryData$nfailed <- summaryData$numModels - summaryData$success

# Open circles sized by the number of models and coloured by success. Filled blue circles indicate # successes
p <- ggplot()
p + geom_point( data=summaryData[ summaryData$success == 0, ],
				aes(x=resolution, y=fastaLength, size=numModels), shape=1, colour=fcolour ) +
	# Success as open circles coloured blue
	geom_point( data=summaryData[ summaryData$success > 0, ],
			aes(x=resolution, y=fastaLength, size=numModels), shape=1, colour=scolour ) +
	geom_point( data=summaryData[ summaryData$success > 0, ],
			aes(x=resolution, y=fastaLength, size=success), colour=scolour)	+
	xlab("Resolution (\uc5)") +
	ylab("Target chain length in residues") +
#	theme(legend.position="none",
#			axis.line=element_line(size = 0.5),
#			axis.title=element_text(size=15),
#			axis.text=element_text(size=12) )
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_resolutionVsChainLength.eps",scale=1.5)
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_resolutionVsChainLength.png")
		
# Open circles sized by the number of models and coloured by success. Filled blue circles indicate # successes
p <- ggplot()
p + geom_point( data=summaryData[ summaryData$success == 0, ],
				aes(x=resolution, y=numResidues, size=numModels), shape=1, colour=fcolour ) +
		# Success as open circles coloured blue
		geom_point( data=summaryData[ summaryData$success > 0, ],
				aes(x=resolution, y=numResidues, size=numModels), shape=1, colour=scolour ) +
		geom_point( data=summaryData[ summaryData$success > 0, ],
				aes(x=resolution, y=numResidues, size=success), colour=scolour)	+
		xlab("Resolution (\uc5)") +
		ylab("Num. residues in ASU") 
		#theme(legend.position="none")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_resolutionVsNumResidues.eps",scale=1.5)
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_resolutionVsNumResidues.png")


#p <- ggplot()
## All failures
#p + geom_point( data=summaryData[ summaryData$success == 0, ],
#		        aes(x=resolution, y=fastaLength, size=numModels), colour=fcolour ) +
#	# Successes with successes > failures
#    geom_point( data=summaryData[ summaryData$success-summaryData$nfailed > 0, ],
#		aes(x=resolution, y=fastaLength, size=success), colour=scolour) +
#    geom_point( data=summaryData[ summaryData$success-summaryData$nfailed > 0, ],
#		aes(x=resolution, y=fastaLength, size=nfailed), colour=fcolour) +
#	# successes with successes < failures
#	geom_point( data=summaryData[ summaryData$success > 0 & summaryData$success-summaryData$nfailed < 0, ],
#			aes(x=resolution, y=fastaLength, size=success), colour=scolour) +
#	geom_point( data=summaryData[ summaryData$success > 0 & summaryData$success-summaryData$nfailed < 0, ],
#			aes(x=resolution, y=fastaLength, size=nfailed), shape=1, colour=scolour) +
#	theme(legend.position="none")
#ggsave("resolutionVsFastaLengthLayered.eps",scale=1.5)



# simple
#p <-ggplot(data=summaryData, aes(x=resolution, y=fastaLength, colour=factor(worked) ) )
#p + geom_point( aes(size=numModels), guide=FALSE  ) +
#		scale_colour_manual( values=c(fcolour,scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success"),
#				guide=FALSE) +
#		xlab("Resolution (\uc5)") +
#		ylab("Chain length (residues)") +
#		ggtitle("Resolution vs chain length")
#ggsave("resolutionVsLength2.eps",scale=1.5)
#
#p <-ggplot(data=summaryData, aes(x=resolution, y=numResidues, colour=factor(worked) ) )
#p + geom_point( aes(size=numModels), guide=FALSE  ) +
#		scale_colour_manual( values=c(fcolour,scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success"),
#				guide=FALSE) +
#		xlab("Resolution (\uc5)") +
#		ylab("Residues in ASU") +
#		ggtitle("Resolution vs residues in ASU")
#ggsave("resolutionVsOverallLength.eps",scale=1.5)
#
#p <-ggplot(data=summaryData, aes(x=resolution, y=numResidues, colour=factor(worked) ) )
#p + geom_point(guide=FALSE) +
#		scale_colour_manual( values=c(fcolour,scolour),
#				name="Success/Failure",
#				labels=c("Failure", "Success"),
#				guide=FALSE) +
#		xlab("Resolution (\uc5)") +
#		ylab("Residues in ASU") +
#		ggtitle("Resolution vs residues in ASU")
#ggsave("resolutionVsOverallLength2.eps",scale=1.5)


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
ggsave("targetsByLength.eps",scale=1.5)

p <-ggplot( data=summaryData, aes( numResidues, fill=factor(worked) ) )
p + geom_histogram() +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		#scale_y_continuous(limits=c(0.10),breaks=c(0,10,1))+
		ylab("Number of targets") +
		xlab("Residues in ASU") +
		ggtitle("Histogram of target size for successful and failing cases")
ggsave("targetsByResASU.eps",scale=1.5)

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
ggsave("targetByLength.eps",scale=1.5)


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
#ggsave("shelxePropVsCC.eps",scale=1.5)


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
		xlab("TM-score of model to crystal structure") +
		#ggtitle("Histogram of TM-score for models for successful and failing cases")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_modelTM.png")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_modelTM.eps",scale=1.5)

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
ggsave("modelRMSD.eps",scale=1.5)

# Plot of Reforigin RMSD for successes and failures
p <-ggplot(data=data, aes(x=reforiginRMSD, fill=factor(success) ) )
#p + geom_histogram(alpha = 0.5, position = 'identity', binwidth = 5 ) +
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("REFORIGIN RMSD of placed model to crystal structure") +
		geom_vline(aes(xintercept=3), colour="#000000", linetype="dashed")
		#ggtitle("Histogram of reforigin RMSD scores by success")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_reforiginRMSD.png")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_reforiginRMSD.eps",scale=1.5)

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
ggsave("rioHistogram.eps",scale=1.5)

p <-ggplot(data=odata, aes(x=ccmtzRioNoCat, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1 ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("Number non-RIO C-alpha") +
		ggtitle("non-RIO contacts")
ggsave("noRioHistogram.eps",scale=1.5)


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
ggsave("rioModelHistogramProp.eps",scale=1.5)


#odata$propRioIn <- replace( odata$propRioIn, is.na(odata$propRioIn), -1 )

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
ggsave("rioNativeHistogramProp.eps",scale=1.5)

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
ggsave("inDensityModelHistogramProp.eps",scale=1.5)

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
ggsave("inDensityNativeHistogramProp.eps",scale=1.5)

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
ggsave("inDensityHistogramNum.eps",scale=1.5)


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
ggsave("outOfDensityHistogramNum.eps",scale=1.5)

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
ggsave("outOfDensityHistogramModel.eps",scale=1.5)

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
ggsave("outOfDensityHistogramNative.eps",scale=1.5)

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
		#xlab("Proportion of RIO_out C") +
p <-ggplot(data=odata[ odata$success==1, ], aes(x=(ccmtzRioOoRegister/ccmtzRioGood)*100) )
p + geom_histogram( binwidth = 1, fill="#3333FF" ) +
		ylab("Number of cases") +
		xlab(expression(paste("% RIO_out C",alpha))) +
		#ggtitle(expression(paste("Percentage of RIO_out C",alpha," for successes")))
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_percentOutRegister.eps",scale=1.5)
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_percentOutRegister.png")

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
ggsave("rioDistribution.eps",scale=1.5)

p <-ggplot(data=odata[ odata$ccmtzRioGood >0, ], aes(x=ccmtzRioGood, fill=factor(success) ) )
p + geom_histogram( position = 'dodge', binwidth = 1  ) +
		scale_fill_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		ylab("Number of cases") +
		xlab("RIO score") +
		ggtitle("Histogram of RIO scores (RIO >0)")
ggsave("rioDistributionNoZero.eps",scale=1.5)


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
ggsave("lengthRioHelices.eps",scale=1.5)


#
# Graphs
#

labeller = function( variable, value ) {
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

#Comparison of the RIO figures for successful search models with the number of in-register 
# residues similarly defined.
p <-ggplot(data=odata, aes(x=ccmtzRioOoRegister, y=ccmtzRioInregister, colour=factor(success) ) ) 
p + geom_point() +
		stat_sum( aes(size=..n..) ) +
		scale_colour_manual( values=c(fcolour,scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("RIO_out") +
		ylab("RIO_in") +
		ggtitle("RIO_in versus RIO_out")
ggsave("rioInVsRioOut.eps",scale=1.5)

p <-ggplot(data=odata[odata$ccmtzRioOoRegister + odata$ccmtzRioInregister > 0, ],
		aes(x=ccmtzRioOoRegister, y=ccmtzRioInregister, colour=factor(success) ) ) 
p + 
	#geom_point() +
	stat_sum( aes(size=..n..) ) +
	scale_colour_manual( values=c(fcolour,scolour),
			guide=FALSE
			#name="Success/Failure",
			#labels=c("Failure", "Success"),
			) +
	xlab("RIO_out") +
	ylab("RIO_in") +
	scale_size(name="Number of\nsearch models",
			breaks=c(1,10,100,300),
			range=c(1,10) ) +
	theme(legend.position="right")
	#ggtitle("RIO_in versus RIO_out (RIO > 0)")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_rioInVsRioOutGtZero.png")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_rioInVsRioOutGtZero.eps",scale=1.5)

#p <-ggplot(data=odata[ odata$ccmtzRioGood > 0, ], aes(x=propRioIn, fill=factor(success)) )
#p + geom_histogram( binwidth = 0.05, position='dodge' ) 
p <-ggplot(data=odata[ odata$ccmtzRioGood > 0, ],
		aes(x=propRioIn,
				y=numResidues,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( .~ success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("In-register proportion of RIO") +
		ylab("Number of residues in ASU") +
		ggtitle("In-register CA as a prop. of RIO vs num. residues for RIO > 0")
ggsave("rioInRegisterPropVsLength.eps",scale=1.5)


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
ggsave("CCVsRIO.eps",scale=1.5)

#
# inDensity vs unmatched
# ronan's idea to divide by number of unmatched
#
# facet_grid( ensembleSideChainTreatment ~ .) +
# label_both


p <-ggplot(data=odata,
		aes(x=ccmtzAaNumContacts,
				y=ccmtzRioGood,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("AIO") +
		ylab("RIO") +
		ggtitle("AIO vs RIO")
ggsave("aioVsRio.eps",scale=1.5)


p <-ggplot(data=odata,
		aes(x=(ccmtzAaNumContacts/numAtoms)*100,
				y=(ccmtzRioGood/numResidues)*100,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab(" % AIO") +
		ylab("% RIO") +
		ggtitle("Percentage AIO vs RIO")
ggsave("percentAioVsRio.eps",scale=1.5)

p <-ggplot(data=odata,
		aes(x=numPlacedAtoms-ccmtzAaNumContacts,
				y=ccmtzRioGood,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("Number of non-AIO model atoms") +
		ylab("RIO score") +
		ggtitle("Non-AIO atoms against RIO score")
ggsave("nonAioVsRIO.eps",scale=1.5)


# Number outside of density as proportion of model shows negative signal as proportion of the input signal
# Number in density shows what actually matched
p <-ggplot(data=odata,
		aes(x=(ccmtzRioGood/numResidues)*100,
				y=((numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms)*100,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("% RIO of native") +
		ylab("% non-AIO atoms of model") +
		ggtitle("RIO/non-AIO as percentages")
ggsave("rioPercentAioPercentFres.eps",scale=1.5)

p <-ggplot(data=odata,
		aes(x=(ccmtzRioGood/numResidues)*100,
				y=((numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms)*100,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid(  ~success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("% RIO of native") +
		ylab("% non-AIO atoms of model") +
		ggtitle("RIO/non-AIO as percentages")
ggsave("rioPercentAioPercent.eps",scale=1.5)

p <-ggplot(data=odata,
		aes(x=(ccmtzRioGood/numResidues)*100,
				y=(numPlacedAtoms-ccmtzAaNumContacts),
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("%RIO of native") +
		ylab("Number non-AIO atoms") +
		ggtitle("RIO as percentage vs number non-AIO")
ggsave("rioPercentAioNumFres.eps",scale=1.5)


p <-ggplot(data=odata,
		aes(x=ccmtzRioGood,
				y=((numPlacedAtoms-ccmtzAaNumContacts)/numPlacedAtoms)*100,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("RIO score") +
		ylab("% non-AIO atoms of model") +
		ggtitle("RIO vs non-AIO as percentage of model")
ggsave("rioNumAioPercentFres.eps",scale=1.5)

p <-ggplot(data=odata,
		aes(x=ccmtzRioGood,
				y=numPlacedAtoms-ccmtzAaNumContacts,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid( resCat ~ success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("RIO score") +
		ylab("Number non-AIO atoms") +
		ggtitle("RIO vs number non-AIO atoms")
ggsave("rioNumAioNumFres.eps",scale=1.5)

p <-ggplot(data=odata,
		aes(x=ccmtzRioGood,
				y=numPlacedAtoms-ccmtzAaNumContacts,
				colour=factor(success) ) )
p + geom_point( size=1 ) +
		stat_sum( aes(size=..n..) ) +
		facet_grid(  ~success, labeller=labeller) +
		scale_colour_manual( values=c(fcolour, scolour),
				name="Success/Failure",
				labels=c("Failure", "Success"),
				guide=FALSE) +
		xlab("RIO score") +
		ylab("Number non-AIO atoms") +
		ggtitle("RIO vs number non-AIO atoms")
ggsave("rioNumAioNum.eps",scale=1.5)



###############################################################################################################################
#
# TFZ/LLG
#
###############################################################################################################################
#
# Split plot so successes in front
p <-ggplot()
p + geom_point(data=data[ ! is.na( data$phaserLLG) & ! is.na( data$phaserTFZ) & data$success==0,], 
		aes(x=phaserTFZ, y=phaserLLG ),
		colour=fcolour
		) + 
		geom_point(data=data[ ! is.na( data$phaserLLG) & ! is.na( data$phaserTFZ) & data$success==1,], 
				aes(x=phaserTFZ, y=phaserLLG ),
				colour=scolour
		) + 	
		xlab("Phaser TFZ") +
		ylab("Phaser LLG") +
		#ggtitle("Phaser LLG vs Phaser TFZ")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_LLGvsTFZ.png")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_LLGvsTFZ.eps",scale=1.5)

# Below for truncating the y-axis
# scale_y_continuous( limits=c(-3000, max(data$phaserLLG, na.rm=TRUE) )  )+

###############################################################################################################################
#
# Additionl stuff
#
###############################################################################################################################


#scolour="#3333FF"
#fcolour="#FF0000"
#library("ggplot2")

# For analysing polyalanine helix results - source updateData from above
hsdata <- read.table(file="/home/jmht/Documents/work/CC/polya_helices/summary.csv",sep=',', header=T)

p <- ggplot()
p + geom_point( data=hsdata,
				aes(x=resolution, y=fastaLength, colour=factor(worked)), shape=16) +
		scale_colour_manual( values=c(fcolour, scolour)) +
		xlab("Resolution (\uc5)") +
		ylab("Target chain length in residues") +
		theme(legend.position="none")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_polyaResultsChain.png")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_polyaResultsChain.eps",scale=1.5)

p <- ggplot()
p + geom_point( data=hsdata,
				aes(x=resolution, y=numResidues, colour=factor(worked)), shape=16) +
		scale_colour_manual( values=c(fcolour, scolour)) +
		xlab("Resolution (\uc5)") +
		ylab("Residues in ASU") +
		theme(legend.position="none")
ggsave("polyaResultsASU.eps",scale=1.5)

#q()
library("ggplot2")
library("VennDiagram")
# Below for comparison of single runs (/home/jmht/Documents/work/CC/all_results) or comparisonResults.csv in google
scolour="#3333FF"
fcolour="#FF0000"
pcolour='#FF00FF'
data <- read.table(file="/home/jmht/Documents/work/CC/all_results/comparisonResults.csv",sep=',', header=T)

rerunSolved <- c("1MI7", "2BEZ", "2Q5U", "2ZZO", "3H00", "3H7Z", "3TYY", "3U1A", "3U1C")
#manual <- c("1M3W", "3BAS", "3CVF") # No manual ones solved with anything else
manual <- c("3BAS", "3CVF") # No manual ones solved with anything else
onlySingle <- c("4DZK")

all <- data$pdbCode
esolved <- data[ data$success >0, ]$pdbCode
ssolved <- data[ data$successSingleStructure >0, ]$pdbCode
hsolved <- data[ data$successHelix >0, ]$pdbCode

# All ones that solved
allsolved <- union(esolved,union(ssolved,hsolved)) # 70

# Ones that didn't solve
failed <- setdiff(all,allsolved) # 24

# ones that solved only with ensembles
onlye <- esolved[ ! esolved %in% union(ssolved,hsolved) ] # 1YBK 2B22 2YKT 1X8Y 1D7M 1KQL

# ones that only solved with helices
onlyh <- hsolved[ ! hsolved %in% union(esolved,ssolved) ] # 3H00 - solved in rerun

# ones that only solved with single structures
onlys <- ssolved[ ! ssolved %in% union(esolved,hsolved) ] # 2Q5U 4DZK 3U1C # 4DZK only solved with the single structure


allFailed <- setdiff(all, union(allsolved,union(rerunSolved,manual) ))
# Ones that solved by everything - excluding reruns (44)
everything <- intersect(intersect(esolved, ssolved), hsolved)

everythingRerun <- union(everything,rerunSolved)

# Bug with postscript or venn?
#postscript("z_finalResultsVenn.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
pdf("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_finalResultsVenn.pdf")
#par(oma=c(2,2,2,2), mar=c(2,2,2,2)) 
venn.plot <- draw.triple.venn(
		area1 = length(esolved),
		area2 = length(ssolved),
		area3 = length(hsolved),
		n12 = length(intersect(esolved,ssolved)),
		n23 = length(intersect(ssolved,hsolved)),
		n13 = length(intersect(esolved,hsolved)),
		n123 = length(everything),
		category = c("Ensembles", "Single\nStructures ", "Ideal\nHelices"),
		euler.d = FALSE,
		scaled=FALSE,
		fill = c("blue", "red", "green"),
		cex = 1.5, # Numbers in cirlces
		#cat.pos = c(330,30,180),
		cat.dist = c(0.05,0.07,0.05),
		cat.cex = 1.1, # labels
		cat.fontface = rep("plain", 3),
		cat.fontfamily = rep("sans", 3)
		#cat.col = c("blue", "red", "green")
		)
#grid.draw(venn.plot)
dev.off()


# Categories
# - solved with something -1
# - failed - 0
# - solved with everything - 1
# - only helix - 2
# - only single - 3
# - only ensemble - 4
# - failed

data$howSolved <- 0 # default - solved with something
data[ data$pdbCode %in%  allFailed, ]$howSolved <- -1 # didn't solve
data[ data$pdbCode %in%  everythingRerun, ]$howSolved <- 1 # solved with everything - including the reruns
data[ data$pdbCode %in%  rerunSolved, ]$howSolved <- 2 # just rerun
data[ data$pdbCode %in%  manual, ]$howSolved <- 3 # manual
data[ data$pdbCode %in%  onlySingle, ]$howSolved <- 4 # only with single-structure

p <- ggplot()
p + geom_point( data=data[ data$howSolved==-1, ],
				aes(x=resolution, y=fastaLength), shape=2, colour=fcolour ) +
		geom_point( data=data[ data$howSolved==0, ],
				aes(x=resolution, y=fastaLength), shape=1, colour=scolour ) + # 1 is empty cirlce
		geom_point( data=data[ data$howSolved==1, ],
				aes(x=resolution, y=fastaLength), shape=16, colour=scolour ) +
		geom_point( data=data[ data$howSolved==2, ],
				aes(x=resolution, y=fastaLength), shape=1, colour=pcolour ) +
		geom_point( data=data[ data$howSolved==3, ],
				aes(x=resolution, y=fastaLength), shape=0, colour=pcolour ) +
		geom_point( data=data[ data$howSolved==4, ],
				aes(x=resolution, y=fastaLength), shape=13, colour=scolour ) +
		xlab("Resolution (\uc5)") +
		ylab("Target chain length in residues") +
		theme(legend.position="none")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_finalSummaryChainLength.eps",scale=1.5)
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/z_finalSummaryChainLength.png")

p <- ggplot()
p + geom_point( data=data[ data$howSolved==-1, ],
				aes(x=resolution, y=numResidues), shape=2, colour=fcolour ) +
		geom_point( data=data[ data$howSolved==0, ],
				aes(x=resolution, y=numResidues), shape=1, colour=scolour ) + # 1 is empty cirlce
		geom_point( data=data[ data$howSolved==1, ],
				aes(x=resolution, y=numResidues), shape=16, colour=scolour ) +
		geom_point( data=data[ data$howSolved==2, ],
				aes(x=resolution, y=numResidues), shape=1, colour=pcolour ) +
		geom_point( data=data[ data$howSolved==3, ],
				aes(x=resolution, y=numResidues), shape=0, colour=pcolour ) +
		geom_point( data=data[ data$howSolved==4, ],
				aes(x=resolution, y=numResidues), shape=13, colour=scolour ) +
		xlab("Resolution (\uc5)") +
		ylab("Num. residues in ASU") +
		theme(legend.position="none")
ggsave("/home/jmht/Dropbox/PHD/CoiledCoilPaper/finalSummaryNumResidues.eps",scale=1.5)

