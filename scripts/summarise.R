

cls <- c(SHELXE_CC="numeric",
		SHELXE_ACL="numeric",
		SXRARP_final_Rfree="numeric",
		SXRBUCC_final_Rfree="numeric",
		final_Rfree="numeric",
		PHASER_TFZ="numeric",
		PHASER_LLG="numeric"
)

f=allresults.csv
data <- read.csv(file=f,sep=',', header=T, stringsAsFactors=FALSE, na.strings="N/A", colClasses=cls)

# Convert to numeric
# Need to replace all N/As with none
data[data=="N/A"]<-NA

updateData <- function(data){
	data$SHELXE_OK <- as.numeric( data$SHELXE_CC >= 25 & data$SHELXE_ACL >= 10 )
	data$SHELXE_OK <- replace( data$SHELXE_OK, is.na(data$SHELXE_OK), 0 )
	
	# Definition of success where at least one rebuild worked 
	data$minRfree <- pmin( data$SXRARP_final_Rfree, data$SXRBUCC_final_Rfree, na.rm=TRUE )
	data$REBUILD_OK <- as.numeric( data$minRfree <= 0.45 )
	data$REBUILD_OK <- replace( data$REBUILD_OK, is.na(data$REBUILD_OK), 0 )
	
	# Data where MR seemingly worked according to refmac rfree
	data$REFMAC_OK <- as.numeric( data$final_Rfree <= 0.45 )
	data$REFMAC_OK <- replace( data$REFMAC_OK, is.na(data$REFMAC_OK), 0 )
	
	# Data where MR seemingly worked according to phaser
	data$PHASER_OK <- as.numeric( data$PHASER_TFZ > 8.0 | data$PHASER_LLG > 120 )
	data$PHASER_OK <- replace( data$PHASER_OK, is.na(data$PHASER_OK), 0 )
	
	# Gold standard
	data$success <- as.numeric( data$SHELXE_OK == 1 & data$REBUILD_OK == 1 )
	data$success <- replace( data$success, is.na(data$success), 0 )
	
	return(data)
}

data <- updateData(data)


# For each case need to get the best success - i.e. success with max CC
# Order by native_pdb_code and CC
summaryData <- data[ order( data$native_pdb_code, data$success, data$SHELXE_CC, decreasing=TRUE ), ]

# Select top by selecting not duplicates on native_pdb_code
summaryData <- summaryData[ !duplicated(summaryData$native_pdb_code),
		c("native_pdb_code","native_pdb_resolution","native_pdb_num_residues","ensemble_name","truncation_num_residues","subcluster_num_models",
				"PHASER_LLG","PHASER_TFZ", "PHASER_RFZ", "SHELXE_CC", "SHELXE_ACL","SXRBUCC_final_Rfree","SXRARP_final_Rfree")  ]

# Now put in alphabetical order
summaryData <- summaryData[ order( summaryData$native_pdb_code ), ]

# Need to get numbers of success
# This gets the stats  - we can join because the by function is the native_pdb_code which is in similar alphabetic order
summaryData["numModels"] <- aggregate( data$native_pdb_code, by=list(data$native_pdb_code), FUN=length )[2]
summaryData["worked"] <- aggregate( data$success, by=list(data$native_pdb_code), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData$worked <- replace( summaryData$worked, summaryData$worked > 0, 1 )


# Different measures of success
summaryData["success"] <- aggregate( data$success, by=list(data$native_pdb_code), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["SHELXE_OK"] <- aggregate( data$SHELXE_OK, by=list(data$native_pdb_code), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["REFMAC_OK"] <- aggregate( data$REFMAC_OK, by=list(data$native_pdb_code), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["PHASER_OK"] <- aggregate( data$PHASER_OK, by=list(data$native_pdb_code), FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["REBUILD_OK"] <- aggregate( data$REBUILD_OK, by=list(data$native_pdb_code), FUN=function(x){ sum( x == 1 ) } )[2]


# side-chain treatment
summaryData["successPolyAla"]  <- aggregate( 
		data[ data$side_chain_treatment == "polya", ]$success,
		by=list( data[ data$side_chain_treatment == "polya", ]$native_pdb_code ),
		FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["successScwrl"]  <- aggregate(
		data[ data$side_chain_treatment == "reliable", ]$success,
		by=list( data[ data$side_chain_treatment == "reliable", ]$native_pdb_code ),
		FUN=function(x){ sum( x == 1 ) } )[2]
summaryData["allatom"]  <- aggregate(
		data[ data$side_chain_treatment == "allatom", ]$success,
		by=list( data[ data$side_chain_treatment == "allatom", ]$native_pdb_code ),
		FUN=function(x){ sum( x == 1 ) } )[2]

# Subclustering radii
native_pdb_codes <- unique( data$native_pdb_code )

x <- aggregate( success~native_pdb_code,
		data=data,
		FUN=function(x){ sum( x == 1 ) },
		subset=data$subcluster_radius_threshold==1 )
missing <- setdiff( native_pdb_codes, x$native_pdb_code)
y<-rep(0,length(missing)) # Make vector of zeros as long as the missing codes
y <- data.frame( missing, y) # Create a data frame and name it
names(y) <- c("native_pdb_code","success")
y<-rbind(x,y)
y <- y[ order(y$native_pdb_code), ]
summaryData$subclustering1A <- y$success

x <- aggregate( success~native_pdb_code,
		data=data,
		FUN=function(x){ sum( x == 1 ) },
		subset=data$subcluster_radius_threshold==2 )
missing <- setdiff( native_pdb_codes, x$native_pdb_code)
y<-rep(0,length(missing)) # Make vector of zeros as long as the missing codes
y <- data.frame( missing, y) # Create a data frame and name it
names(y) <- c("native_pdb_code","success")
y<-rbind(x,y)
y <- y[ order(y$native_pdb_code), ]
summaryData$subclustering2A <- y$success

x <- aggregate( success~native_pdb_code,
		data=data,
		FUN=function(x){ sum( x == 1 ) },
		subset=data$subcluster_radius_threshold==3 )
missing <- setdiff( native_pdb_codes, x$native_pdb_code)
y<-rep(0,length(missing)) # Make vector of zeros as long as the missing codes
y <- data.frame( missing, y) # Create a data frame and name it
names(y) <- c("native_pdb_code","success")
y<-rbind(x,y)
y <- y[ order(y$native_pdb_code), ]
summaryData$subclustering3A <- y$success

# NB LOOK AT REORDER
# Now put in order by success, resolution
#summaryData <- summaryData[ order( -summaryData$worked, summaryData$native_pdb_resolution ), ]
#x <- x[ order( x$success, decreasing=TRUE ), ]
write.csv(summaryData, "summary.csv", row.names=FALSE)
