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


cls <- c(SHELXE_CC="numeric",
		SHELXE_ACL="numeric",
		SXRARP_final_Rfree="numeric",
		SXRBUCC_final_Rfree="numeric",
		final_Rfree="numeric",
		PHASER_TFZ="numeric",
		PHASER_LLG="numeric"
)

f="/media/data/shared/testset/truncate_single/results_truncate_single.csv"
dts <- read.csv(file=f,sep=',', header=T, stringsAsFactors=FALSE, na.strings="N/A", colClasses=cls)
dts[dts=="N/A"]<-NA
f="/media/data/shared/testset/percent_chadwick/results_percent.csv"
dp <- read.csv(file=f,sep=',', header=T, stringsAsFactors=FALSE, na.strings="N/A", colClasses=cls)
dp[dp=="N/A"]<-NA
f="/media/data/shared/testset/thresh_morrigan/results_thresh.csv"
dt <- read.csv(file=f,sep=',', header=T, stringsAsFactors=FALSE, na.strings="N/A", colClasses=cls)
dt[dt=="N/A"]<-NA


dt <- updateData(dt)
dp <- updateData(dp)
dts <- updateData(dts)

# Get list of ensemble_percent_model values
epr <- c(dp$ensemble_percent_model,dt$ensemble_percent_model, dts$ensemble_percent_model)
epr <- sort(unique(epr))

ids <- sort(unique(dp$native_pdb_code))
column_labels <- c("native_pdb_code","treatment",epr)

d <- data.frame(matrix(NA, nrow = length(ids)*3, ncol = length(column_labels)))
colnames(d) <- column_labels

x=1
treatments <- c("threshold","percent","truncate_single")
for (p in ids) {
    for (i in c(1,2,3)){
        d$native_pdb_code[x] <- p
        d$treatment[x] <- treatments[i]
        x <- x+1
    }
}


for ( p in ids ){
    for ( i in epr ) {
        #d[ d$native_pdb_code == p, match(i,column_labels) ] <- sum( dp$native_pdb_code == p & dp$ensemble_percent_model == i & dp$success == 1, na.rm=TRUE)
        d[ d$native_pdb_code == p & d$treatment == "threshold", match(i,column_labels) ] <- sum( dt$native_pdb_code == p & dt$ensemble_percent_model == i & dt$success == 1, na.rm=TRUE)
        d[ d$native_pdb_code == p & d$treatment == "percent", match(i,column_labels) ] <- sum( dp$native_pdb_code == p & dp$ensemble_percent_model == i & dp$success == 1, na.rm=TRUE)
        d[ d$native_pdb_code == p & d$treatment == "truncate_single", match(i,column_labels) ] <- sum( dts$native_pdb_code == p & dts$ensemble_percent_model == i & dts$success == 1, na.rm=TRUE)
    }
}

write.csv(d, "all_map.csv", row.names=FALSE)