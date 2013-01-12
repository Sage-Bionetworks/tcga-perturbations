qry <- synapseQuery('select id, name, numSamples, study, platform from entity where entity.benefactorId == "syn1450028" and entity.status=="processed"')
qry <- qry[grepl("GSE",qry$entity.name),]
qry$study <- sapply(strsplit(qry$entity.name,'_'), function(x){ x[[1]]})
allCancers <- read.table("~/allCancers.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)

toProcess <- setdiff(allCancers$name, qry$study)
toProcess <- allCancers[allCancers$name %in% toProcess,]

toProcess <- toProcess[grep("^h", toProcess$platform),]
toProcess <- toProcess[grepl("^GSE", toProcess$name),]
toProcess <- toProcess[order(toProcess$numSamples),]
toProcess <- toProcess[c(-1,-45),]
toProcess <- toProcess[!duplicated(toProcess$name),]
toProcess <- cbind(toProcess[,c(1,2,10,11,12)],rep(1:9,times=5)[1:nrow(toProcess)])
colnames(toProcess)[5] <- 'bin'
toProcess <- toProcess[order(toProcess[,5]),]

ids <- which(toProcess$bin != 0)
cmds <- paste('Rscript affy3Prime.R', toProcess$name[ids], ';')
sapply(cmds, function(x){
			cat(x, ' & ', "\n")
		})
write.table(cmds, row.names=FALSE, quote=FALSE)

Rscript affy3Prime.R GSE30258 ;
Rscript affy3Prime.R GSE17818 ;
Rscript affy3Prime.R GSE23980 ;
Rscript affy3Prime.R GSE7390 ;
Rscript affy3Prime.R GSE22470 ;
Rscript affy3Prime.R GSE26971 ;
Rscript affy3Prime.R GSE37382 ;
Rscript affy3Prime.R GSE28497 ;
Rscript affy3Prime.R GSE25055 ;
Rscript affy3Prime.R GSE19784 ;
Rscript affy3Prime.R GSE14520 ;
Rscript affy3Prime.R GSE19301 ;
Rscript affy3Prime.R GSE13159 ;

.getRawEntity <- function(id){
	qry <- synapseQuery(paste('select id, name, numSamples, study,platform from entity where entity.benefactorId == "syn1450028" and entity.status=="raw" and entity.study=="',id,'"', sep=""))
}

.loadEntity <- function(id){ 
	ent <- loadEntity(id); 
	dat <- as.matrix(read.table(file.path(ent$cacheDir, ent$files),sep="\t",stringsAsFactors=FALSE, 
					strip.white=TRUE,
					row.names=1))
	return(dat)
}

qry <- .getRawEntity('GSE21036')
rawEntity <- loadEntity(qry$entity.id[1])
rawDataDir <- file.path(rawEntity$cacheDir)

res <- runWorkflow(rawDataDir, workflow="snm")
# where is the data going?
parentId <- synapseQuery(paste('select id, name, numSamples, study,platform from folder where folder.benefactorId == "syn1450028" and folder.name=="', annotValue(rawEntity,'study'),'"',sep=""))$folder.id
res <- res[which(sapply(res,class) == "list")]
numPlatforms <- length(res)
for(i in 1:numPlatforms){
	ficEnt <- .addFIC(res[i], rawEntity, parentId)
	dataEnt <- .addGeneExpression(res[i], rawEntity, parentId); 
	wtsEnt <- .addProbeWeights(res[i], rawEntity, parentId)
}

.addFIC <- function(workflowObj, rawEntity, parentId){
	# FIC
	fic <- workflowObj[[1]]$statistics$singular.values
	ficName <- annotValue(rawEntity,'study')
	ficName <- paste(ficName, '_',names(workflowObj[1]),'_FIC.tsv', sep="")	
	fic <- round(fic,3)
	write.table(fic, file=ficName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=ficName, parentId=parentId))
	propertyValue(newEnt, 'platform') <- names(workflowObj)[1]
	propertyValue(newEnt, 'species') <- 'Homo sapiens'
	propertyValue(newEnt, 'numSamples') <- ncol(exprs(workflowObj[[1]]$eset))
	newEnt <- addFile(newEnt,ficName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(ficName)$size)
	annotValue(newEnt, 'status') <- 'fic'
	annotValue(newEnt, 'summaryType') <- 'gene'
	annotValue(newEnt, 'study') <- annotValue(rawEntity,'study')
	annotValue(newEnt, 'molecFeatureType') <- 'RNA'
	annotValue(newEnt, 'dataType') <- 'mRNA'
	ficEnt <- storeEntity(newEnt)
}

.addProbeWeights <- function(workflowObj, rawEntity, parentId){
	probeWeights <- workflowObj[[1]]$statistics$probe.weights
	
# Probe weights
	probeWeights <- data.frame(symbols = names(probeWeights),
			probeWeights = round(probeWeights,3))
	wtsName <- annotValue(rawEntity,'study')
	wtsName <- paste(wtsName, '_',names(workflowObj[1]),'_probeWeights.tsv', sep="")	
	write.table(probeWeights, file=wtsName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=wtsName, parentId=parentId))
	propertyValue(newEnt, 'platform') <- names(workflowObj)[1]
	propertyValue(newEnt, 'species') <- 'Homo sapiens'
	propertyValue(newEnt, 'numSamples') <- ncol(exprs(workflowObj[[1]]$eset))
	newEnt <- addFile(newEnt,wtsName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(wtsName)$size)
	annotValue(newEnt, 'status') <- 'probeWeights'
	annotValue(newEnt, 'study') <- annotValue(rawEntity,'study')
	annotValue(newEnt, 'molecFeatureType') <- 'RNA'
	annotValue(newEnt, 'dataType') <- 'mRNA'
	annotValue(newEnt, 'summaryType') <- 'gene'
	wtsEnt <- storeEntity(newEnt)
	
}

.addGeneExpression <- function(workflowObj, rawEntity, parentId){
	dat <- exprs(workflowObj[[1]]$eset)
	datName <- annotValue(rawEntity,'study')
	datName <- paste(datName, '_',names(workflowObj[1]),'_gene.tsv', sep="")
	rownames(dat) <- gsub("_mt", "_eg", rownames(dat))
	dat <- round(dat,3)
	write.table(dat, file=datName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=datName, parentId=parentId))
	propertyValue(newEnt, 'platform') <- names(workflowObj)[1]
	propertyValue(newEnt, 'species') <- 'Homo sapiens'
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	newEnt <- addFile(newEnt,datName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(datName)$size)
	annotValue(newEnt, 'status') <- 'processed'
	annotValue(newEnt, 'summaryType') <- 'gene'
	annotValue(newEnt, 'study') <- annotValue(rawEntity,'study')
	annotValue(newEnt, 'molecFeatureType') <- 'RNA'
	annotValue(newEnt, 'dataType') <- 'mRNA'
	dataEnt <- storeEntity(newEnt)	
}

.prettifyFileSize <- function(fileSize){
	if(fileSize > 1073741824){
		# Its bigger than a gigabyte, so translate to something GB
		return(paste(round(fileSize / 1073741824,1), 'GB'))
	}else{
		size <- round(fileSize / 1048576,1)
		if(size < 0.1){ 
			size <- 0.1
		}
		size <- paste(size, 'MB',sep=" ")
		return(size)
	}
}
