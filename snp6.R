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

qry <- synapseQuery('select id, name, numSamples, study,platform from entity where entity.benefactorId == "syn1450028" and entity.status=="processed"')
qry <- qry[grepl("GSE",qry$entity.name),]
qry$study <- sapply(strsplit(qry$entity.name,'_'), function(x){ x[[1]]})
allCancers <- read.table("~/allCancers.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)

toProcess <- setdiff(allCancers$name, qry$study)
toProcess <- allCancers[match(toProcess, allCancers$name),]

qry <- .getRawEntity('GSE32462')
rawEntity <- loadEntity(qry$entity.id[1])
rawDataDir <- file.path(ent$cacheDir)

res <- runWorkflow(rawDataDir, workflow="snm")

for(i in 3:length(ents)){
	newName <- gsub('_processed', '', propertyValue(ents[[i]],'name'))
	newName <- paste(newName, '_pd.genomewidesnp.6_probe.tsv', sep="")
	parentId <- toProcess[i,4]
	ent <- ents[[i]]
	dat <- try(exprs(ent$objects$eset),silent=TRUE)
	annot <- try(annotation(ent$objects$eset),silent=TRUE)
	if(class(dat) == "try-error"){
		dat <- try(exprs(ent$objects[[1]]$eset),silent=TRUE)
		annot <- try(annotation(ent$objects[[1]]$eset),silent=TRUE)
	}
	if(class(dat) == "try-error" | class(annot) == "try-error"){
		cat("ERROR!!!")
		return("Error")
	}
	rownames(dat) <- gsub("_mt", "_eg", rownames(dat))
	dat <- round(dat,3)
	write.table(dat,
			file=newName,
			sep="\t",
			quote=FALSE)
	newEnt <- Data(list(name=newName, parentId=parentId))
	propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
	propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
	propertyValue(newEnt, 'platform') <- 'pd.genomewidesnp.6'
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,newName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(newName)$size)
	annotValue(newEnt, 'status') <- 'processed'
	annotValue(newEnt, 'summaryType') <- 'probe'
	annotValue(newEnt, 'molecFeatureType') <- 'DNA'
	annotValue(newEnt, 'dataType') <- 'CNV'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	probeEnt <- storeEntity(newEnt)
	
	fits <- fit.pset(dat, all.gene2row)
	sapply(fits$probe.weights, function(x){
				nms <- names(x)
				if(is.null(nms)){
					rep("Unkown", length(x))
				}else{
					nms
				}
			}) -> probeNames
	
	# Probe weights
	probeWeights <- data.frame(symbols = rep(names(fits$probe.weights), 
					sapply(fits$probe.weights,length)),
			probeNames=unlist(probeNames),
			probeWeights = round(unlist(fits$probe.weights),3))
	rownames(probeWeights) <- NULL
	wtsName <- gsub('_processed', '', propertyValue(ents[[i]],'name'))
	wtsName <- paste(wtsName, '_pd.genomewidesnp.6_probeWeights.tsv', sep="")	
	write.table(probeWeights, file=wtsName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=wtsName, parentId=parentId))
	propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
	propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
	propertyValue(newEnt, 'platform') <- 'pd.genomewidesnp.6'
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,wtsName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(wtsName)$size)
	annotValue(newEnt, 'status') <- 'probeWeights'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	annotValue(newEnt, 'molecFeatureType') <- 'DNA'
	annotValue(newEnt, 'dataType') <- 'CNV'
	annotValue(newEnt, 'summaryType') <- 'gene'
	wtsEnt <- storeEntity(newEnt)
	
	# FIC
	fic <- fits$singular.values
	ficName <- gsub('_processed', '', propertyValue(ents[[i]],'name'))
	ficName <- paste(ficName, '_pd.genomewidesnp.6_FIC.tsv', sep="")	
	fic <- round(fic,3)
	write.table(fic, file=ficName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=ficName, parentId=parentId))
	propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
	propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
	propertyValue(newEnt, 'platform') <- 'pd.genomewidesnp.6'
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,ficName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(ficName)$size)
	annotValue(newEnt, 'status') <- 'fic'
	annotValue(newEnt, 'summaryType') <- 'gene'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	annotValue(newEnt, 'molecFeatureType') <- 'DNA'
	annotValue(newEnt, 'dataType') <- 'CNV'
	ficEnt <- storeEntity(newEnt)
	
	# Summarized Data
	dat <- fits$estimated.rna.concentration
	datName <- gsub('_processed', '', propertyValue(ents[[i]],'name'))
	datName <- paste(datName, '_pd.genomewidesnp.6_gene.tsv', sep="")	
	dat <- round(dat,3)
	write.table(dat, file=datName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=datName, parentId=parentId))
	propertyValue(newEnt, 'disease') <- propertyValue(ent, 'disease')
	propertyValue(newEnt, 'tissueType') <- propertyValue(ent, 'tissueType')
	propertyValue(newEnt, 'platform') <- 'pd.genomewidesnp.6'
	propertyValue(newEnt, 'species') <- propertyValue(ent, 'species')
	propertyValue(newEnt, 'numSamples') <- ncol(dat)
	annotations(newEnt) <- annotations(ent)
	newEnt <- addFile(newEnt,datName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(datName)$size)
	annotValue(newEnt, 'status') <- 'processed'
	annotValue(newEnt, 'summaryType') <- 'gene'
	annotValue(newEnt, 'study') <- toProcess[i,3]
	annotValue(newEnt, 'molecFeatureType') <- 'DNA'
	annotValue(newEnt, 'dataType') <- 'CNV'
	dataEnt <- storeEntity(newEnt)
	
}
