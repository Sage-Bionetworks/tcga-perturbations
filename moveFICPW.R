#1) Figure out what hasn't been moved
fic <- synapseQuery('select id, name, numSamples, platform, study, tissueType from entity where entity.benefactorId=="syn1450028" and entity.status=="fic"');
proc <- synapseQuery('select id, name, numSamples, platform, study, tissueType from entity where entity.benefactorId=="syn1450028" and entity.status=="processed"');
proc <- proc[grep("GSE", proc$entity.name),]
proc <- proc[!grepl('unsupervised.', proc$entity.name),]
proc <- proc[!duplicated(proc$entity.name),]
toMove <- setdiff(proc$entity.study, fic$entity.study)
#2) Read in all the cancer studies we currently have
allCancers <- read.table("~/allCancers.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
#3) Get the folders for the previously curated breast and prostate data
folders <- synapseQuery('select id, name from folder where folder.parentId=="syn1353107"')

qry <- paste('select id, name, platform, numSamples from entity where entity.parentId=="syn1353107"')
qryRes <- bigQuery(qry,limit=250)
qryRes <- qryRes[grep('^GSE\\d+_processed', qryRes$entity.name),]
qryRes$study <- gsub("_processed", "", qryRes$entity.name)
	
ids <- match(toMove, qryRes$study)
ids <- ids[!is.na(ids)]
easy <- qryRes[ids,]
easy <- easy[which(!grepl(";",easy$entity.platform) & easy$entity.platform != ""),]


for(i in 32:nrow(easy)){
	if(easy$entity.platform[j] == 'hg-cgh-244a' | easy$entity.platform[j] == 'pd.genomewidesnp.6'){
		cat("Skipping...",i,"\n")
		next
	}else{
		cat("Moving...", i,"\n")
		ent <- loadEntity(easy$entity.id[i])
		studyQry <- synapseQuery(paste('select id, name from folder where folder.parentId=="syn1491485" and folder.name=="',gsub("_processed", "", propertyValue(ent,'name')),'"',sep=""))
		if(is.null(studyQry)){
			cat("ERROR: No study for", x$entity.id[j],"\n")
			next;
		}
		ficEnt <- .addFIC(ent, studyQry$folder.id)
		pwEnt <- .addProbeWeights(ent, studyQry$folder.id)
		cat(propertyValue(ficEnt,'id'),"\n")
	}
}

#
#ids <- match(toMove, qryRes$study)
#hard <- toMove[which(is.na(ids))]
#ids <- match(hard, folders$folder.name)
#ids <- ids[!is.na(ids)]
#hard <- folders[ids,]
#hits <- lapply(hard$folder.id, function(x){ synapseQuery(paste('select id, name, platform from entity where entity.parentId=="',x,'"',sep=""))})
#
#for(i in 10){
#	x <- hits[[i]]
#	for(j in 1:nrow(x)){
#		if(x$entity.platform[j] == 'hg-cgh-244a' | x$entity.platform[j] == 'pd.genomewidesnp.6'){
#			cat("Skipping...",j,"\n")
#			next
#		}else{
#			cat("Moving...", j,"\n")
#			ent <- loadEntity(x$entity.id[j])
#			studyQry <- synapseQuery(paste('select id, name from folder where folder.parentId=="syn1491485" and folder.name=="',annotValue(ent,'study'),'"',sep=""))
#			if(is.null(studyQry)){
#				cat("ERROR: No study for", x$entity.id[j],"\n")
#				next;
#			}
#			ficEnt <- .addFIC(ent, studyQry$folder.id)
#			pwEnt <- .addProbeWeights(ent, studyQry$folder.id)
#			cat(propertyValue(ficEnt,'id'),"\n")
#		}
#	}
#}

.addFIC <- function(ent, parentId){
	# FIC
	if(!('eset' %in% names(ent$objects))){
		fic <- ent$objects[[1]]$statistics$singular.values
		numSamples = ncol(exprs(ent$objects[[1]]$eset))
	}else{
		fic <- ent$objects$statistics$singular.values
		numSamples <- ncol(exprs(ent$objects$eset))
	}
	ficName <- gsub("_processed", "", propertyValue(ent,'name'))
	ficName <- paste(ficName, '_',propertyValue(ent, 'platform'),'_FIC.tsv', sep="")	
	fic <- round(fic,3)
	write.table(fic, file=ficName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=ficName, parentId=parentId))
	propertyValue(newEnt, 'platform') <- propertyValue(ent, 'platform')
	propertyValue(newEnt, 'species') <- 'Homo sapiens'
	propertyValue(newEnt, 'numSamples') <- numSamples
	newEnt <- addFile(newEnt,ficName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(ficName)$size)
	annotValue(newEnt, 'status') <- 'fic'
	annotValue(newEnt, 'summaryType') <- 'gene'
	annotValue(newEnt, 'study') <- gsub("_processed", "", propertyValue(ent,'name'))
	annotValue(newEnt, 'molecFeatureType') <- 'RNA'
	annotValue(newEnt, 'dataType') <- 'mRNA'
	ficEnt <- storeEntity(newEnt)
}

.addProbeWeights <- function(ent, parentId){
	if(!('eset' %in% names(ent$objects))){
		probeWeights <- ent$objects[[1]]$statistics$probe.weights
		numSamples = ncol(exprs(ent$objects[[1]]$eset))
	}else{
		probeWeights <- ent$objects$statistics$probe.weights
		numSamples <- ncol(exprs(ent$objects$eset))
	}	
# Probe weights
	probeWeights <- data.frame(symbols = names(probeWeights),
			probeWeights = round(probeWeights,3))
	wtsName <- gsub("_processed", "", propertyValue(ent,'name'))
	wtsName <- paste(wtsName, '_',propertyValue(ent, 'platform'),'_probeWeights.tsv', sep="")	
	write.table(probeWeights, file=wtsName, sep="\t", quote=FALSE)
	newEnt <- Data(list(name=wtsName, parentId=parentId))
	propertyValue(newEnt, 'platform') <- propertyValue(ent, 'platform')
	propertyValue(newEnt, 'species') <- 'Homo sapiens'
	propertyValue(newEnt, 'numSamples') <- numSamples
	newEnt <- addFile(newEnt,wtsName)
	annotValue(newEnt, 'fileSize') <- .prettifyFileSize(file.info(wtsName)$size)
	annotValue(newEnt, 'status') <- 'probeWeights'
	annotValue(newEnt, 'study') <- gsub("_processed", "", propertyValue(ent,'name'))
	annotValue(newEnt, 'molecFeatureType') <- 'RNA'
	annotValue(newEnt, 'dataType') <- 'mRNA'
	annotValue(newEnt, 'summaryType') <- 'gene'
	wtsEnt <- storeEntity(newEnt)
}




bigQuery <- function(qry, limit=100) 
{
	allQryRes <- matrix()
	offset <- 1
	for(i in 1:100000){
		cat("\rQuery",i)
		qryRes <- synapseQuery(paste(qry, 'limit', limit, 'offset', offset))
		offset <- offset + limit 
		if(i == 1){
			allQryRes <- qryRes
		}else{
			allQryRes <- merge(allQryRes, qryRes, all=TRUE)
		}
		if(nrow(qryRes) != limit){
			break;
		}
	}
	allQryRes
}
.writeEntity <- function(content, newName){
	write.table(content, sep="\t",quote=FALSE, file=newName)
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