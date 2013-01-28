allCancers <- read.table("~/allCancers.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
metadata <- bigQuery('select id, name, numSamples, study, platform, parentId from entity where entity.benefactorId == "syn1450028" and entity.status=="metadata"')
raw <- bigQuery('select id, name, numSamples, study, platform, parentId from entity where entity.benefactorId == "syn1450028" and entity.status=="raw"')
processed <- bigQuery('select id, name, numSamples, study, platform, parentId from entity where entity.benefactorId == "syn1450028" and entity.status=="processed"')
fic <- bigQuery('select id, name, numSamples, study, platform, parentId from entity where entity.benefactorId == "syn1450028" and entity.status=="fic"')
probeWeights <- bigQuery('select id, name, numSamples, study, platform, parentId from entity where entity.benefactorId == "syn1450028" and entity.status=="probeWeights"')
studies <- bigQuery('select id, name from folder where folder.benefactorId == "syn1450028"', limit=500)


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

# For each row in allCancers
# 1. Determine if we have a study for it
# 2. Add appropriate annotations to it

id <- match(allCancers$study, studies$folder.name)
id <- id[!is.na(id)]
validCancers <- merge(allCancers, studies, by.x="study", by.y="folder.name")
validCancers$metadata.id <- validCancers$probeWeights.id <- validCancers$raw.id <- validCancers$data.id <- validCancers$fic.id <- rep(NA, nrow(validCancers))

for(i in 11:nrow(validCancers)){
	cat("\r", i, "of", nrow(validCancers))
	ent <- getEntity(validCancers$folder.id[i])
	annotValue(ent, 'tissueType') <- validCancers$tissueType[i]
	annotValue(ent, 'disease') <- validCancers$disease[i]
	annotValue(ent, 'cancerType') <- validCancers$cancerType[i]
	annotValue(ent, 'cancerSubType') <- validCancers$cancerSubType[i]
	for(j in 1:17){
		ent <- .addAnnot(ent, names(validCancers)[j], validCancers[i,j])		
	}
	ent <- updateEntity(ent)	
	validCancers[i,] <- .getIds(validCancers,i)
#	annots <- .clean(validCancers,i)
	# get folder	
}

.addAnnot <- function(ent, name, value){
	if(any (value == "",  is.na(value), value == "NULL", is.null(value))){
		cat("Nothing for:", name, "\n")
		return(ent)
	}
	annotValue(ent, name) <- value
	return(ent)
}

.clean <- function(validCancers,i){
	validCancers[i,which(validCancers[i,] != "" & validCancers[i,] != "NA")]
}

.getIds <- function(validCancers,i){
	platform <- validCancers$platform[i]
	# Get the raw data
	res <- synapseQuery(paste('select id, name, platform from entity where entity.benefactorId == "syn1450028" and entity.study == "', validCancers$study[i], '" and entity.status=="raw"',sep=""))
	if(!is.null(res)){
		# Store the rawId
		res <- res[!duplicated(res$entity.id),]
		# Can potentially be many raw entities.  Use array express as the example.
		if(nrow(res) > 0){
			validCancers$raw.id[i] <- paste(res$entity.id,collapse=",")
		}else{
			validCancers$raw.id[i] <- res$entity.id
		}
	}
	# Query for the processed data summarized to the gene level.  If available then store it
	res <- synapseQuery(paste('select id, name from entity where entity.benefactorId == "syn1450028" and entity.platform=="', platform,'" and entity.status=="processed" and entity.summaryType == "gene" and entity.study=="', validCancers$study[i],'"',sep=""))
	if(!is.null(res)){
		res <- res[!duplicated(res$entity.id),]
		# Store the processed Id
		validCancers$data.id[i] <- res$entity.id
		ent <- .annotate(res$entity.id, validCancers, i)
	}
	# Get the probe weights
	res <- synapseQuery(paste('select id, name from entity where entity.benefactorId == "syn1450028" and entity.platform=="', platform,'" and entity.status=="probeWeights" and entity.study=="', validCancers$study[i],'"',sep=""))
	if(!is.null(res)){
		res <- res[!duplicated(res$entity.id),]
		# Store the probe weights Id
		validCancers$probeWeights.id[i] <- res$entity.id
		ent <- .annotate(res$entity.id, validCancers, i)
	}
	# Get the fic
	res <- synapseQuery(paste('select id, name from entity where entity.benefactorId == "syn1450028" and entity.platform=="', platform,'" and entity.status=="fic" and entity.study=="', validCancers$study[i],'"',sep=""))
	if(!is.null(res)){
		res <- res[!duplicated(res$entity.id),]
		# Store the fic Id
		validCancers$fic.id[i] <- res$entity.id
		ent <- .annotate(res$entity.id, validCancers, i)
	}
	# Get the metadata
	res <- synapseQuery(paste('select id, name from entity where entity.benefactorId == "syn1450028" and entity.status=="metadata" and entity.study=="', validCancers$study[i],'"',sep=""))
	if(!is.null(res)){
		res <- res[!duplicated(res$entity.id),]
		# Store the metadata id
		validCancers$metadata.id[i] <- res$entity.id
		ent <- .annotate(res$entity.id, validCancers, i)
	}
	validCancers[i,]
}

.annotate <- function(entId, validCancers, i){
	ent <- getEntity(entId)
	propertyValue(ent, 'tissueType') <- validCancers$tissueType[i]
	propertyValue(ent, 'disease') <- validCancers$disease[i]
	annotValue(ent, 'cancerType') <- validCancers$cancerType[i]
	annotValue(ent, 'cancerSubType') <- validCancers$cancerSubType[i]
	ent <- updateEntity(ent)
}

getChildren <- function(x){ 
	qry <- synapseQuery(paste('select id, name from entity where entity.parentId=="',x,'"',sep=""))
	qry <- qry[!duplicated(qry$entity.id),]
	qry
}

qry <- synapseQuery('select id,name from entity where entity.parentId=="syn1488308"')

 
