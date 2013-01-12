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

id <- match(allCancers$name, studies$folder.name)
id <- id[!is.na(id)]
validCancers <- merge(allCancers, studies, by.x="name", by.y="folder.name")
validCancers$metadata <- validCancers$probeWeights.id <- validCancers$raw.id <- validCancers$data.id <- validCancers$fic.id <- rep(NA, nrow(validCancers))

for(i in 1:nrow(validCancers)){
	cat("\r", i, "of", nrow(validCancers))
	platform <- validCancers$platform[i]
	# Get the raw data
	res <- synapseQuery(paste('select id, name from entity where entity.study == "', validCancers$name[i], '" and entity.platform=="', platform,'" and entity.status=="raw"',sep=""))
	if(!is.null(res)){
		# Store the rawId
		res <- res[!duplicated(res$entity.id),]
		# Can potentially be many raw entities.  Use array express as the example.
		validCancers$raw.id[i] <- res$entity.id
	}
	
	# Query for the processed data.  If available then store it
	res <- synapseQuery(paste('select id, name from entity where entity.parentId=="',validCancers$folder.id[i],'" and entity.platform=="', platform,'" and entity.status=="processed"',sep=""))
	if(!is.null(res)){
		res <- res[!duplicated(res$entity.id),]
		# Store the processed Id
		validCancers$data.id[i] <- res$entity.id
	}
	# Get the probe weights
	res <- synapseQuery(paste('select id, name from entity where entity.parentId=="',validCancers$folder.id[i],'" and entity.platform=="', platform,'" and entity.status=="probeWeights"',sep=""))
	if(!is.null(res)){
		res <- res[!duplicated(res$entity.id),]
		# Store the probe weights Id
		validCancers$probeWeights.id[i] <- res$entity.id
	}
	# Get the fic
	res <- synapseQuery(paste('select id, name from entity where entity.parentId=="',validCancers$folder.id[i],'" and entity.platform=="', platform,'" and entity.status=="fic"',sep=""))
	if(!is.null(res)){
		res <- res[!duplicated(res$entity.id),]
		# Store the fic Id
		validCancers$fic.id[i] <- res$entity.id
	}
	# Get the metadata
	res <- synapseQuery(paste('select id, name from entity where entity.parentId=="',validCancers$folder.id[i],'" and entity.status=="metadata"',sep=""))
	if(!is.null(res)){
		res <- res[!duplicated(res$entity.id),]
		# Store the fic Id
		validCancers$metadata[i] <- res$entity.id
	}
	
}

getChildren <- function(x){ 
	qry <- synapseQuery(paste('select id, name from entity where entity.parentId=="',x,'"',sep=""))
	qry <- qry[!duplicated(qry$entity.id),]
	qry
}

qry <- synapseQuery('select id,name from entity where entity.parentId=="syn1488308"')

 
