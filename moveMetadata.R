#1) Get all the entity ids for the new folders in the SCR.
allFolders <- synapseQuery('select id, name from folder where folder.parentId=="syn1491485"')
#2) Read in all the cancer studies we currently have
allCancers <- read.table("~/allCancers.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
#3) Get the folders for the previously curated breast and prostate data
folders <- synapseQuery('select id, name from folder where folder.parentId=="syn1353107"')
#4) Get all the existing metadata files in the new SCR.
qry <- synapseQuery('select id, name, numSamples, study,platform from entity where entity.benefactorId == "syn1450028" and entity.status=="metadata"')

toMove <- unique(setdiff(allCancers$name, qry$entity.study))

sapply(toMove, function(x){ 
		nm <- paste(x, '_metadata', sep="")
		synapseQuery(paste('select id, name from entity where entity.name=="', nm,'"',sep=""))
	}) -> existing

ids <- which(sapply(existing, function(x){ grepl("GSE\\d+\\_metadata", x$entity.name)}) == "TRUE")
toMove <- toMove[ids]
sapply(existing[ids], function(x){ 
			ent <- loadEntity(x$entity.id)
			file.copy(file.path(ent$cacheDir, ent$files), paste('newMeta/',ent$files,sep=""))
			ent$files
		}) -> files

# For each study, we need to target folder, the existing id,  
ids2 <- which(sapply(files,length) == 1)


sapply(existing[ids[ids2]][c(-1,-2,-3)], function(x){ 
			cat(unlist(x),"\n")
			inputId <- x$entity.id
			study <- gsub("_metadata", "", x[1])
			qry <- synapseQuery(paste('select id, name from folder where folder.benefactorId == "syn1450028" and folder.name=="', study,'"', sep=""))
			parentId <- qry$folder.id
			newMetadata <- .moveMetadata(inputId, study, parentId)
			#file.link(file.path(ent$cacheDir, ent$files), paste('newMeta/',ent$files,sep=""))
		}) -> newEnts

files <- list.files()
for(i in 33:length(files)){
	study <- gsub("_metadata.tsv", "", files[i])
	qry <- synapseQuery(paste('select id, name from folder where folder.benefactorId == "syn1450028" and folder.name=="', study,'"', sep=""))
	parentId <- qry$folder.id
	newEnt <- .moveLocalMetadata(file[i], study, parentId)
}

.moveLocalMetadata <- function(file, study, parentId){
	newName <- paste(study,'_metadata.tsv', sep="")
	parent <- getEntity(parentId)
	ent <- Data(list(name=newName,parentId=parentId))
	annotValue(ent, 'fileSize') <- .prettifyFileSize(file.info(newName)$size)
	ent <- addFile(ent, newName)
	annotValue(ent, 'study') <- study
	annotValue(ent, 'repository') <- annotValue(parent,'repository')
	annotValue(ent, 'status') <- 'metadata'
	ent <- storeEntity(ent)
}


.moveMetadata <- function(inputId, study, parentId){
	parent <- getEntity(parentId)
	ent <- loadEntity(inputId)
	content <- read.table(file.path(ent$cacheDir, ent$files),sep="\t",stringsAsFactors=FALSE,header=TRUE)
	newName <- paste(study,'_metadata.tsv', sep="")
	eme <- .writeEntity(content, newName)
	ent <- Data(list(name=newName,parentId=parentId))
	annotValue(ent, 'fileSize') <- .prettifyFileSize(file.info(newName)$size)
	ent <- addFile(ent, newName)
	annotValue(ent, 'study') <- study
	annotValue(ent, 'repository') <- annotValue(parent,'repository')
	annotValue(ent, 'status') <- 'metadata'
	ent <- storeEntity(ent)
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
