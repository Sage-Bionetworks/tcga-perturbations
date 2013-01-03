

mgQuery <- function(tissue){
	require(xlsx)

	metadata <- synapseQuery('select id, name, numSamples, study, platform, tissueType from entity where entity.benefactorId == "syn1450028" and entity.status=="metadata"')
	processed <- synapseQuery('select id, name, numSamples, study, platform, tissueType from entity where entity.benefactorId == "syn1450028" and entity.status=="processed"')
	
	folders <- c(ncbi = 'syn1353107',
			tcga = 'syn1393205',
			arrayExpress = 'syn1125842',
			other = 'syn1125846')
	lapply(folders, function(x){ 
				synapseQuery(paste('select * from folder where folder.',tissue,'=="TRUE" and folder.cancer=="TRUE" and folder.parentId=="', x,'"',sep=""))
			}) -> tmp
	tmp <- tmp[which(!sapply(tmp, is.null))]	
	tmp2 <- tmp[[1]]
	for(i in 2:length(tmp)){
		tmp2 <- merge(tmp2, tmp[[i]], all=TRUE)
	}
	dataEntities <- synapseQuery(paste('select * from entity where entity.parentId == "', tmp2$folder.id[1], '" and entity.disease=="cancer" and entity.tissueType=="',tissue,'"',sep=""))
	if(nrow(tmp2) > 1){
		for(i in 2:nrow(tmp2)){
			tmp <- synapseQuery(paste('select * from entity where entity.parentId == "', tmp2$folder.id[i], '" and entity.disease=="cancer" and entity.tissueType=="',tissue,'"',sep=""))
			dataEntities <- merge(dataEntities,tmp, all=TRUE)
		}
	}
	countsByType <- sapply(split(as.numeric(dataEntities$entity.numSamples), dataEntities$entity.dataType), sum)
	tmp <- dataEntities[, c('entity.study', 'entity.numSamples', 'entity.dataType')]
	tmp <- tmp[which(rowSums(is.na(tmp)) == 0),]
	qryMat <- matrix(NA, nr=length(unique(tmp[,1])), nc=length(unique(tmp[,3])))
	rownames(qryMat) <- unique(tmp[,1])
	colnames(qryMat) <- unique(tmp[,3])
	for(i in 1:nrow(tmp)){
		qryMat[tmp[i,1], tmp[i,3]] <- tmp[i,2]
	}
	mode(qryMat) <- "numeric"
	
	metadataIds <- 1:(nrow(qryMat)-1)
	for(i in 2:nrow(qryMat)){
		qryRes <- synapseQuery(paste('select id, name from entity where entity.createdByPrincipalId == "273975" and entity.name=="', rownames(qryMat)[i],'_metadata"',sep=""))	
		metadataIds[i-1] <- qryRes$entity.id
	}
	ents <- sapply(metadataIds, function(x) { loadEntity(x)})

	# Log in 
	auth = getGoogleAuth("brig.mecham@sagebase.org", "STL82sage")
	con = getGoogleDocsConnection(auth)
	
	# Create a new spreadsheet
	sh <- addSpreadsheet(con, dim = c(20, 10) , name = names(ents)[1])
	
	# Log in again
	auth = getGoogleAuth("brig.mecham@sagebase.org", "STL82sage", "wise")
	con = getGoogleDocsConnection(auth)

	# Read in the metadata
	md <- read.table(file=file.path(ents[[1]]$cacheDir,ents[[1]]$files), sep="\t", header=FALSE); 
	
	#Create a new worksheet
	ws <- addWorksheet(names(ents)[1], con, dim = c(nrow(md), ncol(md)), title = names(ents)[1], asSheetRef = TRUE)
	
	con2 <- getGoogleDocsConnection('brig.mecham@sagebase.org','STL82sage', "wise")
	sheets = getWorksheets(names(ents)[1], con2)
	sheet1 <- sheets[[2]]
	
	for(i in 1:nrow(md)){
		for(j in 1:ncol(md)){
			a <- RGoogleDocs:::setCell(sheet1, c(i,j), md[i,j], con=sheet1@connection)
		}
	}
	
	for(i in 1:length(ents)){
		cat("\r", i)
		md <- read.table(file=file.path(ents[[i]]$cacheDir,ents[[i]]$files), sep="\t", header=TRUE); 
		if(i == 1){ 			
			write.xlsx(md,"Breast.xlsx", sheetName=names(ents)[i], append=FALSE, row.names=FALSE)
		}else{
			write.xlsx(md,"Breast.xlsx", sheetName=names(ents)[i], append=TRUE, row.names=FALSE)			
		}
	}
	
	list(hits=dataEntities, countsByType=countsByType, studyByType = qryMat)
}




wb <- loadWorkbook("~/Breast_standardized.xlsx")
synapseIds <- names(getSheets(wb))

for(i in 3:length(synapseIds)){
	tmp <- read.xlsx2("~/Breast_standardized.xlsx",sheetIndex=i)
	ent <- loadEntity(synapseIds[i])
	fname <- ent$files[1]
	write.table(tmp, file=fname, sep="\t", quote=FALSE, row.names=FALSE)
	ent <- deleteFile(ent, ent$files)
	ent <- addFile(ent, fname)
	ent <- storeEntity(ent)
}




