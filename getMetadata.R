.loadEntity <- function(id){ 
	ent <- loadEntity(id); 
	dat <- try(read.table(file.path(ent$cacheDir, ent$files),
					sep="\t",
					stringsAsFactors=FALSE,
					quote="",
					header=TRUE,
					strip.white=TRUE), silent=TRUE)
	if(class(dat) == "try-error"){
		dat <- read.xlsx(file.path(ent$cacheDir, ent$files), sheetIndex=1)
	}
	return(dat)
}

	require(xlsx)
	res <- synapseQuery(paste('select id, name, tissueType from entity where entity.benefactorId == "syn1450028" and entity.status=="metadata"',sep=""))
	res <- res[!duplicated(res$entity.id),]
	allRes <- res
	res <- res[which(tolower(res$entity.tissueType)=="liver"),]
	for(i in 1:nrow(res)){
		cat("\n\nsample: ",i)
		md <- .loadEntity(res$entity.id[i])
		if(i == 1){ 			
			write.xlsx(md,"liver.xlsx", sheetName=res$entity.id[i], append=FALSE, row.names=FALSE)
		}else{
			write.xlsx(md,"liver.xlsx", sheetName=res$entity.id[i], append=TRUE, row.names=FALSE)			
		}
	}
	
	
	qryRes <- synapseQuery(paste('select id, name from entity where entity.benefactorId == "syn1450028" and entity.id=="', rownames(res)[i],'_metadata"',sep=""))	
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




