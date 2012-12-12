# Load the annotation
library(org.Hs.eg.db)
eg2sym <- as.character(org.Hs.egSYMBOL)
names(eg2sym) <- paste(names(eg2sym), '_eg', sep="")

# Build an AR signature
rps <- names(eg2sym)[grep("^RPS\\d+", eg2sym, perl=TRUE)]; rps <- rps[!grepl("P\\d+", eg2sym[rps], perl=TRUE)]
rpl <- names(eg2sym)[grep("^RPL\\d+", eg2sym, perl=TRUE)]; rpl <- rpl[!grepl("P\\d+", eg2sym[rpl], perl=TRUE)]
arEffect <- read.table("~/arEffect",sep="\t",header=TRUE,stringsAsFactors=FALSE)
rnms <- rownames(arEffect[1:300,])
arSig <- names(eg2sym)[match(rnms, eg2sym)]

# Build Myc signature
ent <- loadEntity('syn319533')
dat <- exprs(ent$objects$eset)
mycSig <- calcSig(dat, model.matrix(~factor(rep(1:2,each=9))))
mycSigGenes <- names(eg2sym[gsub("_mt", "_eg", names(which(mycSig$cfs[,2] < -0.75)))])

# TCell Signature
load("~/jak2/tcellSignature.Rda")
tcell <- names(eg2sym)[match(tCellGenes, eg2sym)];


# calculate signature enrichment for a single signature
calcSignatureEnrichment <- function (gic, sig) 
{
	sig <- intersect(sig, rownames(gic))
	other <- setdiff(rownames(gic), sig)
	ks.results <- apply(gic, 2, function(x) {
				res <- ks.test(x[other], x[sig], alternative = "greater")
				c(res$p.value, res$statistic)
			})
}

# compare signature enrichment between two signatures
compareSignatureEnrichment <- function (gic, sig1, sigBreast) 
{
	sig1 <- intersect(sig1, rownames(gic))
	sigBreast <- intersect(sigBreast, rownames(gic))
	ks.results <- apply(gic, 2, function(x) {
				res <- ks.test(x[sig1], x[sigBreast], alternative = "greater")
				c(res$p.value, res$statistic)
			})
}

.loadEntity <- function(id){ 
	ent <- loadEntity(id); 
	dat <- as.matrix(read.table(file.path(ent$cacheDir, ent$files),sep="\t",stringsAsFactors=FALSE, 
			strip.white=TRUE,
			row.names=1))
	return(dat)
}

res <- synapseQuery('select id, name, numSamples, platform from entity where entity.benefactorId=="syn1450028" and entity.status=="processed"');
res <- res[grep("GSE", res$entity.name),]
ids <- c(which(res$entity.platform == "hgu133plus2"),which(res$entity.platform == "hgu133a"))

for(i in 1:length(ids)){ 
	cat("\r", i); 
	ents[[i]] <- .loadEntity(res$entity.id[ids[i]])
}

arScores <- list()
rplScores <- list()
rpsScores <- list()
mycScores <- list()
tcellScores <- list()

nms <- rownames(ents[[1]])
for(i in 2:length(ents)){
	nms <- intersect(nms, rownames(ents[[i]]))
}

for(i in 1:length(ents)){
	cat("\r",i)
	dat2 <- as.matrix(ents[[i]])[nms,]
	mode(dat2) <- "numeric"
	rownames(dat2) <- gsub("_mt", "_eg", rownames(dat2))
	#	sig <- calcSignatureEnrichment(dat2, mycSigGenes)
	sig <- calcSignatureEnrichment(dat2, arSig)
	arScores[[i]] <- sig[2,]
	sig <- calcSignatureEnrichment(dat2, rpl)
	rplScores[[i]] <- sig[2,]
	sig <- calcSignatureEnrichment(dat2, rps)
	rpsScores[[i]] <- sig[2,]
	sig <- calcSignatureEnrichment(dat2, mycSigGenes)
	mycScores[[i]] <- sig[2,]
	sig <- calcSignatureEnrichment(dat2, tcell)
	tcellScores[[i]] <- sig[2,]
}

scores <- data.frame(ar = unlist(arScores),
		rpl=unlist(rplScores),
		rps=unlist(rpsScores),
		myc=unlist(mycScores),
		tcell=unlist(tcellScores))


res$study <- sapply(strsplit(res$entity.name , '_'),function(x){ x[1]})
tissues <- sapply(res$entity.id, function(x){ ent <- getEntity(x); propertyValue(ent, 'tissueType')})
ids <- which(tissues=="")
tissues[ids]<- tmp$tissue[ids]

for(i in 1:length(ids)){
	ent <- getEntity(res$entity.id[ids[i]]); 
	propertyValue(ent, 'tissueType') <- tissues[ids[i]]
	ent <- updateEntity(ent)
}


allCancers <- read.table("~/allCancers.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)
tmp <- allCancers[match(res$study, allCancers$name),]
tmp <- tmp[res$study,]
tmp$tissue[is.na(tmp$tisse)] <- "unknown"

d <- list()
for(i in 1:length(a)){
	cat("\r",i)
	dat2 <- a[[i]]
	rownames(dat2) <- gsub("_mt", "_eg", rownames(dat2))
	#	sig <- calcSignatureEnrichment(dat2, mycSigGenes)
	sig <- compareSignatureEnrichment(dat2, rpl, rps)
	b[[i]] <- sig[2,]
}



png(file="arTissueScores.png", width=1440, height=500)
xyplot(scores$ar ~ 1:length(scores$ar), groups=rep(tolower(tissue), sapply(b,length)), 
		auto.key=list(columns=5), 
		xlab="Samples", 
		ylab="Signature Score",
		main="AR-LNCAP Signature Across 20k Tissue Samples" )
dev.off()

png(file="arTissueScoresByTissue.png")
bwplot(rep(tolower(tissue),sapply(b,length)) ~ scores$ar, main='AR Signature Scores', xlab="Signature Scores")
dev.off()

xyplot(unlist(b) ~ 1:length(unlist(b)), groups=rep(tissue, sapply(b,length)), auto.key=list(columns=5))



xyplot(unlist(b) ~ 1:length(unlist(b)), 
		groups=rep(1:length(b), sapply(b, length)))


sapply(synapseIds, function(x){ 
		qry <- paste('select id, name, platform, tissueType, disease from entity where entity.id=="', x,'"')
		synapseQuery(qry)
	}) -> qry

ent <- loadEntity('syn339037')
gse10072_metadata <- read.table(file.path(ent$cacheDir, ent$files), sep="\t", header=TRUE)
ent <- loadEntity('syn365472')
gse10072_data <- exprs(ent$objects$eset)
lab <- ifelse(grepl("Normal", gse10072_metadata$X.Sample_title), 'Normal', 'Cancer')
rownames(gse10072_data) <- gsub("_mt", "_eg", rownames(gse10072_data))
gse10072_sig <- calcSignatureEnrichment(gse10072_data, mycSigGenes)
xyplot(gse10072_sig[2,] ~ 1:ncol(gse10072_sig), groups=lab)

ent <- loadEntity('syn341974')
gse19804_metadata <- read.table(file.path(ent$cacheDir, ent$files), sep="\t", header=TRUE)
ent <- loadEntity('syn370973')
gse19804_data <- exprs(ent$objects$eset)
lab <- ifelse(grepl("Normal", gse19804_metadata$X.Sample_title), 'Normal', 'Cancer')
rownames(gse19804_data) <- gsub("_mt", "_eg", rownames(gse19804_data))
gse19804_sig <- calcSignatureEnrichment(gse19804_data, mycSigGenes)
xyplot(gse19804_sig[2,] ~ 1:ncol(gse19804_sig), groups=lab)


ent  <- loadEntity('syn341738')
gse19188_metadata <- read.table(file.path(ent$cacheDir, ent$files), sep="\t", header=TRUE)

labs <- c(rep("Adenocarcinoma",120), 
		ifelse(grepl("Normal", gse19804_metadata$X.Sample_title), 'Adjacent-Normal', 'Adenocarcinoma'),
		ifelse(grepl("healthy",gse19188_metadata$X.Sample_characteristics_ch1),"Adjacent-Normal", "NSCLC"),
#		ifelse(grepl('\\d+N', gse19188_metadata$X.Sample_title), "Normal", "NSCLC"),
		rep("Adenocarcinoma", 17),
		rep("Epithelial Lining Fluid", 30),
		ifelse(grepl("Normal", gse10072_metadata$X.Sample_title), 'Adjacent-Normal', 'Adenocarcinoma'),		
		rep(c("Adenocarcinoma", "Adjacent-Normal"), c(30,20)))

xyplot(scores$myc[which(tissues == "lung")] ~ 1:sum(tissues=="lung"), groups=labs, auto.key=T)


f <- list()
for(i in 1:length(synapseIds)){
	cat(i, "\n")
	ent <- ents[[i]]
	dat2 <- try(ent$objects$statistics$singular.values[,1], silent=TRUE)
	if(class(dat2) == "try-error" | is.null(dat2)){
		dat2 <- try(ent$objects[[1]]$statistics$singular.values[,1], silent=TRUE)
	}
	if(is.null(dat2)){
		f[[i]] <- NA
	}else{
		names(dat2) <- gsub("_mt", "_eg", names(dat2))
		f[[i]] <- dat2
	}
}
f <- f[which(sapply(f, length) != 1)]
th <- intersect(names(f[[1]]), names(f[[2]]))
fic <- sapply(f, function(x){ x[th]})
X <- model.matrix(~ -1 + tmp5$Tissue)
cfs <- t(solve(t(X) %*% X) %*% t(X) %*% t(fic))
top50 <- names(sort(-rowMeans(cfs))[1:200])

d <- list()
for(i in 1:length(a)){
	corsMyc <- cor(t(a[[i]][top50,]))	
	corsPerm <- cor(t(a[[i]][sample(nrow(a[[i]]), length(top50)),]))
	d[[i]] <- list(perm=unique(as.numeric(corsPerm)), 
			myc=unique(as.numeric(corsMyc)))
}

g <- list()
for(i in 1:length(a)){
	corsMyc <- cor(t(a[[i]][top50,]))	
	corsPerm <- cor(t(a[[i]][sample(nrow(a[[i]]), length(top50)),]))
	g[[i]] <- corsMyc
}

h <- g[[1]]
for(i in 2:length(g)){ h <- h + g[[i]] }
diag(h) <- 0

best <- names(which(apply(h,2,max) > 25))
h2 <- g[[1]][best,best]
for(i in 2:length(g)){ h2 <- h2 + g[[i]][best,best] }
