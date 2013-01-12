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

# Build a HIF1 signature
load("~/hif1Genes.Rda")
hif1 <- names(eg2sym)[match(hif1Genes, eg2sym)];

# Build a stem cell signature
stemCells <- read.table("~/stemCells.txt",stringsAsFactors=FALSE)[,1]
stem <- names(eg2sym)[match(stemCells, eg2sym)];

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
compareSignatureEnrichment <- function (gic, sig1, sig2) 
{
	sig1 <- intersect(sig1, rownames(gic))
	sig2 <- intersect(sig2, rownames(gic))
	ks.results <- apply(gic, 2, function(x) {
				res <- ks.test(x[sig1], x[sig2], alternative = "greater")
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

# Build ER signature
breastIds <- res$entity.id[which(res$entity.tissueType == "breast")]
sapply(breastIds[1:10], function(x){ 
			cat(x,"\n")
			dat <- .loadEntity(x)
			u <- fs(dat)
			cor(u$v[,1], dat['2099_eg',])
		}) -> cors
good <- names(which(abs(cors) > 0.7))
sapply(good, function(x){ 
			cat(x,"\n")
			dat <- .loadEntity(x)
			cor(t(dat), dat['2099_eg',])
		}) -> cors2
rownames(cors2) <- rownames(.loadEntity(good[1]))
erCors <- rowMeans(cors2)
erPos <- names(sort(-erCors))[1:250]
erNeg <- names(sort(erCors))[1:250]

# Get the data sets
res <- synapseQuery('select id, name, numSamples, platform, study, tissueType from entity where entity.benefactorId=="syn1450028" and entity.status=="processed"');
res <- res[grep("GSE", res$entity.name),]
res <- res[!grepl('unsupervised.', res$entity.name),]
res <- res[!duplicated(res$entity.name),]

ids <- c(which(res$entity.platform == "hgu133plus2"),
		which(res$entity.platform == "hgu133a"),
		which(res$entity.platform == "hthgu133a"),
		which(res$entity.platform == "hgu133a2"),
		which(res$entity.platform == "hugene10stv1"),
		which(res$entity.platform == "huex10stv2"))

# Reduce to affymetrix platforms
res <- res[ids,]

# Load the data
ents <- list()
for(i in 1:nrow(res)){ 
	cat("\r", i); 
	ents[[i]] <- .loadEntity(res$entity.id[i])
}

# How do we find metadata
# How do we find studies
# How do we find data
# Why can't he access the account - blogsdon@fhcrc.org

arScores <- list()
rplScores <- list()
rpsScores <- list()
mycScores <- list()
tcellScores <- list()
stemScores <- list()
hif1Scores <- list()
erPosScores <- list()
erNegScores <- list()
erNegvPosScores <- list()
erPosvNegScores <- list()

nms <- rownames(ents[[1]])
for(i in 2:length(ents)){
	nms <- intersect(nms, rownames(ents[[i]]))
}

dat2 <- as.matrix(ents[[1]])[nms,]	
allDat <- matrix(NA, nr=length(nms), nc=sum(sapply(ents, dim)[2,]))

# Scoring each sample and signature of interest
for(i in 1:length(ents)){
	cat("\r",i)
	dat2 <- as.matrix(ents[[i]])[nms,]
	mode(dat2) <- "numeric"
	rownames(dat2) <- gsub("_mt", "_eg", rownames(dat2))
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
	sig <- calcSignatureEnrichment(dat2, hif1)
	hif1Scores[[i]] <- sig[2,]
	sig <- calcSignatureEnrichment(dat2, stem)
	stemScores[[i]] <- sig[2,]
	sig <- calcSignatureEnrichment(dat2, erPos)
	erPosScores[[i]] <- sig[2,]
	sig <- calcSignatureEnrichment(dat2, erNeg)
	erNegScores[[i]] <- sig[2,]
	s1 <- compareSignatureEnrichment(dat2, erPos, erNeg)
	erNegvPosScores[[i]] <- s1[2,]
	s2 <- compareSignatureEnrichment(dat2, erNeg, erPos)
	erPosvNegScores[[i]] <- s2[2,]
}

scores <- data.frame(celFileName =  unlist(sapply(arScores, names)),
		ar = unlist(arScores),
		rpl=unlist(rplScores),
		rps=unlist(rpsScores),
		myc=unlist(mycScores),
		tcell=unlist(tcellScores),
		stemCell=unlist(stemScores),
		hif1 = unlist(hif1Scores),
		erPos = unlist(erPosScores),
		erNeg = unlist(erNegScores),
		erPosvNeg = unlist(erPosvNegScores),
		erNegvPos = unlist(erNegvPosScores))
		
allCancers <- read.table("~/allCancers.txt",sep="\t",stringsAsFactors=FALSE,header=TRUE)

ids2 <- which(res$entity.tissueType=="" | is.na(res$entity.tissueType))
res$entity.tissueType[ids2] <- allCancers$tissue[match(res$entity.study[ids2], allCancers$name)]
res$numSamples <- sapply(arScores, length)

for(i in 1:length(ids2)){
	ent <- getEntity(res$entity.id[ids2[i]]); 
	propertyValue(ent, 'tissueType') <- res$entity.tissueType[ids2[i]]
	ent <- updateEntity(ent)
}

meta <- res
scores$study <- rep(meta$entity.study, sapply(arScores, length))
scores$platform <- rep(meta$entity.platform, sapply(arScores, length))
scores$tissue<- tolower(rep(meta$entity.tissueType, sapply(arScores, length)))
save(meta, scores, file="~/me.Rda")

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
bwplot(reorder(tolower(scores$tissue), scores$ar, median) ~ scores$ar, 
		main='AR Signature Scores', 
		xlab="Signature Scores")
dev.off()

bwplot(reorder(tolower(scores$tissue), scores$erPosvNeg, median) ~ scores$erPosvNeg, 
		main='ER Pos vs Neg Signature Scores', 
		xlab="Signature Scores")

bwplot(reorder(tolower(scores$tissue), scores$hif1, median) ~ scores$hif1, 
		main='HIF1 Signature Scores', 
		xlab="Signature Scores")

bwplot(reorder(tolower(scores$tissue), scores$erPos, median) ~ scores$erPos, 
		main="ER+ Signature Scores",
		xlab="Signature Scores")

bwplot(reorder(tolower(scores$tissue), scores$myc, median) ~ scores$myc, 
		main="MYC Signature Scores",
		xlab="Signature Scores")

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
