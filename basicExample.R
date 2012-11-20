# Log in and load the necessary entities
source("~/login.R")

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

# Load TCGA Breast Cancer Data
ent2 <- loadEntity('syn332584')
dat2 <- exprs(ent2$objects$eset)

# Load TCGA Agilent Breast Cancer Data
ent4 <- loadEntity('syn313753')
dat4 <- exprs(ent4$objects$eset)

# Load TCGA Lung Adenocarcinoma Cancer Data
ent <- loadEntity('syn886797')
dat <- exprs(ent$objects$eset)

# Load the TCGA ovarian hgu133a data
ent3 <- loadEntity('syn313584')
dat3 <- exprs(ent3$objects$eset)

# Load the annotation
library(org.Hs.eg.db)
eg2sym <- as.character(org.Hs.egSYMBOL)
names(eg2sym) <- paste(names(eg2sym), '_eg', sep="")

# Build Myc signature
ent <- loadEntity('syn319533')
dat <- exprs(ent$objects$eset)
mycSig <- calcSig(dat, model.matrix(~factor(rep(1:2,each=9))))
mycSigGenes <- names(eg2sym[gsub("_mt", "_eg", names(which(mycSig$cfs[,2] < -0.75)))])

# Calculate enrichment
sigBreast <- calcSignatureEnrichment(dat2, mycSigGenes)
sigBreast2 <- calcSignatureEnrichment(dat4, mycSigGenes)
sigLung <- calcSignatureEnrichment(dat, mycSigGenes)
sigOvarian <- calcSignatureEnrichment(dat3, mycSigGenes)

# Plot the scores acorss patients
xyplot(sigBreast[2,] ~ 1:ncol(sigBreast), xlab="Samples", ylab="KS Test Statistics",main="Breast")
xyplot(sigLung[2,] ~ 1:ncol(sigLung), xlab="Samples", ylab="KS Test Statistics",main="Lung")
xyplot(sigOvarian[2,] ~ 1:ncol(sigOvarian), xlab="Samples", ylab="KS Test Statistics",main="Ovarian")
xyplot(c(sigLung[2,],sigBreast[2,],sigOvarian[2,],sigBreast2[2,]) ~ 1:(ncol(sigLung)+ncol(sigBreast)+ncol(sigBreast2)+ncol(sigOvarian)),
		groups=rep(c("Lung", "Breast", "Ovarian", "Breast-2"),c(ncol(sigLung), ncol(sigBreast), ncol(sigOvarian), ncol(sigBreast2))),
		auto.key=list(columns=3),
		xlab="Samples", ylab="KS Test Statistics",main="Lung, Breast, and Ovarian")

# Build labels for each gene
lab <- rep("Remainder", nrow(dat2))
names(lab) <- rownames(dat2)
lab[mycSigGenes] <- "Myc Signature"

# Densityplot showing the most significant enrichment
densityplot(log2(dat2[,which.max(sigBreast[2,])]), groups=lab, plot.points=FALSE)

# Densityplot showing the least significant enrichment
densityplot(log2(dat2[,which.min(sigBreast[2,])]), groups=lab, plot.points=FALSE)



