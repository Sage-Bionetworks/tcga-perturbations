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
compareSignatureEnrichment <- function (gic, sig1, sig2) 
{
	sig1 <- intersect(sig1, rownames(gic))
	sig2 <- intersect(sig2, rownames(gic))
	ks.results <- apply(gic, 2, function(x) {
				res <- ks.test(x[sig1], x[sig2], alternative = "greater")
				c(res$p.value, res$statistic)
			})
}

# Load TCGA Breast Cancer Data
ent2 <- loadEntity('syn332584')
esets2 <- exprs(ent2$objects$eset)

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
sig2 <- calcSignatureEnrichment(log2(1+esets2), mycSigGenes)
xyplot(sig2[2,] ~ 1:ncol(sig2), xlab="Samples", ylab="KS Test Statistics")

# Build labels for each gene
lab <- rep("Remainder", nrow(esets2))
names(lab) <- rownames(esets2)
lab[sig] <- "Myc Signature"

# Densityplot showing the most significant enrichment
densityplot(log2(esets2[,which.max(sig2[2,])]), groups=lab, plot.points=FALSE)

# Densityplot showing the least significant enrichment
densityplot(log2(esets2[,which.min(sig2[2,])]), groups=lab, plot.points=FALSE)

