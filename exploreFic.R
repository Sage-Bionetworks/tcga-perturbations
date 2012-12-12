
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

# Explore FIC matrices using Erich's MYC Signature
ficEntity <- loadEntity('syn375592')
fic <- ficEntity$objects$fic
studyEntity <- loadEntity('syn375792')
studies <- studyEntity$objects$gseAnnotation
ficEnrichment <- calcSignatureEnrichment(fic, as.character(eg2sym[mycSigGenes]))
descriptions <- studies[names(which(ficEnrichment[2,] > 0.6)),]

# Do the same for erich's myc factor model
mycFromErich <- loadEntity('syn1513426')
factors <- mycFromErich$objects$MYCannotatedFactors
sapply(1:length(factors), function(i){ 
			ficEnrichment <- calcSignatureEnrichment(fic, as.character(unlist(factors[[i]])))[2,]
		}) -> res

cor(res)


