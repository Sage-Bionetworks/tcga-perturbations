# Build Myc signature
ent <- loadEntity('syn319533')
dat <- exprs(ent$objects$eset)
mycSig <- calcSig(dat, model.matrix(~factor(rep(1:2,each=9))))
mycSigGenes <- names(eg2sym[gsub("_mt", "_eg", names(which(mycSig$cfs[,2] < -0.75)))])
mycSigGenes <- as.character(unlist(eg2sym[mycSigGenes]))
mycEntity <- Data(list(name='MycSignature', parentId='syn1516116'))
mycEntity <- addObject(mycEntity, mycSigGenes)
mycEntity <- storeEntity(mycEntity)

# Build EGFR signature
ent <- loadEntity('syn317799')
dat <- exprs(ent$objects$eset)
rownames(dat) <- gsub("_mt", "_eg", rownames(dat))
egfrSig <- calcSig(dat, model.matrix(~factor(rep(1:2,c(9,10)))))
egfrSigGenes <- names(eg2sym[gsub("_mt", "_eg", names(which(egfrSig$cfs[,2] < -1.3)))])
egfrSigGenes <- as.character(unlist(eg2sym[egfrSigGenes]))
egfrEntity <- Data(list(name='EgfrSignature', parentId='syn1516116'))
egfrEntity <- addObject(egfrEntity, egfrSigGenes)
egfrEntity <- storeEntity(egfrEntity)


