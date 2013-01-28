library(org.Hs.eg.db)
eg2sym <- as.character(org.Hs.egSYMBOL)
names(eg2sym) <- paste(names(eg2sym), '_eg', sep="")
hsym <- as.character(org.Hs.egSYMBOL)
overlapping.genes <- gsub("_eg", "", rownames(fic.p))
hgo <- as.list(org.Hs.egGO2EG)
hgo <- sapply(hgo, function(x){ intersect(x, overlapping.genes)})
hpath <- as.list(org.Hs.egPATH2EG)
hpath <- sapply(hpath, function(x){ intersect(x, overlapping.genes)})

invert.list <- function(il) {
	lens = sapply(il, length);
	v1 = rep(names(il), lens);
	v2 = unlist(il, use.names=FALSE)
	ol = split(v1, v2);
	ol
}

goh <- invert.list(hgo)
pathh <- invert.list(goh)

library(GO.db)
xx <- as.list(GOTERM)

sapply(xx, function(x){ 
			c(GOID(x), Term(x))
		}) -> go2names
go2names <- t(go2names)
go2n <- go2names[,2]
names(go2n) <- go2names[,1]

library(KEGG.db)
kegg2n <- as.character(KEGGPATHID2NAME)
all.genes <- unique(unlist(hgo))

get.enrich <- function(genes, all.genes, hsym, ggo, go2n, string=NULL) {
	ggo <- ggo[which(sapply(ggo,length) > 2)]
	genes <- intersect(genes, unlist(ggo))
	all.genes <- intersect(all.genes, unlist(ggo))
	k <- length(genes)
	n1 <- length(all.genes)
	sapply(ggo, function(go) { 
				c(length(intersect(go, genes)),
						length(unique(go)),
						n1 - length(unique(go)),
						k,
						1-phyper(length(intersect(go, genes)),
								length(unique(go)),
								n1 - length(unique(go)),
								k)
				)
			}) -> enrich
	ggo.i <- sapply(ggo, function(x) { intersect(x, genes)}) 
	enrich <- t(enrich)
	enrich <- enrich[which(enrich[,1] > 0),]
	colnames(enrich) <- c("Num.Sig.Genes.in.GO","Num.Genes.in.GO", "Num.Genes.Not.in.GO","Num.Sig.Genes","P-value")
	colnames(enrich) <- gsub("GO",string,colnames(enrich))
	enrich <- data.frame(GO.name=unlist(go2n[rownames(enrich)]), enrich)
	names(enrich) <- gsub("GO",string,names(enrich))
	enrich <- enrich[ order(enrich$P.value),]
	list(enrich=enrich, ggo.i=ggo.i)
}

genes <- gsub("_eg", "", rownames(fic.p2)[1:1000])

go.enrich <- get.enrich(genes, all.genes, hsym,hgo,go2n,"GO")
kegg.enrich <- get.enrich(genes, all.genes,hsym,hpath,kegg2n,"KEGG")

goe <- go.enrich$enrich
goe <- goe[which(goe[,2] > 1 & goe[,6] < 0.05),]
goe[,6] <- round(goe[,6],3)
write.table(goe,file="gic_go_enrichment.xls",sep="\t",quote=FALSE)

goe <- go.enrich$enrich
goe <- goe[which(goe[,2] > 1 & goe[,6] < 0.05),]
goe[,6] <- round(goe[,6],3)
write.table(goe,file="gic_go_enrichment.xls",sep="\t",quote=FALSE)


koe <- kegg.enrich$enrich
koe <- koe[which(koe[,2] > 1 & koe[,6] < 0.05),]
koe[,6] <- round(koe[,6],3)
write.table(koe,file="gic_kegg_enrichment.xls",sep="\t",quote=FALSE)

