library(ggplot2)
library(org.Mm.eg.db)
library(limma)
library(edgeR)
library(clusterProfiler)

myTheme <- theme(panel.background = element_blank(),
                 axis.line = element_line(color = "black"),
                 text = element_text(color = 'black', size = 12),
                 legend.key=element_blank())

drugData <- read.csv("drugCounts.csv")
crData <- read.csv("crCounts.csv")

rownames(drugData) <- drugData[,1]
drugData <- drugData[,-1]

rownames(crData) <- crData[,1]
crData <- crData[,-1]

drugGroup <- c(rep("Con", 8), rep("RIL", 4))
crGroup <- c(rep("CR", 3), rep("Con", 3), rep("Rap", 3))

xDrug <- DGEList(counts = drugData, genes = rownames(drugData))
xCR <- DGEList(counts = crData, genes = rownames(crData))

xDrug$symbol <- mapIds(org.Mm.eg.db, rownames(xDrug),
                       keytype = "ENSEMBL", column = "SYMBOL")
xDrug$entrez <- mapIds(org.Mm.eg.db, rownames(xDrug),
                       keytype = "ENSEMBL", column = "ENTREZID")

xCR$symbol <- mapIds(org.Mm.eg.db, rownames(xCR),
                     keytype = "ENSEMBL", column = "SYMBOL")
xCR$entrez <- mapIds(org.Mm.eg.db, rownames(xCR),
                     keytype = "ENSEMBL", column = "ENTREZID")

keepD <- rowSums(cpm(xDrug) > 1) >= 4
keepC <- rowSums(cpm(xCR) > 1) >= 4

table(keepD)
table(keepC)

xDrug <- xDrug[keepD, ,keep.lib.sizes = F]
xDrug <- calcNormFactors(xDrug)

xCR <- xCR[keepC, ,keep.lib.sizes = F]
xCR <- calcNormFactors(xCR)

modMatDrug <- model.matrix( ~ drugGroup)
modMatCR <- model.matrix( ~ crGroup)

yDrug <- voom(xDrug, modMatDrug, plot = T)
yCR <- voom(xCR, modMatCR, plot = T)

fitDrug <- lmFit(yDrug, modMatDrug)
fitCR <- lmFit(yCR, modMatCR)

contDrug <- contrasts.fit(fitDrug, coef = 2)
contCR <- contrasts.fit(fitCR, coef = 2)
contRap <- contrasts.fit(fitCR, coef = 3)

contDrug <- eBayes(contDrug)
contCR <- eBayes(contCR)
contRap <- eBayes(contRap)

topDrug <- topTable(contDrug, sort.by = "P", n = Inf)
topCR <- topTable(contCR, sort.by = "P", n = Inf)
topRap <- topTable(contRap, sort.by = "P", n = Inf)

topDrug$sig <- F
topDrug[topDrug$adj.P.Val < 0.05, ]$sig <- T

topCR$sig <- F
topCR[topCR$adj.P.Val < 0.05, ]$sig <- T

topRap$sig <- F
topRap[topRap$adj.P.Val < 0.05 , ]$sig <- T

table(topDrug$adj.P.Val < 0.05)
table(topCR$adj.P.Val < 0.05)
table(topRap$adj.P.Val < 0.99)

