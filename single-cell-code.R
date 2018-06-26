#----------------------------------- Import data ---------------------------------------
# Change the directories for your data
geneData = read.csv("./fluidigmpro/filtergeneRawData.csv", header = TRUE, sep = ",", row.names = 1)
isoformData = read.csv("./fluidigmpro/filterisoformRawData.csv", header = TRUE, sep = ",", row.names = 1)
header = colnames(geneData)
names(geneData) = header
names(isoformData) = header

#----------------------------------- density plot ---------------------------------------
# Load files of different types of cell for density plot (noncsc and csc)
# noncsc and csc are files containing information on reads in each cell
d <- density(noncsc$c4/1000000)
d2 <- density(csc$c4/1000000)
plot(d,col="red",xlim=c(0,6),ylim=c(0,0.6), title = NULL, xlab = "Number of transcripts per cell (millions)")
polygon(d, col="#97C62E97",border = "#ffffff")
par(new=T)
plot(d2,col="blue",xlim=c(0,6),ylim=c(0,0.6),xaxt="n",yaxt="n",bty="n",xlab =" ")
polygon(d2, col="#C9275A97",border = "#ffffff")
legend("topright", bty = "n", legend = c("Non-CSCs\n(13 cells)","CSCs\n(28 cells)"), col = c("#97C62E","#C9275A"), pch = 15, ncol = 1, cex = 1)

#----------------------------------- Differential expression ---------------------------------------
library(ggplot2)
library(methods)
library(edgeR)
options(digits=4, width=190)

# Make DGELists
celllines = as.factor(c(rep('non',17),rep('met',3),'dual',rep('emt',10)))
geneList = DGEList(counts=round(geneData), genes=rownames(geneData), group = celllines)
isoformList = DGEList(counts=round(isoformData), genes=rownames(isoformData), group = celllines)

# filter lowly expressed genes/transcripts and recompute the library sizes
geneList = geneList[rowSums(cpm(geneList) > 1) >= 2, , keep.lib.sizes=FALSE]
isoformList = isoformList[rowSums(cpm(isoformList) > 1) >= 2, , keep.lib.sizes=FALSE]

# Apply TMM normalization
geneList = calcNormFactors(geneList, method="TMM")
isoformList = calcNormFactors(isoformList, method="TMM")

# Create design matrix
design = model.matrix(~celllines)
rownames(design) = header

# Data exploration, MD plot
pdf("mdplot.pdf", width=10, height=10, pointsize=12)
par(mfrow = c(2, 2)) # 2 rows, 2 columns
for (i in 1:ncol(geneList)) {
	plotMD(cpm(geneList, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(geneList)[i])
	abline(h=0, col="red", lty=2, lwd=2)
}

par(mfrow = c(2, 2)) # 2 rows, 2 columns
for (i in 1:ncol(geneList)) {
	plotMD(cpm(isoformList, log=TRUE), column=i, xlab = "Average log-expression", ylab = "Expression log-ratio (this sample vs others)", main = colnames(isoformList)[i])
	abline(h=0, col="red", lty=2, lwd=2)
}
dev.off()

# Data exploration, MDS plot
pdf("mdsplot.pdf", width=12, height=6, pointsize=12)
colors = c("blue", "red")
par(mfrow = c(1, 2))
plotMDS(geneList, cex = 1, labels = c(rep("O",87)), col = colors[celllines], main = "Gene-level")
legend("topleft", legend = levels(celllines), col = colors, pch = c(15,16), ncol = 2)
plotMDS(isoformList, cex = 1, labels = c(rep("O",87)), col = colors[celllines], main = "Transcript-level")
legend("topleft", legend = levels(celllines), col = colors, pch = c(15,16), ncol = 2)
dev.off()

# Estimating dispersion
geneList = estimateDisp(geneList, design, robust=TRUE)
isoformList = estimateDisp(isoformList, design, robust=TRUE)

# Differential expression: quasi-likelihood F-test
geneQLF = glmQLFit(geneList, design, robust=TRUE)
geneQLFT = glmQLFTest(geneQLF)

isoformQLF = glmQLFit(isoformList, design, robust=TRUE)
isoformQLFT = glmQLFTest(isoformQLF)

# Uncomment the code for sanity check
# Print top 10 DE genes
# topTags(geneQLFT)

# Print top 10 DE transcripts
# topTags(isoformQLFT)

is.de.gene = decideTestsDGE(geneQLFT, p.value=0.05)
is.de.isoform = decideTestsDGE(isoformQLFT, p.value=0.05)

# Print number of DE genes
# summary(is.de.gene)

# Print number of DE transcripts
# summary(is.de.isoform)

# Print number of up-regulated and down-regulated genes or isoforms
# summary(dt<-decideTestsDGE(geneQLFT))
# summary(dt<-decideTestsDGE(isoformQLFT))

# Average logCPM vs. logFC
pdf("plotsmear.pdf", width=12, height=6, pointsize=12)
par(mfrow = c(1, 2))
plotSmear(geneQLFT, de.tags=rownames(geneQLFT)[is.de.gene != 0], main = "Gene-level")
abline(h=c(-1,1),col="blue")
plotSmear(isoformQLFT, de.tags=rownames(isoformQLFT)[is.de.isoform != 0], main = "Transcript-level")
abline(h=c(-1,1),col="blue")
dev.off()

# Heatmap with top 1000 variable genes
library("gplots")
library("RColorBrewer")
logcounts_gene = cpm(geneList,log=TRUE)
var_genes = apply(logcounts_gene,1,var)
select_genes = names(sort(var_genes, decreasing=TRUE))[1:1000]
variable_lcpm_genes = logcounts_gene[select_genes,]

logcounts_isoform = cpm(isoformList,log=TRUE)
var_isoform = apply(logcounts_isoform,1,var)
select_isoform = names(sort(var_isoform, decreasing=TRUE))[1:1000]
variable_lcpm_isoform = logcounts_isoform[select_isoform,]

mypalette = brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","orange","red","blue")[celllines]

pdf('Heatmapmarker.pdf', width=12, height=12)
par(mfrow = c(1, 2))
geneheat <- heatmap.2(variable_lcpm_genes,col=rev(morecols(20)),trace="none", main="Top 150 most variable genes across samples",ColSideColors=col.cell,scale="row",dendrogram = "none")
isoheat <- heatmap.2(logcounts_isoform,col=rev(morecols(20)),trace="none", main="Top 150 most variable isoforms across samples",ColSideColors=col.cell,scale="row",dendrogram = "none",Rowv = FALSE, Colv = FALSE)
dev.off()

# Gnerate tabular output for differential expression
geneDE = as.data.frame(topTags(geneQLFT, n = nrow(geneQLFT)))
genede = geneDE[geneDE$FDR<0.05,]
write.csv(genede,'Desktop/csciso_nonvsemt.csv')
write.table(geneDE, file="DE_analysis_cscgene.txt", sep = "\t", quote = F, row.names = F, col.names = T)

isoformDE = as.data.frame(topTags(isoformQLFT, n = nrow(isoformQLFT)))
isoformde = isoformDE[isoformDE$FDR<0.05,]
write.csv(isoformde,'Desktop/cscisode.csv')
write.table(isoformDE, file="DE_analysis_csciso.txt", sep = "\t", quote = F, row.names = F, col.names = T)

#----------------------------------- Volcano plot ---------------------------------------
# Label upregulated and downregulated genes
res = NULL
for (i in 1:nrow(f)) {
  if (f[i,2] < -0.25 && f[i,5]<0.05) {
    res <- c(res, 'Up')
  }
  else if (f[i,2] > 0.25 && f[i,5]<0.05) {
    res <- c(res, 'Down')
  }
  else {
    res <- c(res, 'No Change')
  }
}
Significance = res

r03 <- ggplot(f, aes(f$logFC,-1*log10(f$FDR))) + geom_point()
# Modify color according to significance
r03xy = r03 +geom_point(aes(color =Significance))
# Title, labels
r03xy = r03xy + labs(title="Volcano plot for differential expression",x="log2(Fold Change)",y="-log10(FDR)") + theme(plot.title = element_text(hjust = 0.5))
# Add threshold lines
r03xy <- r03xy + geom_hline(yintercept=1.4,linetype=4, color="#686377")+geom_vline(xintercept=c(-1,1),linetype=4, color="#686377")

# Filter out upregulated and downregulated genes
# down <- f[f$logFC<0,]
# down <- down[order(down$logFC, decreasing = FALSE),]
# down <- down[1:10,]
# up <- f[f$logFC>0,]
# up <- up[order(up$logFC, decreasing = TRUE),]
# up <- up[1:10,]

# Label the top 12 genes 
target <- f[order(f$FDR),]
target <- target[1:12,]

#Add text labels 
r03xy <- r03xy + geom_label_repel(data=target, aes(target$logFC, -1*log10(target$FDR),label=target$genes),fontface = 'bold')

#Save the plot
ggsave("~/Desktop/Volcanoplot.png", device="png", width = 8, height = 8)

#----------------------------------- GO analysis ---------------------------------------
#select first 800 genes for GO analysis
genede <- geneDE[order(-abs(geneDE$logFC),geneDE$PValue)[1:800],]

library(clusterProfiler)
library(DO.db)
library(org.Hs.eg.db)

gene_ensembl <- names(genede)
transid <- select(org.Hs.eg.db, keys=gene_ensembl,columns = c("GENENAME","SYMBOL","ENTREZID"),keytype = "ENSEMBL")
geneid <- transid[,4] 

goenrich_BP <- enrichGO(transid$ENTREZID,OrgDb = org.Hs.eg.db,ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
goenrich_CC <- enrichGO(transid$ENTREZID,OrgDb = org.Hs.eg.db,ont = "CC", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
goenrich_MF <- enrichGO(transid$ENTREZID,OrgDb = org.Hs.eg.db,ont = "MF", pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)

# Plotting GO analysis result
library(topGO)
library(Rgraphviz)
# Dot plot
pdf("GO_dotPlot.pdf", width=12, height=6)
par(mfrow = c(3, 1))
dotplot(goenrich_BP,font.size=8, showCategory = 30, title = "GO analysis on Cellular Component")
dotplot(goenrich_CC,font.size=8, showCategory = 30, title = "GO analysis on Cellular Component")
dotplot(goenrich_MF,font.size=8, showCategory = 30, title = "GO analysis on Molecular Function")
dev.off()

# Bar plot
pdf("GO_barPlot.pdf", width=12, height=6)
par(mfrow = c(3, 1))
barplot(goenrich_BP,font.size = 8, showCategory = 30, title = "GO analysis on Biological Process")
barplot(goenrich_CC,font.size = 8, showCategory = 30, title = "GO analysis on Cellular Component")
barplot(goenrich_MF,font.size = 8, showCategory = 30, title = "GO analysis on Molecular Function")
dev.off()

# Net plot
enrichMap(goenrich_BP, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
enrichMap(goenrich_CC, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)
enrichMap(goenrich_MF, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)

pdf("GO_netGraph.pdf", width=6, height=24)
par(mfrow = c(3, 1))
plotGOgraph(goenrich_BP)
plotGOgraph(goenrich_CC)
plotGOgraph(goenrich_MF)
dev.off()

# KEGG
keggenrich <- enrichKEGG(gene = geneid, organism = "hsa",pvalueCutoff = 0.05)

# View first 6 enriched terms
# head(kk)

# Plotting result
dotplot(keggenrich, showCategory = 30, title = "KEGG analysis")
barplot(keggenrich, showCategory = 30, title = "KEGG analysis")
enrichMap(keggenrich, showCategory = 30, title = "KEGG analysis")

# GSEA analysis
gse_MF <- gseGO(genego, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "MF")
gse_BP <- gseGO(genego, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "BP")
gse_CC <- gseGO(genego, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "CC")

# view first 6 items
# head(gse_MF)
# head(gse_CC)
# head(gse_BP)

gseaplot(gse_CC, geneSetID = "GO:0016020")

#----------------------------------- PCA ---------------------------------------
library(FactoMineR)
library(factoextra)
PCA.allgenes = PCA(t(geneData, ncp=4, graph=T)

# Pretty PCA plot
fviz_pca_ind(PCA.allgenes,label='none', pointsize=8,pointshape=21,fill.ind = types$types,mean.point = F)+
  scale_color_jama()+
  scale_fill_jama()+
  theme(axis.line.x = element_line(size = 1, colour = "black"),
        axis.line.y = element_line(size = 1, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
# Variable correlation plot
fviz_pca_var(PCA.allgenes,pointsize=100,label='none')+
  theme(axis.line.x = element_line(size = 0.5, colour = "black"),
        axis.line.y = element_line(size = 0.5, colour = "black"),
        axis.line = element_line(size=1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 10),
        axis.text.y=element_text(colour="black", size = 10),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        panel.border = element_rect(colour = "black", fill=NA, size=2))


#----------------------------------- Visualization ---------------------------------------

# Barplot for genes with P-values
ggplot(data=genede, aes(x=reorder(as.character(f1$gene),rev(1:30)),y=f$FDR)) + geom_bar(fill="grey35", colour="white",stat="identity") + theme_bw() +
       ylab(expression(paste("corrected p-value (", log[10],")",sep=""))) +
       xlab(expression(paste("marker genes",sep=""))) +
       scale_y_log10() +
       coord_flip()
ggsave(filename = "Barplot.pdf", plot = last_plot(),  path = NULL, scale = 1, width = 100, height = 200, units ="mm", dpi = 300)

# Boxplot for a particular genes
kk$group <- as.factor(kk$group)
p <- ggplot(vio,aes(x=vio$cell,y=vio$ZDHHC23, fill=vio$cell))+geom_boxplot()
p <- p+scale_fill_brewer(palette="Dark2") + theme_bw()
p <- p+scale_y_continuous(name = "log10(Read count)",
                          breaks = seq(0, 175, 25),
                          limits=c(0,2)) +scale_x_discrete(name = "CD44-209")
p <- p+theme(axis.line.x = element_line(size = 0.5, colour = "black"),
             axis.line.y = element_line(size = 0.5, colour = "black"),
             axis.line = element_line(size=1, colour = "black"),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_blank(),
             panel.background = element_blank(),
             axis.text.x=element_text(colour="black", size = 10),
             axis.text.y=element_text(colour="black", size = 10),
             axis.title.y = element_text(size = 12),
             axis.title.x = element_text(size = 12))
ggsave(filename = "cd209.pdf", plot = last_plot(),  path = NULL, scale = 1, width = 100, height = 100, units ="mm", dpi = 300)

# multiple box plot
# each column (gene) generates a box plot
# label the cells with types
# tmp is a file with cell name as row and gene name as column
library(reshape2)
library(ggpubr)
fin=melt(tmp,id.vars = "types")
ggplot(fin,aes(types,value,fill=types))+
  geom_boxplot()+stat_smooth()+facet_wrap(~variable)+
  scale_fill_brewer(palette="Dark2")+
  stat_compare_means(method='wilcox.test')+
  theme_bw()

# multiple dot plot
ggplot(fin,aes(types,value,fill=types))+
  geom_dotplot(binaxis = 'y',dotsize = 4)+
  facet_wrap(~variable)+scale_fill_brewer(palette="Dark2")+
  stat_smooth()+theme_bw()+theme(text = element_text(size=10))+
  ylim(-1,15)

# mutiple violin plot
ggplot(fin,aes(types,value,fill=types))
  +geom_violin()+stat_smooth()+facet_wrap(~variable)+
  scale_fill_brewer(palette="Dark2")+
  stat_summary(fun.y=mean, geom="point",shape=15, size=1, color="black")+
  coord_flip()+theme_bw()

# Saving images compatible with Adobe Illustrator
# ggsave(plot=p,height=6,width=6,dpi=200, filename="~/example.pdf", useDingbats=FALSE)
