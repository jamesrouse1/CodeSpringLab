library("DESeq2")

library(dplyr)
library(ggplot2)
#library(EnhancedVolcano)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)

outpath <- args[3]
refcond <- args[4]
compared <- args[5]
redundant <- args[6]

#count <- read.delim(args[1],sep="",row.names = 1,check.names=FALSE)
count <- read.delim(args[1],sep="",row.names = NULL,check.names=FALSE)
modsum <- as.formula(paste(".~",colnames(count)[1],sep=""))
count <- aggregate(modsum,count,sum)
row.names(count) <- count[,1]
count[,1] <- NULL

design <- read.delim(args[2],header=T,sep="\t",row.names = 1,check.names=FALSE)
if (redundant %in% colnames(design)){
    design <- design[,-match(redundant, colnames(design))]
}

design_all <- design
if (ncol(design_all) > 1) {
    design_all <- design_all[1:(ncol(design_all)-1)]
}
overlap_all <- intersect(row.names(design_all), colnames(count))
count_all <- count[, overlap_all, drop=FALSE]
design_all_header <- colnames(design_all)
design_all <- as.data.frame(design_all[overlap_all,,drop=FALSE])
row.names(design_all) <- overlap_all
colnames(design_all) <- design_all_header
for (i in 1:length(colnames(design_all))){
    design_all[,i] <- factor(design_all[,i])
}

if (ncol(count_all) >= 2 && nrow(design_all) >= 2 && length(colnames(design_all)) >= 1) {
    for (i in 1:length(colnames(design_all))){
      if (i == 1){
        model_all <- paste("~",colnames(design_all)[i])
      } else {
        model_all <- paste(model_all,"+",colnames(design_all)[i])
      }
    }
    dds_all <- DESeqDataSetFromMatrix(countData = count_all, colData = design_all, design = as.formula(model_all))
    keep_all <- rowSums(counts(dds_all)) >= 10
    dds_all <- dds_all[keep_all,]
    if (nrow(dds_all) >= 2) {
        vsdata_all <- vst(dds_all, blind=FALSE)
        for (i in 1:length(colnames(design_all))){
            group_col <- colnames(design_all)[i]
            png(paste(outpath,'/pca_all_samples_',group_col,'.png',sep=""),units="in", width=6.8, height=5.8, res=300,type="cairo")
            pca_data <- plotPCA(vsdata_all, intgroup=group_col, returnData=TRUE)
            percent_var <- round(100 * attr(pca_data, "percentVar"))
            pca_data$sample <- rownames(pca_data)
            pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
              geom_point(size = 3.4, alpha = 0.92) +
              geom_text_repel(aes(label = sample), color = "black", size = 3.0,
                              box.padding = 0.45, point.padding = 0.35,
                              segment.color = "grey45", segment.size = 0.3,
                              max.overlaps = Inf, seed = 8) +
              labs(
                title = "PCA: all samples",
                subtitle = paste("Colored by", group_col),
                x = paste0("PC1: ", percent_var[1], "% variance"),
                y = paste0("PC2: ", percent_var[2], "% variance"),
                color = group_col
              ) +
              coord_fixed() +
              theme_classic(base_size = 12) +
              theme(
                plot.title = element_text(face = "bold", size = 15),
                plot.subtitle = element_text(size = 11, color = "grey35"),
                axis.title = element_text(face = "bold"),
                legend.position = "right",
                legend.title = element_text(face = "bold"),
                panel.border = element_rect(color = "grey25", fill = NA, size = 0.5)
              )
            print(pca)
            dev.off()
        }
    }
}

design <- design[1:(ncol(design)-1)]
comparison_candidates <- colnames(design)[vapply(
    design,
    function(column) all(c(refcond, compared) %in% as.character(column)),
    logical(1)
)]
if (length(comparison_candidates) != 1) {
    stop(
        "Could not uniquely identify the comparison column containing both '",
        refcond, "' and '", compared, "'. Candidates: ",
        paste(comparison_candidates, collapse = ", ")
    )
}
pheno_col <- comparison_candidates[[1]]
design <- design[as.character(design[[pheno_col]]) %in% c(refcond, compared), , drop=FALSE]

overlap<-intersect(row.names(design),colnames(count))

count <- count[,overlap]
design_header<-colnames(design)
design <- as.data.frame(design[overlap,])
row.names(design)<-overlap
colnames(design)<-design_header

for (i in 1:length(colnames(design))){
    design[,i] <-factor(design[,i])
}

for (i in 1:length(colnames(design))){
  if (i == 1){
    model <- paste("~",colnames(design)[i])
  } else {
    model <- paste(model,"+",colnames(design)[i])
  }
}

formula <- as.formula(model)

dds <- DESeqDataSetFromMatrix(countData = count,colData = design,design = formula)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

### PCA Plot #####
vsdata <- vst(dds, blind=FALSE)
for (i in 1:length(colnames(design))){

png(paste(outpath,'/pca_',colnames(design)[i],'_',compared,"_vs_",refcond,'(ref).png',sep=""),units="in", width=6.2, height=5.6, res=300,type="cairo")
pca_data <- plotPCA(vsdata, intgroup=colnames(design)[i], returnData=TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
group_col <- colnames(design)[i]
pca_data$sample <- rownames(pca_data)

pca <- ggplot(pca_data, aes(x = PC1, y = PC2, color = .data[[group_col]])) +
  geom_point(size = 3.2, alpha = 0.92) +
  geom_text_repel(aes(label = sample), color = "black", size = 3.0,
                  box.padding = 0.45, point.padding = 0.35,
                  segment.color = "grey45", segment.size = 0.3,
                  max.overlaps = Inf, seed = 8) +
  labs(
    title = paste("PCA:", compared, "vs", refcond),
    subtitle = paste("Colored by", group_col),
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    color = group_col
  ) +
  coord_fixed() +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11, color = "grey35"),
    axis.title = element_text(face = "bold"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    panel.border = element_rect(color = "grey25", fill = NA, size = 0.5)
  )
print(pca)
dev.off()

}
##################

dds@colData@listData[[pheno_col]] <- relevel(dds@colData@listData[[pheno_col]], ref = refcond)

dds <- DESeq(dds)

normalized_counts<-as.data.frame(counts(dds, normalized=TRUE))
normalized_counts$DESCRIPTION <- 'na'
normalized_counts<-normalized_counts[,c(ncol(normalized_counts),1:ncol(normalized_counts)-1)]

#write.table(normalized_counts, file=paste(outpath,"/normalized_counts.txt",sep=""), sep="\t", quote=F, col.names=NA)
write.table(normalized_counts, file=paste(outpath,"/normalized_counts_",compared,"_vs_",refcond,"(ref).txt",sep=""), sep="\t", quote=F, col.names=NA)

resultsNames(dds)

res <- results(dds,name=paste(pheno_col,compared,"vs",refcond,sep="_"))
resLFC <- lfcShrink(dds, coef=paste(pheno_col,compared,"vs",refcond,sep="_"), type="apeglm")

### Volcano Plot (enhanced) ###
png(paste(outpath,'/volcano_',compared,"_vs_",refcond,'(ref).png',sep=""),units="in", width=7, height=5, res=300, ,type="cairo")

##EnhancedVolcano(resLFC,
##                lab = rownames(resLFC),
##                x = 'log2FoldChange',
##                y = 'pvalue',
##                pCutoff = 1e-4,
##                FCcutoff = 1,
##                pointSize = 2.0,
##                labSize = 3.0)

resLFC_data<-data.frame(resLFC)
resLFC_data$diffexpressed <- "NO"
resLFC_data$diffexpressed[resLFC_data$log2FoldChange > 0 & resLFC_data$padj < 0.05] <- "UP"
resLFC_data$diffexpressed[resLFC_data$log2FoldChange < 0 & resLFC_data$padj < 0.05] <- "DOWN"
resLFC_data$delabel <- NA
resLFC_data <- resLFC_data[order(resLFC_data$padj),]
resLFC_data$delabel[1:50]<-rownames(resLFC_data)[1:50]

ggplot(data = resLFC_data, aes(x = log2FoldChange, y = -log10(pvalue),col = diffexpressed,label=delabel))+geom_vline(xintercept = c(-1, 1), col = "red", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.0001), col = "red", linetype = 'dashed') + geom_point(size=0.5)+
  scale_color_manual(values = c("darkblue", "grey", "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  geom_text_repel(max.overlaps = Inf,color="black",size=1)

dev.off()
###########

######## MA Plot ###########
png(paste(outpath,'/MAplot_',compared,"_vs_",refcond,'(ref).png',sep=""),units="in", width=5, height=5, res=300,type="cairo")
plotMA(resLFC, ylim=c(-2,2))
dev.off()
############################

resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resOrdered)

resOrdered <- as.data.frame(resOrdered)
resOrdered <- na.omit(resOrdered)

res_sig <- resOrdered[resOrdered$padj<0.05,]

write.table(res_sig,file=paste(outpath,"/DEG_fdr0.05_",compared,"_vs_",refcond,"(ref).txt",sep=""),sep="\t",quote=F) 
write.table(resOrdered,file=paste(outpath,"/DEG_",compared,"_vs_",refcond,"(ref).txt",sep=""),sep="\t",quote=F)
