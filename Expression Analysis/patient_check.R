library('plyr')
library('ggplot2')
library('plotrix')
library('cowplot')
library('reshape')
library('gplots')
library('gridExtra')

M <- read.delim('patients.txt', header = TRUE, sep = "\t", dec = ".")
#select <- c('ENSG00000105221','ENSG00000173402','ENSG00000111229','ENSG00000143878', 'ENSG00000114251',
#            'ENSG00000001631','ENSG00000136280','ENSG00000114209')
#ccms <- c('ENSG00000173402','ENSG00000072415','ENSG00000164692')
ccms <- c('ENSG00000173402','ENSG00000072415','ENSG00000164692','ENSG00000001631','ENSG00000136280','ENSG00000114209')
M <- M[M$Geneid %in% ccms,]
#M$Gene <- c('CCM1','AKT2','ARPC3','CCM3','WNT5A','CCM2','RHOB','DAG1')
#M$Geneid <- c('MPP5','COL1A2','DAG1')
M$Geneid <- c('CCM1','MPP5','CCM2','CCM3','COL1A2','DAG1')
ord_genes <- c('MPP5','COL1A2','DAG1','CCM1','CCM2','CCM3')
M <- M[match(ord_genes, M$Geneid),]

x <- as.matrix(M[,c(2:15)])
row.names(x) <- M$Geneid
heatmap.2(x, cellnote=round(x,2), notecol = 'black',
          scale = "row", col = bluered(100), dendogram = 'none', Colv=TRUE, Rowv = FALSE,
          trace = "none", density.info = "none")
heatmap.2(x, cellnote=round(x,2), notecol = 'black',
                scale = "col", col = bluered(100), dendogram = 'none', Colv=TRUE, Rowv = FALSE,
                trace = "none", density.info = "none")

S <- melt(M, by='Geneid')
names(S) <- c('Gene','Type','Value')
ggplot(data = S, aes(x = Type, y = Value, color = Type)) +
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE, lwd = .1) +
  facet_wrap(~Gene, scales = "free_y") 

