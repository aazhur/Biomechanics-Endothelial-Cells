library('plyr')
library('ggplot2')
library('plotrix')
library('cowplot')
library('reshape')
library('gplots')
library('gridExtra')

remove <- c("")
M <- read.delim('FA.txt', header = FALSE, sep = "\t", dec = ".")
M <- c(as.matrix(M))
FA <- M[! M %in% remove]

M <- read.delim('AJ.txt', header = FALSE, sep = "\t", dec = ".")
M <- c(as.matrix(M))
AJ <- M[! M %in% remove]

M <- read.delim('TJ.txt', header = FALSE, sep = "\t", dec = ".")
M <- c(as.matrix(M))
TJ <- M[! M %in% remove]

M <- read.delim('GJ.txt', header = FALSE, sep = "\t", dec = ".")
M <- c(as.matrix(M))
GJ <- M[! M %in% remove]

M <- read.delim('RAC.txt', header = FALSE, sep = "\t", dec = ".")
M <- c(as.matrix(M))
RAC <- M[! M %in% remove]

M <- read.delim('ECM.txt', header = FALSE, sep = "\t", dec = ".")
M <- c(as.matrix(M))
ECM <- M[! M %in% remove]

#INT <- Reduce(intersect, list(ECM,AJ))
#INT <- c(INT, 'IQGAP1')
#INT <- c('CDC42','RAC1','IQGAP1','ARHGDIA','CTNNB1','APC','GSK3B','CTNNA1','DVL2')
#INT <- c('CDC42','RAC1','VCL','MYH9','RHOA','PTK2','SRC','PAK1','PXN')
#INT <- c('RHOB','RAC1','RHOC','RHOG','PTK2','IL8','DKK1','MAPK9','WNT3')
#INT <- c('WNT2','WNT2B','WNT4','WNT9A','WNT8B','WNT11','WNT10B','WNT9B','WNT3','WNT7B','WNT5A')
#INT <- c('ARPC1A','ARPC1A','ARPC2','ARPC3','ARPC4,ARPC4-TTLL3,TTLL3','ARPC5','ACTR2')
#INT <- FA
#INT <- c('DAG1','RHOB','ARPC3','RAC1','AKT2','WNT5A','ROR1')
#INT <- c('WNT3','WNT5A')
#INT <- c('WNT3','WNT5A','FZD2','ROR1','RAC1','RHOA','CTNNB1','DVL1','DVL2','DVL3','DKK1','PTPB1','LRP5','LRP6','GSK3B','CDC42','PXN','JNK','VCL')
#INT <- c('ABCB1','CTNNB1')
INT <- c('C3AR1','C5AR1')

M1 <- as.data.frame(read.table('genes0.txt', header = TRUE))
M1 <- M1[M1$gene_name %in% INT,]
#M1 <- M1[,c(2,4:7)]
M1 <- M1[,c(2,4,6,7)]

M2 <- as.data.frame(read.table('genes2.txt', header = TRUE))
M2 <- M2[M2$gene_name %in% INT,]
#M2 <- M2[,c(2,4:7)]
M2 <- M2[,c(2,4,6,7)]

M3 <- as.data.frame(read.table('genes6.txt', header = TRUE))
M3 <- M3[M3$gene_name %in% INT,]
#M3 <- M3[,c(2,4:7)]
M3 <- M3[,c(2,4,6,7)]

M <- merge(M1, M2, by='gene_name')
M <- merge(M, M3, by='gene_name')

#M <- merge(M1, M2, by="gene_name", all = T)
#M <- merge(M, M3, by="gene_name", all = T)
#M[is.na(M)] <- 0
#names(M) <- c('Gene',
#             'CCM1_0','CCM2_0','CCM3_0','WT_0',
#             'CCM1_2','CCM2_2','CCM3_2','WT_2',
#             'CCM1_6','CCM2_6','CCM3_6','WT_6')
names(M) <- c('Gene',
              'CCM1_0','CCM3_0','WT_0',
              'CCM1_2','CCM3_2','WT_2',
              'CCM1_6','CCM3_6','WT_6')

#M[,c(2:10)] <- M[,c(2:10)] + 0.000001
#M[,2:4] <- M[,2:4]/M[,4]-1
#M[,5:7] <- M[,5:7]/M[,7]-1
#M[,8:10] <- M[,8:10]/M[,10]-1

#M <- MN[,c(1,2,5,8,3,6,9,4,7,10)]
#M <- M[,c(1,2,5,8,3,6,9,4,7,10)]

#M[,2:5] <- M[,2:5] - M[,5]
#M[,6:9] <- M[,6:9] - M[,9]
#M[,10:13] <- M[,10:13] - M[,13]

#x <- as.matrix(M[,c(2:13)])
x <- as.matrix(M[,c(2:10)])
row.names(x) <- M$Gene
#x <- x[,c(1,5,9,2,6,10,3,7,11,4,8,12)]

heatmap.2(x, cellnote=round(x,2), notecol = 'black',
          scale = "row", col = bluered(100), dendogram = 'none', Colv=FALSE, #Rowv = FALSE,
          trace = "none", density.info = "none")

heat <- heatmap.2(x, scale = "row", col = bluered(100), Colv=FALSE, #Rowv = FALSE,
          trace = "none", density.info = "none")
heat

ord_genes <- rownames(x)[heat$rowInd]
M <- M[match(ord_genes, M$Gene),]
ML <- melt(M, by='Gene')
names(ML) <- c('Gene','Type','Percent_Change')
ML$Gene <- factor(ML$Gene, levels = M$Gene)

colors <- c("#FFCCCC", "#FF3333", "#990000", "#99CCFF", "#3399FF", "#004C99", "#E0E0E0", "#A0A0A0", "#404040")
p <- ggplot() + 
  geom_point(data = ML, aes(x=Percent_Change, y=Gene, color=Type)) + 
  scale_color_manual(values=colors) #+ ylim(-1.2,2)
p

p <- list()
d <- list()
i <- 1
for(name in M$Gene){
  d[[i]] <- ML[ML$Gene == name, ]
  p[[i]] <- ggplot(data = d[[i]]) + 
    geom_point(aes(x=Percent_Change, y=Gene, color=Type)) + 
    scale_color_manual(values=colors)
  i <- i+1  
}
do.call("grid.arrange", c(p))

write(ord_genes, 'selected.txt')

