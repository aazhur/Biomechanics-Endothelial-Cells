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

M1 <- as.data.frame(read.table('genes0.txt', header = TRUE))
#M1 <- M1[M1$gene_name %in% INT,]
#M1 <- M1[,c(2,4:7)]
M1 <- M1[,c(2,4,6,7)]
# M1 <- M1[M1$NT_fpkm > 0,]
M1$dec <- sign(M1$ccm1_fpkm-M1$NT_fpkm) + sign(M1$ccm3_fpkm-M1$NT_fpkm)
M1 <- M1[M1$dec == 0,]

M2 <- as.data.frame(read.table('genes2.txt', header = TRUE))
#M2 <- M2[M2$gene_name %in% INT,]
#M2 <- M2[,c(2,4:7)]
M2 <- M2[,c(2,4,6,7)]
# M2 <- M2[M2$NT_fpkm > 0,]
M2$dec <- sign(M2$ccm1_fpkm-M2$NT_fpkm) + sign(M2$ccm3_fpkm-M2$NT_fpkm)
M2 <- M2[M2$dec == 0,]

M3 <- as.data.frame(read.table('genes6.txt', header = TRUE))
#M3 <- M3[M3$gene_name %in% INT,]
#M3 <- M3[,c(2,4:7)]
M3 <- M3[,c(2,4,6,7)]
# M3 <- M3[M3$NT_fpkm > 0,]
M3$dec <- sign(M3$ccm1_fpkm-M3$NT_fpkm) + sign(M3$ccm3_fpkm-M3$NT_fpkm)
M3 <- M3[M3$dec == 0,]

M <- merge(M1, M2, by='gene_name')
M <- merge(M, M3, by='gene_name')

#M <- merge(M1, M2, by="gene_name", all = T)
#M <- merge(M, M3, by="gene_name", all = T)
#M[is.na(M)] <- 0
M <- M[,c(1:4,6:8,10:12)]
#names(M) <- c('Gene',
#             'CCM1_0','CCM2_0','CCM3_0','WT_0',
#             'CCM1_2','CCM2_2','CCM3_2','WT_2',
#             'CCM1_6','CCM2_6','CCM3_6','WT_6')
names(M) <- c('Gene',
              'CCM1_0','CCM3_0','WT_0',
              'CCM1_2','CCM3_2','WT_2',
              'CCM1_6','CCM3_6','WT_6')
#M[,2:5] <- M[,2:5]/(M[,5]+10^(-6))-1
#M[,6:9] <- M[,6:9]/(M[,9]+10^(-6))-1
#M[,10:13] <- M[,10:13]/(M[,13]+10^(-6))-1

M$CCM1 <- sign(M$CCM1_0-M$WT_0) + sign(M$CCM1_2-M$WT_2) + sign(M$CCM1_6-M$WT_6)
M$CCM3 <- sign(M$CCM3_0-M$WT_0) + sign(M$CCM3_2-M$WT_2) + sign(M$CCM3_6-M$WT_6)
M <- M[(M$CCM1 == 3 & M$CCM3 == -3) | (M$CCM3 == 3 & M$CCM1 == -3),]
# 
# M$CCM1_02 <- sign(M$CCM1_2-M$CCM1_0)
# M$CCM1_26 <- sign(M$CCM1_6-M$CCM1_2)
# M$CCM3_02 <- sign(M$CCM3_2-M$CCM3_0)
# M$CCM3_26 <- sign(M$CCM3_6-M$CCM3_2)
# MO <- M[((M$CCM1_02 + M$CCM1_26) == 2 & (M$CCM3_02 + M$CCM3_26) == -2 & M$CCM1 == 3) |
#          ((M$CCM3_02 + M$CCM3_26) == 2 & (M$CCM1_02 + M$CCM1_26) == -2 & M$CCM1 == -3),]
# #
# M$CCM1_02 <- M$CCM1_0/M$CCM1_2 - 1
# M$CCM1_26 <- M$CCM1_2/M$CCM1_6 - 1
# M$CCM3_02 <- M$CCM3_0/M$CCM1_2 - 1
# M$CCM3_26 <- M$CCM3_2/M$CCM1_6 - 1
# MM <- M[((M$CCM1_02 <= 0.3 & M$CCM1_02 >= - 0.3) & (M$CCM1_26 <= 0.3 & M$CCM1_26 >= - 0.3)) &
#           ((M$CCM3_02 <= 0.3 & M$CCM3_02 >= - 0.3) & (M$CCM3_26 <= 0.3 & M$CCM3_26 >= - 0.3)),]
# 
# M$CCM1_02 <- sign(M$CCM1_2-M$CCM1_0)
# M$CCM1_26 <- sign(M$CCM1_6-M$CCM1_2)
# M$CCM3_02 <- M$CCM3_0/M$CCM1_2 - 1
# M$CCM3_26 <- M$CCM3_2/M$CCM1_6 - 1
# MX1 <- M[(((M$CCM1_02 + M$CCM1_26) == 2 & M$CCM1 == 3) | ((M$CCM1_02 + M$CCM1_26) == -2  & M$CCM1 == -3)) &
#           ((M$CCM3_02 <= 0.3 & M$CCM3_02 >= - 0.3) & (M$CCM3_26 <= 0.3 & M$CCM3_26 >= - 0.3)),]
# 
# M$CCM1_02 <- M$CCM1_0/M$CCM1_2 - 1
# M$CCM1_26 <- M$CCM1_2/M$CCM1_6 - 1
# M$CCM3_02 <- sign(M$CCM3_2-M$CCM3_0)
# M$CCM3_26 <- sign(M$CCM3_6-M$CCM3_2)
# MX2 <- M[(((M$CCM3_02 + M$CCM3_26) == 2 & M$CCM3 == 3) | ((M$CCM3_02 + M$CCM3_26) == -2  & M$CCM3 == -3)) &
#            ((M$CCM1_02 <= 0.3 & M$CCM1_02 >= - 0.3) & (M$CCM1_26 <= 0.3 & M$CCM1_26 >= - 0.3)),]
# #
# MN <- rbind(MO, MM, MX1, MX2)
# M <- rbind(MO, MM, MX1, MX2)
# 
MN <- M
M[,c(2:10)] <- M[,c(2:10)] + 0.000001
M[,2:4] <- M[,2:4]/M[,4]-1
M[,5:7] <- M[,5:7]/M[,7]-1
M[,8:10] <- M[,8:10]/M[,10]-1

lower <- -0.5
upper <-  0.5
MN <- MN[((M$CCM1_0 >= upper | M$CCM1_0 <= lower) & (M$CCM1_2 >= upper | M$CCM1_2 <= lower) & (M$CCM1_6 >= upper | M$CCM1_6 <= lower)) |
           ((M$CCM3_0 >= upper | M$CCM3_0 <= lower) & (M$CCM3_2 >= upper | M$CCM3_2 <= lower) & (M$CCM3_6 >= upper | M$CCM3_6 <= lower)),]
#M <- M[(M$CCM1_6 >= upper | M$CCM1_6 <= lower),]
#M <- M[(M$CCM3_6 >= upper | M$CCM3_6 <= lower),]
M <- MN[,c(1,2,5,8,3,6,9,4,7,10)]
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