library('plyr')
library('ggplot2')
library('plotrix')
library('cowplot')
library('reshape')
library('gplots')
library('gridExtra')

M1 <- as.data.frame(read.table('genes0.txt', header = TRUE))
#M1 <- M1[,c(2,4:7)]
M1 <- M1[,c(2,4,6,7)]

M2 <- as.data.frame(read.table('genes2.txt', header = TRUE))
#M2 <- M2[,c(2,4:7)]
M2 <- M2[,c(2,4,6,7)]

M3 <- as.data.frame(read.table('genes6.txt', header = TRUE))
#M3 <- M3[,c(2,4:7)]
M3 <- M3[,c(2,4,6,7)]

M <- merge(M1, M2, by='gene_name')
M <- merge(M, M3, by='gene_name')
names(M) <- c('Gene','CCM1_0','CCM3_0','WT_0','CCM1_2','CCM3_2','WT_2','CCM1_6','CCM3_6','WT_6')
M <- M[,c(1,2,5,8,3,6,9,4,7,10)]
write.csv(M, 'data.csv')

MA <- M[,c(2:10)]
names(MA) <- NULL
write.csv(MA, 'data_values.csv', row.names=FALSE)

IND2 <- as.data.frame(read.table('GeneIND2.txt', header = FALSE))
IND1 <- as.data.frame(read.table('GeneIND1.txt', header = FALSE))
#IND <- intersect(IND1$V1,IND2$V1)
IND <- IND1$V1
S <- M[IND,]

x <- as.matrix(S[,c(2:10)])
row.names(x) <- S$Gene
heatmap.2(x, cellnote=round(x,2), notecol = 'black',
          scale = "row", col = bluered(100), dendogram = 'none', Colv=FALSE, #Rowv = FALSE,
          trace = "none", density.info = "none")

heat <- heatmap.2(x, scale = "row", col = bluered(100), Colv=FALSE, #Rowv = FALSE,
                  trace = "none", density.info = "none")
heat

ord_genes <- rownames(x)[heat$rowInd]
S <- S[match(ord_genes, S$Gene),]
SV <- melt(S, by='Gene')
names(SV) <- c('Gene','Type','Value')
SV$Gene <- factor(SV$Gene, levels = S$Gene)

colors <- c("#FFCCCC", "#FF3333", "#990000", "#99CCFF", "#3399FF", "#004C99", "#E0E0E0", "#A0A0A0", "#404040")
p1 <- ggplot() + 
  geom_point(data = SV, aes(x=Value, y=Gene, color=Type)) + 
  scale_color_manual(values=colors) #+ ylim(-1.2,2)

S <- M[IND,]
S[,2:4] <- (S[,2:4]/S[,8:10]-1)*100
S[,5:7] <- (S[,5:7]/S[,8:10]-1)*100
S[,8:10] <- (S[,8:10]/S[,8:10]-1)*100

x <- as.matrix(S[,c(2:10)])
row.names(x) <- S$Gene
heatmap.2(x, cellnote=round(x,2), notecol = 'black',
          scale = "row", col = bluered(100), dendogram = 'none', Colv=FALSE, #Rowv = FALSE,
          trace = "none", density.info = "none")

ord_genes <- rownames(x)[heat$rowInd]
S <- S[match(ord_genes, S$Gene),]
SP <- melt(S, by='Gene')
names(SP) <- c('Gene','Type','Percent')
SP$Gene <- factor(SP$Gene, levels = S$Gene)

colors <- c("#FFCCCC", "#FF3333", "#990000", "#99CCFF", "#3399FF", "#004C99", "#E0E0E0", "#A0A0A0", "#404040")
p2 <- ggplot() + 
  geom_point(data = SP, aes(x=Percent, y=Gene, color=Type)) + 
  scale_color_manual(values=colors) #+ ylim(-1.2,2)

plot_grid(p1,p2,
          labels = c("Values", "Percent"),
          ncol = 2, nrow = 1)

write(S$Gene, 'SG.txt')

l <- 3*length(S$Gene)
SV$Pheno <- c(rep('CCM1',l),rep('CCM3',l),rep('WT',l))
SP$Pheno <- c(rep('CCM1',l),rep('CCM3',l),rep('WT',l))
  
ggplot(data = SV, aes(x = Pheno, y = Value, color = Type)) +
  geom_point() + 
  scale_color_manual(values=colors) +
  geom_smooth(method = "lm", se = FALSE, lwd = .1) +
  facet_wrap(~Gene, scales = "free_y") 

ggplot(data = SP, aes(x = Pheno, y = Percent, color = Type)) +
  geom_point() + 
  scale_color_manual(values=colors) +
  geom_smooth(method = "lm", se = FALSE, lwd = .1) +
  facet_wrap(~Gene, scales = "free_y")

