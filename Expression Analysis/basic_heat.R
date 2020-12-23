library('ggplot2')
library('plotrix')
library('cowplot')
library('reshape')
library('gplots')
library('gridExtra')
library('scatterplot3d')

M1 <- as.data.frame(read.table('Contacts-0hrs.txt'))
LM1 <- M1[,c(2,4,5,6,7)]
M1 <- M1[,c(2,4,5,6,7)]
names(M1) <- c('Gene',	'CCM1_0',	'CCM2_0',	'CCM3_0',	'WT_0')
names(LM1) <- c('Gene',	'CCM1',	'CCM2',	'CCM3',	'WT')
LM1$Hour <- rep("0", 42)

M2 <- as.data.frame(read.table('Contacts-2hrs.txt'))
LM2 <- M2[,c(2,4,5,6,7)]
M2 <- M2[,c(4,5,6,7)]
names(M2) <- c('CCM1_2',	'CCM2_2',	'CCM3_2',	'WT_2')
names(LM2) <- c('Gene',	'CCM1',	'CCM2',	'CCM3',	'WT')
LM2$Hour <- rep("2", 42)

M3 <- as.data.frame(read.table('Contacts-6hrs.txt'))
LM3 <- M3[,c(2,4,5,6,7)]
M3 <- M3[,c(4,5,6,7)]
names(M3) <- c('CCM1_3',	'CCM2_3',	'CCM3_3',	'WT_3')
names(LM3) <- c('Gene', 'CCM1',	'CCM2',	'CCM3',	'WT')
LM3$Hour <- rep("6", 42);

M <- cbind(M1,M2,M3)
for (i in 2:13){
  M[,i] <- as.double(M[,i])
}

LM1 <- melt(LM1, id=c("Gene","Hour"))
LM1 <- LM1[,c(1,4)]
names(LM1) <- c('Gene','Exp_0')

LM2 <- melt(LM2, id=c("Gene","Hour"))
LM2 <- LM2[,c(1,4)]
names(LM2) <- c('Gene','Exp_2')

LM3 <- melt(LM3, id=c("Gene","Hour"))
LM3 <- LM3[,c(4,3)]
names(LM3) <- c('Exp_6','Type')

LM <- cbind(LM1, LM2, LM3)
LM <- LM[,c(1,2,4,5,6)]

Genes <- M1$Gene

p <- list()
d <- list()
i <- 1
for(name in Genes){
  d[[i]] <- LM[LM$Gene == name, ]
  p[[i]] <- ggplot(data = d[[i]]) + 
    geom_point(aes(x = rep(0,4), y = Exp_0, color = Type), size = 1) + 
    geom_point(aes(x = rep(2,4), y = Exp_2, color = Type), size = 1) + 
    geom_point(aes(x = rep(6,4), y = Exp_6, color = Type), size = 1) + 
    ggtitle(name) + theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())
  i <- i+1  
}
do.call("grid.arrange", c(p, ncol=7))

#p <- ggplot(data = LM, aes(x = Exp_0, y = Exp_6, label = Gene)) + 
#  geom_point(aes(color = Type)) + 
#  geom_text() # +
  #xlim(0,100) + ylim(0,100)
#p

colors <- c("#FF6666", "#FFFF33", "#00CC66", "#99CCFF")
colors <- colors[as.numeric(LM$Type)]
scatterplot3d(LM[,2:4], pch = 16, color=colors)

x <- as.matrix(M[,c(3,7,11,2,6,10,4,8,12,5,9,13)])
row.names(x) <- M$Gene
#x <- x[order(x[,'WT_3']),]

substract = x[,10]
for (i in 1:12){
  x[,i] <- (x[,i]-substract)
}
x <- x/(substract+0.0001)

#catenin <- c("APC","GSK3B","CSNK1A1","DACT1","CTNNB1","WNT2B")
#actinin <- c("ACTN1","ACTN2","ACTN3","ACTN4")
heatmap.2(x[,4:12], cellnote=round(x[,4:12], 2), notecol = 'black',
          scale = "row", col = bluered(100), dendogram = 'none', Colv=FALSE, Rowv = FALSE,
          trace = "none", density.info = "none")

d <- melt(M, id=c('Gene'))
names(d) <- c("Gene","Type","Expression")
ggplot(d, aes(x = Type, y = Expression, colour = Type)) + geom_boxplot() + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)

ggplot(d, aes(x = Gene, y = Expression, colour = Type)) + geom_point()
