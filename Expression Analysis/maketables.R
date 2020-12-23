library('plyr')

M1 <- as.data.frame(read.table('genes0.txt'))
M1 <- M1[-c(1),]
row.names(M1) <- M1[,1]
M1i <- t(M1[,c(4:7)])

M2 <- as.data.frame(read.table('genes2.txt'))
M2 <- M2[-c(1),]
row.names(M2) <- M2[,1]
M2i <- t(M2[,c(4:7)])

M3 <- as.data.frame(read.table('genes6.txt'))
M3 <- M3[-c(1),]
row.names(M3) <- M3[,1]
M3i <- t(M3[,c(4:7)])

M <- rbind.fill.matrix(M1i,M2i,M3i)
M[is.na(M)] <- 0
M <- unique(M)
write.csv(M, 'clust_table.csv')

Ids <- as.data.frame(read.csv('groups.csv'))
ind <- c(Ids[,1])
names <- colnames(M)
Clust <- as.data.frame(cbind(names[ind],Ids[,2:3]))
#row.names(Clust) <- Clust[,1]
MG <- rbind(M1[,1:3],M2[,1:3],M3[,1:3])
MG <- unique(MG)
MG <- MG[Clust[,1],]
Clust <- cbind(Clust, MG[,2:3])
write.csv(Clust, 'genes_grouped.csv')
