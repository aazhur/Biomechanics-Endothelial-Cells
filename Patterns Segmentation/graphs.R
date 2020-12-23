library('plyr')
library('ggplot2')
library('plotrix')
library('cowplot')
library('reshape')
library('gplots')
library('gridExtra')

SA <- as.data.frame(read.delim('sim_area.txt', header = FALSE))
SH <- as.data.frame(read.delim('sim_holes.txt', header = FALSE))
S <- cbind(SA,SH)
names(S) <- c('Area', 'Holes')

EA <- as.data.frame(read.delim('exp_area.txt', header = FALSE))
EH <- as.data.frame(read.delim('exp_holes.txt', header = FALSE))
ET <- as.data.frame(read.delim('exp_types.txt', header = FALSE))
E <- cbind(EA,EH,ET)
names(E) <- c('Area', 'Holes','Type')
EN <- E[(E$Type == '0' | E$Type == '1' |  E$Type == '2' | E$Type == '3'), ]
ER <- E[(E$Type == '0R' | E$Type == '1R' |  E$Type == '2R' | E$Type == '3R'), ]   

p <- ggplot(data = EN, aes(x = Area, y = Holes, color = Type, shape = Type)) + 
  geom_point() + stat_ellipse(type = "norm")
p

p <- ggplot(data = ER, aes(x = Area, y = Holes, color = Type, shape = Type)) + 
  geom_point() + stat_ellipse(type = "norm")
p

p <- ggplot(data = E, aes(x = Area, y = Holes, color = Type)) + 
  geom_point()
p

p <- ggplot(data = S, aes(x = Area, y = Holes)) + 
  geom_point()
p
