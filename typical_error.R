library(tidyverse)
library(rstatix)
library(reshape2)
setwd('D:/GLP1RA/pheno_data/')

glu = read.csv('GLU.csv', row.names = 2)
ins = read.csv('INS.csv', row.names = 2)
homaIR = (ins[, 3:13] * glu[,3:13]/18) / 22.5   # /18 for glulcose unit transter
f = cbind(Animal_ID=glu$Animal_ID, Sample_ID=rownames(glu), Group=glu$Group, homaIR)
df = f[f$Group=='VEH', c("Animal_ID", "Day.1.",  "Day.8.", "Day.15.", "Day.22.", "Day.29.",   
                         "Day.36.", "Day.43.", "Day.50.", "Day.57." )]
colnames(df)[2:10] = paste0('t', 1:9)
df = melt(df, id.vars = "Animal_ID")

aa = lm(value~as.factor(variable)+as.factor(Animal_ID), data=df)
summary(aa) # extract !!Residual standard error!!  16.24
rsd = 21.54
TE = rsd/sqrt(2)
SI = (f$Day..7. - f$Day.57.)/TE
Effect_this = ifelse(SI>2, 'R', 'NR')
Effect_this[which(f$Group=='VEH')] = 'nan'
f = cbind(f, TE, SI, Effect_this)
write.csv(f, 'HOMA_IR.csv', row.names=F)

