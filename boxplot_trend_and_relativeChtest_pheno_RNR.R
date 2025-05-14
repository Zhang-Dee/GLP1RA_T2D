setwd('D:/GLP1RA/')
rnr = read.csv('response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR based on homa-IR improvement
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 
library(reshape2)
library(ggplot2)
library(stringr)
library(ggpubr)

##-------------------------------------- HOMA-IR grouping visual (SI~)----------------------------
f = read.csv('HOMA_IR.csv', row.names = 1)
colnames(f)[3:(ncol(f)-3)] = paste0('d', str_extract(colnames(f)[3:(ncol(f)-3)], '[0-9]+'))
f = f[-which(rownames(f) %in% test_id),]#; f = f[f$Group!='VEH', ]

# 1. homa-IR trend
df = f[f$Group!='VEH', c('Sample_ID', 'Effect_this', 'd7', 'd29', 'd57')]  # Day-1 -> Day57 
colnames(df)=c('Sample_ID', 'Response', 'T1','T2','T3') 
dat = melt(df, id.vars = c('Sample_ID', 'Response'), variable.name = 'time', value.name = 'HOMAIR')
ggplot(dat, aes(x=time, y=HOMAIR)) + 
  geom_boxplot(width=.5, size=0.8, color= c('#D5766D','#84A677','#356A8C')) +  
  geom_jitter(aes(color=time), size=1.8, pch=1, width = 0.1, stroke=1.2) +  
  geom_line(aes(group=Sample_ID, color=Response), size = 0.3, alpha = 1, lty=2) + 
  scale_color_manual(values = c('#D5766D','#84A677','#356A8C', '#E68460', '#8DCDD5')) + 
  xlab('') + ylab('HOMA-IR') +   # Glu(mg/dL)
  theme_classic() + 
  theme( aspect.ratio = 5/4, 
         axis.title = element_text(size=12), 
         axis.text = element_text(size = 12),
         legend.position = 'none') + 
  stat_compare_means(label = "p.signif", method = "t.test", ref.group = "T1") 
#+facet_wrap(~Response, scales = 'free')

df['relative_change'] = 100*(df$T3-df$T1)/df$T1
ggplot(df, aes(x=Response, y=relative_change, color=Response))+
  geom_boxplot(width=0.5, size=0.8) +
  geom_jitter(aes(color=Response), size=2, pch=1, width = 0.2, stroke=1.2) +  
  scale_color_manual(values = c('#E68460', '#8DCDD5')) +
  theme_classic() + ylab('ΔHOMA-IR') + xlab('') + 
  ylim(-110,85) +
  theme( aspect.ratio = 2/1, 
         axis.title = element_text(size=12), 
         axis.text = element_text(size = 12),
         legend.position = 'none') + 
  stat_compare_means(label = 'p.signif', method = 't.test', ref.group = 'NR')

# 2. sensitive-index
si = f[f$Group!='VEH', c('Group', 'Effect_this', 'SI')] 
ggplot(si, aes(1:nrow(si), SI, color=Effect_this)) +
  geom_jitter(size=3, alpha=1) + 
  scale_color_manual(values =c('#E68460', '#8DCDD5')) + 
  theme_classic() + 
  scale_x_discrete(breaks = NULL) + 
  theme(axis.text.y = element_text(size=12), 
        axis.title = element_text(size=12)) +
  xlab('') + ylab('Sensitive Index') +
  theme(legend.position = 'top') + 
  geom_hline(yintercept = 2,lty=2) 

ggplot(si, aes(1:nrow(si), SI, color=Group)) +
  geom_jitter(size=3, alpha=1) + 
  scale_color_manual(values =c('#969696', '#6498b9', '#e3a268', '#8c705b')) + 
  theme_classic() + 
  scale_x_discrete(breaks = NULL) + 
  theme(axis.text.y = element_text(size=12), 
        axis.title = element_text(size=12)) +
  xlab('') + ylab('Sensitive Index') +
  theme(legend.position = 'top') + 
  geom_hline(yintercept = 2, lty=2, color='red') 

#-------------------------------------------- Glu,INs,HbA1c,bw,TEI-Rvs.NR----------------------
##-------------------------------------------- boxplot + relative_change_test ----------------
# 3. GLU,Ins,HbA1c,bw,TEi
#-------------------------------------------- bw,TEI ---------------------
fn = 'BW' # 'TEI'
f = read.csv(paste0(fn,'.csv'), row.names = 1)
f['Response'] = c(rnr[rownames(f)[1:30], 'Effect_this'], rep('VEH', 10)); colnames(f)
f$Response = factor(f$Response, levels = c('VEH', 'NR','R'))
f = f[f$Group!='VEH', c(1,ncol(f),4:21)]  # Day-7 -> Day57; TEI -> f = f[f$Group!='VEH', c(1,ncol(f), 22:79)]
table(f$Response)
colnames(f)[3:ncol(f)] = paste0('d', c(0, str_extract(colnames(f)[4:ncol(f)], '[0-9]+')))
f = f[-which(rownames(f) %in% test_id),]

library(ggalt)
dat = melt(f, id.vars = c('Sample_ID', 'Response'), variable.name = 'time', value.name = 'bw')
ggplot(dat, aes(x=time, y=bw, group=Sample_ID)) +
  geom_boxplot(aes(group=time, color=Response), alpha=1, show.legend = F) +
  geom_jitter(aes(col=Response), show.legend = F, position = position_jitter(0.2), size=1.2, alpha=0.3)+   
  #geom_xspline(aes(y=bw, col=Response), show.legend = F, size=0.8, alpha=0.5) +  
  #geom_line(aes(y=bw, col=Response), show.legend = F, size=0.8, alpha=0.5) + 
  stat_summary(fun = 'median', aes(group=Response, color=Response), geom='line',
               color=rep(c('red','blue'), each=18)) + 
  stat_summary(fun = 'median', aes(group=Response), geom='point', 
               pch=18, cex=3, color=rep(c('red','blue'), each=18)) + 
  theme_bw() + xlab('') + ylab('BW (kg)') +
  scale_color_manual(values = c( '#E68460', '#8DCDD5')) +  #'#9891A3',
  facet_wrap(~Response, ncol = 1, scales = 'free_y')

f['relative_change'] = 100*(f$d57-f$d0)/f$d0
ggplot(f, aes(x=Response, y=relative_change, color=Response))+
  geom_boxplot(width=0.5, size=0.8) +
  geom_jitter(aes(color=Response), size=2, pch=1, width = 0.2, stroke=1.2) +  
  scale_color_manual(values = c('#E68460', '#8DCDD5')) +
  theme_classic() + ylab('ΔBW (%)') + xlab('') + 
  ylim(-25, 5) +
  theme( aspect.ratio = 2/1, 
         axis.title = element_text(size=12), 
         axis.text = element_text(size = 12),
         legend.position = 'none') + 
  stat_compare_means(label = 'p.signif', method = 'wilcox.test', ref.group = 'NR') #t.test

df = f[, c(c('Sample_ID', 'Response', 'd0', 'd29', 'd57'))]
colnames(f)[3:5] = c('T1', 'T2', 'T3')
dat = melt(df, id.vars = c('Sample_ID', 'Response'), variable.name = 'time', value.name = 'bw')
