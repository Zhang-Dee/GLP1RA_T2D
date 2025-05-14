######################################### ANOVA for pheno data difference #############
setwd('D:/GLP1RA/')
rnr = read.csv('response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR based on homa-IR improvement
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

library(reshape2)
library(ggplot2)
f = read.csv('GLU.csv') # HbA1c, Ins
f = f[-which(f$Animal_ID %in% test_id), ] 
df = f[c("Animal_ID", "Group", "Day_.7", "Day_29", "Day_57")]
colnames(df)[3:5] = c('T1','T2','T3')

dat = melt(df, id.vars = c('Animal_ID', 'Group'), value.name = 'FBG', variable.name = 'time')
dat$Animal_ID = factor(dat$Animal_ID)
dat$Group = factor(dat$Group, levels=c('VEH','CpdA', 'CpdC', 'SAR17')) #  c('#9891A3', '#E4937C', '#E2B6AE', '#8BB39A')
dat$time = factor(dat$time)

# repeat measured two way ANOVA
library(rstatix)
f2 = anova_test(data = dat, dv = FBG, wid = Animal_ID, within = time, between = Group) 
get_anova_table(f2)
write.table(get_anova_table(f2), 'CpdABC_compare/two_way_repeated_anova_glu.txt', row.names = F)

## group-time：Simple main effect + Simple pairwise comparisons
"""
mainEff = aov(FBG ~ Group, data=dat)  
summary(mainEff)
TukeyHSD(mainEff, p.adjust.methods="bonferroni") 
"""
one.way = group_by(dat, time) %>% anova_test(dv = FBG, wid = Animal_ID, between = Group) %>% get_anova_table() 
pwc <- group_by(dat, time) %>% pairwise_t_test(FBG ~ Group) # Pairwise comparisons between group levels
## it can be seen that the simple main effect of group was significant at t3 (p = 0.032) but not at t1 (p = 0.672) or t2 (0.188)
## Pairwise comparisons show that the mean FBG was significantly different in grpVEH vs grpCpdA comparison at t3 (p = 0.027)..
write.table(one.way, 'CpdABC_compare/anova_main_effect_glu.txt', row.names = F)
write.table(pwc, 'CpdABC_compare/anova_pairwise_glu.txt', row.names = F)
#write.table(pwc, 'CpdABC_compare/except_test_sample/anova_pairwise_glu.txt', row.names = F)

# Visualization: boxplots with p-values
ggplot(dat, aes(x=time, y=FBG)) +
  geom_boxplot(aes(color=Group), width=0.6, position = position_dodge(0.8),
               outliers = FALSE, lwd=0.8) +
  scale_color_manual(values = c('#9891A3', '#E4937C', '#E2B6AE', '#8BB39A'))+ 
  geom_point(aes(color=Group), position = position_dodge(0.8), alpha=0.5, size=1.2) +
  ylim(20, 200) + ylab('FBG (mg/dL)') + xlab(NULL) +     # HbA1c (%) 3-8   INS (µU/mL) -30-650  HOMA-IR -20-170
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 10),
        legend.direction = 'horizontal',
        legend.position = 'top',
        legend.text = element_text(size = 12,margin = margin(r=8)),
        axis.title = element_text(size = 12))  



######################################### ANOSIM for microbe data difference #############
## ADONIS 
## ANOSIM 

# step1. data prepare
f = read.csv('raw_data/count/species/Unigenes.relative.s.xls', row.names = 1) # psudo path
sort(colnames(f)[grepl('H[0-9]+.T1', colnames(f))])
t1 = f[c(paste0('F',1:10,'.T1'), paste0('G',1:10,'.T1'), 
         paste0('H',1:10,'.T1'), paste0('E',c(1,2,3,5,6,7),'.T1'), 'class')]  

t2 = f[c(paste0('F',1:10,'.T2'), paste0('G',c(1:10)[-7],'.T2'),     
         paste0('H',1:10,'.T2'), paste0('E',1:7,'.T2'), 'class')]           

t3 = f[c(paste0('F',1:10,'.T3'), paste0('G',c(1:10),'.T3'),     
         paste0('H',1:10,'.T3'), paste0('E',c(4,7),'.T2'), 'class')]           

## remove Others -> save >=80% -> >=1e-6
i = 1
for (t in list(t1,t2,t3)) {
  cl = t['class']
  t = t[-nrow(t), -ncol(t)] 
  t = t[rowSums(t>0) > (0.8*ncol(t)), ]
  
  take = c()                            
  for (g in c('F','G','H','E')){   
    sub_t = t[,grepl(g, colnames(t))]
    sub_t = rownames(sub_t)[rowMeans(sub_t) > 1e-6]
    take = union(take, sub_t)
  }
  
  t = t[take, ]          
  cl = cl[take, 1]
  t['annotation'] = cl
  write.csv(t, paste0('./species_data/relative_abundance_T', i, '.csv'))
  i = i+1
} # t1:4843; t2:5011; t3:5315


# step2. Multivariate nonparametric test
library(vegan)
library(stringr)
f = read.csv('./species_data/relative_abundance_T1.csv', row.names = 1)
df = f[grepl('[FGH]', colnames(f))]  
df = data.frame(t(df), row.names = str_sub(colnames(df), 1, -4))  
df = df[ -which(rownames(df) %in% rnr[rnr$type=='test', 'Sample_ID']), ] 

group = ifelse(grepl('F', rownames(df)), 'CpdA', 
               ifelse(grepl('G', rownames(df)), 'CpdC', 'SAR17'))
group = data.frame(Group=group, row.names = rownames(df))

## ADONIS   
bray = adonis2(df ~ Group, data = group, permutations = 999, method="bray") # 用group分组df样本，并进行Adonis分析
bray 

## ANOSIM    
anosim_bray <- anosim(df, group$Group, permutations = 999, distance = "bray")
plot(anosim_bray, col =c('gray',  '#E4937C', '#E2B6AE', '#8BB39A'), 
     boxwex=1, xlab=NULL, ylab='Dissimilarity rank distribution', 
     title = 'ANOSIM (Distance function: Bray-Curtis)')


## PCov+ADOSIM
##PCoA
dist_bc = vegdist(df, 'bray',  na.rm = T)
pcoa = cmdscale(dist_bc, k=4, eig = T)
points = as.data.frame(pcoa$points)
eig = pcoa$eig
points = cbind(points, group)
colnames(points) = c('PC1', 'PC2', 'PC3', 'PC4', 'Group') 
points$Group = factor(points$Group)
ggplot(points, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point(aes(shape = Group), alpha=0.9, size=2.5) + 
  stat_ellipse(aes(fill = Group), linetype = 3, geom = 'polygon', alpha = 0.05, 
               level = 0.95, size = 0.1, show.legend = F) + 
  scale_color_manual(values = c('#E4937C', '#E2B6AE', '#8BB39A')) +
  scale_fill_manual(values = c('#E4937C', '#E2B6AE', '#8BB39A')) + 
  scale_shape_manual(values = c(15,16,17)) + 
  theme_classic() + 
  labs(x=paste('PCoA 1 (', format(100*eig[1]/sum(eig), digits = 4), '%)', sep = ''), 
                         y=paste('PCoA 2 (', format(100*eig[2]/sum(eig), digits = 4), '%)', sep = '')) + 
  theme(legend.title = element_text(size=12), legend.text = element_text(size = 12)) +
  labs(title = paste('ADOSIM: P=', signif(bray$`Pr(>F)`[1], 4)))  # adosim


