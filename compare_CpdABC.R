######################################### ANOVA for pheno data difference #############
setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Animal_ID']
library(reshape2)
library(ggplot2)
f = read.csv('GLU.csv') # HbA1c, Ins
f = f[-which(f$Animal_ID %in% test_id), ] # ./CpdABC_compare/except_test_sample/
df = f[c("Animal_ID", "Group", "Day_.7", "Day_29", "Day_57")]
colnames(df)[3:5] = c('T1','T2','T3')

dat = melt(df, id.vars = c('Animal_ID', 'Group'), value.name = 'FBG', variable.name = 'time')
dat$Animal_ID = factor(dat$Animal_ID)
dat$Group = factor(dat$Group, levels=c('VEH','CpdA', 'CpdC', 'SAR17')) #  c('#9891A3', '#E4937C', '#E2B6AE', '#8BB39A')
dat$time = factor(dat$time)

# 重复测量两因素多水平方差分析
"""
f1 = aov(FBG ~ time * Group + Error(Animal_ID/(time)), data = dat)
summary(f1) # 2table, 1st 组件比较方差分析表,2nd处理前后比较和交互作用

with(dat,
     interaction.plot(time, Group, FBG, type = "b", col = c("red","blue"), 
                      pch = c(12,16), main = "两因素两水平重复测量方差分析"))
boxplot(FBG ~ Group*time, data = dat, col = c("gold","green"),
        main = "两因素两水平重复测量方差分析")
"""

library(rstatix)
f2 = anova_test(data = dat, dv = FBG, wid = Animal_ID, within = time, between = Group) 
get_anova_table(f2)
write.table(get_anova_table(f2), 'CpdABC_compare/two_way_repeated_anova_glu.txt', row.names = F)

## group-time交互作用显著时需要做事后检验：Simple main effect + Simple pairwise comparisons
"""
mainEff = aov(FBG ~ Group, data=dat)  # 将group不同处理作为主效应,检验不同组间差异
summary(mainEff)
TukeyHSD(mainEff, p.adjust.methods="bonferroni") #事后两两比较
"""
one.way = group_by(dat, time) %>% anova_test(dv = FBG, wid = Animal_ID, between = Group) %>% get_anova_table() # main effect的单因素
pwc <- group_by(dat, time) %>% pairwise_t_test(FBG ~ Group) # Pairwise comparisons between group levels
## it can be seen that the simple main effect of group was significant at t3 (p = 0.032) but not at t1 (p = 0.672) or t2 (0.188)
## Pairwise comparisons show that the mean FBG was significantly different in grpVEH vs grpCpdA comparison at t3 (p = 0.027)..
write.table(one.way, 'CpdABC_compare/anova_main_effect_glu.txt', row.names = F)
write.table(pwc, 'CpdABC_compare/anova_pairwise_glu.txt', row.names = F)
#write.table(pwc, 'CpdABC_compare/except_test_sample/anova_pairwise_glu.txt', row.names = F)

# 可视化  Visualization: boxplots with p-values
ggplot(dat, aes(x=time, y=FBG)) +
  geom_boxplot(aes(color=Group), width=0.6, position = position_dodge(0.8),
               outliers = FALSE, lwd=0.8) +
  scale_color_manual(values = c('#9891A3', '#E4937C', '#E2B6AE', '#8BB39A'))+ 
  geom_point(aes(color=Group), position = position_dodge(0.8), alpha=0.5, size=1.2) +
  ylim(20, 200) + ylab('FBG (mg/dL)') + xlab(NULL) +     # HbA1c (%) 3-8   INS (µU/mL) -30-650  HOMA-IR -20-170
  theme_classic() +
  theme(axis.text = element_text(face = "bold",size = 10), #坐标轴刻度标签加粗
        legend.direction = 'horizontal',
        legend.position = 'top',
        legend.text = element_text(size = 12,margin = margin(r=8)),
        axis.title = element_text(size = 12))  
# 输出pdf后用AI加pval, 颜色调整


######################################### ANOSIM for microbe data difference #############
## ADONIS 多元非参数方差分析, 通过线性模型分析不同分组因素或环境因子(如土壤理化性质等)对样品差异的解释度，并使用置换检验进行显著性分析。
## ANOSIM 组间差别分析, 用于比较组间的差异是否显著大于组内差异的非参数检验，从而判断分组是否有意义

# step1. 数据准备
f = read.csv('raw_data/count/species/Unigenes.relative.s.xls', row.names = 1) # psudo path
sort(colnames(f)[grepl('H[0-9]+.T1', colnames(f))])
t1 = f[c(paste0('F',1:10,'.T1'), paste0('G',1:10,'.T1'), 
         paste0('H',1:10,'.T1'), paste0('E',c(1,2,3,5,6,7),'.T1'), 'class')]  # control在T1的宏基因组数据有6个样本，少4个

t2 = f[c(paste0('F',1:10,'.T2'), paste0('G',c(1:10)[-7],'.T2'),      # G7在T2缺失
         paste0('H',1:10,'.T2'), paste0('E',1:7,'.T2'), 'class')]            # H对照组在T2有7个样本

t3 = f[c(paste0('F',1:10,'.T3'), paste0('G',c(1:10),'.T3'),      # F9在T3有两个数据，不用F9.T3.2
         paste0('H',1:10,'.T3'), paste0('E',c(4,7),'.T2'), 'class')]            # control在T3只两个样本数据

## 清洗数据: 去除Others -> 保留80%样本中存在的物种 -> 在样本中均值高于1e-6的物种
i = 1
for (t in list(t1,t2,t3)) {
  cl = t['class']
  t = t[-nrow(t), -ncol(t)] # 最后一行是未识别的:Others; 最后一列是物种注释:type
  t = t[rowSums(t>0) > (0.8*ncol(t)), ] # 80%样本中都测到的保留
  
  take = c()                            # 在每组样本中检测到均值高于1e-6的物种保留
  for (g in c('F','G','H','E')){   
    sub_t = t[,grepl(g, colnames(t))]
    sub_t = rownames(sub_t)[rowMeans(sub_t) > 1e-6]
    take = union(take, sub_t)
  }
  
  t = t[take, ]             # 各组样本检测均值高于1e-6的物种保留
  cl = cl[take, 1]
  t['annotation'] = cl
  write.csv(t, paste0('./species_data/relative_abundance_T', i, '.csv'))
  i = i+1
} # t1:4843; t2:5011; t3:5315


# step2. 多元非参数检验
library(vegan)
library(stringr)
f = read.csv('./species_data/relative_abundance_T1.csv', row.names = 1)
df = f[grepl('[FGH]', colnames(f))]  # 只取3trt group   
df = data.frame(t(df), row.names = str_sub(colnames(df), 1, -4))   # 样本*物种
df = df[ -which(rownames(df) %in% rnr[rnr$type=='test', 'Sample_ID']), ]  # ./CpdABC_compare/except_test_sample/

group = ifelse(grepl('F', rownames(df)), 'CpdA', 
               ifelse(grepl('G', rownames(df)), 'CpdC', 'SAR17'))
group = data.frame(Group=group, row.names = rownames(df))

## ADONIS    置换多元方差分析或非参数多元方差分析，利用各种组间距离指数对总方差进行分解，可以分析不同分类因子对群落差异的解释度，并使用置换检验进行统计学检验。
bray = adonis2(df ~ Group, data = group, permutations = 999, method="bray") # 用group分组df样本，并进行Adonis分析
bray # 只能看分类因子整体上是否有差异
      # p>0.05表明基于给药分组对群落结构差异不具统计学意义


## ANOSIM    检验不同群落相似性,判断分组是否有意义
anosim_bray <- anosim(df, group$Group, permutations = 999, distance = "bray")
            ## R值表示群落相似性，R值越接近1表示组间差异大于组内差异，值越大表示组间差异越大。
plot(anosim_bray, col =c('gray',  '#E4937C', '#E2B6AE', '#8BB39A'), 
     boxwex=1, xlab=NULL, ylab='Dissimilarity rank distribution', 
     title = 'ANOSIM (Distance function: Bray-Curtis)')


## PCov可视化+ADOSIM结果
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
  labs(title = paste('ADOSIM: P=', signif(bray$`Pr(>F)`[1], 4)))  # 加adosim结果


