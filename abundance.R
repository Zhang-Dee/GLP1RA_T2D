<1>.alpha-abundance 统计箱线图
#组间比较（先正态性检验,方差齐性检验后）
setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院/RNR_compare/abundance_compare/')
alpha_diversity = read.csv('alpha_diversity_T1.csv', row.names = 1)  # T2-T3同
test_id = c("F3",  "F5",  "F6",  "G2",  "G3" , "H5",  "H8",  "F9" , "H2",  "H10")
alpha_diversity = alpha_diversity[!rownames(alpha_diversity) %in% test_id, ] # train set
alpha_diversity$effect = factor(alpha_diversity$effect)

r <- alpha_diversity[alpha_diversity$effect == 'R', 1:6]
nr <- alpha_diversity[alpha_diversity$effect == 'NR', 1:6]
pr <- apply(r, MARGIN = 2, shapiro.test); pr <- as.numeric(sapply(pr, '[', 2)) #也可以直接unlist(pr)转化为向量
pnr <-apply(nr, MARGIN = 2, shapiro.test); pnr <- as.numeric(sapply(pnr, '[', 2))
pr;pnr  # 大都p>0.05
#shapiro.test检验正态性，p>0.05认为正态分布
#bartlett.test检验方差齐性，p>0.05认为齐
var_equal <- function(x) {
  bartlett.test(x = x, g = alpha_diversity$effect)$p.value
}
apply(alpha_diversity[ ,1:6], MARGIN = 2, var_equal) #不算coverage
#正态且方差齐：t检验
ttest <- function(x){
  t.test(x[alpha_diversity$effect == 'R'], x[alpha_diversity$effect == 'NR'], paired = FALSE, var.equal = TRUE)$p.value
}
apply(alpha_diversity[,1:6], MARGIN = 2, ttest)  #都不显著, beta多样性也在T1-3都不显著

#可视化
library(ggplot2)
library(ggpubr)
library(reshape2)
dat = alpha_diversity[grepl('[FGH]', rownames(alpha_diversity)), -7]
dat[1:6] = apply(dat[1:6], 2, scale) 
dat <- melt(dat, id.vars = 'effect')
dat$variable = apply(dat['variable'], 1, function(x){paste0(toupper(substr(x,1,1)), str_sub(x,2))}) # 首字母大写
dat$variable = factor(dat$variable, levels=c('Shannon', 'Chao1', 'Richness', 'Ace', 'Simpson', 'Pielou'))

ggplot(dat, aes(x=variable, y=value, color=effect))+
  geom_boxplot(width=0.6, size=0.8) +    #, outlier.size = 0.3
  scale_color_manual(values = c('#E68460', '#8DCDD5')) +
  theme_bw() + labs(x='',y='Scaled α-Diversity') +
  theme( #aspect.ratio = 2/1, 
    axis.title = element_text(size=12), 
    axis.text = element_text(size = 12),
    legend.position = 'right') + 
  stat_compare_means(label = 'p.signif', method = 't.test') #t.test

## T2-T3时刻的多样性需要重复测量方差分析
library(rstatix)
ab_t1 = read.csv('alpha_diversity_T1.csv', row.names = 1)
ab_t2 = read.csv('alpha_diversity_T2.csv', row.names = 1)
ab_t3 = read.csv('alpha_diversity_T3.csv', row.names = 1)
co_s = intersect(intersect(rownames(ab_t1), rownames(ab_t2)), rownames(ab_t3))
co_s = co_s[grepl('[FGH]', co_s)]
co_s = co_s[!co_s %in% test_id]

anova_result = data.frame(row.names = c('p_response', 'p_time', 'p_interact'))
for (ab in 1:6) {
  ab = colnames(ab_t1[ab])
  dat = cbind(ab_t1[co_s, c('effect', ab)], ab_t2[co_s, ab], ab_t3[co_s, ab])
  dat['Animal_ID'] = co_s
  colnames(dat)[1:4] = c('Response', 'T1','T2','T3')
  
  dat = melt(dat, id.vars = c('Animal_ID', 'Response'), value.name = 'FBG', variable.name = 'time')
  dat$Animal_ID = factor(dat$Animal_ID)
  dat$Response = factor(dat$Response, levels = c('NR','R')); table(dat$Response)
  dat$time = factor(dat$time)
  
  aov = anova_test(data = dat, dv = FBG, wid = Animal_ID, within = time, between = Response) 
  anova_result[ab] = get_anova_table(aov)$p # 交互效应如果显著做后面simple main effect 不显著算到这就可以
}
anova_result
write.csv(anova_result, 'two_way_repeat_anova_T3_abundance.csv')

</1> alpha-abnundace 可视化


<2> Shannon单独展示
...从- dat = alpha_diversity[grepl('[FGH]', rownames(alpha_diversity)), -7]
dat = dat[c('shannon', 'effect')]
dat = melt(dat, id.vars = 'effect')
ggplot(dat, aes(effect, value, color = effect)) + 
  geom_boxplot(width=0.6, size=1) +
  geom_jitter(size=3, width = 0.25, alpha=0.5) +
  scale_color_manual(values = c('#E67460','#8DCDD5')) +
  theme_classic() + ylab('Shannon Index') + xlab('') +
  theme(legend.position = 'none', 
        aspect.ratio = 2/1) + 
  stat_compare_means(label='p.signif', method='t.test')
  
  
</2> Shannon单独展示





<3> DEG- heatmap
setwd("D:/Zhang Di/linshibangong/linshiMonkey/杭高院/direct_data")
library(pheatmap)
f = read.csv('spe_7level_T1.csv', row.names = 1, encoding = 'utf-8')
f = f[grepl('[FGH]', colnames(f))]
colnames(f) = str_sub(colnames(f), 1, -4)
df = data.frame(t(f[!colnames(f) %in% test_id]))  # 只取train set
colnames(df) = rownames(f)
dim(df)

deg = read.csv('../DEGanalysis_202406/DEspe_7level_T1.csv', row.names = 1)
table(deg$tax)
# deg = deg[deg$tax=='s', ]; p = deg$label %>% str_extract('p__[\\S ]+;c__') %>% str_sub(1, -5); table(p) %>% sort()
# deg = deg[(deg$meanNR > 1e-5) & (deg$meanR > 1e-5),]
df = df[rownames(deg)[abs(deg$log2FC) > 1]]

df['Effect'] = rnr_[rownames(df)]
table(df$Effect)

dat_R = df[which(df$Effect == 'R'),]
dat_NR = df[which(df$Effect == 'NR'),]
dat = rbind(dat_R, dat_NR)

annotation_col = data.frame(SG = dat$Effect, row.names=rownames(dat))  #annotation_row = data.frame(SG = p, row.names=rownames(deg))
annotation_col$SG = factor(annotation_col$SG, levels = c('R', 'NR'))   # annotation_row$SG = factor(annotation_row$SG, levels = names(sort(table(p))))
annotation_color = list(SG = c('R' = '#8DCDD5', 'NR' = '#E6846D'))  
dat = t(dat[-ncol(dat)])

dat_ = dat
rownames(dat_) = str_extract(rownames(dat_), 's__[\\S ]+') %>% str_sub(4)
for (i in 1:nrow(dat_)) {
  dat_[i,] = (dat_[i,]-min(dat_[i,]))/(max(dat_[i,])-min(dat_[i,]))
}


pheatmap(dat_, scale = "none", clustering_method = "complete", 
         clustering_distance_rows = "correlation", 
         border = T, # border_color='darkgrey',
         color = colorRampPalette(c('#FDF5E6', '#FF4500'))(100), 
         show_rownames = T, show_colnames = F, cluster_cols = F, treeheight_row = 70, # show_rownames = F
         annotation_col = annotation_col, annotation_legend = T, 
         annotation_colors = annotation_color, legend_labels = NA,
         annotation_names_col = FALSE)

</3> DEG- heatmap


<4> T1vsT3
ggplot(dat, aes(x=effect, y=value, fill = effect, alpha=time))+
  geom_boxplot(width=0.6, size=0.5, show.legend = F, outliers = F) +   
  scale_fill_manual(values = c('#DC8070', '#6188B5')) + 
  geom_point( position = position_dodge(0.6), size=1, show.legend = F) +
  theme_bw() + labs(x='',y='') +
  theme( aspect.ratio = 3/2, 
         axis.title = element_text(size=12), 
         axis.text = element_text(size = 12),
         legend.position = 'right') + 
  stat_compare_means(label = 'p.signif', method = 't.test', paired = T) + facet_wrap(~variable, nrow = 1, scales = 'free')
</4> T1vsT3
