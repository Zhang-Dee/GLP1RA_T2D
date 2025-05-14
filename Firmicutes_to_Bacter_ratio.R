setwd('D:/GLP1RA/')
rnr = read.csv('response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR based on homa-IR improvement
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

<1>. Firmicutes/Bacteroidetes Ratio

library(stringr)
library(ggplot2)
library(reshape2)
library(ggpubr)
f = read.csv('spe_7level_T1.csv', row.names = 1)

f = f[rownames(f) %in% c("k__Bacteria;p__Firmicutes", "k__Bacteria;p__Bacteroidetes"), ]  
rownames(f) = str_sub(str_extract(rownames(f), ';[a-z]__.*'), 5) 
colnames(f) = str_sub(colnames(f), 1, -4)

f = f[!colnames(f) %in% test_id]  

r = f[colnames(f) %in% names(rnr_)[rnr_=='R']]
nr = f[colnames(f) %in% names(rnr_)[rnr_=='NR']]

rm_outlier = function(vec){
  lower_bound <- quantile(vec, 0.25) - 1.5 * IQR(vec)
  upper_bound <- quantile(vec, 0.75) + 1.5 * IQR(vec)
  filtered_data <- vec[vec >= lower_bound & vec <= upper_bound]
  if (length(filtered_data)==length(vec)){print('0 be filtered.')} else {print(which(!vec %in% filtered_data )); print('st be filtered')} 
  filtered_data
}
a = rm_outlier(as.numeric(r[1,]/r[2,]))
b = rm_outlier(as.numeric(nr[1,]/nr[2,]))  
wilcox.test(a,b)


dat = data.frame(t(cbind(r[1,]/r[2,], nr[1,]/nr[2,])), 
                 Response = c(rep('R', ncol(r)), rep('NR', ncol(nr))))
dat = dat[rownames(dat)!='G9', ]    
ggplot(dat, aes(Response, Firmicutes, color=Response)) +
  geom_boxplot(width=0.6, size=0.9, alpha=0.5) +
  geom_jitter(size=1.3, pch=1, width = 0.2, stroke=1.5, alpha=0.6) +  
  ylab('Firmicutes/Bacteroidets ratio') + xlab('') +
  theme_classic() +
  scale_color_manual(values = c('#E68460', '#8DCDD5')) +
  theme(aspect.ratio = 3/2, 
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        legend.position = 'none') +
  stat_compare_means(aes(group=Response)) 

dat_T1 = dat


f = read.csv('spe_7level_T3.csv', row.names = 1)


sid = intersect(rownames(dat_t1), rownames(dat))
dat_new = cbind(dat_T1[sid,], dat[sid, ])[-2]
colnames(dat_new)[1:2] = c('Before Treatment', 'After Treatment')

df = melt(dat_new, 'Response')
df$Response = factor(df$Response); df$variable = factor(df$variable)
ggplot(df, aes(variable, value, fill=Response)) +
  geom_boxplot(width=0.6, alpha=0.8) + facet_grid(.~Response)+
  theme_bw()+
  scale_fill_manual(values = c('#E68460', '#8DCDD5')) + 
  scale_y_continuous(name = 'Firmicutes/Bacteroidets ratio')+
  theme(strip.text.x = element_text(size = 18, colour = "black")) +#分面字体和背景
  theme(axis.title.y = element_text(size=20), axis.title.x = element_text(size=0),
        axis.text = element_text(size=10),
        legend.position = 'none',
        axis.line = element_line(color = 'black'), panel.background = element_blank()) + 
  stat_compare_means(aes(group=variable))

</1>. Firmicutes/Bacteroidetes Ratio



<2> DEG- heatmap
library(pheatmap)
f = read.csv('spe_7level_T1.csv', row.names = 1, encoding = 'utf-8')
f = f[grepl('[FGH]', colnames(f))]
colnames(f) = str_sub(colnames(f), 1, -4)
df = data.frame(t(f[!colnames(f) %in% test_id]))  # 只取train set
colnames(df) = rownames(f)
dim(df)

deg = read.csv('../DEGanalysis_202406/DEspe_7level_T1.csv', row.names = 1)
table(deg$tax)
# deg = deg[deg$tax=='s', ]
df = df[rownames(deg)[abs(deg$log2FC) > 1]]

df['Effect'] = rnr_[rownames(df)]
table(df$Effect)

dat_R = df[which(df$Effect == 'R'),]
dat_NR = df[which(df$Effect == 'NR'),]
dat = rbind(dat_R, dat_NR)

annotation_col = data.frame(SG = dat$Effect, row.names=rownames(dat))
annotation_col$SG = factor(annotation_col$SG, levels = c('R', 'NR'))
annotation_color = list(SG = c('R' = '#8DCDD5', 'NR' = '#E6846D'))  
dat = t(dat[-ncol(dat)])

dat_ = dat
for (i in 1:nrow(dat_)) {
  dat_[i,] = (dat_[i,]-min(dat_[i,]))/(max(dat_[i,])-min(dat_[i,]))
}


pheatmap(dat_, scale = "none", clustering_method = "complete", 
         clustering_distance_rows = "correlation", 
         border = T, # border_color='darkgrey',
         color = colorRampPalette(c('#FDF5E6', '#FF4500'))(100), 
         show_rownames = F, show_colnames = T, cluster_cols = F, treeheight_row = 70, 
         annotation_col = annotation_col, annotation_legend = T, 
         annotation_colors = annotation_color, legend_labels = NA,
         annotation_names_col = FALSE)
