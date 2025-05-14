setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 
#test_id = c("F3",  "F5",  "F6",  "G2",  "G3",  "H5",  "H8",  "F9",  "H2",  "H10")

<1.DEG-RF feature可视化>  feature的箱线图
setwd("D:/Zhang Di/linshibangong/linshiMonkey/杭高院/DEGanalysis_202406")
library(reshape2)
library(ggplot2)
library(ggpubr)
# gene = c("K03471", "K09985", "K16788", "K13479", "K00226", "K01436", "K01457", "K10122", "K16875", "K13059") 这是直接根据DEG重要性前十，最终不用
gene = c("K16788", "K01457", "K10122", "K13677", "K01692", "K16958", "K10121", "K05878", "K18303", "K10120", "K03562", "K01941") # 用T2D通路相关12gene
f = read.csv('../direct_data/gene_4level_T1.csv', row.names = 1)
f = f[gene, grepl('[FGH]', colnames(f))]  # DEG的样本特征矩阵, feature*sample
colnames(f) = str_sub(colnames(f), 1, -4)
f = data.frame(t(f[gene, !colnames(f) %in% test_id]))
f['Effect'] = factor(rnr_[rownames(f)])
deg = read.csv("DEgene_4level_T1.csv")


dat = melt(f, id.vars = 'Effect')
dat$variable = factor(dat$variable, levels =  gene[order(deg[match(gene, deg$name), 'log2FC'])]) # 按FC顺序作图
dat$value = -log10(dat$value + 1e-8)
ggplot(data=dat, aes(x=variable, y=value, fill=Effect))+
  geom_boxplot(width=0.8, outlier.size = 0.3, alpha=0.8)+
  scale_fill_manual(values = c('#E6846D', '#8DCDD5')) +
  stat_compare_means(aes(group=Effect), method = 'wilcox.test', label="p.signif", hide.ns = F)+
  labs(x='',y='Relative Abundance (-Log10)')+
  theme_classic() + theme(axis.text.y = element_text(size=10), 
                          axis.text.x = element_text(size=10, angle = 45, hjust = 1),
                          axis.title = element_text(size=12),
                          legend.title = element_blank(),
                          plot.margin = unit(c(1,1,1,1),'cm'),
                          legend.text = element_text(size=10))

</1.DEG-RF feature可视化>  feature的箱线图



<2.DEspe-RF feature可视化>  feature的箱线图
setwd("D:/Zhang Di/linshibangong/linshiMonkey/杭高院/DEGanalysis_202406")
library(reshape2)
library(ggplot2)
library(ggpubr)
primary_f = c('g__Prevotella;s__Prevotella sp. P3-120',
              'g__Bacteroides;s__Bacteroides sp. 3_1_40A',
              'g__Roseburia;s__Roseburia intestinalis CAG:13', 
              'g__Eubacterium;s__Eubacterium sp. CAG:86',
              'g__Eubacterium;s__Eubacterium sp. 36_13',
              'g__Clostridium;s__Clostridium sp. 7_2_43FAA',
              'g__Clostridium;s__Clostridium sp. HMSC19A10',
              'g__Clostridium;s__Clostridium disporicum',
              'g__Clostridium;s__Clostridium sp. 29_15',
              'g__Clostridium;s__Clostridium sp. C8',
              'f__Acidaminococcaceae;g__Phascolarctobacterium', 
              'f__Peptoniphilaceae;g__Finegoldia',
              'f__Sphingobacteriaceae;g__Solitalea')


f = read.csv('../direct_data/spe_7level_T1.csv', row.names = 1)
f = f[deg$name, grepl('[FGH]', colnames(f))]  # DEG的样本特征矩阵, feature*sample
colnames(f) = str_sub(colnames(f), 1, -4)
df = data.frame(t(f[primary_f, !colnames(f) %in% test_id]))
colnames(df) = primary_f

df['Effect'] = factor(rnr_[rownames(df)])
table(df$Effect)
deg = read.csv("DEspe_7level_T1.csv")


dat = melt(df, id.vars = 'Effect')
dat$variable = factor(dat$variable, levels =  primary_f[order(deg[match(primary_f, deg$name), 'log2FC'])]) # 按FC顺序作图
dat$value = -log10(dat$value + 1e-8)
xlabel = matrix(unlist(strsplit(levels(dat$variable), ';')), nrow=2)[2,]
levels(dat$variable) = str_sub(xlabel, 4)
ggplot(data=dat, aes(x=variable, y=value, fill=Effect))+
  geom_boxplot(width=0.8, outlier.size = 0.3, alpha=0.8)+
  scale_fill_manual(values = c('#E6846D', '#8DCDD5')) +
  stat_compare_means(aes(group=Effect), method = 'wilcox.test', label="p.signif", hide.ns = F)+
  labs(x='',y='Relative Abundance (-Log10)')+
  theme_classic() + theme(axis.text.y = element_text(size=10), 
                          axis.text.x = element_text(size=10, angle = 45, hjust = 1),
                          axis.title = element_text(size=12),
                          legend.title = element_blank(),
                          plot.margin = unit(c(1,1,1,1),'cm'),
                          legend.text = element_text(size=10),
                          legend.position = 'leftbottom')

</2.DEspe-RF feature可视化>  feature的箱线图


<3> gene的importance条形图
# score = read.csv('randomForest_DEG/geneFC2Feature_importance_score.csv')
score = read.csv('randomForest_DEG/DEgeneFeature_T2Dpathway.csv', row.names = 1)[c(2,1)]
score = score[grepl('K',score$name),]
colnames(score) = c('name', 'MeanDecreaseGini')

ggplot(score, aes(x=MeanDecreaseGini, y=reorder(name, MeanDecreaseGini))) +
  geom_segment(aes(yend=name), xend=0, colour="grey50")+   #geom_segment()函数用"以数据点为端点的线段"代替网格线
  geom_point(size=1) +
  theme_classic() + 
  theme(plot.margin = unit(c(1,4,1,1), 'lines'),
        axis.text.y = element_text(size=6), 
        axis.title = element_text(size=12))+
  xlab('Importance score') + ylab('')
  
</3> gene的importance条形图


<4> gene-Edge的importance条形图
setwd("D:/Zhang Di/linshibangong/linshiMonkey/杭高院")
score = read.csv('new_ssn_202406/geneFeature_importance_score.csv', row.names = 1)
colnames(score) = c('name', 'MeanDecreaseGini')
... ggplot...

</4> gene-Edge的importance条形图

