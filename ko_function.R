setwd("D:/Zhang Di/linshibangong/linshiMonkey/杭高院/DEGanalysis_202406")
rnr = read.csv('../new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

<0>. 整理富集reference文件
#不用食蟹猴的一因为匹配不到二因为就不是根据食蟹猴来的ko. 
ko_brite = read.table('KEGG_brite_KO.txt', sep = '\t', header = T, quote = '')
ko_brite = ko_brite[!ko_brite$l3_name %in% c('Others', 'Function unknown'), ]
Nn = unique(unlist(strsplit(ko_brite$l4_ids,'/')))   # 背景基因集- 参考集中所有的基因
N = length(Nn)        # 背景基因集大小- 所有参与KEGG通路的gene数目
enrich = function(refRow){
  nn = unlist(strsplit(refRow[7], '/'))  #参与某通路的cpd
  n = length(nn)     # 特定通路基因集大小
  m = length(intersect(my_gene, nn)) # 交集基因集大小- 富集到某通路的差异基因数目
  pvalue = phyper(m-1, M, N-M, n, lower.tail = FALSE) 
  GeneRatio = paste(m,M,sep='/')
  BgRatio = paste(n,N,sep='/')
  fold = (m/M)/(n/N)
  Count = m
  r = c(GeneRatio,BgRatio,fold,Count,pvalue)
  return(r)
}
</0>
  
<1>. 富集
f = read.csv('../direct_data/gene_4level_T1.csv', row.names = 1)
my_gene = intersect(rownames(f)[f$label=='ko.xls'], Nn)  # 特定基因集- 兴趣基因子集(包含在背景基因集的)
M = length(my_gene) 

enrich_result = data.frame(t(apply(ko_brite, 1, enrich)))
colnames(enrich_result) = c('GeneRatio','BgRatio','fold','Count','pvalue')
enrich_result[, c('fold','Count','pvalue')] = apply(enrich_result[, c('fold','Count','pvalue')], 2, as.numeric)
enrich_result$p.adjust = p.adjust(enrich_result$pvalue, method = 'BH')

enrich_result = cbind(ko_brite[c('l1_name', 'l2_name', 'l3_name')], enrich_result)
enrich_result = enrich_result[enrich_result$pvalue < 0.05, ]
enrich_result = enrich_result[order(enrich_result$pvalue), ]
enrich_result = enrich_result[enrich_result$l1_name!='Not Included in Pathway or Brite', ]
write.csv(enrich_result, 'enrichment/ALLgene_enrichment_microbe.csv', row.names = F)
</1>


<2> 可视化 所有测到的基因富集的桑吉图
a = data.frame(table(enrich_result$l2_name))  # 按level3通路条目计数
colnames(a) = c('level2', 'count'); a$level2=as.character(a$level2)
aa = na.omit(unique(enrich_result[1:2])); rownames(aa)=aa$l2_name
a['level1'] = aa[a$level2, 'l1_name']

library(RColorBrewer)
library(randomcoloR)
library(networkD3)
data = arrange(a[c('level1', 'level2', 'count')], desc(a$count))
nodes = data.frame(name = unique(c(data$level1, data$level2)), nodeGroup = paste0('c', 1:30))
data$IDsource = match(data$level1, nodes) - 1
data$IDtarget = match(data$level2, nodes) - 1
data$linkGroup = paste0('c',1:nrow(data))
## c( brewer.pal(6, 'Set1'), distinctColorPalette(nrow(data))) 生成颜色手动输
color <- 'd3.scaleOrdinal().range(["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
            "#DCA896", "#CEBBDB", "#97AE6F", "#DBE4A0", "#5EE1C9",  "#6E4FCD",
            "#DF98D0", "#82E849", "#CD77E0", "#797FE0", "#D8DC50", "#6E9ADB",
            "#6F6E6A", "#885084", "#B0E9D5", "#E6E3DA", "#D7944A", "#B843F0",
            "#80E991", "#83CAE1", "#DD43D4", "#DF6069", "#DE4A9D", "#6536E8"])'

sankeyNetwork(Links = data, Nodes = nodes, Source = 'IDsource', Target = 'IDtarget',
              Value = 'count', NodeGroup = 'nodeGroup', LinkGroup = 'linkGroup',
              sinksRight = F, fontSize = 11, nodeWidth = 30, nodePadding = 10,
              colourScale = color, height = 550, width = 550)
              # to.pdf: 202406fig/ko_function/ko_detect_all.pdf
## 如果要AI编辑需要输出html->pdf,字符很难改因此另外输出相同字符->
test = data
test$level2 = factor(test$level2, levels = test$level2[c(nrow(test):1)])
ggplot(test, aes(count, level2, fill=level1)) +
  geom_bar(stat = 'identity') +
  scale_y_discrete(position = 'right')

</2> 可视化

<3> 说明 + DEG根据富集参与跟T2D有关的通路的基因提取作为RF feature
如果分别做RNR->在1-富集处重新提取R/NR基因集
Rf = f[, names(rnr_)[rnr_=='R']] -> filter(gene)
但是没有必要,这里只需查看肠道菌群参与的功能,RNR差异会在DEG部分比较

ls = c('Transporters', 'Glycerolipid metabolism', 'Pyruvate metabolism', 'Arginine biosynthesis',
       'Fatty acid degradation')
a = ko_brite[ko_brite$l3_name %in% ls, 'l4_ids']
ref = c()
for (i in a){ref = unique(c(ref, unlist(strsplit(i, '/'))))}

f = read.csv('DEgene_4level_T1.csv', row.names = 1)
f = f[abs(f$log2FC)>1, ]
feature = intersect(rownames(f), ref)
write.csv(feature, 'randomForest_DEG/geneFeature_T2Dpathway.csv')
</3> 说明+ DEG根据富集参与跟T2D有关的通路的基因提取作为RF feature


<4>. PCOV图ko的
library(vegan)
library(stringr)
f = read.csv('../direct_data/gene_4level_T1.csv', row.names = 1)  # T3
f = f[f$label=='ko.xls', ]
df = f[grepl('[FGH]', colnames(f))]  # 只取3trt group   
df = data.frame(t(df), row.names = str_sub(colnames(df), 1, -4))   # 样本*基因
df = df[ !rownames(df) %in% test_id, ]  
dim(df)

## ADONIS
group = data.frame(Group=rnr_[rownames(df)], row.names = rownames(df))
bray = adonis2(df ~ Group, data = group, permutations = 999, method="bray") # 用group分组df样本，并进行Adonis分析
bray   # 只能看分类因子整体上是否有差异, p>0.05表明基于给药分组对群落结构差异不具统计学意义
## ANOSIM    检验不同群落相似性,判断分组是否有意义
anosim_bray <- anosim(df, group$Group, permutations = 999, distance = "bray")
## PCoA
dist_bc = vegdist(df, 'bray',  na.rm = T)
pcoa = cmdscale(dist_bc, k=4, eig = T)
points = as.data.frame(pcoa$points)
eig = pcoa$eig
points = cbind(points, group)
colnames(points) = c('PC1', 'PC2', 'PC3', 'PC4', 'Group') 
points$Group = factor(points$Group)

r <- points[points$Group == 'R', c(1,2,5)]
nr <- points[points$Group == 'NR', c(1,2,5)]
border1<- r[chull(r),]
border2<- nr[chull(nr),]
points$PC1.mean = structure( c(rep(mean(r$PC1), nrow(r)), 
                               rep(mean(nr$PC1), nrow(nr))),
                             names=c(rownames(r), rownames(nr)))[rownames(points)]
points$PC2.mean = structure( c(rep(mean(r$PC2), nrow(r)), 
                               rep(mean(nr$PC2), nrow(nr))),
                             names=c(rownames(r), rownames(nr)))[rownames(points)]
ggplot(points, aes(x=PC1, y=PC2, color=Group)) + 
  geom_point(aes(shape = Group), alpha=1, size=3) + 
  geom_polygon(data=border2, aes(color = Group), fill='#E6846D', alpha=0.1, show.legend = FALSE) +
  geom_polygon(data=border1, aes(color = Group), fill='#8DCDD5', alpha=0.1, show.legend = FALSE) + 
  geom_segment(aes(x = PC1.mean, y = PC2.mean, xend= PC1, yend = PC2, color = Group),
               alpha= 0.5, show.legend = FALSE, linetype='dashed' )+
  scale_color_manual(values = c( '#E6846D', '#8DCDD5')) +
  scale_fill_manual(values = c( '#E6846D', '#8DCDD5')) + 
  scale_shape_manual(values = c(19,17)) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) +
  labs(x=paste('PCoA 1 (', format(100*eig[1]/sum(eig), digits = 4), '%)', sep = ''), 
       y=paste('PCoA 2 (', format(100*eig[2]/sum(eig), digits = 4), '%)', sep = '')) +
  labs(title = paste('ANOSIM: P=', signif(anosim_bray$signif, 3))) 
</4>. PCOV图ko的


<5>. DEG heatmap
library(pheatmap)
deg = read.csv('DEgene_4level_T1.csv')
deg = deg[grepl('K', deg$name), ]
df = df[deg[abs(deg$log2FC) > 1, 'name']]
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

</5>. DEG heatmap


<6>. level3 pathway abundance 放弃！！！！！！！！！！
signif_level3 ='RIG-I-like receptor signaling pathway' # 'ko04622'
f = read.csv('../direct_data/gene_4level_T1.csv', check.names = F, row.names = 1)
f = f[f$label == 'level3', ]
label = f['tax']
f = f[grepl('[FGH]', colnames(f))]; colnames(f) = str_sub(colnames(f), 1, -4)
f = f[!colnames(f) %in% test_id]
f = f[rowMeans(f) > 1e-5, ]; label = data.frame(label[rownames(f), 'tax'], row.names = rownames(f))

r = f[colnames(f) %in% names(rnr_)[rnr_=='R']]
nr = f[colnames(f) %in% names(rnr_)[rnr_=='NR']] # 就6个统计差异显著
logFC = log2(rowMeans(nr) / rowMeans(r))
logFC = logFC[order(abs(logFC))]
label['log2FC'] = logFC[rownames(label)]
tail(label)
length(logFC[abs(logFC) > 0.585]) # 1.5倍FC

dat = label[abs(label$log2FC) > 0.585, ]
dat['name'] = unlist(strsplit(dat[,1], ';'))[seq(3, nrow(dat)*3, 3)]
dat = dat[- which(rownames(dat) %in% c('ko05145', 'ko05169', 'ko05143', 'ko05166')), -1]
xxx
代码没写完，没必要做这部分
dat_level3 = dat_level3[dat_level3$pvals < 0.05 & abs(dat_level3$log2FC) > 0.1, ]
column = c('level_3', 'level_1', 'log2FC')
dat_bar = dat_level3[,column];dat_bar=dat_bar[order(dat_bar$log2FC)[59:93],]
library(ggplot2)
dat_bar$level_1 = factor(dat_bar$level_1, levels = c('Metabolism', 'Environmental Information Processing', 'Cellular Processes', 'Organismal Systems', 'Human Diseases','Genetic Information Processing'))
dat_bar$level_3 = factor(dat_bar$level_3, levels = unique(dat_bar$level_3))
ggplot(dat_bar, aes(x = log2FC, y = level_3, color = level_1, fill = level_1)) + geom_bar(stat = 'identity', position = 'identity', alpha=0.7, width = 0.7) + 
  geom_vline(xintercept = c(0), linetype = 1, colour = "black", size = 1) +  
  theme(panel.grid=element_blank(), panel.background=element_rect(color='black', fill='transparent')) + 
  theme(axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 8), legend.text = element_text(size = 8, face='italic')) + 
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.title = element_text(size = 10, face='italic')) + theme(legend.position = c(0.2,0.2))+
  guides(col = guide_legend(nrow = 2))
</6>. level3 pathway abundance






#######################################################
# <废-桑吉图>
library(ggalluvial)
library(ggplot2)
library(randomcoloR)
library(RColorBrewer)
sank = arrange(a[c('level1', 'level2', 'count')], desc(a$count))
sank$level1 = factor(sank$level1, levels = names(sort(table(sank$level1), decreasing = T)))

# 绘制桑基图
color = c(brewer.pal(6,'Set1'), distinctColorPalette(nrow(sank))) # 为level2生成差异明显的颜色
ggplot(data, aes(axis1 = level1, axis2 = level2, value = count)) +
  geom_alluvium(aes(fill = level2), width = 2/5) + # , width = 2/5  控制色带
  geom_stratum(aes(fill = level1), width = 2/5, alpha=1) + # width = 2/5  控制柱子
  scale_fill_manual(values = c(color)) + 
  theme_minimal() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        legend.position = 'none') 


# </废-桑吉图>




# <废> 弦状图,因标签太长不好看放弃
a = data.frame(table(enrich_result$l2_name))  # 按level3通路条目计数
colnames(a) = c('level2', 'count'); a$level2=as.character(a$level2)
aa = na.omit(unique(enrich_result[1:2])); rownames(aa)=aa$l2_name
a['level1'] = aa[a$level2, 'l1_name']

library(circlize)
library(RColorBrewer)
sank = arrange(a[c('level1', 'level2', 'count')], desc(a$count))
level1 = sort(table(sank$level1), decreasing = T)
node = c(rep(names(level1), level1), sank$level2) # 按数目排序的

grid.col <- c( brewer.pal(length(unique(sank$level1)), 'Set1'), rainbow(nrow(sank)))
names(grid.col) = unique(node)
group <- structure(rep(c('level1', 'level2'), each=nrow(sank)), names = node)
#grid.col <- structure(grid.col[node])

par(oma=c(0.01,0.01,0.01,0.01), mar=c(1,1,1,1) + 0.1)
circos.par(start.degree = 250)
chordDiagram(sank, group = group, grid.col = grid.col, 
             annotationTrackHeight = c(0.05, 0.1), #外面一圈的宽度,改第一个值
             transparency = 0.5, # 越大越透明
             annotationTrack = "grid",
             directional = 1,  #表示线条的方向，0代表没有方向，1代表正向，-1代表反向，2代表双向
             big.gap = 20) # small.gap = 1,
#direction.type = c("diffHeight","arrows"), #线条是否带有箭头
#link.arr.type = "big.arrow",#另一个选择是巨丑无比的尖头
#diffHeight  = -0.04#外圈和中间连线的间隔
# order = #顺序
abline(a=1.2, b=15, lty = 2, col = "#fdb46280", lwd = 2) # a,b截距斜率;hv水平垂直线刻度
#add labels          
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.85)
}, bg.border = NA)

circos.clear()
# </废>



