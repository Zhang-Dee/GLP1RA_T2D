setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 


library(igraph)
library(stringr)
setwd("D:/Zhang Di/linshibangong/linshiMonkey/杭高院/direct_data")
f = read.csv('spe_7level_T1.csv', row.names = 1, check.names = F)   # gene同  T3时要跟T1做co-spe以保证网络大小同
tax = f[c('tax', 'label')]
f = f[grepl('[FGH]', colnames(f))]; colnames(f) = str_sub(colnames(f),1,-4)
f = f[!colnames(f) %in% test_id]

f3 = read.csv('spe_7level_T3.csv', row.names = 1, check.names = F)   # gene同  T3时要跟T1做co-spe以保证网络大小同
f3 = f3[grepl('[FGH]', colnames(f3))]; colnames(f3) = str_sub(colnames(f3),1,-4)
f3 = f3[!colnames(f3) %in% test_id]

m1 = f[rnr_[colnames(f)]=='R'] %>% rowMeans() %>% sort(); m1 = names(m1)[m1 > 1e-4]
m2 = f[rnr_[colnames(f)]=='NR'] %>% rowMeans() %>% sort(); m2 = names(m2)[m2 > 1e-4]
m = intersect(m1, m2)  # 只用相对丰度>1e-4的物种作共存网络

m1 = f3[rnr_[colnames(f3)]=='R'] %>% rowMeans() %>% sort(); m1 = names(m1)[m1 > 1e-4]
m2 = f3[rnr_[colnames(f3)]=='NR'] %>% rowMeans() %>% sort(); m2 = names(m2)[m2 > 1e-4]
m = intersect(m1, m2) %>% intersect(m) # 只用相对丰度>1e-4的物种作共存网络
m = m[tax[m, 'label'] == 's'] # m = m[tax[m, 'label'] == 'ko.xls']


dat_R = f[m, rnr_[colnames(f)]=='R']
dat_NR = f[m, rnr_[colnames(f)]=='NR']

dat_R_t3 = f3[m, rnr_[colnames(f3)]=='R']
dat_NR_t3 = f3[m, rnr_[colnames(f3)]=='NR']


# 将dat_NR样本数缩减至7和R一样
a = corr.test(dat_NR, method = 'spearman')$r
rst = 0.63 # T1-spe: 0.63; gene: 0.88; T3-spe: 0.65   gene:0.93
s = which(rowSums(a > rst) == 7) %>% names() # 此时和H7相关性高于0.88的正好6个，根据样本间相关性只保留7个NR样本，因为网络大小和样本量有关，需要一致
dat_NR = dat_NR[colnames(a)[a[s, ] > rst]]  # T1-species是0.62-G7; 1e-6-0.77  
# dat_NR_t3 = dat_NR_t3[colnames(dat_NR)]
dim(dat_NR)

# 计算spearman相似性-调python
library(reticulate)
library(psych)
use_condaenv('base')  # 指定Python所在conda虚拟环境
py_run_string('from scipy.stats import spearmanr')

corr_R = py$spearmanr(t(dat_R))
corp_R = corr_R[1]; corr_R = corr_R[0]

corr_NR = py$spearmanr(t(dat_NR))
corp_NR = corr_NR[1]; corr_NR = corr_NR[0]

corr_R_t3 = py$spearmanr(t(dat_R_t3))
corp_R_t3 = corr_R_t3[1]; corr_R_t3 = corr_R_t3[0]

corr_NR_t3 = py$spearmanr(t(dat_NR_t3))
corp_NR_t3 = corr_NR_t3[1]; corr_NR_t3 = corr_NR_t3[0]

taxonomy = tax[m, 'tax'] %>% str_split(';', 3) %>% unlist() %>% matrix(nrow = 3)
taxonomy = taxonomy[2, ]; table(taxonomy) %>% sort()

color = structure(c("#f57c6e", "#71b7ed", "#f2b56f", "#b8aeeb", "#fae69e",
                    "#84c3b7", "#cccccc"),  # "#f8f8f8" 
                  names = c('p__Firmicutes','p__Bacteroidetes', 'p__Proteobacteria',
                            'p__Spirochaetes', 'p__Actinobacteria', 'p__Lentisphaerae', 
                            'Others'))
taxonomy = data.frame(phylum = taxonomy, color = color[taxonomy], 
                      relative_abundance_R = rowMeans(dat_R),
                      relative_abundance_NR = rowMeans(dat_NR),
                      ra_R_t3 = rowMeans(dat_R_t3),
                      ra_NR_t3 = rowMeans(dat_NR_t3),
                      row.names = m)
taxonomy$color[is.na(taxonomy$color)] = '#cccccc'

{KO的taxonomy表制作-->
  # 为每个ko节点确定功能大类："Metabolism" > "Genetic Information Processing" > 'Environ' > cellular > orga > human disease
  tax = read.delim('../DEGanalysis_202406/KEGG_brite_KO.txt')[c(2,4,6,7)]
  a = c()
  for (i in m){
    j = c(1:nrow(tax))[grepl(i, tax$l4_ids)]
    j = unique(tax[j, 1])[1] # 如果ko参与多种功能，只取其一
    j = ifelse(length(j)==1, j, paste0(j, collapse = '/') )
    a = c(a, j)
  }
  # a[a %in% c('Brite Hierarchies', 'Not Included in Pathway or Brite')] = 'Others'
  color = structure(c('#f8f8f8', '#cccccc', '#84c3b7', '#fae69e', '#b8aeeb', '#f2b56f', '#71b7ed', '#f57c6e'),
                    names = names(sort(table(a))))
  ### c('#F5F5F5', '#BAB682', '#4169E1', '#b3de69', '#fdb462', '#008B00', '#EE0000'
  taxonomy = data.frame(phylum = a, color = color[a], 
                        relative_abundance_R = rowMeans(dat_R),
                        relative_abundance_NR = rowMeans(dat_NR),
                        ra_R_t3 = rowMeans(dat_R_t3),
                        ra_NR_t3 = rowMeans(dat_NR_t3),
                        row.names = m)
}



dat_R_pval = ifelse(corp_R <= 0.001, 1, 0)
dat_R_cor = ifelse(corr_R >= 0.8, 1, ifelse(corr_R<=-0.8, -1, 0))

adj = dat_R_cor * dat_R_pval
rownames(adj) = m; colnames(adj) = m
diag(adj) = 0 
table(adj)
adj_R = adj

dat_NR_pval = ifelse(corp_NR <= 0.001, 1, 0)
dat_NR_cor = ifelse(corr_NR >= 0.8, 1, ifelse(corr_NR<=-0.8, -1, 0))
adj = dat_NR_cor * dat_NR_pval
rownames(adj) = m; colnames(adj) = m
diag(adj) = 0 
table(adj)
adj_NR = adj

adj_R_T3 =...; adj_NR_T3 =... 同

# 绘图
adj = adj_R  # = adj_NR
library(igraph)

igraph_matrix = graph_from_adjacency_matrix(adj, mode = 'undirected', weighted = TRUE, diag = FALSE)
igraph_matrix    #igraph 的邻接矩
# min(degree(igraph_matrix))
degree(igraph_matrix)[order(degree(igraph_matrix))] %>% table()
bad.vs = V(igraph_matrix)[degree(igraph_matrix) == 0] # 这里先算网络属性，再删低度节点作图
# bad.vs = V(igraph_matrix)[degree(igraph_matrix) <= 1] 
igraph_matrix = delete_vertices(igraph_matrix, bad.vs)
igraph.weight = E(igraph_matrix)$weight
E(igraph_matrix)$weight = NA

E(igraph_matrix)$color = as.character(ifelse(igraph.weight>0, "red", "blue")) %>% adjustcolor(alpha.f = 0.5) # "#00A9E8" 蓝
# E(igraph_matrix)$lty = as.character(ifelse(igraph.weight>0, 'solid','dashed'))
V(igraph_matrix)$color = taxonomy[V(igraph_matrix)$name, ]$color %>% adjustcolor(alpha.f = 0.8)
V(igraph_matrix)$size = log(taxonomy[V(igraph_matrix)$name, ]$relative_abundance_R * 1e6)/1.2   # NR

#V(igraph_matrix)$label = taxonomy_R[V(igraph_matrix)$name,]$ID
V(igraph_matrix)$label = NA
c=c()
for (i in V(igraph_matrix)$name){n=strsplit(i,';s__')[[1]][2]; c=c(c,n)}
V(igraph_matrix)$label = c
set.seed(233)
plot(igraph_matrix,
     main="Co-occurrence network",
     vertex.frame.color= 'black',  # 'white', #
     vertex.label=V(igraph_matrix)$label, vertex.label.dist=1, vertex.label.cex=0.45, vertex.label.font=2,vertex.label.color ='black',
     edge.width = 1.8, edge.lty = 1, edge.curved = F, margin=c(0,0,0,0), 
     layout=layout.fruchterman.reingold)

tiff('../new_ssn_202406/co_abundance_network/co_speR_T1.tiff', width = 5, units = 'in', height = 5, res = 300)
set.seed(233)
plot(igraph_matrix,
     vertex.frame.color= 'black',  # 'white', #
     vertex.label=V(igraph_matrix)$label, vertex.label.dist=1, vertex.label.cex=0.45, vertex.label.font=2,vertex.label.color ='black',
     edge.width = 1.8, edge.lty = 1, edge.curved = F, margin=c(0,0,0,0), 
     layout=layout.fruchterman.reingold)
dev.off()

# 图例legend
color = rev(color)
plot(1)
legend(x=0.6, y=1.2, names(color), fill = color, cex = 0.8, ncol = 2)
...

tuli_ = taxonomy[V(igraph_matrix)$name, 'color']
tuli_= tuli_[!duplicated(tuli_)]
c=c()
for (i in tuli_){n=as.character(strsplit(i,'__')[[1]][2]);c=c(c,n)}
tuli=c
tuli_c=c()
for (i in tuli_){c=taxonomy[which(taxonomy$Phylum==i)[1],9]; tuli_c=c(tuli_c,c)}
legend(x=-2,y=-0,tuli,fill=tuli_c)

#网络属性
num_node = length(V(igraph_matrix))
num_edge = length(E(igraph_matrix))
num_average_degree = mean(igraph::degree(igraph_matrix))
num_average_path = mean_distance(igraph_matrix, directed = F)

# spe_R: 422, 1430, 6.777251, 8.433923; speNR: 415, 613, 2.954217, 11.72065
# ko_R: 963, 5072, 10.53375, 8.20776;   ko_NR: 978, 2985, 6.1043, 9.3569

# T1vs.T3 -> co-spe = 483 -> 
                          T1: spe_R: 405, 1328, 6.558025, 10.75423; spe_NR: 376, 705, 3.75, 7.436652
                          T3: spe_R: 395, 933, 4.724051, 8.410017;  spe_NR: 396, 896, 4.525253, 9.118901
         -> co-gene = ->
                          T1: spe_R: 860, 4100, 9.534884, 8.151264; spe_NR: 875, 2320, 5.302857, 9.638213
                          T3: spe_R: 898, 8153, 18.15813, 7.497237; spe_NR: 874, 5176, 11.84439, 7.768497
# T3-> ko_R: 963, 5072, 10.53375, 8.20776;   ko_NR: 978, 2985, 6.1043, 9.3569

## 社区
g = igraph_matrix
communities <- cluster_walktrap(g) # 默认steps = 4
community_ids <- membership(communities)  # 获取社区编号
length(community_ids %>% unique())
# layout_sphere <- layout_on_sphere(igraph_matrix)  # 使用 sphere 布局
plot(communities, g, 
     col = V(g)$color, # membership(communities), 
     edge.color = E(g)$color,
     mark.groups = communities(communities),
     vertex.label= '', edge.width=1.5, edge.curved = F, #vertex.frame.color = 'darkgrey',
     main = 'mcf-before',
     mark.border = NA, # 去除社区边框
     mark.col = NA,#community_colors[community_ids], # 设置社区填充颜色
     #mark.expand = 5, # 扩大社区填充区域
     #mark.groups = unique(community_ids), # 为不同社区设置不同填充
     layout = layout.fruchterman.reingold) # 设置透明度


color = rev(color)
plot(1)
legend(x=0.6, y=1.2, names(color), fill = color, cex = 0.8, ncol = 2)

# save(g_t1r, g_t1nr, g_t3r, g_t3nr, taxonomy, file = '../new_ssn_202406/co_abundance_network/co_ko_igraph.RData')

# 查看中心节点
d = degree(g) %>% sort(decreasing = T)
quantile(d, 0.95)
d[d >= quantile(d, 0.95)]
# c_node = 'g__Unclassified;s__Ruminococcaceae bacterium D5'
c_nodes = d[d >= quantile(d, 0.95)] %>% names()
sub_g = induced_subgraph(g, vid = c(which(V(g)$name %in% c_nodes), # 提取某节点的一阶子图
                                    unlist(adjacent_vertices(g, v=which(V(g)$name %in% c_nodes)) )))
genera = V(sub_g)$name %>% str_match('g__[\\S]+;') %>% table()
plot(sub_g, main = "Subgraph of node-", vertex.label='')
# 制表
a = d[d >= quantile(d, 0.95)]
a = data.frame(species = str_match(names(a), 's__[\\S ]+'), 
               phylum = taxonomy[names(a), 'phylum'], 
               degree = a)
a$species = str_sub(a$species, 4)
a$phylum = str_sub(a$phylum, 4)
library(gridExtra)
library(grid)
table_plot <- tableGrob(a, rows = NULL, theme = ttheme_default(padding = unit(c(9,3), 'mm'),
                                                               core=list(fg_params=list(hjust=0, x=0.05)),  # hjust0-0.5-1 左中右对齐
                                                               rowhead=list(fg_params=list(hjust=0, x=0)))) # 创建表格图形对象

grid.newpage()  # 绘制表格图形
grid.draw(table_plot)


{
  <1>. 中心节点-基因富集分析
  deg = V(sub_g)$name
  
  {# </0> reference 数据, 只第一次做一次即可
  ko_brite = read.table('../DEGanalysis_202406/KEGG_brite_KO.txt', sep = '\t', header = T, quote = '')
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
  # </0>
  }
  
  my_gene = intersect(deg, Nn)
  M = length(my_gene) 
  
  enrich_result = data.frame(t(apply(ko_brite, 1, enrich)))
  colnames(enrich_result) = c('GeneRatio','BgRatio','fold','Count','pvalue')
  enrich_result[, c('fold','Count','pvalue')] = apply(enrich_result[, c('fold','Count','pvalue')], 2, as.numeric)
  enrich_result$p.adjust = p.adjust(enrich_result$pvalue, method = 'BH')
  
  enrich_result = cbind(ko_brite[c('l1_name', 'l2_name', 'l3_name')], enrich_result)
  enrich_result = enrich_result[enrich_result$pvalue < 0.05, ]
  enrich_result = enrich_result[order(enrich_result$pvalue), ]
  enrich_result = enrich_result[enrich_result$l1_name!='Not Included in Pathway or Brite', ]
  enrich_result = enrich_result[enrich_result$l1_name!='Brite Hierarchies', ]
  
  ## 可视化
  library(ggplot2)
  library(dplyr)
  data = arrange(enrich_result, enrich_result$Count)
  data$l3_name = factor(data$l3_name, levels = data$l3_name)
  # data = data[data$pvalue<0.01, ]
  dim(data)
  
  color = structure(c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33"),
                    names=c('Metabolism','Brite Hierarchies','Cellular Processes', 
                            'Genetic Information Processing', 
                            'Environmental Information Processing', 'Human Diseases'))
  ggplot(data, aes(Count, l3_name, fill=pvalue)) +
    geom_bar(stat = 'identity', show.legend = T) +
    scale_fill_gradient(low = '#e22834', high = '#f8e3d3') +
    theme_classic() +
    ylab(NULL) + xlab('Count - T1-NR') +
    theme(axis.text.y = element_text(color = color[data$l1_name]))
  </1>. 中心节点-基因富集分析
}

# 属性柱状图-T1-RNR
a = read.csv('../new_ssn_202406/co_abundance_network/spe_Co_abundance_network_property.csv', row.names = 1)[5:6]
a = t(a) %>% data.frame()
a['response'] = str_sub(rownames(a), 1, -4) %>%
  str_replace('NR', 'LR')

a = melt(a)
ggplot(a, aes(response, value, fill=response)) + 
  geom_bar(stat = 'identity', show.legend = F, width = 0.7) + 
  scale_fill_manual(values =  c('#DC8070', '#6188B5')) +
  theme_bw() + theme(aspect.ratio = 3/2) +
  ylab(NULL) + xlab(NULL) +
  facet_wrap(~variable, scales = 'free', nrow = 1)


# 属性柱状图-T1T3  co_spe_T1T3_attr.pdf
a = read.csv('../new_ssn_202406/co_abundance_network/spe_Co_abundance_network_property.csv', row.names = 1)
a = t(a[3:6]) %>% data.frame()
a['response'] = str_extract(rownames(a), '[NR]+')
a['time'] = str_extract(rownames(a), 'T[13]+')
a = melt(a)
a$response = str_replace(a$response, 'NR', 'LR')
ggplot(a, aes(time, value, fill=response)) +
  geom_bar(stat = 'identity', 
           position = 'dodge', #使用position=position_dodge(width=0.9),可使组内柱子间隔,自行试一下。
           width = 0.75) +
  ylab(NULL) + xlab(NULL) +
  scale_fill_manual(values = c('#DC8070', '#6188B5')) +
  theme_bw() + theme(aspect.ratio = 3/2) +
  facet_wrap(~variable, scales = 'free', nrow =1)

ggplot(a, aes(response, value, fill=response, alpha=time)) +
  geom_bar(stat = 'identity', 
           position = position_dodge(width = 0.85), #使用position=position_dodge(width=0.9),可使组内柱子间隔,自行试一下。
           width = 0.75) +
  ylab(NULL) + xlab(NULL) +
  scale_fill_manual(values = c('#DC8070', '#6188B5')) +
  theme_bw() + theme(aspect.ratio = 3/2) +
  facet_wrap(~variable, scales = 'free', nrow =2)

