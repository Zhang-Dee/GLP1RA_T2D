setwd('D:/GLP1RA/')
rnr = read.csv('response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR based on homa-IR improvement
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

library(igraph)
library(stringr)
setwd("direct_data")
f = read.csv('spe_7level_T1.csv', row.names = 1, check.names = F)  
tax = f[c('tax', 'label')]
f = f[grepl('[FGH]', colnames(f))]; colnames(f) = str_sub(colnames(f),1,-4)
f = f[!colnames(f) %in% test_id]

f3 = read.csv('spe_7level_T3.csv', row.names = 1, check.names = F)  
f3 = f3[grepl('[FGH]', colnames(f3))]; colnames(f3) = str_sub(colnames(f3),1,-4)
f3 = f3[!colnames(f3) %in% test_id]

m1 = f[rnr_[colnames(f)]=='R'] %>% rowMeans() %>% sort(); m1 = names(m1)[m1 > 1e-4]
m2 = f[rnr_[colnames(f)]=='NR'] %>% rowMeans() %>% sort(); m2 = names(m2)[m2 > 1e-4]
m = intersect(m1, m2) 

m1 = f3[rnr_[colnames(f3)]=='R'] %>% rowMeans() %>% sort(); m1 = names(m1)[m1 > 1e-4]
m2 = f3[rnr_[colnames(f3)]=='NR'] %>% rowMeans() %>% sort(); m2 = names(m2)[m2 > 1e-4]
m = intersect(m1, m2) %>% intersect(m) 
m = m[tax[m, 'label'] == 's'] # m = m[tax[m, 'label'] == 'ko.xls']


dat_R = f[m, rnr_[colnames(f)]=='R']
dat_NR = f[m, rnr_[colnames(f)]=='NR']

dat_R_t3 = f3[m, rnr_[colnames(f3)]=='R']
dat_NR_t3 = f3[m, rnr_[colnames(f3)]=='NR']


# sample size balance
a = corr.test(dat_NR, method = 'spearman')$r
rst = 0.63 # T1-spe: 0.63; gene: 0.88; T3-spe: 0.65   gene:0.93
s = which(rowSums(a > rst) == 7) %>% names() 
dat_NR = dat_NR[colnames(a)[a[s, ] > rst]] 
# dat_NR_t3 = dat_NR_t3[colnames(dat_NR)]
dim(dat_NR)

# python-spearman
library(reticulate)
library(psych)
use_condaenv('base')  
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

{KO-taxonomy table-->
  # "Metabolism" > "Genetic Information Processing" > 'Environ' > cellular > orga > human disease
  tax = read.delim('../DEGanalysis_202406/KEGG_brite_KO.txt')[c(2,4,6,7)]
  a = c()
  for (i in m){
    j = c(1:nrow(tax))[grepl(i, tax$l4_ids)]
    j = unique(tax[j, 1])[1] 
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

# visualization
adj = adj_R  # = adj_NR
library(igraph)

igraph_matrix = graph_from_adjacency_matrix(adj, mode = 'undirected', weighted = TRUE, diag = FALSE)
igraph_matrix   
# min(degree(igraph_matrix))
degree(igraph_matrix)[order(degree(igraph_matrix))] %>% table()
bad.vs = V(igraph_matrix)[degree(igraph_matrix) == 0] 
# bad.vs = V(igraph_matrix)[degree(igraph_matrix) <= 1] 
igraph_matrix = delete_vertices(igraph_matrix, bad.vs)
igraph.weight = E(igraph_matrix)$weight
E(igraph_matrix)$weight = NA

E(igraph_matrix)$color = as.character(ifelse(igraph.weight>0, "red", "blue")) %>% adjustcolor(alpha.f = 0.5)
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
     vertex.frame.color= 'black',  # 'white', 
     vertex.label=V(igraph_matrix)$label, vertex.label.dist=1, vertex.label.cex=0.45, vertex.label.font=2,vertex.label.color ='black',
     edge.width = 1.8, edge.lty = 1, edge.curved = F, margin=c(0,0,0,0), 
     layout=layout.fruchterman.reingold)

tiff('../new_ssn_202406/co_abundance_network/co_speR_T1.tiff', width = 5, units = 'in', height = 5, res = 300)
set.seed(233)
plot(igraph_matrix,
     vertex.frame.color= 'black',  # 'white', 
     vertex.label=V(igraph_matrix)$label, vertex.label.dist=1, vertex.label.cex=0.45, vertex.label.font=2,vertex.label.color ='black',
     edge.width = 1.8, edge.lty = 1, edge.curved = F, margin=c(0,0,0,0), 
     layout=layout.fruchterman.reingold)
dev.off()

# legend
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

# network attr
num_node = length(V(igraph_matrix))
num_edge = length(E(igraph_matrix))
num_average_degree = mean(igraph::degree(igraph_matrix))
num_average_path = mean_distance(igraph_matrix, directed = F)


# hubs
d = degree(g) %>% sort(decreasing = T)
quantile(d, 0.95)
d[d >= quantile(d, 0.95)]
# c_node = 'g__Unclassified;s__Ruminococcaceae bacterium D5'
c_nodes = d[d >= quantile(d, 0.95)] %>% names()
sub_g = induced_subgraph(g, vid = c(which(V(g)$name %in% c_nodes), 
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
                                                               core=list(fg_params=list(hjust=0, x=0.05)),  
                                                               rowhead=list(fg_params=list(hjust=0, x=0)))) 

grid.newpage() 
grid.draw(table_plot)


# attr barplot-T1-RNR
a = read.csv('spe_Co_abundance_network_property.csv', row.names = 1)[5:6]
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


# attr barplot-T1T3  
a = read.csv('spe_Co_abundance_network_property.csv', row.names = 1)
a = t(a[3:6]) %>% data.frame()
a['response'] = str_extract(rownames(a), '[NR]+')
a['time'] = str_extract(rownames(a), 'T[13]+')
a = melt(a)
a$response = str_replace(a$response, 'NR', 'LR')
ggplot(a, aes(time, value, fill=response)) +
  geom_bar(stat = 'identity', 
           position = 'dodge',
           width = 0.75) +
  ylab(NULL) + xlab(NULL) +
  scale_fill_manual(values = c('#DC8070', '#6188B5')) +
  theme_bw() + theme(aspect.ratio = 3/2) +
  facet_wrap(~variable, scales = 'free', nrow =1)

ggplot(a, aes(response, value, fill=response, alpha=time)) +
  geom_bar(stat = 'identity', 
           position = position_dodge(width = 0.85), 
           width = 0.75) +
  ylab(NULL) + xlab(NULL) +
  scale_fill_manual(values = c('#DC8070', '#6188B5')) +
  theme_bw() + theme(aspect.ratio = 3/2) +
  facet_wrap(~variable, scales = 'free', nrow =2)

