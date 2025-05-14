setwd('D:/GLP1RA/')
rnr = read.csv('response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR based on homa-IR improvement
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

<1> global network
library(stringr)
library(psych)
f = read.csv('gene_4level_T1.csv', row.names = 1)
f = f[grepl('[FGH]', colnames(f))]; colnames(f) = str_sub(colnames(f),1,-4)
f = f[!colnames(f) %in% test_id]

f = f[grepl('K', rownames(f)), ]
dat_R = f[colnames(f) %in% names(rnr_)[rnr_=='R']]
dat_NR = f[colnames(f) %in% names(rnr_)[rnr_=='NR']]


filt_low = function(df, cutOff = 1e-4){
  # feature * sample
  a = rowMeans(df)
  a = names(a)[a > cutOff]
  df[a,]
}
dat_R = filt_low(dat_R)
dat_NR = filt_low(dat_NR)
co = intersect(rownames(dat_R), rownames(dat_NR)) 
dat_R = t(dat_R[co, ]); dat_NR = t(dat_NR[co, ]);  dat_R0 = dat_R; dat_NR0 = dat_NR

dim(dat_R); dim(dat_NR)  

a = corr.test(t(dat_NR), method = 'spearman')$r
which(rowSums(a > 0.88) == 7) 
dat_NR = dat_NR0[colnames(a)[a['H7', ] > 0.88], ]  # species-1e-4

# python for correlation
library(reticulate)
use_python('E:/miniconda/install/python.exe')
py_available()
py_module_available('numpy')

np = import('numpy')
stats = import('scipy.stats')

dat_R_spearmanr = stats$spearmanr(dat_R)
dat_R_cor = dat_R_spearmanr[0]
dat_R_pval = dat_R_spearmanr[1]

dat_NR_spearmanr = stats$spearmanr(dat_NR)
dat_NR_cor = dat_NR_spearmanr[0]
dat_NR_pval = dat_NR_spearmanr[1]

BH_adj = function(p_mat){
  p_mat[upper.tri(p_mat)] = p.adjust(p_mat[upper.tri(p_mat)], method = 'BH')
  p_mat[lower.tri(p_mat)] = p.adjust(p_mat[lower.tri(p_mat)], method = 'BH')
  diag(p_mat) = 1
  p_mat
}
dat_R_pval = BH_adj(dat_R_pval); length(dat_R_pval[dat_R_pval < 0.05])/2 
dat_NR_pval = BH_adj(dat_NR_pval); length(dat_NR_pval[dat_NR_pval < 0.05])/2
min(abs(dat_R_cor[dat_R_pval < 0.05])); min(abs(dat_NR_cor[dat_NR_pval < 0.05])) # cor ~= 1 when pval fits

# visualization
# R - co-abundance network
library(igraph)
pst = 0.001  
dat_R_pval = ifelse(dat_R_pval < pst, 1, 0)
dat_R_cor[abs(dat_R_cor) < 0.8] = 0 
adj_R = as.matrix(dat_R_cor) * as.matrix(dat_R_pval)
diag(adj_R) = 0 
rownames(adj_R) = co; colnames(adj_R) = co
adj_R = adj_R[-c(1:nrow(adj_R))[rowSums(abs(adj_R))==0], -c(1:ncol(adj_R))[rowSums(abs(adj_R))==0]] # 有些节点没有边,删去
length(adj_R[adj_R!=0])

# "Metabolism" > "Genetic Information Processing" > 'Environ' > cellular > orga > human disease
tax = read.delim('../DEGanalysis_202406/KEGG_brite_KO.txt')[c(2,4,6,7)]
a = c()
for (i in co){
  j = c(1:nrow(tax))[grepl(i, tax$l4_ids)]
  j = unique(tax[j, 1])[1] 
  j = ifelse(length(j)==1, j, paste0(j, collapse = '/') )
  a = c(a, j)
}
# a[a %in% c('Brite Hierarchies', 'Not Included in Pathway or Brite')] = 'Others'
color = structure(c('#f8f8f8', '#cccccc', '#84c3b7', '#fae69e', '#b8aeeb', '#f2b56f', '#71b7ed', '#f57c6e'), names = names(sort(table(a))))
### c('#F5F5F5', '#BAB682', '#4169E1', '#b3de69', '#fdb462', '#008B00', '#EE0000'
taxonomy = data.frame(row.names = co, level1 = a, color = color[a], 
                      relative_abundance_R = colMeans(dat_R),
                      relative_abundance_NR = colMeans(dat_NR))

igraph_matrix = graph_from_adjacency_matrix(adj_R, mode = 'undirected', weighted = TRUE, diag = FALSE)
igraph_matrix    

table(degree(igraph_matrix))
igraph.weight = E(igraph_matrix)$weight
E(igraph_matrix)$weight = NA

E(igraph_matrix)$color = adjustcolor(as.character(ifelse(igraph.weight>0, "#FF0000", "#0000FF")), alpha.f = 0.5) ##D74E96红#00A9E8蓝##E1E0E1灰
E(igraph_matrix)$lty=as.character(ifelse(igraph.weight>0, 'solid','dashed'))
V(igraph_matrix)$color = adjustcolor(taxonomy[V(igraph_matrix)$name,]$color, alpha.f = 1)
V(igraph_matrix)$size = (taxonomy[V(igraph_matrix)$name,]$relative_abundance_R*1e5)**0.5
V(igraph_matrix)$label = V(igraph_matrix)$name
set.seed(2222)
plot(igraph_matrix, 
     main="Co-occurrence network of functional genes",
     vertex.frame.color=  'black', #NA, #
     vertex.label= '',                  # vertex.label.dist=1, vertex.label.cex=0.25, vertex.label.font=2, vertex.label.color ='black',
     edge.width=1.1,                      #edge.lty=1, vertex.label=V(igraph_matrix)$label
     edge.curved = F,
     margin=c(0,0,0,0), 
     layout=layout_as_tree(igraph_matrix, circular = T)) 
#legend('topright',c('p__Firmicutes','p__Proteobacteria','p__Verrucomicrobia','p__Bacteroidetes','p__Actinobacteria','p__Spirochaetes','p__Chlamydiae','p__Synergistetes','p__Elusimicrobia','p__Fibrobacteres','p__Ascomycota','p__Euryarchaeota','p__Tenericutes','p__Unclassified','p__Mucoromycota','p__Candidatus_Melainabacteria','p__Lentisphaerae'),fill = c('#EE0000','#008B00','#FFF8DC','#fb8072','#fdb462','#87CEFF','#b3de69','#FFEFD5','#ffffb3','#E9967A','#a4c2f4','#bb7dbc','#FFEC8B','#FFA54F','#FFD39B','#8dd3c7','#DEB887'),col=c('#EE0000','#008B00','#FFF8DC','#fb8072','#fdb462','#87CEFF','#b3de69','#FFEFD5','#ffffb3','#E9967A','#a4c2f4','#bb7dbc','#FFEC8B','#FFA54F','#FFD39B','#8dd3c7','#DEB887'),pch = 1)
color = rev(color)
plot(1)
legend(x=0.6, y=1.2, names(color), fill = color, cex = 0.8, ncol = 2)
# attr
num_node = length(V(igraph_matrix))
num_edge = length(E(igraph_matrix))
num_average_degree = mean(igraph::degree(igraph_matrix))
num_average_path = mean_distance(igraph_matrix)

a = data.frame(property = c('num_node', 'num_edge', 'average_degree', 'average_path'),
               R = c(963, 5072, 10.53375, 8.20776), NR = c(978, 2985, 6.1043, 9.3569))
write.csv(a, 'co_abundance_network/ko_Co_abundance_network_property.csv', row.names = F)
</1> global network



<3> co-network species
f = read.csv('../direct_data/spe_7level_T1.csv', row.names = 1)
tax = f[c('tax', 'label')]
f = f[grepl('[FGH]', colnames(f))]; colnames(f) = str_sub(colnames(f),1,-4)
f = f[!colnames(f) %in% test_id]

f = f = f[tax$label=='s',]  # 只取level7-species
dat_R = f[colnames(f) %in% names(rnr_)[rnr_=='R']]
dat_NR = f[colnames(f) %in% names(rnr_)[rnr_=='NR']]

tax = tax[co, ]
tax['label'] = str_extract(tax$tax, 'p__[a-zA-Z0-9]+') 
color = structure(c('#cccccc', '#84c3b7', '#fae69e', '#b8aeeb', '#f2b56f', '#71b7ed', '#f57c6e'),
                  names = c('Others', 'p__Euryarchaeota', 'p__Lentisphaerae', 'p__Spirochaetes', 'p__Proteobacteria', 
                            'p__Bacteroidetes', 'p__Firmicutes'))

tax['color'] = color[tax$label];  tax$color[is.na(tax$color)] = '#cccccc' 
taxonomy = data.frame(row.names = co, level1 = tax$label, color = tax$color, 
                      relative_abundance_R = colMeans(dat_R),
                      relative_abundance_NR = colMeans(dat_NR))

############## same...

a = data.frame(property = c('num_node', 'num_edge', 'average_degree', 'average_path'),
               R = c(422, 1430, 6.777251, 8.433923), NR = c(415, 613, 2.954217, 11.72065))
write.csv(a, 'co_abundance_network/spe_Co_abundance_network_property.csv', row.names = F)

</3> co-network species




<4> ssn-spe
library(Hmisc)
ctl = read.csv('../direct_data/spe_7level_T1.csv', row.names = 1)[co, ]  #  '.../gene_4level_T1.csv'
colnames(ctl) = str_sub(colnames(ctl), 1, -4)
test = ctl[colnames(ctl) %in% test_id]
ctl = ctl[grepl('E', colnames(ctl))]; 
SSN = t(rbind(dat_R, dat_NR0, t(test)))

for (i in 1:ncol(SSN)) {
  dat = as.data.frame(SSN)[i]
  colnames(dat) = paste0('SSN_', colnames(SSN)[i])
  dat_exp = cbind(ctl, dat)
  
  corr_ctl = rcorr(t(ctl), type = 'spearman')
  ctl_cor = corr_ctl$r; diag(ctl_cor) = 0; ctl_cor[ctl_cor==1] = 0.99999; ctl_cor[ctl_cor==-1] = -0.99999
  ctl_pval = corr_ctl$P; diag(ctl_pval) = 1
  
  corr_s = rcorr(t(dat_exp), type = 'spearman')
  s_cor = corr_s$r; diag(s_cor) = 0; s_cor[s_cor==1] = 0.99999; s_cor[s_cor==-1] = -0.99999
  s_pval = corr_s$P; diag(s_pval) = 1
  
  pst = 0.01
  length(ctl_cor[(abs(ctl_cor) > 0.8) * (ctl_pval < pst)])/2  
  
  ctl_cor[!(abs(ctl_cor) > 0.8) * (ctl_pval < pst)] = NA
  s_cor[!(abs(s_cor) > 0.8) * (s_pval < pst)] = NA
  
  delta = s_cor - ctl_cor
  z = delta / ((1 - ctl_cor ** 2) / (6 - 1))   # 6 ctl
  p = data.frame(stats$norm$sf(abs(z)), row.names = rownames(ctl_cor))
  colnames(p) = colnames(ctl_cor)
  
  for (k in 1:ncol(p)) { a = p[k]; a[!is.na(a)] = p.adjust(a[!is.na(a)], method = 'BH'); p[k]=a }
  
  p[is.na(p)] = 1; length(p[p<0.05])
  p = ifelse(p<0.05, 1, 0)
  p = p*delta; p[is.na(p)] = 0; p = ifelse(p<0, -1, ifelse(p>0, 1, 0))
  
  if (length(which(rowSums(p)==0))>0) { p = p[-which(rowSums(p)==0), -which(colSums(p)==0)]}  # 删 孤点
  
  print(i); print(length(p[p!=0]))
  write.csv(p, paste0('ssn_1e4_for_fig/SSN_', colnames(SSN)[i], '.csv')) # 'ssn_1e4_for_fig/SSN_ko/SSN_'
}
</4> ssn-spe/ko  


<5> ssn-RF-feature visualization
library(igraph)
df = read.csv('gene_allMightFeatures_need_deeper_filter_r5_nr11.csv', row.names = 1, check.names = F) 
dat = read.csv('geneRF_allMightFeatures_need_deeper_filter_r5_nr11.csv', row.names = 2, check.names = F) 
colnames(dat) = str_replace(colnames(dat), '--', '++')
score = read.csv('geneFeature_importance_score.csv', row.names = 1)
need = str_replace_all(score$edge, '--', '++')

need_ = unique(unlist(strsplit(need, '[+]{2}')))
f = na.omit(read.delim('../raw_data/count result/KEGG/Relative/Unigenes.relative.ko.xls',  # spe:'../raw_data/count result/species/Relative/Unigenes.relative.s.xls'
                     quote = '', row.names = 1)[need_, ])
colnames(f) = c(colnames(f)[2:ncol(f)], 'tax')
f = f[grepl('[FGH][0-9]+.T1', colnames(f))]
colnames(f) = str_sub(colnames(f), 1, -4)

tax = read.delim('../DEGanalysis_202406/KEGG_brite_KO.txt')[c(2,4,6,7)]
a = c()
for (i in need_){
  j = c(1:nrow(tax))[grepl(i, tax$l4_ids)]
  j = unique(tax[j, 1])[1] 
  j = ifelse(length(j)==1, j, paste0(j, collapse = '/') )
  a = c(a, j)
}
color = structure(c('#f8f8f8', '#cccccc', '#84c3b7', '#fae69e', '#b8aeeb', '#f2b56f', '#71b7ed', '#f57c6e'),
                  names = c("Human Diseases", "Not Included in Pathway or Brite", "Organismal Systems",
                            "Genetic Information Processing", "Cellular Processes", "Brite Hierarchies",
                            "Environmental Information Processing", "Metabolism" ))
taxonomy = data.frame(row.names = need_, level1 = a, color = color[a])

dat_ = dat[dat$type=='train', c('group', need)]
rownames(df) = str_replace(rownames(df), '--','++')
par(mfrow=c(2,5), mar=c(0,0,0,0))
for (i in 1:7) {    # 8:20  par(mfrow=c(3,5), mar=c(0,0,0,0))
  sid = rownames(dat_)[i]; effect = dat_[i, 'group'] 
  a = dat_[i, 2:ncol(dat_)]
  a_ = matrix(unlist(strsplit(names(a), '[+]{2}')), nrow=2)
  b = data.frame(n1 = a_[1,], n2 = a_[2, ], e = as.numeric(a))
  zero_node = unique(as.vector(as.matrix(b[b$e==0, 1:2])))
  b = b[b$e == 1, 1:2]
  zero_node = zero_node[! zero_node %in% c(b$n1, b$n2)]
  
  g = graph_from_edgelist(as.matrix(b), directed = F)
  # if (length(zero_node)>0) {g = add_vertices(g, length(zero_node), name = zero_node)}  # 是否画孤点
  g.weight = df[paste0(b$n1,'++',b$n2), 'enrich_group']
  
  E(g)$color = adjustcolor(ifelse(g.weight=='R', 'blue', 'red'), alpha.f = 0.5)  ##D74E96红#00A9E8蓝
  E(g)$lty = 'solid' # ifelse(wt>0, 'solid', 'dashed')
  V(g)$color = adjustcolor(taxonomy[V(g)$name,]$color, alpha.f = 0.8)
  V(g)$size = log2((f[V(g)$name, sid]+1e-8) * 1e8)  # 8-20NR时 V(g)$size = log2((f[V(g)$name, sid]+1e-6) * 1e8)*2
  plot(g,
       vertex.frame.color='black',
       vertex.label= NA, #V(g)$label,
       edge.width = 2,
       #vertex.label.dist=0,vertex.label.degree=pi,vertex.label.cex=0.2,vertex.label.font=2,vertex.label.color ='black',
       edge.curved=F, margin=c(0,0,0,0),
       layout= layout.sphere) # layout_as_tree(g, circular = T))   #

}


</5> ssn-RF-feature visualization
