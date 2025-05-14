setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 


library(igraph)
library(stringr)
setwd("D:/Zhang Di/linshibangong/linshiMonkey/杭高院/new_ssn_202406/")
SSN = c( "F1", "F8", "F10", "G4", "G5", "G10", "H3", # 7 R sampleID - train
         "F2", "F4", "F7", "G1", "G6", "G7",  "G8", "G9", "H1", "H4", "H6", "H7", "H9",# 13 NR sampleID - train
         "F3", "F5", "F6", "F9", "G2", "G3", "H2", "H5", "H8", "H10" ) # 10 test sampleID



{<species SSN> --------------------------------------------------------------
# 因为每个SSN中的节点并不一致，在网络结构上不好直接比较
# 因此只提取每个SSN都有的节点绘制网络图
a = c()
for (i in 1:20){
  fname = paste0('ssn_1e4_for_fig/SSN_', SSN[i], '.csv')
  sp = read.csv(fname, row.names = 1, encoding = 'utf-8') %>% rownames()
  a = c(a, sp)
}
a = table(a)
co = names(a)[a==20] # 在20个SSN每个SSN都存在的节点


## 将每个SSN (提co节点), 剔除在75%另一组样本中出现的边后, 画SSN
par(mfrow=c(2,5), mar=c(0,0,0,0))
st = 0.75; st_r = round(13*st)-1; st_nr = round(7*st)-1
for (i in 1:7){
  fname = paste0('ssn_1e4_for_fig/SSN_', SSN[i], '.csv')
  f = read.csv(fname, row.names = 1, encoding = 'utf-8'); colnames(f)=rownames(f)
  f = f[co, co] #%>% abs() #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #tax = taxonomy; tax['this_abundace'] = SSN[i]  # 不需要对节点设置颜色太花，想设置参考ssn_fig.R中的taxonomy
  t = f
  
  for (j in 8:20){
    fname = paste0('ssn_1e4_for_fig/SSN_', SSN[j], '.csv')
    ft = read.csv(fname, row.names = 1, encoding = 'utf-8'); colnames(ft)=rownames(ft)
    ft = ft[co, co]  #%>% abs() #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    f[t==ft] = (f[t==ft] + t[t==ft])
    #print(table(as.matrix(f)))
  }

  t[f < (-st_r)] =0; t[f > st_r] = 0  # when R
  #t[f < (-st_nr) ] =0; t[f > st_nr] = 0  # when NR  i,j 反一下
  f = t
  print(table(as.matrix(f)))

  g = graph_from_adjacency_matrix(as.matrix(f), mode = 'undirected', weighted = T, diag = FALSE)
  g   #igraph 的邻接矩
  
  #bad.vs = V(g)$name[degree(g) < 2]
  #g = delete_vertices(g, bad.vs)
  #bad.vs = V(g)$name[degree(g) == 0]
  #g = delete_vertices(g, bad.vs)
  wt = E(g)$weight; E(g)$weight = NA
  V(g)$label = str_sub(str_extract(V(g)$name, 's__[\\S ]+'), 4)
  E(g)$color = adjustcolor(ifelse(wt>0, 'red', '#3b7cb7'), alpha.f = 0.5)  ##D74E96红#00A9E8蓝
  E(g)$lty = 'solid' # ifelse(wt>0, 'solid', 'dashed')
  #V(g)$color = adjustcolor(tax[V(g)$name,]$color, alpha.f = 0.8)
  V(g)$color = 'white'
  V(g)$size = 15 #*(tax[V(g)$name, 'this_abundace']*2000)
  plot(g,
       vertex.frame.color='black',vertex.label= NA, #V(g)$label,
       edge.width = 2,
       #vertex.label.dist=0,vertex.label.degree=pi,vertex.label.cex=0.2,vertex.label.font=2,vertex.label.color ='black',
       edge.curved=F, margin=c(0,0,0,0),
       layout= layout.sphere)
}


for (i in 8:20){
  fname = paste0('ssn_1e4_for_fig/SSN_', colnames(SSN)[i], '.csv')
  f = read.csv(fname, row.names = 1, encoding = 'utf-8'); colnames(f)=rownames(f)
  f = f[co, co] #%>% abs() #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tax = taxonomy; tax['this_abundace'] = SSN[i]  # 继承上面co-abundance的tax信息
  t = f
  
  for (j in 1:7){
    fname = paste0('ssn_1e4_for_fig/SSN_', colnames(SSN)[j], '.csv')
    ft = read.csv(fname, row.names = 1, encoding = 'utf-8'); colnames(ft)=rownames(ft)
    ft = ft[co, co]  #%>% abs() #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    f[t==ft] = (f[t==ft] + t[t==ft])
    
    #print(table(as.matrix(f)))
  }
  #t[f < (-st_r)] =0; t[f > st_r] = 0  # when R
  t[f < (-st_nr) ] =0; t[f > st_nr] = 0  # when NR  i,j 反一下
  f = t
  print(table(as.matrix(f)))
  
  g = graph_from_adjacency_matrix(as.matrix(f), mode = 'undirected', weighted = T, diag = FALSE)
  g   #igraph 的邻接矩
  
  #bad.vs = V(g)$name[degree(g) < 2]
  #g = delete_vertices(g, bad.vs)
  bad.vs = V(g)$name[degree(g) == 0]
  g = delete_vertices(g, bad.vs)
  wt = E(g)$weight; E(g)$weight = NA
  V(g)$label = str_sub(str_extract(V(g)$name, 's__[\\S ]+'), 4)
  E(g)$color = adjustcolor(ifelse(wt>0, 'red', '#3b7cb7'), alpha.f = 0.5)  ##D74E96红#00A9E8蓝
  E(g)$lty = 'solid' # ifelse(wt>0, 'solid', 'dashed')
  #V(g)$color = adjustcolor(tax[V(g)$name,]$color, alpha.f = 0.8)
  V(g)$color = 'white'
  V(g)$size = 15 #*(tax[V(g)$name, 'this_abundace']*2000)
  plot(g,
       vertex.frame.color='black',vertex.label= NA, #V(g)$label,
       edge.width = 2,
       #vertex.label.dist=0,vertex.label.degree=pi,vertex.label.cex=0.2,vertex.label.font=2,vertex.label.color ='black',
       edge.curved=F, margin=c(0,0,0,0),
       layout= layout.sphere)
}
</species SSN> --------------------------------------------------------------}
  
<species SSN signature>
s = read.csv('speFeature_importance_score.csv')$edge[1:50] 
ss = str_split(s, '[+-]{2}') %>% unlist() %>% matrix(nrow = 2) %>% t() %>% data.frame()
ss = cbind(ss, ifelse(grepl('\\++', s), 1, 0))
colnames(ss) = c('from', 'to', 'weight')
g = graph_from_data_frame(data.frame(ss), directed = F)
weight = E(g)$weight  # 1:R feature 0: NR feature
E(g)$weight = NA


E(g)$color = adjustcolor(as.character(ifelse(weight>0, "#0000FF", "#FF0000" )), alpha.f = 0.5) ##D74E96红#00A9E8蓝##E1E0E1灰
E(g)$lty=as.character(ifelse(igraph.weight>0, 'solid','dashed'))
V(g)$color = adjustcolor(taxonomy[V(g)$name,]$color, alpha.f = 1)
V(g)$size =15 #(taxonomy[V(g)$name,]$relative_abundance_R*1e5)**0.5
V(g)$label = V(g)$name
plot(g, vertex.label= '', 
     edge.width=1.5,                      #edge.lty=1, vertex.label=V(igraph_matrix)$label
     edge.curved = F,
     layout = layout.sphere)
  
  
  
  
  
{<gene SSN> --------------------------------------------------------------
f = read.csv('../direct_data/gene_4level_T1.csv', row.names = 1)
f = f[grepl('[FGH]', colnames(f))]; colnames(f) = str_sub(colnames(f),1,-4)
f = f[!colnames(f) %in% test_id]

m1 = f[SSN[1:7]] %>% rowMeans() %>% sort(); m1 = names(m1)[m1 > 4e-4]
m2 = f[SSN[8:20]] %>% rowMeans() %>% sort(); m2 = names(m2)[m2 > 4e-4]

m = intersect(m1, m2)
# 因为每个SSN中的节点并不一致，在网络结构上不好直接比较
# 因此只提取每个SSN都有的节点绘制网络图
a = c()
for (i in 1:20){
  fname = paste0('ssn_1e4_for_fig/SSN_ko/SSN_', SSN[i], '.csv')
  sp = read.csv(fname, row.names = 1, encoding = 'utf-8') %>% rownames()
  a = c(a, sp)
}
a = table(a)
co = names(a)[a==20] # 在20个SSN每个SSN都存在的节点
co = co[co %in% m]


## 将每个SSN (提co节点), 剔除在75%另一组样本中出现的边后, 画SSN
par(mfrow=c(2,5), mar=c(0,0,0,0))
st = 0.75; st_r = round(13*st) -1; st_nr = round(7*st)-1
#st_r = 10; st_nr = 4
for (i in 1:7){
  fname = paste0('ssn_1e4_for_fig/SSN_ko/SSN_', SSN[i], '.csv')
  f = read.csv(fname, row.names = 1, encoding = 'utf-8'); colnames(f)=rownames(f)
  f = f[co, co] #%>% abs() #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  #tax = taxonomy; tax['this_abundace'] = SSN[i]  # 不需要对节点设置颜色太花，想设置参考ssn_fig.R中的taxonomy
  t = f
  
  for (j in 8:20){
    fname = paste0('ssn_1e4_for_fig/SSN_ko/SSN_', SSN[j], '.csv')
    ft = read.csv(fname, row.names = 1, encoding = 'utf-8'); colnames(ft)=rownames(ft)
    ft = ft[co, co]  #%>% abs() #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    f[t==ft] = (f[t==ft] + t[t==ft])
    #print(table(as.matrix(f)))
  }
  
  t[f < (-st_r)] =0; t[f > st_r] = 0  # when R
  #t[f < (-st_nr) ] =0; t[f > st_nr] = 0  # when NR  i,j 反一下
  f = t
  print(table(as.matrix(f)))
  
  g = graph_from_adjacency_matrix(as.matrix(f), mode = 'undirected', weighted = T, diag = FALSE)
  g   #igraph 的邻接矩
  
  #bad.vs = V(g)$name[degree(g) < 2]
  #g = delete_vertices(g, bad.vs)
  bad.vs = V(g)$name[degree(g) == 0]
  g = delete_vertices(g, bad.vs)
  wt = E(g)$weight; E(g)$weight = NA
  E(g)$color = adjustcolor(ifelse(wt>0, 'red', '#3b7cb7'), alpha.f = 0.5)  ##D74E96红#00A9E8蓝
  E(g)$lty = 'solid' # ifelse(wt>0, 'solid', 'dashed')
  #V(g)$color = adjustcolor(tax[V(g)$name,]$color, alpha.f = 0.8)
  V(g)$color = 'white'
  V(g)$size = 15 #*(tax[V(g)$name, 'this_abundace']*2000)
  plot(g,
       vertex.frame.color='black',vertex.label= NA, #V(g)$label,
       edge.width = 2,
       #vertex.label.dist=0,vertex.label.degree=pi,vertex.label.cex=0.2,vertex.label.font=2,vertex.label.color ='black',
       edge.curved=F, margin=c(0,0,0,0),
       layout= layout.sphere)
}


for (i in 8:20){
  fname = paste0('ssn_1e4_for_fig/SSN_ko/SSN_', SSN[i], '.csv')
  f = read.csv(fname, row.names = 1, encoding = 'utf-8'); colnames(f)=rownames(f)
  f = f[co, co] #%>% abs() #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tax = taxonomy; tax['this_abundace'] = SSN[i]  # 继承上面co-abundance的tax信息
  t = f
  
  for (j in 1:7){
    fname = paste0('ssn_1e4_for_fig/SSN_ko/SSN_', SSN[j], '.csv')
    ft = read.csv(fname, row.names = 1, encoding = 'utf-8'); colnames(ft)=rownames(ft)
    ft = ft[co, co]  #%>% abs() #!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    f[t==ft] = (f[t==ft] + t[t==ft])
    
    #print(table(as.matrix(f)))
  }
  #t[f < (-st_r)] =0; t[f > st_r] = 0  # when R
  t[f < (-st_nr) ] =0; t[f > st_nr] = 0  # when NR  i,j 反一下
  f = t
  print(table(as.matrix(f)))
  
  g = graph_from_adjacency_matrix(as.matrix(f), mode = 'undirected', weighted = T, diag = FALSE)
  g   #igraph 的邻接矩
  
  #bad.vs = V(g)$name[degree(g) < 2]
  #g = delete_vertices(g, bad.vs)
  bad.vs = V(g)$name[degree(g) == 0]
  g = delete_vertices(g, bad.vs)
  wt = E(g)$weight; E(g)$weight = NA
  E(g)$color = adjustcolor(ifelse(wt>0, 'red', '#3b7cb7'), alpha.f = 0.5)  ##D74E96红#00A9E8蓝
  E(g)$lty = 'solid' # ifelse(wt>0, 'solid', 'dashed')
  #V(g)$color = adjustcolor(tax[V(g)$name,]$color, alpha.f = 0.8)
  V(g)$color = 'white'
  V(g)$size = 15 #*(tax[V(g)$name, 'this_abundace']*2000)
  plot(g,
       vertex.frame.color='black',vertex.label= NA, #V(g)$label,
       edge.width = 2,
       #vertex.label.dist=0,vertex.label.degree=pi,vertex.label.cex=0.2,vertex.label.font=2,vertex.label.color ='black',
       edge.curved=F, margin=c(0,0,0,0),
       layout= layout.sphere)
}
</gene SSN> --------------------------------------------------------------}