setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

<1>. 门水平差异丰度
library(stringr)
f = read.csv('direct_data/spe_7level_T1.csv', row.names = 1)
f = f[f$label == 'p', 1:36] 
rownames(f) = str_extract(rownames(f), 'p__.*')
colnames(f) = str_sub(colnames(f), 1, -4)
f = f[grepl('[FGH]', colnames(f))]
f = f[!colnames(f) %in% test_id]
dim(f)

r = rowMeans(f[colnames(f) %in% names(rnr_)[rnr_=='R']])
nr = rowMeans(f[colnames(f) %in% names(rnr_)[rnr_=='NR']])
df = data.frame(Phylum = rownames(f),
                mean_R = r, Proportion_R = r/sum(r),
                mean_NR = nr, Proportion_NR = nr/sum(nr),
                Overall = rowMeans(f), Proprotion = rowMeans(f)/sum(rowMeans(f)))

df = df[order(df$Proprotion, decreasing = T), ]
write.csv(df, '../202406fig/tables/phylum_mean-abundance.csv', row.names = F)
</1>. 门水平差异丰度

<2>. 差异菌种
f = read.csv('DEGanalysis_202406/DEspe_7level_T1.csv', row.names = 1)
f = f[f$tax=='s', ]
tax = t(matrix(unlist(strsplit(f$label, ';')), nrow=7))

df = cbind(f[c("meanR", "meanNR", "log2FC", "pwil")], data.frame(tax[,1:6]))
rownames(df) = str_extract(rownames(f), 's__.*')
colnames(df) = c('avg_R', 'avg_NR', 'log2FC', 'p-value', 'Kingdom',
                 'Phylum', 'Class', 'Order', 'Family', 'Genus')
write.csv(df, '../202406fig/tables/DEspecies.csv')
</2>. 差异菌种


<3>. 差异KO
f = read.csv('DEGanalysis_202406/DEgene_4level_T1.csv', row.names = 1)
f = f[f$tax=='ko.xls', 1:4]
colnames(f) = c('avg_R', 'avg_NR', 'log2FC', 'p-value')
write.csv(f, '../202406fig/tables/DE_KO.csv')
</3>. 差异KO

<4> SPEcies features in RF
primary = read.csv('DEGanalysis_202406/randomForest_DEG/DEspeFC2Feature_13primary_disbiome.csv')[,1]
Level = ifelse(substr(primary,1,1)=='g', 'species', 'genus')
f = read.csv('DEGanalysis_202406/DEspe_7level_T1.csv', row.names = 1)[primary,]
Enrichment = ifelse(f$log2FC > 0, 0, 1) # 1: R中富集； 0：在NR中富集
avg_R = f$meanR; avg_NR = f$meanNR; pval = f$pwil
Name = primary

df = data.frame(Name, Level, Enrichment, avg_R, avg_NR, pval)
write.csv(df, '../202406fig/tables/DEspe13_RFfeature.csv', row.names = F)
</4> SPEcies features in RF

<5> DEG features in RF
deg = read.csv('DEGanalysis_202406/randomForest_DEG/DEgeneFeature_T2Dpathway.csv')[,3][1:12]
f = read.csv('DEGanalysis_202406/DEgene_4level_T1.csv', row.names = 1)[deg,]
Enrichment = ifelse(f$log2FC > 0, 0, 1) # 1: R中富集； 0：在NR中富集
avg_R = f$meanR; avg_NR = f$meanNR; pval = f$pwil
Name = deg
Description = f$label

df = data.frame(Name, Description, Enrichment, avg_R, avg_NR, pval)
write.csv(df, '../202406fig/tables/DEG12_RFfeature.csv', row.names = F)
</5> DEG features in RF


<6> speEdge features in RF
edge = read.csv('new_ssn_202406/speFeature_importance_score.csv')[,2][1:50]
nodes = t(matrix(unlist(str_split(edge, '[-+]{2}')), nrow=2))
Enrichment = ifelse(grepl('[+]{2}', edge), 1, 0)  #  1: R中富集； 0：在NR中富集

df = data.frame(Edge_node1 = nodes[,1], Edge_node2 = nodes[,2], Enrichment)
write.csv(df, '../202406fig/tables/speEdge50_RFfeature.csv', row.names = F)
</6> speEdge features in RF



<7> koEdge features in RF
edge = read.csv('new_ssn_202406/geneFeature_importance_score.csv')[,2][1:50]
nodes = t(matrix(unlist(str_split(edge, '[-+]{2}')), nrow=2))
Enrichment = ifelse(grepl('[+]{2}', edge), 1, 0)  #  1: R中富集； 0：在NR中富集

annot = read.delim2('raw_data/count result/KEGG/Relative/Unigenes.relative.ko.xls', 
                    row.names = 1, quote = '')[, 207]
n1_desc = annot[nodes[, 1], 1]
n2_desc = annot[nodes[, 2], 1]

df = data.frame(Edge_node1 = nodes[,1], Description_node1 = n1_desc,
                Edge_node2 = nodes[,2], Description_node2 = n2_desc, Enrichment)
write.csv(df, '../202406fig/tables/geneEdge50_RFfeature.csv', row.names = F)

</7> koEdge features in RF