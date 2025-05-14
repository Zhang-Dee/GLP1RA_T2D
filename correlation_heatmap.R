setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

<1>. spe_edge里面包含的物种和表型的相关性
(CHO,HDL,LDL,VLDL,NEFA,HBA1C,FBG,INS,TG,免疫细胞？)
library(stringr)
library(corrplot)
library(psych)

ls = c("GLU.csv", "INS.csv", "HOMA_IR.csv",  "HbA1c.csv", "TEI.csv", "CHO.csv", "HDL.csv", 
       "LDL.csv", "NEFA.csv", "TG.csv", "VLDL.csv")
p = data.frame(row.names = c(paste0('F', 1:10), paste0('G', 1:10), paste0('H', 1:10),
                             paste0('C', 1:10)))
for (fn in ls){
  pheno = read.csv(paste0('pheno_data/', fn), row.names = 2, check.names = F)
  if (sum(rownames(pheno)==rownames(p)) < 40) {return(fn)}
  colnames(pheno) = str_replace(colnames(pheno), '-', '.')
  pheno = pheno[c('Day_.7', "Day_29", "Day_57")]
  colnames(pheno) = paste0(str_sub(fn, 1, -5), c('_T1', '_T2', '_T3'))
  p = cbind(p, pheno)
}
dim(p); p[1:3,1:4]
p = p[grepl('T1', colnames(p))]


feat = read.csv('new_ssn_202406/speFeature_importance_score.csv')[, 2][1:50] %>%
  strsplit('[-+]{2}') %>% unlist() %>% unique()
f = read.delim('raw_data/count result/species/Relative/Unigenes.relative.s.xls', 
               row.names = 1, check.names = F)[feat, ]
colnames(f) = colnames(f)[c(2:ncol(f), 1)]
f = f[grepl('[FGH][0-9]+.T1', colnames(f))]
colnames(f) = str_sub(colnames(f), 1, -4)
rownames(f) = str_extract(rownames(f), 's__.*')

！step1-所有样本做相关
samples = colnames(f)
f[1:3, 1:4]; p[1:3, 1:4]
corr = corr.test(t(f), p[samples, ], method = 'spearman')$p

a = 1 - as.matrix(corr)

color = rev(COL2('RdBu', 20)); color = c(color[1], rep(color[2:10], each=2),  "#B31B2C") # 换这个颜色
corrplot(a,  col.lim = c(0,1), col = color,  #tl.pos = 'n',
         is.corr = F, method = 'circle', tl.cex=0.4, cl.ratio=0.4)   # legend 表示1-p颜色手动换一下


# 如果要设置circle大小相同，需要更改源代码成这段
if (method == "circle" && plotCI == "n") {
  symbols(Pos, add = TRUE, inches = FALSE, circles = asp_rescale_factor * 
            0.9 * rep(0.8, nrow(a)*ncol(a)) ^0.5/2, fg = col.border, bg = col.fill)
}

</1>. spe_edge里面包含的物种和表型的相关性

<2>. 换变化的相关
f = read.delim('raw_data/count result/species/Relative/Unigenes.relative.s.xls', 
               row.names = 1, check.names = F)[feat, ]
colnames(f) = colnames(f)[c(2:ncol(f), 1)]
t1 = f[grepl('[FGH][0-9]+.T1', colnames(f))]
t3 = f[grepl('[FGH][0-9]+.T3', colnames(f))][-1]
colnames(t1) = str_sub(colnames(t1), 1, -4); colnames(t3) = str_sub(colnames(t3), 1, -4)
rownames(t1) = str_extract(rownames(t1), 's__.*'); rownames(t3) = str_extract(rownames(t3), 's__.*')
assert colnames(t3) == colnames(t1)
f = (t3-t1)/t1

...p1 = p[grepl('T1', colnames(p))]; p3 = p[grepl('T3', colnames(p))]
p = (p3-p1)/p1

... samples = ... row38

这是换RNR的相关
f0 = f
f = f0[!colnames(f0) %in% test_id]
f = f[colnames(f) %in% names(rnr_)[rnr_=='R']]  
fnr = f0[c("F2", "F4", "F6", "G7", "G8", "H4", "H6", "H7", "H8", "H9")]
</2>. 换变化的相关


<3> RvsNRf = f0[!colnames(f0) %in% test_id]
fr = f0[colnames(f0) %in% names(rnr_)[rnr_=='R']]  
fnr = f0[c("F2", "F4", "F6", "G7", "G8", "H4", "H6", "H7", "H8", "H9")]

samples = colnames(fr)
corr_r = corr.test(t(fr), p[samples, ], method = 'spearman')
a = rowSums(corr_r$p <= 0.05); a = names(a)[a > 0]
corp_r = corr_r$p[a, ]; corr_r = corr_r$r[a, ]


samples = colnames(fnr)
corr_nr = corr.test(t(fnr), p[samples, ], method = 'spearman')
a = rowSums(corr_nr$p <= 0.05); a = names(a)[a > 0]
corp_nr = corr_nr$p[a, ]; corr_nr = corr_nr$r[a, ]

a = rownames(corp_r) %>% intersect(rownames(corp_nr))
corr_r = rbind(corr_r[a, ], corr_r[!rownames(corr_r) %in% a, ])
corp_r = corp_r[rownames(corr_r),]

corr_nr = rbind(corr_nr[a, ], corr_nr[!rownames(corr_nr) %in% a, ])
corp_nr = corp_nr[rownames(corr_nr),]

squre = function(df){df**2 * df/abs(df)}

作图
corrplot(squre(corr_r), is.corr = T, p.mat = corp_r, sig.level = c(0.001, 0.01, 0.05), 
         method = 'circle', insig = 'label_sig', # 区分显著不显著，可选项 p-value, blank, n, label_sig, pch
         number.cex = 0.8, pch.col = 'black', pch.cex = 1, 
         tl.cex = 0.6, title = 'T1-R'
         )

corrplot(corr_nr, is.corr = T, p.mat = corp_nr, sig.level = c(0.001, 0.01, 0.05), 
         method = 'circle', insig = 'label_sig', # 区分显著不显著，可选项 p-value, blank, n, label_sig, pch
         number.cex = 0.8, pch.col = 'black', pch.cex = 1, 
         tl.cex = 0.6, title = 'T1-NR',
)


</3> RvsNR