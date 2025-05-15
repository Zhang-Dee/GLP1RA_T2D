setwd('D:/GLP1RA/')
rnr = read.csv('response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR based on homa-IR improvement
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 


library(stringr)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsci)
library(ggalluvial)
library(RColorBrewer)
f3 = read.csv('spe_7level_T3.csv', row.names = 1)
f3 = f3[grepl('[FGH]', colnames(f3))]
colnames(f3) = str_sub(colnames(f3), 1, -4)
f3 = f3[!colnames(f3) %in% test_id ] 


f1 = read.csv('spe_7level_T1.csv', row.names = 1)  # spe_7level_T1.csv
label = f1[c('tax', 'label')]
f1 = f1[grepl('[FGH]', colnames(f1))]
colnames(f1) = str_sub(colnames(f1), 1, -4)
f1 = f1[!colnames(f1) %in% test_id ]   

coT1T3 = intersect(rownames(f1), rownames(f3))
f1 = f1[coT1T3, ]; f3 = f3[coT1T3, ]

top = function(df){ df = df[order(rowMeans(df), decreasing = T)[1:9], ]; df['Others',] = 1-colSums(df); df['Taxonomy']=rownames(df); df}


f = f3  # f3
f = f[grepl(';p__', rownames(f)), ] 
rownames(f) = str_sub(str_extract(rownames(f), ';[a-z]__.*'), 2)
top_f = top(f)
dat = melt(top_f, id.vars = 'Taxonomy', variable.name = 'Sample')
dat$Sample = as.character(dat$Sample); dat['Response'] = rnr_[dat$Sample]
dat$Sample = factor(dat$Sample, levels = names(sort(rnr_[colnames(f)])))
dat$Taxonomy = factor(dat$Taxonomy, levels=rev(rownames(top_f)))

color = rev(brewer.pal(10, 'Paired'))
ggplot(dat, aes(x=Sample, y=100 * value, alluvium = Taxonomy, stratum = Taxonomy)) + 
  geom_alluvium(aes(fill = Taxonomy), width = 0) + 
  geom_stratum(aes(fill = Taxonomy), width = 0) + 
  scale_fill_manual(values = color) + 
  theme_classic() + 
  scale_y_continuous(expand = c(0, 0)) + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_wrap(~Response, scales = 'free_x', ncol = 5) + theme(strip.text = element_text(size = 14)) + 
  labs(x = '', y = 'Relative Abundance(%)')


# average-value R vs NR
dat = data.frame(cbind(rowMeans(top_f[colnames(top_f) %in% names(rnr_)[rnr_=='R']]),
                       rowMeans(top_f[colnames(top_f) %in% names(rnr_)[rnr_=='NR']])))
colnames(dat) = c('R', 'NR')
dat['Taxonomy'] = rownames(dat)

dat = melt(dat, id.vars = 'Taxonomy', variable.name = 'Response')
dat$Taxonomy = factor(dat$Taxonomy, levels=rev(rownames(top_f)))
ggplot(dat, aes(x=Response, y=100 * value, alluvium = Taxonomy, stratum = Taxonomy)) + 
  geom_alluvium(aes(fill = Taxonomy), width = 0.4, alpha=0.4) + 
  geom_stratum(aes(fill = Taxonomy), width = 0.4, alpha=0.6) + 
  scale_fill_manual(values =color) +
  labs(x = '', y = 'Relative Abundance(%)') + 
  coord_fixed(ratio = 0.04) +
  theme_classic() 



############################################# human 16S top10 genera ########################

hsa = read.csv('PRJNA665123/filter_mito/ra_level_6.csv', row.names = 1, check.names = F)
rownames(hsa) = str_match(rownames(hsa), 'g__[\\S]+')
hsa = hsa[rownames(hsa) != 'g__Unclassified', ]
sid = read.csv('../../PRJNA665123/filter_mito/pheno.csv', row.names=1)
hsa_t0 = hsa[, colnames(hsa) %in% sid$L0_SRR]
hsa_t3 = hsa[, colnames(hsa) %in% sid$L4_SRR]

# average-value of top10
hsa_f = hsa_t3
top_f = top(hsa_f)


dat = data.frame(rowMeans(top_f[-ncol(top_f)]))
colnames(dat) = 'value'
dat['Taxonomy'] = rownames(dat)
dat['time'] = 'baseline'

dat$Taxonomy = factor(dat$Taxonomy, levels=rev(dat$Taxonomy))
ggplot(dat, aes(x=time, y=100 * value, alluvium = Taxonomy, stratum = Taxonomy)) + 
  geom_alluvium(aes(fill = Taxonomy), width = 0.4, alpha=0.4) + 
  geom_stratum(aes(fill = Taxonomy), width = 0.4, alpha=0.6) + 
  scale_fill_manual(values =color) +
  labs(x = '', y = 'Relative Abundance(%)') + 
  coord_fixed(ratio = 0.04) +
  theme_classic() 



