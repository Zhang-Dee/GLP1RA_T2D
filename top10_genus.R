setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

<phylum-top10> 门属
setwd("D:/Zhang Di/linshibangong/linshiMonkey/杭高院/direct_data/")
library(reshape2)
library(ggplot2)
library(ggsci)
library(ggalluvial)
library(RColorBrewer)
#Samples R vs NR
f = read.csv('spe_7level_T1.csv', row.names = 1)
f = f[f$label=='p', grepl('[FGH]', colnames(f))]  # g; 只取给药组
colnames(f) = str_sub(colnames(f), 1, -4)
rownames(f) = str_sub(str_extract(rownames(f), ';[a-z]__.*'), 5) # 删掉前缀
f = f[!colnames(f) %in% test_id]  # 只取train set
top = function(df){ df = df[order(rowMeans(df), decreasing = T)[1:9], ]; df['Others',] = 1-colSums(df); df['Taxonomy']=rownames(df); df}

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
  theme(strip.text = element_text(size = 14),
        legend.title = element_text(size=14, face="bold"), 
        legend.text = element_text(size = 14), axis.title = element_text(size=12)) 

</phylum-top10> 门属
  
