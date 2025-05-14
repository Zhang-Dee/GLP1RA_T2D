setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Animal_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

library(reshape2)
library(ggplot2)

# f = read.csv('BW.csv', row.names = 1); f = f[f$Group!='VEH', c(1,2,4:21)]  # Day-7 -> Day57 每周两次
# colnames(f)[3:20] = paste0('d', c(0, str_extract(colnames(f)[4:20], '[0-9]+')))
f = read.csv('HOMA_IR.csv', row.names = 1)
colnames(f)[3:ncol(f)] = paste0('d', str_extract(colnames(f)[3:ncol(f)], '[0-9]+'))
f = f[f$Group!='VEH', c(1:2,4:13)]  # Day-1 -> Day57 每周一次数据
f = f[-which(rownames(f) %in% test_id),]

library(ggalt)
dat = melt(f, id.vars = c('Sample_ID', 'Group'), variable.name = 'time', value.name = 'HOMA-IR')
ggplot(dat, aes(x=time, group=Sample_ID)) +
  geom_point(aes(y=bw, col=Group), show.legend = F)+   
  geom_xspline(aes(y=bw, col=Group), show.legend = F, size=0.8, alpha=0.5) +  # 平滑曲线
  #geom_line(aes(y=bw, col=Group), show.legend = F, size=0.8, alpha=0.5) +  # 直线
  theme_bw() + 
  scale_color_manual(values = c('#E4937C', '#E2B6AE', '#8BB39A')) + 
  facet_wrap(~Group, ncol = 1, scales = 'free_y')
