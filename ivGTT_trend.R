setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
library(ggplot2)
library(reshape2)

f = read.csv('ivGTT_T0.csv', row.names = 2)
f$sample = paste(f$group, 1:40, sep='_')
a = f[c("group", "BS", "X1min", "X3min", "X5min", "X10min", "X20min", "X40min", "X60min")]

mu = aggregate(a[-1], by=list(a$group), FUN = 'mean')
sem = aggregate(a[-1], by=list(a$group), FUN = function(x){sd(x)/sqrt(length(x))})

aa = melt(mu, id.vars = 'Group.1', value.name = 'mean')
aa['sem'] = melt(sem, id.vars = 'Group.1')$value

aa['min'] = structure(c(0,0.5,1.8,3,5,8,12,18), names=levels(aa$variable))[aa$variable]
ggplot(aa, aes(min, mean, color=Group.1)) +
    geom_line(linewidth=0.9, alpha=0.3) +
    geom_point(cex=2, alpha=0.5) +
    theme_classic() +
    xlab(NULL) + ylab('Plasma Glucose (mg/dL)') + ggtitle('ivGTT before treatment') +
    geom_errorbar(aes(ymin=mean-sem, ymax=mean+sem),
                  width=.8, position=position_dodge(0.05)) +
    scale_x_continuous(breaks = c(0.5,1.8,3,5,8,12,18),
                       labels = c('1min','3min','5min','10min','20min','40min','60min')) #+
    #geom_area(aes(fill=Group.1), linewidth=0.01, position = 'identity',color='white', alpha=0.3)


########################### R vs. NR
rnr = read.csv('R_NR.csv', row.names = 1)
a['response'] = rnr[rownames(a), 'Effect']

mu = aggregate(a[-c(1,10)], by=list(a$response), FUN = 'mean')
sem = aggregate(a[-c(1,10)], by=list(a$response), FUN = function(x){sd(x)/sqrt(length(x))})
...


