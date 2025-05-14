setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Animal_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 

setwd('pheno_data/')
library(ggpubr)
library(effects) # 消除协变量影响
library(multcomp)
library(rstatix)
library(reshape2)
ls = c('BW', 'GLU', 'HbA1c', 'HOMA_IR', 'INS', 'LDL', 'NEFA', 'TEI', 'TG', 'VLDL', 'CHO', 'HDL')

muSEM = c()
for (fn in ls){
  f =  read.csv(paste0(fn,'.csv'))
  f = f[f$Group != 'VEH', ]
  f['Response'] = rnr_[f$Sample_ID]
  f = f[- which(f$Animal_ID %in% test_id), ] # 是否纳入测试集的区别
  print(table(f$Response))
  
  df = f[c("Animal_ID", "Response", "Day_.7", "Day_29", "Day_57")]
  colnames(df)[3:5] = c('T1', 'T2', 'T3')
  df['relative_change'] = (df$T3 - df$T1) / df$T1
  df$Response = factor(df$Response)
  
  # R-NR 前后paired-ttest
  r_t13 = t.test(df[df$Response=='R', 'T1'], df[df$Response=='R', 'T3'], paired = T)$p.value
  nr_t13 = t.test(df[df$Response=='NR', 'T1'], df[df$Response=='NR', 'T3'], paired = T)$p.value
  overall_t13 = t.test(df[, 'T1'], df[, 'T3'], paired = T)$p.value
  
  rm_covarP = function(mat){
    p = aov(T3~T1+Response, data = mat) %>%  #T1_Glu是协变量,写在分组前  summary(fit)
      glht(linfct = mcp(Response = 'Tukey')) %>%
      summary()
    p$test$pvalues[1]
  }
  
  way2repeated = function(mat){
    mat = melt(mat[-6], id.vars = c('Animal_ID', 'Response'), value.name = 'FBG', variable.name = 'time')
    aov_ = aov(FBG ~ time * Response + Error(Animal_ID/time), data = mat) %>% summary()
    aov_ = aov_[[3]][[1]]
    if (aov_['time:Response', 'Pr(>F)'] > 0.05){
      rp = aov_['Response ', 'Pr(>F)'] } 
    else{
      simpleEff = group_by(mat, time) %>% anova_test(dv = FBG, wid = Animal_ID, between = Response) %>% get_anova_table()
      rp = simpleEff[simpleEff$time=='T3', ]$p }
    rp
  }
  
  
  mu = aggregate(df[,3:5], by=list(df$Response), FUN = 'mean')
  sem = aggregate(df[,3:5], by=list(df$Response), FUN = function(x){sd(x)/sqrt(length(x))})
  musem = data.frame(NR = paste(round(mu[1,2:4],2), round(sem[1,2:4],2), sep=' ± '),
                     R = paste(round(mu[2,2:4],2), round(sem[2,2:4],2), sep=' ± '),
                     p_t1 = compare_means(T1~Response, data =df, method = 't.test', paired = F)$p,
                     p_t3_effect = rm_covarP(df),
                     p_t3_anova = way2repeated(df),
                     p_t3_relaCh_ttest = compare_means(relative_change~Response, data = df, method = 't.test', paired = F)$p,
                     p_t3_relaCh_wilc = compare_means(relative_change~Response, data = df, method = 'wilcox.test', paired = F)$p,
                     pheno = fn,
                     row.names = c('T1','T2','T3'))
  musem['R_T1vsT3'] = r_t13; musem['NR_T1vsT3'] = nr_t13
  
  musem['all_T1vsT3'] = overall_t13
  musem['all'] = paste(round(colMeans(df[c('T1','T2','T3')]), 2), 
                          apply(df[c('T1','T2','T3')], MARGIN = 2, function(x){round(sd(x)/sqrt(length(x)), 2)}),
                          sep=' ± ')
  
  
  muSEM = rbind(muSEM, musem)
}

muSEM['time'] = substr(rownames(muSEM),1,2)
muSEM = muSEM[c(ncol(muSEM), 1:(ncol(muSEM)-1))]
write.csv(muSEM, 'pheno_table_train.csv', row.names = F, fileEncoding = 'GBK') # pheno_table_train




