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
library(emmeans)

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
  
  # R-NR-all 前后(T2) paired-ttest
  r_t12 = t.test(df[df$Response=='R', 'T1'], df[df$Response=='R', 'T2'], paired = T)$p.value
  nr_t12 = t.test(df[df$Response=='NR', 'T1'], df[df$Response=='NR', 'T2'], paired = T)$p.value
  overall_t12 = t.test(df[, 'T1'], df[, 'T2'], paired = T)$p.value
  
  # R-NR-all 前后(T3) Repeated Measures ANOVA
  dat = melt(df, id.vars = c('Animal_ID', 'Response'), value.name = 'FBG', variable.name = 'time')
  dat_r = dat[dat$Response == 'R', ]; dat_nr = dat[dat$Response == 'NR', ]
  repeat_r <- aov(FBG ~ time + Error(Animal_ID/time), data = dat_r) %>% summary(); repeat_r_p <- repeat_r[[length(repeat_r)]][[1]]$`Pr(>F)`[1]
  repeat_nr <- aov(FBG ~ time + Error(Animal_ID/time), data = dat_nr) %>% summary(); repeat_nr_p <- repeat_nr[[length(repeat_nr)]][[1]]$`Pr(>F)`[1]
  repeat_all <- aov(FBG ~ time + Error(Animal_ID/time), data = dat) %>% summary(); repeat_all_p <- repeat_all[[length(repeat_all)]][[1]]$`Pr(>F)`[1]
  
  # RvsNR-T1: independentStudent’sttest.
  vs_t1 = t.test(df[df$Response=='R', 'T1'], df[df$Response=='NR', 'T1'], paired = F)$p.value
  
  # RvsNR-T2: ANCOVA controlling for baseline measurements 协方差分析
  vs_t2 = aov(T2 ~ Response + T1, data = df) %>% summary(); vs_t2 = vs_t2[[1]]['Response', 'Pr(>F)']
  
  #RvsNR-T3-a: ANCOVA controlling for baseline measurements 不必考虑T2的数据(除非你希望分析它们的影响——GPT4o)
  vs_t3 = aov(T3 ~ Response + T1, data = df) %>% summary(); vs_t3 = vs_t3[[1]]['Response', 'Pr(>F)']
  
  #RvsNR-T3-b: 双因素方差分析,T3比较最终用这个,更好地处理每个样本的重复测量特性。
  data_long <- melt(df, id.vars = c("Animal_ID", "Response", "T1"),
                    measure.vars = c("T2", "T3"),
                    variable.name = "time", value.name = "FBG")
  
  ## rm_anova_model <- aov(FBG ~ Response * time + Error(Animal_ID/time), data = data_long)  # 不控制基线
  ## aov_ = summary(rm_anova_model)[[3]][[1]]     # 不控制基线
  rm_ancova_model <- aov(FBG ~ Response * time + T1 + Error(Animal_ID/time), data = data_long)  # 同时控制基线数据
  aov_ = summary(rm_ancova_model)[[3]][[1]]
  if (aov_['Response:time', 'Pr(>F)'] > 0.05){rp = aov_['Response ', 'Pr(>F)'] } else {
    emmeans_results <- emmeans(rm_ancova_model, ~ Response | time)     # rm_anova_model
    rp =  pairs(emmeans_results, adjust = "tukey") %>% summary() # 查看两周后的组间比较
    rp = rp[rp$time=='T3', 'p.value']
  }
  vs_t3 = rp
  
  # R-NR-all的 mu±sem
  mu = aggregate(df[,3:5], by=list(df$Response), FUN = 'mean')
  sem = aggregate(df[,3:5], by=list(df$Response), FUN = function(x){sd(x)/sqrt(length(x))})
  musem = data.frame(NR = paste(round(mu[1,2:4],2), round(sem[1,2:4],2), sep=' ± '),
                     R = paste(round(mu[2,2:4],2), round(sem[2,2:4],2), sep=' ± '),
                     all = paste(round(colMeans(df[3:5]), 2), 
                                 round(apply(df[3:5], 2, function(x){sd(x)/sqrt(length(x))}), 2), 
                                 sep=' ± '),
                     p_nr = c(NaN, nr_t12, repeat_nr_p),
                     p_r = c(NaN, r_t12, repeat_r_p),
                     p_all = c(NaN, overall_t12, repeat_all_p),
                     pRNR = c(vs_t1, vs_t2, vs_t3),
                     pheno = fn,
                     row.names = c('T1','T2','T3'))
  
  muSEM = rbind(muSEM, musem)
}

muSEM['time'] = substr(rownames(muSEM),1,2)
muSEM = muSEM[c(ncol(muSEM), 1:(ncol(muSEM)-1))]
muSEM['note'] = c('p_r,p_nr,p_all use pairedTtest for T2; repeated ANOVA for T3', 
                  'pRNR use ANCOVA controlling baseline for T2; two-way ANOVA for T3',
                  rep(NULL, ncol(muSEM)-2))
muSEM['note'] = c('p_r,p_nr,p_all use pairedTtest for T2; repeated ANOVA for T3', 
                  'pRNR use ANCOVA controlling baseline for T2; two-way ANOVA controlling baseline for T3',
                  rep(NULL, ncol(muSEM)-2))
write.csv(muSEM, 'pheno_table_train_new.csv', row.names = F, fileEncoding = 'GBK') # pheno_table_train





# 制作三种Cpd的表
muSEM = c()
for (fn in ls){
  f =  read.csv(paste0(fn,'.csv'))
  f = f[- which(f$Animal_ID %in% test_id), ] # 是否纳入测试集的区别
  print(table(f$Group))
  
  df = f[c("Animal_ID", "Group", "Day_.7", "Day_29", "Day_57")]
  colnames(df)[3:5] = c('T1', 'T2', 'T3')
  df['relative_change'] = (df$T3 - df$T1) / df$T1
  df$Group = factor(df$Group)
  
  # cpdA-cpdC-sar17-veh 前后(T2) paired-ttest
  cpdA_12 = t.test(df[df$Group=='CpdA', 'T1'], df[df$Group=='CpdA', 'T2'], paired = T)$p.value
  cpdC_12 = t.test(df[df$Group=='CpdC', 'T1'], df[df$Group=='CpdC', 'T2'], paired = T)$p.value
  sar17_12 = t.test(df[df$Group=='SAR17', 'T1'], df[df$Group=='SAR17', 'T2'], paired = T)$p.value
  veh_12 = t.test(df[df$Group=='VEH', 'T1'], df[df$Group=='VEH', 'T2'], paired = T)$p.value
  
  
  # cpdA-cpdC-sar17-veh 前后(T3) Repeated Measures ANOVA
  cpd_13 = c(); cpd_veh_t123 = list()
  dat = melt(df[-ncol(df)], id.vars = c('Animal_ID', 'Group'), value.name = 'FBG', variable.name = 'time')
  for (g in c('CpdA', 'CpdC', 'SAR17', 'VEH')){
    dat_cpd = dat[dat$Group == g, ]
    repeat_cpd <- aov(FBG ~ time + Error(Animal_ID/time), data = dat_cpd) %>% summary()
    repeat_cpd_p <- repeat_cpd[[length(repeat_cpd)]][[1]]$`Pr(>F)`[1]
    cpd_13 = structure(c(cpd_13, repeat_cpd_p), names=c(names(cpd_13), g))
    
    if (g=='VEH'){next}
    # CpdvsVEH-T1: independentStudent’sttest.
    vs_t1 = t.test(df[df$Group == g, 'T1'], df[df$Group == 'VEH', 'T1'], paired = F)$p.value
    
    # CpdvsVEH-T2: ANCOVA controlling for baseline measurements 协方差分析
    vs_t2 = aov(T2 ~ Group + T1, data = df[df$Group %in% c(g, 'VEH'), ]) %>% summary(); vs_t2 = vs_t2[[1]]['Group', 'Pr(>F)']
    
    # pdvsVEH-T3-a: ANCOVA controlling for baseline measurements 不必考虑T2的数据(除非你希望分析它们的影响——GPT4o)
    vs_t3 = aov(T3 ~ Group + T1, data = df[df$Group %in% c(g, 'VEH'), ]) %>% summary(); vs_t3 = vs_t3[[1]]['Group', 'Pr(>F)']
    
    #RvsNR-T3-b: 双因素方差分析,T3比较最终用这个,更好地处理每个样本的重复测量特性。不用控制T1因为T1无显著差异。
    data_long <- melt(df[df$Group %in% c(g, 'VEH'), -ncol(df)], id.vars = c("Animal_ID", "Group", "T1"),
                      measure.vars = c("T2", "T3"),
                      variable.name = "time", value.name = "FBG")
    
    rm_anova_model <- aov(FBG ~ Group * time + Error(Animal_ID/time), data = data_long) 
    aov_ = summary(rm_anova_model)[[3]][[1]]
    if (aov_['Group:time', 'Pr(>F)'] > 0.05){rp = aov_['Group ', 'Pr(>F)'] } else {
      emmeans_results <- emmeans(rm_anova_model, ~ Group | time)
      rp =  pairs(emmeans_results, adjust = "tukey") %>% summary() # 查看两周后的组间比较
      rp = rp[rp$time=='T3', 'p.value']
    }
    vs_t3 = rp
    
    cpd_veh_t123[[g]] = c(vs_t1, vs_t2, vs_t3)
  }
  
  # cpd-VEH的 mu±sem
  mu = aggregate(df[,3:5], by=list(df$Group), FUN = 'mean')
  sem = aggregate(df[,3:5], by=list(df$Group), FUN = function(x){sd(x)/sqrt(length(x))})
  musem = data.frame(CpdA = paste(round(mu[1,2:4],2), round(sem[1,2:4],2), sep=' ± '),
                     CpdC = paste(round(mu[2,2:4],2), round(sem[2,2:4],2), sep=' ± '),
                     SAR17 = paste(round(mu[3,2:4],2), round(sem[3,2:4],2), sep=' ± '),
                     VEH = paste(round(mu[4,2:4],2), round(sem[4,2:4],2), sep=' ± '),
                     
                     p_cpdA_t1_23 = c(NaN, cpdA_12, cpd_13['CpdA']),
                     p_cpdC_t1_23 = c(NaN, cpdC_12, cpd_13['CpdC']),
                     p_sar_t1_23 = c(NaN, sar17_12, cpd_13['SAR17']),
                     p_veh_t1_23 = c(NaN, veh_12, cpd_13['VEH']),
                     
                     p_veh_cpdA_t123 = cpd_veh_t123[['CpdA']],
                     p_veh_cpdC_t123 = cpd_veh_t123[['CpdC']],
                     p_veh_sar17_t123 = cpd_veh_t123[['SAR17']],
                     
                     pheno = fn,
                     row.names = c('T1','T2','T3'))
  
  muSEM = rbind(muSEM, musem)
}

muSEM['time'] = substr(rownames(muSEM),1,2)
muSEM = muSEM[c(ncol(muSEM), 1:(ncol(muSEM)-1))]
muSEM['note'] = c('p_t1_23 use pairedTtest for T2; repeated ANOVA for T3', 
                  'p_veh_cpd use ANCOVA controlling baseline for T2; two-way ANOVA for T3',
                  'n_CpdA=6; n_CpdC=8; n_SAR17=6; n_VEH=10',
                  'CpdA= liraglutide; CpdC= Sanofi dual GLP-1R/GCGR agonist; SAR17= semaglutide',
                  rep(NA, nrow(muSEM)-4))
write.csv(muSEM, 'Cpd_table.csv', row.names = F, fileEncoding = 'GBK') # pheno_table_train




# 整理所有样本的原始数据："Day_.7", "Day_29", "Day_57"
library(stringr)
for (fn in ls){
  f =  read.csv(paste0(fn,'.csv'))
  f['Response'] = rnr_[f$Sample_ID]
  f['Batch'] = ifelse(f$Animal_ID %in% test_id, 'validation', 'discovery')
  
  df = f[c("Animal_ID", "Group", "Batch", "Response", "Day_.7", "Day_29", "Day_57")]
  colnames(df)[5:7] = paste(c('T1', 'T2', 'T3'), fn, sep='_')
  
  if (fn==ls[1]){data = df; next}
  data = merge(data, df, by = c('Animal_ID', 'Group','Batch', 'Response'))
}

colnames(data) = str_replace(colnames(data), 'GLU', 'FBG (mg/dL)') %>% 
  str_replace('BW', 'BW (kg)') %>%
  str_replace('INS', 'FINS (µU/mL)') %>%
  str_replace('HOMA_IR', 'HOMA-IR') %>%
  str_replace('HbA1c', 'HbA1c (%)') %>%
  str_replace('LDL', 'LDL (mmol/L)') %>%
  str_replace('HDL', 'HDL (mmol/L)') %>%
  str_replace('NEFA', 'NEFA (mmol/L)') %>%
  str_replace('TG', 'TG (mmol/L)') %>%
  str_replace('CHO', 'CHO (mmol/L)') %>%
  str_replace('TEI', 'TEI (Kcal/d)')
  
  
data$Response = gsub('NR', 'LR', data$Response)
data = data[order(data$Batch, data$Group), ]
write.csv(data, 'data_tidy.csv', row.names = F) # pheno_table_train


















## 把三种药物作为协变量做混合线性模型的话
library(lme4)
a = melt(df, id.vars = c('Animal_ID', 'Response'), variable.name = 'Time', value.name = 'FBG')
a$Animal_ID = as.character(a$Animal_ID)
a['Drug_Type'] = f[a$Animal_ID, 'Group']
lmm_model <- lmer(FBG ~ Time * (Response + Drug_Type) + (1 | Animal_ID), data = a)
summary(lmm_model)

r = Anova(lmm_model, type = 3)
print(r)

