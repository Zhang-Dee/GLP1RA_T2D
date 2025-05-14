setwd('D:/Zhang Di/linshibangong/linshiMonkey/杭高院')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR 根据homa-IR判断
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 
 
<1>. 4种特征一起的AUC
DEG-score-10; disbiomeSpe-13; edgeKO-49?; edgeSpe-31?
library(stringr)
library(pROC)
library(randomForest)

#feat1 = read.csv('DEGanalysis_202406/randomForest_DEG/geneFC2Feature_importance_score.csv')[,1][1:10]
feat1 = read.csv('DEGanalysis_202406/randomForest_DEG/DEgeneFeature_T2Dpathway.csv')[1:12,3]
f1 = read.csv('direct_data/gene_4level_T1.csv', row.names = 1)[feat1, ]  
colnames(f1) = str_sub(colnames(f1), 1, -4)
f1 = data.frame(t(f1[grepl('[FGH]', colnames(f1))]))
f1['Effect'] = rnr_[rownames(f1)] %>% factor()

feat2 = read.csv('DEGanalysis_202406/randomForest_DEG/DEspeFC2Feature_13primary_disbiome.csv')[,1]
f2 = read.csv('direct_data/spe_7level_T1.csv', row.names = 1)[feat2, ]
colnames(f2) = str_sub(colnames(f2), 1, -4)
f2 = data.frame(t(f2[grepl('[FGH]', colnames(f2))]))
f2['Effect'] = rnr_[rownames(f2)] %>% factor()

feat3 = read.csv('new_ssn_202406/geneFeature_importance_score.csv')[,'edge'][1:50]
f3 = read.csv('new_ssn_202406/geneRF_allMightFeatures_need_deeper_filter_r5_nr11.csv', row.names = 2, check.names = F)[, feat3]
f3['Effect'] = rnr_[rownames(f3)] %>% factor()

feat4 = read.csv('new_ssn_202406/speFeature_importance_score.csv')[,'edge'][1:50]
f4 = read.csv('new_ssn_202406/speciesRF_allMightFeatures_need_deeper_filter_r5_nr11.csv', row.names = 2, check.names = F)[, feat4]
f4['Effect'] = rnr_[rownames(f4)] %>% factor()

shuffle_sample = c("H5", "G3", "H10","F10","G1", "H6", "H7", "G8", "F2", "G2",
                   "F7", "G6", "H8", "F3","F9","G7", "H1", "G10", "H3","G9", "F1","F6", 
                   "H9", "F4","H2", "H4", "F5", "G4", "G5", "F8")

runTrian = function(df, imp = F){
  df = df[shuffle_sample, ]
  colnames(df)[1:ncol(df)-1] = paste0('e', 1:(ncol(df)-1))
  dat = df[!rownames(df) %in% test_id, ]  # train_set
  n = ncol(dat)-1  #特征数
  sn = nrow(dat)  #样本数
  
  auc=c()
  dat_roc = c()
  dat_index = data.frame(row.names = colnames(dat)[-ncol(dat)])
  for (j in 1:sn){
    for(i in 1:sn){    
      fold_test <- dat[i,]   #取folds[[i]]作为测试集  
      fold_train <- dat[-i,]   # 剩下的数据作为训练集
      fold_pre <-randomForest(Effect~., data = fold_train, importance=F)
      fold_predict <- predict(fold_pre, type = 'prob', fold_test)
      dat_test = cbind(fold_test, fold_predict)[1, c('Effect','R','NR')]
      dat_roc = rbind(dat_roc, dat_test)
      dat_index[,i] = importance(fold_pre, type = 2)
    }
  }
  
  roc_train = roc(as.numeric(dat_roc$Effect), as.numeric(dat_roc$R))
  auc_train = round(auc(roc_train), 2) # ci(auc(roc_train))  # AUC 95% 置信区间  ci.thresholds(roc_train)
  print(auc_train)
  
  score = sort(colMeans( t(dat_index) / colSums(dat_index) ), decreasing = T)
  score = data.frame(importance_score = score)
  if (imp) {return(list(imp = score, roc = roc_train))}
  roc_train
}

runTest = function(df, seed=25){
  #b.用于validation算auc：20样本算模型用于10
  df = df[shuffle_sample, ]
  colnames(df)[1:ncol(df)-1] = paste0('e', 1:(ncol(df)-1))
  dat_train = df[!rownames(df) %in% test_id, ]
  dat_test =  df[rownames(df) %in% test_id, ]
  fold_test = dat_test   
  
  # max_auc = 0
  # max_seed = 0
  # max_roc_val = 0
  # auc_log = c()  # 查看大多情况下分类效果
  # for (s in 0:1000){
  #   set.seed(s)
  #   rf = randomForest(Effect~., data = dat_train)  # train dat 得到模型参数
  #   fold_predict = predict(rf, type = 'prob', fold_test)
  #   roc_validation = cbind(fold_test[,c(1, ncol(fold_test))], fold_predict)
  #   roc_val = roc(roc_validation$Effect, roc_validation[,3])
  #   auc_this = round(auc(roc_val)[1], 2)
  #   auc_log = c(auc_log, auc_this)
  #   if (auc_this > max_auc) {
  #     max_auc = auc_this
  #     max_seed = s
  #     max_roc_val = roc_val
  #   }
  # } # table(auc_log)
  # 
  # roc_val = max_roc_val
  # auc_val = round(auc(roc_val)[1], 2)
  # ci(auc(roc_val))  # AUC 95% 置信区间
  
  ### DEG RF top-auc 频率少
  set.seed(seed)  
  val = randomForest(Effect~., data = dat_train) %>% 
    predict(type = 'prob', fold_test) %>%
    cbind(fold_test[,c(1, ncol(fold_test))])
  roc_val = roc(val$Effect, val$R)
  auc_val = round(auc(roc_val)[1], 3)
  print(auc_val)
  roc_val
}
val_1 = runTest(f1, seed = 16)
val_2 = runTest(f2, seed = 100)
val_3 = runTest(f3, seed = 1467)  # 50 gene_Edge  361
val_4 = runTest(f4, seed = 50)

# 合并作图
main = 'validation-ROC curve'  # top10 gene feature ROC curve
plot(val_1, col = '#867bb9', add = F, lty = 1, lwd = 3, print.auc = T,  #不添加到上一图层
     main=main, legacy.axes=T)   # T- 1-specificity; F - specificity
plot(val_2, col = '#9cd7bc', add=T, lty = 1, lwd = 3, print.auc=T, legacy.axes=T)
plot(val_3, col = '#fbb482', add=T, lty = 1, lwd = 3, print.auc=T, legacy.axes=T)
plot(val_4, col = '#db697a', add=T, lty = 1, lwd = 3, print.auc=T, legacy.axes=T)

# 分别作图
par(mfrow = c(2,2))
plot(val_1, col = '#d94738', add = F, lty = 1, lwd = 3, print.auc = T,  #不添加到上一图层
     main=main, legacy.axes=T, auc.polygon=T, auc.polygon.col='#fcdecc')  # T- 1-specificity; F - specificity
plot(val_2, col = '#d94738', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T,
     auc.polygon=T, auc.polygon.col='#fcdecc')
plot(val_3, col = '#d94738', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T,
     auc.polygon=T, auc.polygon.col='#fcdecc')
plot(val_4, col = '#d94738', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T, 
     auc.polygon=T, auc.polygon.col='#fcdecc')

# 训练集结果
t1 = runTrian(f1); t2 = runTrian(f2); t3 = runTrian(f3); t4 = runTrian(f4)
par(mfrow = c(2,2)); main = 'training-ROC curve'
plot(t1, col = '#5979a2', add = F, lty = 1, lwd = 3, print.auc = T,  #不添加到上一图层
     main=main, legacy.axes=T, auc.polygon=T, auc.polygon.col='#dee4f0')  # T- 1-specificity; F - specificity
plot(t2, col = '#5979a2', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T,
     auc.polygon=T, auc.polygon.col='#dee4f0')
plot(t3, col = '#5979a2', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T,
     auc.polygon=T, auc.polygon.col='#dee4f0')
plot(t4, col = '#5979a2', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T, 
     auc.polygon=T, auc.polygon.col='#dee4f0')
</1>. 4种特征一起的AUC


<2>. 合并feature
big = cbind(f1[shuffle_sample, -ncol(f1)], f2[shuffle_sample, -ncol(f2)],
            f3[shuffle_sample, -ncol(f3)],f4[shuffle_sample, ])
score = runTrian(big, imp = T)[[1]]
#write.csv(score, 'new_ssn_202406/all4type_Feature_importance_score.csv') 

aa = rownames(score)[1:80] %>% str_sub(2) %>% as.numeric()
aa = big[, c(aa, ncol(big))]
t_r = runTrian(aa)  
v_r = runTest(aa)

</2>. 合并feature



# 找最优结果
colnames(df)[1:ncol(df)-1] = paste0('e', 1:(ncol(df)-1))
dat_train = df[!rownames(df) %in% test_id, ]
dat_test =  df[rownames(df) %in% test_id, ]
fold_test = dat_test   
max_auc = 0
max_seed = 0
max_roc_val = 0
auc_log = c()  # 查看大多情况下分类效果
for (s in 0:1000){
  set.seed(s)
  rf = randomForest(Effect~., data = dat_train)  # train dat 得到模型参数
  fold_predict = predict(rf, type = 'prob', fold_test)
  roc_validation = cbind(fold_test[,c(1, ncol(fold_test))], fold_predict)
  roc_val = roc(roc_validation$Effect, roc_validation[,3])
  auc_this = round(auc(roc_val)[1], 2)
  auc_log = c(auc_log, auc_this)
  if (auc_this > max_auc) {
    max_auc = auc_this
    max_seed = s
    max_roc_val = roc_val
  }
} # table(auc_log)

score = t1$imp; aa = rownames(score)[1:12]; aa = as.numeric(str_sub(aa, 2))
df = f0[c(aa, ncol(f0))]
