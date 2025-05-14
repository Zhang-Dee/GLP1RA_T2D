setwd('D:/GLP1RA/')
rnr = read.csv('new_ssn_202406/response_by_homaIR.csv')[c("Animal_ID", "Sample_ID", "Effect_this", "type")]  # R/NR based on HOMA-IR
test_id = rnr[rnr$type=='test', 'Sample_ID']
rnr_ = structure(rnr$Effect_this, names=rnr$Sample_ID) 
 
<1>. feature AUC
library(stringr)
library(pROC)
library(randomForest)

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
  n = ncol(dat)-1  
  sn = nrow(dat)  
  
  auc=c()
  dat_roc = c()
  dat_index = data.frame(row.names = colnames(dat)[-ncol(dat)])
  for (j in 1:sn){
    for(i in 1:sn){    
      fold_test <- dat[i,]   
      fold_train <- dat[-i,]   
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
  df = df[shuffle_sample, ]
  colnames(df)[1:ncol(df)-1] = paste0('e', 1:(ncol(df)-1))
  dat_train = df[!rownames(df) %in% test_id, ]
  dat_test =  df[rownames(df) %in% test_id, ]
  fold_test = dat_test     
 
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
val_3 = runTest(f3, seed = 1467)  # 50 gene_Edge  
val_4 = runTest(f4, seed = 50)

# visualization
main = 'validation-ROC curve'  # top10 gene feature ROC curve
plot(val_1, col = '#867bb9', add = F, lty = 1, lwd = 3, print.auc = T,  
     main=main, legacy.axes=T)   # T- 1-specificity; F - specificity
plot(val_2, col = '#9cd7bc', add=T, lty = 1, lwd = 3, print.auc=T, legacy.axes=T)
plot(val_3, col = '#fbb482', add=T, lty = 1, lwd = 3, print.auc=T, legacy.axes=T)
plot(val_4, col = '#db697a', add=T, lty = 1, lwd = 3, print.auc=T, legacy.axes=T)


par(mfrow = c(2,2))
plot(val_1, col = '#d94738', add = F, lty = 1, lwd = 3, print.auc = T,  
     main=main, legacy.axes=T, auc.polygon=T, auc.polygon.col='#fcdecc')  # T- 1-specificity; F - specificity
plot(val_2, col = '#d94738', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T,
     auc.polygon=T, auc.polygon.col='#fcdecc')
plot(val_3, col = '#d94738', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T,
     auc.polygon=T, auc.polygon.col='#fcdecc')
plot(val_4, col = '#d94738', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T, 
     auc.polygon=T, auc.polygon.col='#fcdecc')


t1 = runTrian(f1); t2 = runTrian(f2); t3 = runTrian(f3); t4 = runTrian(f4)
par(mfrow = c(2,2)); main = 'training-ROC curve'
plot(t1, col = '#5979a2', add = F, lty = 1, lwd = 3, print.auc = T, 
     main=main, legacy.axes=T, auc.polygon=T, auc.polygon.col='#dee4f0')  # T- 1-specificity; F - specificity
plot(t2, col = '#5979a2', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T,
     auc.polygon=T, auc.polygon.col='#dee4f0')
plot(t3, col = '#5979a2', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T,
     auc.polygon=T, auc.polygon.col='#dee4f0')
plot(t4, col = '#5979a2', add=F, lty = 1, lwd = 3, print.auc=T, legacy.axes=T, 
     auc.polygon=T, auc.polygon.col='#dee4f0')
</1>. feature AUC
