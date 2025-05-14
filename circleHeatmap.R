setwd('D:/Zhang Di/linshibangong/linshiMonkey')
f = read.csv('liverFunc.csv')

# 加载包

library(devtools)
library(circlize)
library(tidyverse)
library(ComplexHeatmap)
library(gridBase)
#生成KEGG数据矩阵（矩阵1）
data1<-matrix(rnorm(670,mean=0.5),nrow=67)
rownames(data1)<-c("K01446","K01971","K01142","K01151","K01246","K00784","K02031","K01644","K02037","K02065","K01448",
                   "K01890","K00266","K01725","K00806","K00231","K01737","K00858","K00019","K01715","K01692","K00249",
                   "K00023","K00626","K00101","K00803","K01710","K01791","K01176","K00799","K00800","K01667","K01668",
                   "K01712","K00053","K01696","K01697","K00108","K00639","K01489","K00226","K01488","K02339","K01428",
                   "K01438","K02124","K02275","K01796","K00632","K00648","K00849","K01805","K01685","K00065","K00090",
                   "K01619","K01834","K00121","K02182","K02082","K02005","K01266","K01990","K01463","K02217","K01174","K02003")
colnames(data1) <- paste0('S_',seq(1:10))
#生成EC数据矩阵（矩阵2）
data2<-matrix(rnorm(160,mean=1),nrow = 16)
rownames(data2)<-c('Ribosome','DNA replication',  # Genetic Information Processing
                   'Notch signaling pathway', 'Ras signaling pathway', # Environmental Information Processing
                   'Purine metabolism','Biotin metabolism','Pyrimidine metabolism','Starch and sucrose metabolism','Glycolysis / Gluconeogenesis', # Metabolism
                   'Focal adhesion', 'Motor proteins', #Cellular Processes
                   'Type II diabetes mellitus', 'Type I diabetes mellitus', #Human Diseases
                   'Circadian rhythm', 'Insulin secretion','PPAR signaling pathway') # Organismal Systems
colnames(data2) <- paste0('S_',seq(1:10))
#生成细菌数据矩阵（矩阵3）
data3<-matrix(rnorm(60,mean=1),nrow = 6)
rownames(data3)<-c("Metabolism", 'Genetic Information Processing', 'Environmental Information Processing',
                   'Cellular Processes', 'Human Diseases', 'Organismal Systems')
#将三个矩阵按行合并
mat_data<-rbind(data1,data2,data3)


#设置热图颜色范围：
colpattern = colorRamp2(c(-1, 0, 2), c('#1E466E','white','#E65F54'))# c("#2574AA", "white", "#ED7B79"))
#设置扇区，这里划分了三个扇区，KEGG,EC和细菌种类。
level_test<-c(rep("level4-KO",67),rep("level3-pathway",16),rep("level1",6)) %>% factor()

# 扇区，单元格边框
pdf("plot2.pdf",width = 8, height = 6)
circos.par(gap.after = c(10,10,12))  # 扇区间隔
circos.heatmap(mat_data, split = level_test,
               col = colpattern, rownames.side = "outside",
               cluster = T,cell.lwd=0.8,
               cell.border="white",track.height = 0.2)

# 添加列名 panel.fun
circos.track(track.index = get.current.track.index(),
             panel.fun = function(x, y) {
                 if(CELL_META$sector.numeric.index < 3) { # the last sector
                     cn = colnames(mat_data)
                     n = length(cn)
                     circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
                                 1:n - 0.5, cn, cex = 0.3, adj = c(0, 0.5), facing = "inside")
                 }
            }, bg.border = NA)


# 添加连接线
#灰色连接线数据
df_link = data.frame(
    from_index = c(1:67, sample(67, 60, replace = F)),
    to_index = sample(68:83, 127, replace = T)
)
#红色连接线数据
red_df_link<-data.frame(from_index = 68:83,
                        to_index = c(85,85,86,86,84,84,84,84,84,87,87,88,88,89,89,89))

for(i in seq_len(nrow(df_link))) {
    circos.heatmap.link(df_link$from_index[i],
                        df_link$to_index[i],
                        col = "grey")
}


for(i in 1:2) {
    circos.heatmap.link(red_df_link$from_index[i],
                        red_df_link$to_index[i],
                        col = "#f2a7a2",
                        lwd = 3)
}


for(i in 3:4) {
    circos.heatmap.link(red_df_link$from_index[i],
                        red_df_link$to_index[i],
                        col = "#fce694",
                        lwd = 3)
}

for(i in 5:9) {
    circos.heatmap.link(red_df_link$from_index[i],
                        red_df_link$to_index[i],
                        col = '#eb6468',
                        lwd = 3)
}

for(i in 10:11) {
    circos.heatmap.link(red_df_link$from_index[i],
                        red_df_link$to_index[i],
                        col = '#c3e4f5',
                        lwd = 3)
}

for(i in 12:13) {
    circos.heatmap.link(red_df_link$from_index[i],
                        red_df_link$to_index[i],
                        col = '#f4cfd6',
                        lwd = 3)
}

for(i in 14:16) {
    circos.heatmap.link(red_df_link$from_index[i],
                        red_df_link$to_index[i],
                        col = '#6a88c2',
                        lwd = 3)
}


circos.clear()
dev.off()
