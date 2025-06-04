# 一、导入R包 ----------------------------------------------------------------
#读取数据用，默认已安装
library(openxlsx)
library(readr)
library(hdf5r)
library(patchwork)
#单细胞分析的包
library(glmGamPoi) #SeuratV5之后的SCT标准化依赖的包
library(sctransform) #SeuratV5之后的SCT标准化用的包
library(glmGamPoi) #SCT标准化加速包
library(presto) #findmarkers时的加速工具，下载需用devtools::install_github('immunogenomics/presto')
library(Seurat) #5.2.1，核心包，与4版本的代码有差别，官网https://satijalab.org/seurat/articles/get_started_v5_new
library(monocle3) #拟时序分析
library(SingleR) #自动注释包
library(celldex) #自动注释包的参考数据集单独成包
library(harmony) #批次效应
library(DoubletFinder) #找双细胞
#下游分析
library(rjson)
library(GEOquery)
library(limma) #微阵列和RNA-seq数据的差异表达分析
library(clusterProfiler) #功能富集分析
library(glmnet) #正则计算回归
library(GSVA) #基因集变异分析,单样本富集方法
#数据清洗，默认已安装
library(reshape2)
library(future) #并行计算
library(stringr) #字符处理
library(dplyr)
library(data.table)
library(Matrix)
library(tidyverse)
#绘图包，默认已安装
library(grid)
library(gridExtra)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(patchwork) #图形排列对其用的，针对ggplot2
library(RColorBrewer)
library(viridis)
library(tibble) #排版
#生存分析
library(survival)
library(timeROC)
library(survminer) #生存分析绘图





# 二、读取数据 --------------------------------------------------------------------
#1.基础准备
samples <- list.files("samples250529/")
sample250529_list <- list()

#2.读取数据并创建Seurat对象：删除“小于200个基因表达的细胞”与“小于3个细胞表达的基因”
for(sample in samples){
  sample.path <- paste0("samples250529/", sample)
  sample.data <- Read10X(data.dir=sample.path) #读取x10数据
  sample.SeuratObject <- CreateSeuratObject(counts = sample.data, 
                         project = sample, #样本名，后续要改循环的地方要用
                         min.cells = 3, 
                         min.features = 200) #创建Seurat对象，删除表达基因数小于200个的细胞，删除小于3个细胞表达的基因（质控）
  sample250529_list <- append(sample250529_list, sample.SeuratObject)
}

#3.合并并保存RData
sample250529_combined <- merge(sample250529_list[[1]], y=sample250529_list[-1], add.cell.ids=samples) #在meta.data的行名中添加样品id，
sample250529_combined@meta.data$label.main <- ifelse(sample250529_combined@meta.data$orig.ident %in% c("sample1", "sample2", "sample3", "sample4", "sample5"),
                                                     "B1-5","B6-10")
sample250529_combined_new <- JoinLayers(sample250529_combined) #sample250529_combined_2_10@assays[["RNA"]]@layers有多个层次，不利于后续分析，本步合并
table(sample250529_combined_new@meta.data$orig.ident) #查看样品名
table(sample250529_combined_new@meta.data$label.main) #查看分组名
save(sample250529_combined_new, file="samples250529_combined.RData")






# 三、质控1：按特征值质控 --------------------------------------------------------------
#1.清空环境，导入数据
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="samples250529_combined.RData")

#2.获取线粒体基因表达情况，MT或者mt都是需要考虑的
sample250529_combined_new[["percent.mt"]] <- PercentageFeatureSet(sample250529_combined_new, pattern="^mt-") #有些会以大写MT开头，这里正则表达需要改一下
table(sample250529_combined_new@meta.data$percent.mt) #检查线粒体情况的步骤；疑惑，质量这么好吗，还是cell ranger自带去除mt的效果？

#3.创建标准矩阵
sample250529_combined_new <- NormalizeData(sample250529_combined_new)

#4.质控前绘图
##其实个人认为用ggsave储存更方便（方便后期改格式和进PS调整位置）
pdf("samples250529_all_quality_control_before.pdf", width=8, height=10)
VlnPlot(sample250529_combined_new, features =c("nFeature_RNA","nCount_RNA","percent.mt"), ncol =3)
FeatureScatter(sample250529_combined_new, "nCount_RNA", "percent.mt", group.by="orig.ident" )
FeatureScatter(sample250529_combined_new, "nCount_RNA", "nFeature_RNA", group.by="orig.ident" )
dev.off() #这次没有明显的双峰，不需要直接去

#5.数值查看
##凭感觉去筛除极端值？这个质控范围怎么定的还需要读读文献
quantile(sample250529_combined_new$nFeature_RNA, seq(0.01,0.1,0.01)) #显示前10%
quantile(sample250529_combined_new$nFeature_RNA, seq(0.9,1,0.01))#显示90%-100% #这个样品nFeature设99%
quantile(sample250529_combined_new$nCount_RNA, seq(0.01,0.1,0.01))
quantile(sample250529_combined_new$nCount_RNA, seq(0.9,1,0.01)) #这里nCount_RNA设98%
quantile(sample250529_combined_new$percent.mt, seq(0.01,0.1,0.01)) 
quantile(sample250529_combined_new$percent.mt, seq(0.9,1,0.01)) #这个样品Mt设95%

#6.1传递质控1后的矩阵
sample250529_qc <- subset(sample250529_combined_new, subset = nFeature_RNA > 900 & nFeature_RNA < 6000 &
                            nCount_RNA < 23000 & percent.mt < 4)
#6.2质控后绘图
pdf("samples250529_all_quality_control_after.pdf", width=8, height=10)
VlnPlot(sample250529_qc, features =c("nFeature_RNA","nCount_RNA","percent.mt"), ncol =3)
FeatureScatter(sample250529_qc, "nCount_RNA", "percent.mt", group.by="orig.ident" )
FeatureScatter(sample250529_qc, "nCount_RNA", "nFeature_RNA", group.by="orig.ident" )
dev.off()

#7储存生成的文件
save(sample250529_qc, file="samples250529_all_qc.RData")





# 三、质控2：SCT标准化[替代标准化、去心化与高变检索三步] --------------------------------------------
#1.清空环境，导入数据
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="samples250529_all_qc.RData")
plan(multisession, workers=8) #使用8核
#nbrOfWorkers() #确认核心数
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存

#2.分步运行版本
#2.1 标准化normalization
sample250529_normalization <- NormalizeData(sample250529_qc, 
                                            normalization.method = "LogNormalize") #标准化方法，包括LogNormalize、CLR、RC三种，默认为LogNormalize
#2.2 识别高变特征（特征选择）Identification of highly variable features (feature selection)
sample250529_variable <- FindVariableFeatures(sample250529_normalization, 
                                              selection.method = "vst", 
                                              nfeatures = 2000) #识别变化比较大的基因，默认值是2000（据说这样比较准），但SCT优势就是能找3000个
sample250529_variable_top100 <- head(VariableFeatures(sample250529_variable), 100)
write.table(top100, file="sample250529_variable_top100.txt") #写出前100个突变基因并保存
sample250529_variable_top20 <- head(VariableFeatures(sample250529_variable), 20) #写出前10个突变基因
#2.2补充 绘图：平均表达量-SV变异关系图
pdf("average_VS_SV_top20.pdf", width=8, height=10)
VariableFeaturePlot(sample250529_variable)
LabelPoints(plot = plot1, points = sample250529_variable_top20, repel = TRUE)
dev.off()
#2.3归一化scale
all.genes <- rownames(sample250529_variable)
sample250529_scale <- ScaleData(sample250529_variable, features = all.genes)
save(sample250529_scale, file="samples250529_all_scale_all.RData")

#2.一步版本：SCT标准化
##[实际上，如果不需要绘图的话，10万细胞以内直接使用SCT标准化会更快]
sample250529_sct <- SCTransform(sample250529_qc, 
                                vars.to.regress="percent.mt", #去除percent.mt的列（因为质控过可以不用了）
                                verbose=T, #显示进度条，可以设置成false
                                method="glmGamPoi") #SCT运行加速用

#3.储存生成的文件
save(sample250529_sct, file="samples250529_all_sct.RData")





# 四、线性降维分析（粗降） --------------------------------------------------------------
#1.清空环境，导入数据
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="samples250529_all_sct.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存

#2.PCA降维并绘图
sample250529_PCA <- RunPCA(object = sample250529_sct, 
                           features = VariableFeatures(object=sample250529_sct)) #设定特征对象，可以换用其他层次，会对结果产生相应的干扰
pdf("samples250529_all_PCAPlot.pdf", width=8, height=10)
DimPlot(sample250529_PCA, reduction = "pca", group.by="label.main")
DimPlot(sample250529_PCA, reduction = "pca")
dev.off()

#3.确定维度数
#3.1 指标法，这一指标没有后面两种重要常用
#3.1.1 确定每个PC 贡献所占百分比（一般到阈值的时候应小于5%），用于过滤技术噪音/测序误差
stdev <- sample250529_PCA[["pca"]]@stdev
pct <- stdev^2 / sum(stdev^2) * 100
pct
#3.1.2 主成分累积贡献（即计算每添加一个PC，累计贡献的百分比变为多少，一般应大于90%），保留主要生物学变异
cumu <-cumsum(pct)
cumu
#3.1.3 确定信息增益饱和点（一般应为0.1%），确定信息增益饱和点（肘点）
diff <- c(0, abs(diff(pct)))
diff
#3.2 生物学信息法
##在特定的维度中如果发现熟悉的marker，则这个维度通常是应该被纳入考量的（有利于下游分析）
##具体而言有三种可视化方法VizDimReduction() 、 DimPlot() 和 DimHeatmap()，主要是展示形式不同，最推荐最后一个，信息量最大
#VizDimReduction()
pdf("samples250529_all_PCA_VizDimReduction.pdf", width=16, height=20)
VizDimLoadings(sample250529_PCA, 
               dims=1:16, #维度数是人为设置的，最大为50，这里为了避免重复下面没有弄相同范围的维度数了
               reduction="pca")
dev.off()
#DimHeatmap()
pdf("samples250529_all_PCA_DimHeatmap.pdf", width=16, height=20)
DimHeatmap(sample250529_PCA, dims = 1:16, cells = 500, balanced = TRUE)
dev.off()
#3.3 肘部图,本图是最重要的降维标准
##确定维度数时，可以选择拐点数，但一般不应小于6
pdf("samples250529_all_ElbowPlot.pdf", width=8, height=10)
ElbowPlot(sample250529_PCA, ndims=50) #感觉其实在10左右开始变平滑
dev.off()

#4.储存生成的文件
save(sample250529_PCA, file="samples250529_all_PCA.RData")

##其实seuratV5官网还说到了一种监督下的方法，不过这里先跳过了





# 五、Harmony去除批次效应 -----------------------------------------------------------
#options(future.globals.maxSize = 8 * 1024^3) #设置运行内存
# sample250529_PCA <- RunHarmony(sample250529_PCA, group.by.vars="orig.ident", assay.use="SCT", 
#                                max_iter = 20)
#table(sample250529_PCA@meta.datasorig.ident)

##因为这里不确定批次效应的去处是否符合差异分析的需求，所以我没有去除；
##此外，从UMAP的结果来看，批次差异并不明显，反而是分组带来的数据差异问题值得关注





# 六、UMAP降维聚类分析（细降） ----------------------------------------------------------
#1.清空环境，导入数据
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="samples250529_all_PCA.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存

#2.降维聚类
sample250529_cluster <- FindNeighbors(sample250529_PCA, dims = 1:16)  #dims决定分群情况 ##如果用了上面的去处批次效应的数据的话，需要添加参数reduction="harmony"
sample250529_cluster <- FindClusters(sample250529_cluster, resolution = 0.5)

#3.UMAP降维并绘图
sample250529_UMAP <- RunUMAP(sample250529_cluster, dims = 1:16) #dims决定绘图形状
pdf("sample250529_UMAPPlot.pdf", width=8, height=10)
DimPlot(sample250529_UMAP, reduction = "umap", label=T) #调色参数cols=a，其中a是预设好的颜色vector
DimPlot(sample250529_UMAP, reduction = "umap", label=F, group.by="orig.ident")#样品来源
DimPlot(sample250529_UMAP, reduction = "umap", label=F, group.by="label.main")#组别
FeaturePlot(sample250529_UMAP, 
            features = c("Ptprc", "Sdc1", "Ighd", "Fas" #细胞marker，MBC
            ))
FeaturePlot(sample250529_UMAP, 
            features = c("Il9r","Ptpn22","Tlr7","Zbtb18","Zbtb32")) #重要基因
# FeaturePlot(sample250529_UMAP, 
#             features = c("Il9r","Cyb561a3","Cfp","Ptpn22","Mpeg1","Zfp385a")) #31个差异基因，第1-6个
# FeaturePlot(sample250529_UMAP, 
#             features = c("Sdc4","Zmynd11","St8sia6","Tspan32","Tlr7","Naip5")) #31个差异基因，第7-12个
# FeaturePlot(sample250529_UMAP, 
#             features = c("Pdia4","Ache",#"AhnaK",
#                          "Cdkn1b","Smad3","9930111J21Rik1","Gpr174")) #31个差异基因，第13-19个，其中第15个基因"AhnaK"消失了
# FeaturePlot(sample250529_UMAP, 
#             features = c("Cdh17","Fcgrt","Ptpre","Aff1","Lgals3bp","BC147527")) #31个差异基因，第20-25个
# FeaturePlot(sample250529_UMAP, 
#             features = c("Serpina3f","Pglyrp1","Rnasel","Tnfrsf18","Ifit3","Pira2")) #31个差异基因，第26-31个
dev.off()

#4.储存生成的文件
saveRDS(sample250529_UMAP, file = "samples250529_UMAP.rds") #注意这里变格式了，后续读取格式也应改变





# 七、自动注释（粗定义群） --------------------------------------------------------------
##自动注释，此步可跳过，但如果觉得样本不够纯净可以做做去除

#1.导入数据，设置环境
rm(list=ls())
samples250529_UMAP <- readRDS("samples250529_UMAP.rds")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存
##个人感觉主要是看样本纯不纯净，除非我们能找到亚群高度重合的论文数据集作为参考，否则意义不大

#2.传入参考数据集并注释
ref1_mus_immune <- celldex::ImmGenData() #传入参考基因集
sample250529_SingleR <- GetAssayData(object=sample250529_UMAP@assays[["SCT"]], layer="counts") #获取标准化矩阵层
sample250529_SingleR.ref1 <- SingleR(test = sample250529_SingleR, ref = ref1_mus_immune, 
                                     labels = ref1_mus_immune$label.main) #这里ref和Labels可以用多个数据集综合注释，将传入参数变为list()即可
table(sample250529_SingleR.ref1$labels) #查看细胞集情况
sample250529_UMAP@meta.data$labels <- sample250529_SingleR.ref1$labels #回传labels到UMAP中（不过懒得单独存为一个文件了，毕竟意义不大）

#3.绘图
pdf("sample250529_UMAP_SingleR.auto.annotation.pdf", width=8, height=10)
DimPlot(sample250529_UMAP, group.by = "labels",reduction = "umap") #生成UMAP的自动注释绘图
dev.off()

#4.保存生成的文件
save(sample250529_UMAP, file="sample250529_UMAP_SingleR.RData")



# 八、基于UMAP的精细质控 -------------------------------------------------------------------
## （一）双细胞筛除 -------------------------------------------------------------------
#1.设置环境 
# rm(list = setdiff(ls(), "sample250529_UMAP")) 
# gc()
# Sys.setenv(R_MAX_NUM_DLLS=999)
# load(file="sample250529_UMAP_SingleR.RData")
# plan(multisession, workers=1) #前面设多核的话，这里运行会超内存（除非用服务器），所以调回来
# options(future.globals.maxSize = 25 * 1024^3) #10万个细胞差不多需要22GB，这里需要重新设置运行内存

#2.进行计算
#2.1 首先获得最佳的K值
##pK表示领域大小
# sweep.res.list <- paramSweep(sample250529_UMAP, PCs =1:16, sct=F)
# sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
# bcmvn <- find.pK(sweep.stats)
# pk_best = bcmvn %>%
#   dplyr::arrange(desc(BCmetric))%>%
#   dplyr::pull(pK)%>%
#   .[1]%>% as.character()%>% as.numeric()
#2.2 估算出双细胞群中，homotypic doublets的比例
# annotations <- sample250529_UMAP$seurat_clusters
# homotypic.prop <- modelHomotypic(annotations)
# print(homotypic.prop)
##双细胞占比为7%左右
# nExp_poi <- round( 0.07 * nrow( sample250529_UMAP@meta.data) )
# nExp_poi.adj <- round( nExp_poi * (1-homotypic.prop) )
#2.3 模拟出artificial doublet数量。不同取值对识别结果影响不大，默认为0.25
# sample250529_UMAP_test <- sample250529_UMAP
# DefaultAssay(sample250529_UMAP_test) <- "RNA"
# sample250529_UMAP_test[["SCT"]] <- NULL  # 如果只用 SCT
# gc()  # 回收内存
# sample250529_UMAP_test <- doubletFinder(sample250529_UMAP_test, PCs=1:16,
#                                    pN=0.25, pK=pk_best, nExp=nExp_poi.adj, sct=F) #把默认值reuse.pANN=F加上反而会报错，不知道是为什么
# colnames(sample250529_UMAP@meta.data) #内存会爆，先不跑了
#2.4将列名改为"Double_score"和"Is_Double"
# colnames(sample250529_UMAP@meta.data)
# colnames(sample250529_UMAP@meta.data)[length(colnames(sample250529_UMAP@meta.data))-1] <- "Double_score"
# colnames(sample250529_UMAP@meta.data)[length(colnames(sample250529_UMAP@meta.data))] <- "Is_Double"
#2.5查看Doubletrinder分析结果
# head(sample250529_UMAP@meta.data[,c("Double_score","Is_Double")])
# table(sample250529_UMAP@meta.data["Is_Double"])

#3.绘制DoubletFinder分类的UMAP图并保存
# pdf("UMAP_DoubletFinder.pdf", width=8, height=10)
# DimPlot(sample250529_UMAP, reduction="umap", group.by="Is_Double")
# dev.off()





## （二）针对细胞周期进行了解 --------------------------------------------------------------
#1.获取基因
#1.1 获取G2M期相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <-CaseMatch(search=g2m_genes, match=rownames(sample250529_UMAP))
#1.2获取S期相关基因基因
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search=s_genes,match=rownames(sample250529_UMAP)) #主要目的是匹配

#2.细胞周期阶段评分
sample250529_UMAP <- CellCycleScoring(sample250529_UMAP, g2m.features=g2m_genes, s.features=s_genes)
colnames(sample250529_UMAP@meta.data)
table(sample250529_UMAP$Phase)

#3.画图
pdf("sample250529_UMAP_CellCycle.pdf", width=8, height=10)
DimPlot(sample250529_UMAP, reduction= "umap", group.by="Phase")
dev.off()

#4.保存生成的文件
#save(sample250529_UMAP, file="sample250529_UMAP_CellCycle.RData")




# 九、提取亚群再分析 -----------------------------------------------------------------
## （一）根据SingleR提取细胞类型进行细分亚群定义 -----------------------------------------------------------------
#1.加载数据设置环境
rm(list=ls())
load(file="sample250529_UMAP_CellCycle.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存

#2.提取亚群并绘图查看
sample250529_UMAP.B.more = sample250529_UMAP[,sample250529_UMAP@meta.data$labels %in% c("B cells","B cells, pro")]
sample250529_UMAP.B.less = sample250529_UMAP[,sample250529_UMAP@meta.data$seurat_clusters %in% c(0,1,2,3,4,5,6,7,9,14)] #-8,10,11,12,13,15,16,17
DimPlot(sample250529_UMAP.B.more, reduction = "umap", label=T)
DimPlot(sample250529_UMAP.B.less, reduction = "umap", label=T)##查看去除干不干净

#3.程序重运行（即四、六）
#3.1 PCA降维
sample250529_UMAP.B.more <- RunPCA(object = sample250529_UMAP.B.more, 
                           features = VariableFeatures(object=sample250529_UMAP.B.more))
sample250529_UMAP.B.less <- RunPCA(object = sample250529_UMAP.B.less, 
                              features = VariableFeatures(object=sample250529_UMAP.B.less))
#3.2 降维聚类
sample250529_cluster.B.more <- FindNeighbors(sample250529_UMAP.B.more, dims = 1:16)  #dims决定分群情况 ##如果用了上面的去处批次效应的数据的话，需要添加参数reduction="harmony"
sample250529_cluster.B.more <- FindClusters(sample250529_cluster.B.more, resolution = 0.5)
sample250529_cluster.B.less <- FindNeighbors(sample250529_UMAP.B.less, dims = 1:16)  #dims决定分群情况 ##如果用了上面的去处批次效应的数据的话，需要添加参数reduction="harmony"
sample250529_cluster.B.less <- FindClusters(sample250529_cluster.B.less, resolution = 0.5)
#3.3 UMAP降维并绘图
sample250529_UMAP.B.more <- RunUMAP(sample250529_cluster.B.more, dims = 1:16) #dims决定绘图形状
pdf("sample250529_UMAP.B.more.pdf", width=8, height=10)
DimPlot(sample250529_UMAP.B.more, reduction = "umap", label=T)
dev.off()
sample250529_UMAP.B.less <- RunUMAP(sample250529_cluster2, dims = 1:16) #dims决定绘图形状
pdf("sample250529_UMAP.B.less.pdf", width=8, height=10)
DimPlot(sample250529_UMAP.B.less, reduction = "umap", label=T)
dev.off()

#4.保存数据
save(sample250529_UMAP.B.more, file="sample250529_UMAP.B.more.RData")
save(sample250529_UMAP.B.less, file="sample250529_UMAP.B.less.RData")



## （二）手动注释（细定义群） --------------------------------------------------------------
###有监督注释有两种方法：看文献 + 查数据库，下面是两个推荐较多的数据库
# CellMarker, 网址是http://xteam.xbio.top/CellMarker/index.jsp
# PanglaoDB, 网址是https://panglaodb.se/index.html#google_vignette
###注释的时候要注意，很多基因是有别名的，比如什么CD45别名Ptprc之类的（打PTPRC都搜不出来）
# 推荐使用genecard来确认注释的基因的名称
# genecard，网址是https://www.genecards.org/

### 设置环境 --------------------------------------------------------------
rm(list=ls())
load(file="sample250529_UMAP.B.less.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存



### 有监督注释 --------------------------------------------------------------
#1 参考ref1进行注释并绘图
#参考文献reef1为《Single-cell BCR and transcriptome analysis after influenza infection reveals spatiotemporal dynamics of antigen-specific B cells》f1d
# sample250529_markers_ref1 <- c("Cd1d1","Cd9", #边缘区（MZ）B 细胞【C6簇】——文献明确提及的经典标准
#                                "Cr2","Plac8","Lars2", #边缘区（MZ）B 细胞【C6簇】——文献高表达的其他基因
#                                "Ccr7","Ebf1","Cd74","Nr4a1", #Naive细胞【C1.2.4簇】——文献明确提及的经典标准
#                                "Btg1","Cd79a","Ighm","Fcmr","Ighd","Sell","Fcer2a","Junb","Cd69",#Naive细胞【C1.2.4簇】——文献高表达的其他基因（3个簇，所以略有差异）
#                                "Cd83", #似乎是一种比较交界状态的表达？
#                                "Npm1","Eif5a","Ran","Mif","Eif4a1", #preGC【C3簇】——文献高表达的基因
#                                "Ppia","Top1","Hmgb1","Hmgb2", #EarlyGC【C14簇】——包括上面的preGC的高表达基因，文献还有这些高表达的基因
#                                "Hmces","Ptma", #GCDZ【G2M期C15簇、S期C10簇】——文献高表达的基因
#                                "Tuba1b", "Top2a","Tubb5", #GCDZ【G2M期C16簇、S期C8簇】，文献除上一行外，还高表达的基因
#                                "Aicda","Tnfrsf13c", "Actb", #GCLZ【S期C13簇、不明期C12簇】——文献高表达的基因
#                                "Bcl6","Bach2","Zbtb20", "Mki67", "Foxo1","Ptprc","Apoe",#GCLZ/G2M【C11簇】和preMem【C9簇】——文献高表达的基因
#                                "Vim", "Cd38", #Bmem【C5簇】——文献高表达的基因
#                                "Slpi","Sdc1","Prdm1", #浆细胞【C7簇】——文献提及的经典标准
#                                "Xbp1","Jchain" #浆细胞——文献高表达的其他基因
# )
sample250529_markers_ref1 <- c("Cr2","Cd1d1","Cd9","Lars2","Ccr7",
                               "Ebf1","Btg1","Cd79a","Ighm","Fcmr",
                               "Ighd","Cd74","Sell","Fcer2a","Junb",
                               "Nr4a1","Cd69","Cd83","Npm1","Eif5a",
                               "Ran","Mif","Eif4a1","Hmces","Ppia",
                               "Top1","Hmgb1","Hmgb2","Tuba1b", "Top2a",
                               "Ptma","Tubb5","Aicda","Tnfrsf13c", "Actb",
                               "Bcl6","Bach2","Zbtb20", "Mki67", "Foxo1",
                               "Ptprc","Apoe","Plac8","Vim", "Cd38",
                               "Slpi","Sdc1","Prdm1","Xbp1","Jchain") #按文献原有顺序重新排列
dotplot_sample250529_UMAP_ref1 <- DotPlot(sample250529_UMAP, features=sample250529_markers_ref1, cols = c("lightblue","orange"),col.min=-1
) + RotatedAxis()
ggsave("sample250529_UMAP_ref1_dotplot.pdf", dotplot_sample250529_UMAP_ref1, width=14, height =6)
vlnplot_sample250529_UMAP_ref1<- VlnPlot(sample250529_UMAP, features=sample250529_markers_ref1, stack =T, flip =T) #+ NoLegend()
ggsave("sample250529_UMAP_ref1_vlnplot.pdf", vlnplot_sample250529_UMAP_ref1, width=14, height=30)
dotplot_sample250529_UMAP.B.less_ref1 <- DotPlot(sample250529_UMAP.B.less, features=sample250529_markers_ref1, cols = c("lightblue","orange"),col.min=-1
) + RotatedAxis()
ggsave("sample250529_UMAP.B.less_ref1_dotplot.pdf", dotplot_sample250529_UMAP.B.less_ref1, width=14, height =6)
vlnplot_sample250529_UMAP.B.less_ref1<- VlnPlot(sample250529_UMAP.B.less, features=sample250529_markers_ref1, stack =T, flip =T) #+ NoLegend()
ggsave("sample250529_UMAP.B.less_ref1_vlnplot.pdf", vlnplot_sample250529_UMAP.B.less_ref1, width=14, height=30)

#2 换用更全面的marker
#2.1 Marker
genes_to_plot <- c(
  # B细胞标志物
  "Cr2", "Cd1d1", "Cd9", "Lars2", "Ccr7", "Ebf1", "Btg1", "Cd79a",
  "Ighm", "Fcmr", "Ighd", "Cd74", "Sell", "Fcer2a", "Junb", "Nr4a1",
  "Cd69", "Cd83", "Npm1", "Eif5a", "Ran", "Mif", "Eif4a1", "Hmces",
  "Ppia", "Top1", "Hmgb1", "Hmgb2", "Tuba1b", "Top2a", "Ptma", "Tubb5",
  "Aicda", "Tnfrsf13c", "Actb", "Bcl6", "Bach2", "Zbtb20", "Mki67",
  "Foxo1", "Ptprc", "Apoe", "Plac8", "Vim", "Cd38", "Xbp1", "Jchain",
  "Slpi", "Sdc1",
  "Prdm1", "Zbtb18", "Klf2", "Il9r", "Ccr6", "S1pr2", #剩下8个是值得关注的标志物
  "Ccnb1", "Ccnb2", "Zbtb32",
  # 非B细胞标志物
  "Cd3e", "Cd4", "Cd8a", "Cd8b1",  # T细胞标志物
  "Nkg7", "Gzmb", "Prf1",         # NK细胞标志物
  "Cd14", "Cd68", "Adgre1",       # 单核/巨噬细胞标志物
  "Itgam",                        # 髓系标志物 (已移除Cd11b)
  "Siglech", "Ly6c2", "Ly6g",     # 粒细胞标志物
  "Col1a1", "Col1a2", "Acta2",    # 基质/成纤维细胞标志物
  "Pecam1", "Cdh5", "Vwf",       # 内皮细胞标志物
  "Epcam", "Krt19", "Krt8",      # 上皮细胞标志物
  "H2-Aa", "H2-Eb1",             # MHC II标志物
  "Foxp3", "Il2ra",              # Treg标志物
  "Il7r", "S100a4", "S100a8", "S100a9" # 其他标志物
)
#2.2 生成 DotPlot 数据（不直接绘制）
#2.2.1 导入数据
dot_data <- DotPlot(sample250529_UMAP.B.less,features = genes_to_plot,
                    scale = TRUE,cluster.idents = FALSE
                    )$data
#2.2.2 组别数据 #似乎有点小bug？
gene_groups <- data.frame(
  gene = genes_to_plot,
  group = c(
    rep("B cells", 57),rep("T cells", 4),rep("NK cells", 3),rep("Myeloid cells", 4),
    rep("granulocytes", 3),rep("stromal cells", 3),rep("endothelial cells", 3),rep("epithelial cells", 3),
    rep("MHC II", 2),rep("Treg cells", 2),rep("other", 5)
  )
)
#2.2.3 合并分组信息
dot_data <- merge(dot_data, gene_groups, by.x = "features.plot", by.y = "gene")
#2.2.4固定基因顺序（避免字母排序）
dot_data$features.plot <- factor(dot_data$features.plot, levels = rev(genes_to_plot))
#2.3 绘图
#2.3.1 主图（点图）
main_plot <- ggplot(dot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(
    low = "blue",mid = "grey90",high = "red",
    midpoint = 0,limits = c(-2.5, 2.5)
  ) +
  scale_size(range = c(0, 6)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  labs(color = "Expression\n(Z-score)",size = "Percent\nExpressed")
#2.3.2 分组注释条（颜色块）
annotation_plot <- ggplot(dot_data, aes(x = 1, y = features.plot, fill = group)) +
  geom_tile() + # 用色块表示分组
  scale_fill_manual(values = c(
    "B cells" = "#1F78B4","T cells" = "#33A02C","NK cells" = "#E31A1C","Myeloid cells" = "#FF7F00",
    "granulocytes" = "#6A3D9A","stromal cells" = "#B15928","endothelial cells" = "#A6CEE3","epithelial cells" = "#B2DF8A",
    "MHC II" = "#FB9A99","Treg cells" = "#FDBF6F","other" = "#CAB2D6"
  )) +
  theme_void() +
  theme(legend.position = "left",plot.margin = margin(0, 0, 0, 0)) +
  labs(fill = "Cell Type")
#2.3.3 组合点图和注释条
final_plot <- annotation_plot + main_plot + 
  plot_layout(widths = c(0.05, 1)) # 注释条宽度5%，点图宽度95%
#2.4 保存PDF
pdf("sample250529_UMAP.B.less_ref1_gene.expression.dotplot.pdf", width = 14, height = 30, useDingbats = FALSE)
print(final_plot)
dev.off()



### 无监督注释  -----------------------------------------------------------------
#0.找出各簇细胞的标志marker的其他函数
#sample250529_markers.cluster0 <- FindMarkers(sample250529_UMAP, ident.1=0, ident.2=c(1:3)) #只看某一簇细胞（0）与其他指定簇（1到3）细胞的markers；若不指定ident.2，则为其余全部细胞
#sample250529_markers <- FindAllMarkers(sample250529_UMAP, only.pos=TRUE, logfc.threshold=1) #只看上调，一簇细胞与其他所有细胞簇的log2fc阈值为1

#1.严格
sample250529_markers_all_strict <- FindAllMarkers(sample250529_UMAP.B.less, only.pos=F, 
                                                  logfc.threshold=1, 
                                                  min.pct=0.5, #min_pct是两组细胞中至少应有一组的表达比例超过这一值
                                                  min.diff.pct=0.3) #较为严格的阈值，min_diff_pct为两组表达比例差异
write.csv(sample250529_markers_all_strict, file="sample250529_all_SCT_markers_strict.csv")

#2.宽松
sample250529_markers_all_lenient <- FindAllMarkers(sample250529_UMAP.B.less, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) #宽松一点的阈值
sample250529_markers_all_lenient$label <- ifelse(sample250529_markers_all_lenient$avg_log2FC>0, "positive", "negative") #赋值
sample250529_markers_all_lenient$label <- paste0(sample250529_markers_all_lenient$cluster, "_", sample250529_markers_all_lenient$label) #贴标签，方便找各簇上下调基因top
write.csv(sample250529_markers_all_lenient, file="sample250529_all_SCT_markers_lenient.csv")

#2.1 宽松+按p值排序
#上下调各top30
sample250529_markers_all_lenient <- read.csv(file="sample250529_all_SCT_markers_lenient.csv") #导入数据
sample250529_markers_all_lenient <- arrange(sample250529_markers_all_lenient, p_val_adj, by_group="label") #分label按p.adj排序，从小到大
sample250529_markers_all_lenient_ptop30 <- sample250529_markers_all_lenient %>%
  group_by(label) %>%
  filter(row_number() <= 30 ) %>%
  ungroup() #按p值，提取每簇细胞上下调基因各top30
write.csv(sample250529_markers_all_lenient_ptop30, file="sample250529_all_SCT_markers_lenient_ptop30.csv") #保存
#上下调各top5
sample250529_markers_all_lenient_ptop5 <- sample250529_markers_all_lenient %>%
  group_by(label) %>%
  filter(row_number() <= 5 ) %>%
  ungroup() #按p值，提取每簇细胞上下调基因各top5
write.csv(sample250529_markers_all_lenient_ptop5, file="sample250529_all_SCT_markers_lenient_ptop5.csv") #保存
#top5绘图
sample250529_markers_all_lenient_ptop5 <- arrange(sample250529_markers_all_lenient_ptop5, cluster, avg_log2FC)
genes_to_plot <- rev(unique(sample250529_markers_all_lenient_ptop5$gene)) #反序保证下面翻转坐标之后与csv文件顺序一致
pdf("sample250529_UMAP.B.less_gene.expression.ptop5.dotplot.pdf", width = 14, height = 30, useDingbats = FALSE)
DotPlot(sample250529_UMAP.B.less, features = genes_to_plot,) + 
  ggplot2::scale_color_gradient2(
    low = "#2E9FDF",mid = "white",high = "#FC4E07",
    midpoint = 0,limits = c(-2.5, 2.5)
  ) + #按需求填色
  theme_minimal() + #画上格线方便阅读
  coord_flip() #坐标翻转
dev.off()

#2.2 宽松+按log2FC值排序
#上下调各top30
sample250529_markers_all_lenient <- arrange(sample250529_markers_all_lenient, desc(abs(avg_log2FC)), by_group="label") #绝对值降序
sample250529_markers_all_lenient_log2FCtop30 <- sample250529_markers_all_lenient %>%
  group_by(label) %>%
  filter(row_number() <= 30) %>%
  ungroup() #按log2FC，提取每簇细胞上下调基因各top30
write.csv(sample250529_markers_all_lenient_log2FCtop30, file="sample250529_all_SCT_markers_lenient_log2FCtop30.csv")
#上下调各top5
sample250529_markers_all_lenient_log2FCtop5 <- sample250529_markers_all_lenient %>%
  group_by(label) %>%
  filter(row_number() <= 5) %>%
  ungroup() #按log2FC，提取每簇细胞上下调基因各top5
write.csv(sample250529_markers_all_lenient_log2FCtop5, file="sample250529_all_SCT_markers_lenient_log2FCtop5.csv")
#宽松log2FC top5绘图
sample250529_markers_all_lenient_log2FCtop5 <- arrange(sample250529_markers_all_lenient_log2FCtop5, cluster)
genes_to_plot <- rev(unique(sample250529_markers_all_lenient_log2FCtop5$gene)) #反序保证下面翻转坐标之后与csv文件顺序一致
pdf("sample250529_UMAP.B.less_gene.expression.log2FCtop5.dotplot.pdf", width = 14, height = 30, useDingbats = FALSE)
DotPlot(sample250529_UMAP.B.less, features = genes_to_plot) + 
  ggplot2::scale_color_gradient2(
    low = "#2E9FDF",mid = "white",high = "#FC4E07",
    midpoint = 0,limits = c(-2.5, 2.5)
  ) + #按需求填色
  theme_minimal() + #画上格线方便阅读
  coord_flip() #坐标翻转
dev.off()

#3.保存数据
save(sample250529_markers_all_strict, sample250529_markers_all_lenient, file="sample250529_markers_all.RData")



## （三）再次对周期、细胞类型绘图 ---------------------------------------------------------------
###这里沿用上一步的数据，不额外储存和读取了
#1.1定义细胞名
sample250529_UMAP.B.less$celltype.main <- recode(sample250529_UMAP.B.less@meta.data[["seurat_clusters"]],
                                                 "13"="pre GC",
                                                 "8"="GC DZ", "9"="GC DZ", "17"="GC DZ",
                                                 "1"="GC LZ", "2"="GC LZ", "11"="GC LZ",
                                                 "4"="GC MP",
                                                 "5"="Naive", "10"="Naive",
                                                 "3"="MBC","12"="MBC","6"="maybe MBC ?",
                                                 "0"="MZ","7"="MZ","15"="maybe MZ ?",
                                                 "14"="PC","16"="PC"
)#后续依次类推，不同的簇可以根据需求定义为相同的细胞
#1.2按celltype绘图
pdf("sample250529_UMAP.B.less_celltype.pdf", width = 8, height = 10)
DimPlot(sample250529_UMAP.B.less, reduction = "umap", label=T, group.by="celltype.main") #绘图，这里写celltype.main是因为还有可能要细化分型
dev.off()

#2.再次查看细胞周期是否符合预期
g2m_genes <- cc.genes$g2m.genes
g2m_genes <-CaseMatch(search=g2m_genes, match=rownames(sample250529_UMAP.B.less))## 获取G2M期相关基因
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search=s_genes,match=rownames(sample250529_UMAP.B.less))## 获取S期相关基因基因
sample250529_UMAP.B.less <- CellCycleScoring(sample250529_UMAP.B.less, g2m.features=g2m_genes, s.features=s_genes)
colnames(sample250529_UMAP.B.less@meta.data)
table(sample250529_UMAP.B.less$Phase)## 细胞周期阶段评分
pdf("sample250529_UMAP.B.less_CellCycle.pdf", width=8, height=10)
DimPlot(sample250529_UMAP.B.less, reduction= "umap", group.by="Phase")
FeaturePlot(sample250529_UMAP.B.less, 
            features = c("Cdkn1c","Cdkn1a","Cdkn1b","Hes1","Btg1" #G0细胞表达marker
            ))
FeaturePlot(sample250529_UMAP.B.less, 
            features = c("Rb1","Id2","Pcna","Mcm2","Mki67" #第一个为G0细胞表达marker，后三个为G0期细胞不应表达
            ))
dev.off()## 画图




# 十、查看WT与KO分组是否符合预期  -----------------------------------------------------------------
## （零）导入数据，设置环境 --------------------------------------------------------------
rm(list=ls())
load(file="sample250529_UMAP.B.less.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存



## （一）总体分析与无监督注释 -----------------------------------------------------------------
#1.添加WT、KO标签与celltype/clusters的综合标签
sample250529_UMAP.B.less@meta.data$labels_complete <- paste0(sample250529_UMAP.B.less@meta.data$label.main, "_", sample250529_UMAP.B.less@meta.data$celltype.main)
sample250529_UMAP.B.less@meta.data$labels_complete2 <- paste0(sample250529_UMAP.B.less@meta.data$label.main, "_", sample250529_UMAP.B.less@meta.data$seurat_clusters)
table(sample250529_UMAP.B.less@meta.data$labels_complete)
table(sample250529_UMAP.B.less@meta.data$labels_complete2)

#2.提取总体的WT、KO、only.MP.MBC群并保存
sample250529_UMAP.B.less.WT <- sample250529_UMAP.B.less[,sample250529_UMAP.B.less@meta.data$label.main %in% "B1-5"]
sample250529_UMAP.B.less.KO <- sample250529_UMAP.B.less[,sample250529_UMAP.B.less@meta.data$label.main %in% "B6-10"]
sample250529_UMAP.B.less_only.MP.MBC <- sample250529_UMAP.B.less[,sample250529_UMAP.B.less@meta.data$celltype.main %in% c("maybe MBC ?", "MBC", "pre GC", "GC LZ", "GC MP")]
save(sample250529_UMAP.B.less, file="sample250529_UMAP.B.less2.RData")
save(sample250529_UMAP.B.less.WT, sample250529_UMAP.B.less.KO, file="sample250529_UMAP.B.less2_WT.and.KO.RData")
save(sample250529_UMAP.B.less_only.MP.MBC, file="sample250529_UMAP.B.less2_only.MP.MBC.RData")

#3.绘制UMAP图
# rm(list=ls())
# load(file="sample250529_UMAP.B.less2.RData")
# load(file="sample250529_UMAP.B.less2_WT.and.KO.RData")
# load(file="sample250529_UMAP.B.less2_only.MP.MBC.RData")
# plan(multisession, workers=8) #使用8核
# options(future.globals.maxSize = 20 * 1024^3) #设置运行内存 #加载数据，设置环境
# color_for_labels.complete.20 <- c(  "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
#                                  "#b15928", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f",
#                                  "#cab2d6", "#ffff99", "#8dd3c7", "#ff33cc",  "#FFB36D",
#                                  "#6DFF79","#FF8F6D", "#7D6DFF","#DBFF6D", "#6D85FF")
# color_for_labels.complete.10.1 <- c(  "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
#                                     "#b15928", "#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f")
# color_for_labels.complete.10.2 <- c("#cab2d6", "#ffff99", "#8dd3c7", "#ff33cc",  "#FFB36D",
#                                     "#6DFF79","#FF8F6D", "#7D6DFF","#DBFF6D", "#6D85FF")
# color_for_labels.complete.10 <- c(  "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a",
#                                     "#b15928","#a6cee3","#ff33cc","#b2df8a", "#fb9a99") #准备色板
# pdf("sample250529_UMAP.B.less_labels.complete.pdf", width=8, height=10)
# DimPlot(sample250529_UMAP.B.less, reduction = "umap", label=T, group.by="labels_complete", 
#         cols=color_for_labels.complete.20, 
#         alpha = 0.5)
# DimPlot(sample250529_UMAP.B.less.WT, reduction = "umap", label=T, group.by="labels_complete", 
#         cols=color_for_labels.complete.10.1,
#         alpha = 0.5)
# DimPlot(sample250529_UMAP.B.less.KO, reduction = "umap", label=T, group.by="labels_complete", 
#         cols=color_for_labels.complete.10.2,
#         alpha = 0.5)
# DimPlot(sample250529_UMAP.B.less_only.MP.MBC, reduction = "umap", label=T, group.by="labels_complete", 
#         cols=color_for_labels.complete.10,
#         alpha = 0.5)
# dev.off()
pdf("sample250529_UMAP.B.less_labels.complete2.pdf", width=8, height=10)
DimPlot(sample250529_UMAP.B.less, reduction = "umap", label=T, group.by="labels_complete2", alpha = 0.5)
DimPlot(sample250529_UMAP.B.less.WT, reduction = "umap", label=T, group.by="labels_complete2", alpha = 0.5)
DimPlot(sample250529_UMAP.B.less.KO, reduction = "umap", label=T, group.by="labels_complete2", alpha = 0.5)
DimPlot(sample250529_UMAP.B.less_only.MP.MBC, reduction = "umap", label=T, group.by="labels_complete2", alpha = 0.5)
dev.off()

#4.查看关键基因的UMAP图与小提琴图
#labels <- c("B1-5_GC LZ", "B1-5_GC MP", "B1-5_maybe MBC ?","B1-5_MBC",
#            "B6-10_GC LZ", "B6-10_GC MP", "B6-10_maybe MBC ?", "B6-10_MBC")
#combination_list <- combn(labels, 2, simplify = FALSE) #生成比对list
six_of_important_genes <-  c("Il9r","Fas","Ptpn22","Tlr7","Zbtb18","Zbtb32")
pdf("sample250529_UMAP.B.less_symbol.genes.pdf", width=8, height=10)
FeaturePlot(sample250529_UMAP.B.less, 
            features = six_of_important_genes) #重要基因
VlnPlot(sample250529_UMAP.B.less_only.MP.MBC, features=c("Il9r","Fas","Ptpn22","Tlr7","Zbtb18","Zbtb32"),
        group.by="labels_complete2",
        stack =T, flip =T
        )
  # stat_compare_means(method = "wilcox.test", # 检验方法：wilcox.test, t.test, anova等
  #       comparisons = combination_list, # 指定要比较的组
  #       label = "p.signif")+
  # scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
dev.off()

#5.导出小提琴图数据单独绘制方差分析结果
data1 <- VlnPlot(sample250529_UMAP.B.less_only.MP.MBC, features=c("Il9r","Fas","Ptpn22","Tlr7","Zbtb18","Zbtb32"),
        group.by="labels_complete2",stack =T, flip =T)$data
write_csv(data1, file="six_of_important_genes_expression2.csv")

#6.按簇分WT与KO（即按label_complete2分群），再提取相近的细胞亚群做差异markers分析
#6.1 划分1.2.4.11等可疑GC LZ作进一步分析
sample250529_UMAP.B.less_only.1.2.4.11 <- sample250529_UMAP.B.less_only.MP.MBC[,sample250529_UMAP.B.less_only.MP.MBC$seurat_clusters %in% c(1,2,4,11)]
sample250529_UMAP.B.less_only.2.4.11 <- sample250529_UMAP.B.less_only.MP.MBC[,sample250529_UMAP.B.less_only.MP.MBC$seurat_clusters %in% c(2,4,11)]
sample250529_UMAP.B.less_only.2.4.11_wt <- sample250529_UMAP.B.less_only.2.4.11[,sample250529_UMAP.B.less_only.2.4.11$label.main %in% "B1-5"]
sample250529_UMAP.B.less_only.2.4.11_ko <- sample250529_UMAP.B.less_only.2.4.11[,sample250529_UMAP.B.less_only.2.4.11$label.main %in% "B6-10"]
#6.2 find markers启动：
sample250529_markers_wt.Vs.ko_strict <- FindAllMarkers(sample250529_UMAP.B.less, group.by="labels_complete2",only.pos=F, logfc.threshold=1, min.pct=0.5, min.diff.pct=0.3) 
sample250529_markers_only.MP.MBC_wt.Vs.ko_lenient <- FindAllMarkers(sample250529_UMAP.B.less_only.MP.MBC, group.by="labels_complete2", only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2)
sample250529_markers_only.1.2.4.11_lenient <- FindAllMarkers(sample250529_UMAP.B.less_only.1.2.4.11,  only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2)
sample250529_markers_only.1.2.4.11_wt.Vs.ko_lenient <- FindAllMarkers(sample250529_UMAP.B.less_only.1.2.4.11, group.by="labels_complete2", only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2)
sample250529_markers_only.2.4.11_wt.Vs.ko_lenient <- FindAllMarkers(sample250529_UMAP.B.less_only.2.4.11, group.by="labels_complete2", only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2)
sample250529_markers_only.2.4.11_wt_lenient <- FindAllMarkers(sample250529_UMAP.B.less_only.2.4.11_wt, group.by="labels_complete2", only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2)
sample250529_markers_only.2.4.11_ko_lenient <- FindAllMarkers(sample250529_UMAP.B.less_only.2.4.11_ko, group.by="labels_complete2", only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2)
sample250529_markers_cluster2.Vs.cluster1 <- FindMarkers(sample250529_UMAP.B.less_only.MP.MBC, ident.1=2, ident.2=1, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 
sample250529_markers_cluster2.Vs.cluster4 <- FindMarkers(sample250529_UMAP.B.less_only.MP.MBC, ident.1=2, ident.2=4, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 
sample250529_markers_cluster2.Vs.cluster11 <- FindMarkers(sample250529_UMAP.B.less_only.MP.MBC, ident.1=2, ident.2=11, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 
sample250529_markers_cluster4.Vs.cluster1 <- FindMarkers(sample250529_UMAP.B.less_only.MP.MBC, ident.1=4, ident.2=1, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 
sample250529_markers_cluster4.Vs.cluster11 <- FindMarkers(sample250529_UMAP.B.less_only.MP.MBC, ident.1=4, ident.2=11, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 
#6.3 保存结果
lists_markers <- list("all_wt.Vs.ko_strict" = sample250529_markers_wt.Vs.ko_strict ,
                      "LZ+MBC_wt.Vs.ko_lenient" = sample250529_markers_only.MP.MBC_wt.Vs.ko_lenient,
                      "1.2.4.11" = sample250529_markers_only.1.2.4.11_lenient,
                      "1.2.4.11_wt.Vs.ko" = sample250529_markers_only.1.2.4.11_wt.Vs.ko_lenient,
                      "2.4.11_wt.Vs.ko" = sample250529_markers_only.2.4.11_wt.Vs.ko_lenient,
                      "2.4.11_wt" = sample250529_markers_only.2.4.11_wt_lenient,
                      "2.4.11_ko" = sample250529_markers_only.2.4.11_ko_lenient,
                      "2.Vs.1" = sample250529_markers_cluster2.Vs.cluster1,
                      "2.Vs.4" = sample250529_markers_cluster2.Vs.cluster4,
                      "2.Vs.11" = sample250529_markers_cluster2.Vs.cluster11,
                      "4.Vs.1" = sample250529_markers_cluster4.Vs.cluster1,
                      "4.Vs.11" = sample250529_markers_cluster4.Vs.cluster11)
write.xlsx(lists_markers, file="sample250529_wt.Vs.ko_SCT_markers.xlsx", rowNames=T)#添加sheet名并保留行名



## （二）提取MBC与LZ的可疑群体，进行有监督注释 -----------------------------------------------------------------
### 加载数据，设置环境  -----------------------------------------------------------------
rm(list=ls())
load(file="sample250529_UMAP.B.less2.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存
### 有监督注释 -----------------------------------------------------------------
#1.提取亚群
sample250529_UMAP.B.less_sub.MBC <- sample250529_UMAP.B.less[,sample250529_UMAP.B.less@meta.data$seurat_clusters %in% c("3", "6", "13")]
sample250529_UMAP.B.less_sub.LZ <- sample250529_UMAP.B.less[,sample250529_UMAP.B.less@meta.data$seurat_clusters %in% c("1", "2", "4", "11")]
sample250529_UMAP.B.less_sub.MBC.LZ <- sample250529_UMAP.B.less[,sample250529_UMAP.B.less@meta.data$seurat_clusters %in% c("3", "6", "1", "2", "4", "11")]
sample250529_UMAP.B.less_sub.LZ.and.WT <- sample250529_UMAP.B.less_sub.LZ[,sample250529_UMAP.B.less_sub.LZ@meta.data$label.main %in% c("B1-5")]
sample250529_UMAP.B.less_sub.LZ.and.KO <- sample250529_UMAP.B.less_sub.LZ[,sample250529_UMAP.B.less_sub.LZ@meta.data$label.main %in% c("B6-10")]

#2.有监督注释,参考文献ref2
#《Germinal-center development of memory B cells driven by IL-9 from follicular helper T cells》f3a
#2.1 导入参考markers，区分GC MP与LZ用
sample250529_markers_ref2 <- c("Bcl2","Bmf","Ccnb2","Cdk2ap1","Fas","Timeless", #Apoptosis cell cycle
                               "Ccr6","Ccr7","S1pr1","S1pr2","S1pr4","Sell","Gpr183", #Localization
                               "Cd38","Cd55","Cd97","Ly6a","Ly6d", #Surface markers
                               "Il9r","Nlrx1","Mx1","Tlr1","Tlr7", #Immune response
                               "Bach2","Klf2","Klf3","Runx3","Stat1","Zbtb32" #Transcription factors
                               )
#2.2 绘图
#2.2.1点图
pdf("sample250529_UMAP.B.less_ref2_gene.expression.dotplot.pdf", width = 14, height = 30, useDingbats = FALSE)
DotPlot(sample250529_UMAP.B.less, features = rev(sample250529_markers_ref2), group.by = "labels_complete2") + 
  ggplot2::scale_color_gradient2(
    low = "#2E9FDF",mid = "white",high = "#FC4E07",
    midpoint = 0,limits = c(-2.5, 2.5)
  ) + #按需求填色
  theme_minimal() + #画上格线方便阅读
  coord_flip() #坐标翻转
DotPlot(sample250529_UMAP.B.less_sub.LZ, features = rev(sample250529_markers_ref2), group.by = "labels_complete2") + 
  ggplot2::scale_color_gradient2(
    low = "#2E9FDF",mid = "white",high = "#FC4E07",
    midpoint = 0,limits = c(-2.5, 2.5)
  ) + #按需求填色
  theme_minimal() + #画上格线方便阅读
  coord_flip() #坐标翻转
dev.off()
#2.2.2热图
#2.2.2.1 总体
#按细胞绘制的热图
DoHeatmap(sample250529_UMAP.B.less_sub.LZ, features = sample250529_markers_ref2, group.by = "labels_complete2")+
  scale_fill_gradientn(colors = c("navy","white","firebrick3"))
#按样品绘制的热图
sample250529_UMAP.B.less_sub.LZ$orig.ident_clusters <- paste0(sample250529_UMAP.B.less_sub.LZ$seurat_clusters, "_", sample250529_UMAP.B.less_sub.LZ$orig.ident)
data1 <- DotPlot(sample250529_UMAP.B.less_sub.LZ, features = sample250529_markers_ref2, group.by = "orig.ident_clusters")$data
data1_wide <- dcast(data1, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data1_wide) <- data1_wide$features.plot
data1_wide <- data1_wide[,-1]
p1 <- pheatmap(data1_wide,cluster_rows = FALSE,cluster_cols = FALSE)
#按组别绘制的热图
data2 <- DotPlot(sample250529_UMAP.B.less_sub.LZ, features = sample250529_markers_ref2, group.by = "labels_complete2")$data
data2_wide <- dcast(data2, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data2_wide) <- data2_wide$features.plot
data2_wide <- data2_wide[,-1]
p2 <- pheatmap(data2_wide,cluster_rows = FALSE,cluster_cols = FALSE)
#2.2.2.2 WT
#按样品绘制的热图
sample250529_UMAP.B.less_sub.LZ.and.WT$orig.ident_clusters <- paste0(sample250529_UMAP.B.less_sub.LZ.and.WT$seurat_clusters, "_", sample250529_UMAP.B.less_sub.LZ.and.WT$orig.ident)
data1 <- DotPlot(sample250529_UMAP.B.less_sub.LZ.and.WT, features = sample250529_markers_ref2, group.by = "orig.ident_clusters")$data
data1_wide <- dcast(data1, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data1_wide) <- data1_wide$features.plot
data1_wide <- data1_wide[,-1]
p1.WT <- pheatmap(data1_wide,cluster_rows = FALSE,cluster_cols = FALSE)
#按组别绘制的热图
data2 <- DotPlot(sample250529_UMAP.B.less_sub.LZ.and.WT, features = sample250529_markers_ref2, group.by = "labels_complete2")$data
data2_wide <- dcast(data2, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data2_wide) <- data2_wide$features.plot
data2_wide <- data2_wide[,-1]
p2.WT <- pheatmap(data2_wide,cluster_rows = FALSE)
#2.2.2.3 KO
#按样品绘制的热图
sample250529_UMAP.B.less_sub.LZ.and.KO$orig.ident_clusters <- paste0(sample250529_UMAP.B.less_sub.LZ.and.KO$seurat_clusters, "_", sample250529_UMAP.B.less_sub.LZ.and.KO$orig.ident)
data1 <- DotPlot(sample250529_UMAP.B.less_sub.LZ.and.KO, features = sample250529_markers_ref2, group.by = "orig.ident_clusters")$data
data1_wide <- dcast(data1, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data1_wide) <- data1_wide$features.plot
data1_wide <- data1_wide[,-1]
p1.KO <- pheatmap(data1_wide,cluster_rows = FALSE,cluster_cols = FALSE)
#按组别绘制的热图
data2 <- DotPlot(sample250529_UMAP.B.less_sub.LZ.and.KO, features = sample250529_markers_ref2, group.by = "labels_complete2")$data
data2_wide <- dcast(data2, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data2_wide) <- data2_wide$features.plot
data2_wide <- data2_wide[,-1]
p2.KO <- pheatmap(data2_wide,cluster_rows = FALSE)
#2.2.2.4 汇总保存
#将热图转换为gtable对象
g1 <- p1$gtable
g2 <- p2$gtable
g1.WT <- p1.WT$gtable
g2.WT <- p2.WT$gtable
g1.KO <- p1.KO$gtable
g2.KO <- p2.KO$gtable
layout <- rbind(c(1),c(2, 3))
pdf("sample250529_UMAP.B.less_sub.LZ_ref2_gene.expression.pheatmap.pdf", width = 30, height = 30, useDingbats = FALSE)
grid.arrange(g1, g1.WT, g1.KO, ncol = 3, layout_matrix = layout)
grid.arrange(g2, g2.WT, g2.KO, ncol = 3, layout_matrix = layout)
dev.off()

#3.有监督注释,参考文献ref3
#《CR6 Defines Memory B Cell Precursors in Mouse and Human Germinal Centers, 
# Revealing Light-Zone Location and Predominant Low Antigen Affinity》
#3.1 导入参考markers
sample250529_markers_ref3.1 <- c("Ccr6", #关键marker
                                 "Plp2", "Plac8", "Maml2", "Plxnc1", "Tbc1d9",
                                 "Il23a", "Rnase6", "P2ry10", "Faim3",
                                 "Celf2", "Tmem66", #上面应为MBC偏阳性的情况
                                 "Abca6", "Ndfip2", "Prrg4", "Ccdc144cp", "Rapgef5",
                                 "Aff2", "Cd38", "Mybl1", "Spred2", "Nuf2",
                                 "Chek1", "Exo1", "Mlf1ip", "Ctnnal1", "Cdk1",
                                 "Ska3", "Gmnn", "Dscc1", "Mnd1", "Ccna2",
                                 "Kif11", "Ndc80", "Mki67", "Nusap1", "Cd10",
                                 "Casc5", "Clspn", "Tyms") #f7a的marker（Ccr6-LZ与MBC差异最大）
#3.1.1点图
pdf("sample250529_UMAP.B.less_sub.MBC.LZ_ref3.1_gene.expression.dotplot[Ccr6-LZ.Vs.MBC].pdf", width = 14, height = 30, useDingbats = FALSE)
DotPlot(sample250529_UMAP.B.less_sub.MBC.LZ, features = rev(sample250529_markers_ref3.1), group.by = "labels_complete2") + 
  ggplot2::scale_color_gradient2(
    low = "#2E9FDF",mid = "white",high = "#FC4E07",
    midpoint = 0,limits = c(-2.5, 2.5)
  ) + #按需求填色
  theme_minimal() + #画上格线方便阅读
  coord_flip() #坐标翻转
dev.off()
#3.1.2按细胞簇+样品绘制的热图
sample250529_UMAP.B.less_sub.MBC.LZ$orig.ident_clusters <- paste0(sample250529_UMAP.B.less_sub.MBC.LZ$seurat_clusters, "_", sample250529_UMAP.B.less_sub.MBC.LZ$orig.ident)
data1 <- DotPlot(sample250529_UMAP.B.less_sub.MBC.LZ, features = sample250529_markers_ref3.1, group.by = "orig.ident_clusters")$data
data1_wide <- dcast(data1, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data1_wide) <- data1_wide$features.plot
data1_wide <- data1_wide[,-1]
p1_new <- pheatmap(data1_wide,cluster_rows = FALSE,cluster_cols = FALSE)
#3.1.3按细胞簇+组别绘制的热图
data2 <- DotPlot(sample250529_UMAP.B.less_sub.MBC.LZ, features = sample250529_markers_ref3.1, group.by = "labels_complete2")$data
data2_wide <- dcast(data2, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data2_wide) <- data2_wide$features.plot
data2_wide <- data2_wide[,-1]
p2_new <- pheatmap(data2_wide,cluster_rows = FALSE)
#3.1.4将两个热图绘图保存
pdf("sample250529_UMAP.B.less_sub.MBC.LZ_ref3.1_gene.expression.pheatmap[Ccr6-LZ.Vs.MBC].pdf", width = 15, height = 10)
p1_new
grid.newpage()
p2_new
dev.off()

#3.2 导入参考markers
sample250529_markers_ref3.2 <- c("Ccr6", #关键marker
                                 "Cdk1", "Cysltr1", "Camk1d", "Zmynd11", "S1pr1",
                                 "Col19a1", "Slc12a6", "Rnase6", "P2ry10", "Fam177b",
                                 "Zeb2", "Tbc1d9", "Maml2", "Plac8", "Tmem66","Faim3",
                                 "Chi3l2", "Celf2", "Capn3", "Mtss1", "Mir3139",
                                 "Aff1", "Utrn", "App", "Rp9p") #f7b的marker（Ccr6-LZ与Ccr6+LZ差异最大）
#3.2.1 点图
pdf("sample250529_UMAP.B.less_sub.MBC.LZ_ref3.2_gene.expression.dotplot[Ccr6-LZ.Vs.[Ccr6+LZ].pdf", width = 14, height = 30, useDingbats = FALSE)
DotPlot(sample250529_UMAP.B.less_sub.MBC.LZ, features = rev(sample250529_markers_ref3.2), group.by = "labels_complete2") + 
  ggplot2::scale_color_gradient2(
    low = "#2E9FDF",mid = "white",high = "#FC4E07",
    midpoint = 0,limits = c(-2.5, 2.5)
  ) + #按需求填色
  theme_minimal() + #画上格线方便阅读
  coord_flip() #坐标翻转
dev.off()
#3.2.2 按样品绘制的热图
sample250529_UMAP.B.less_sub.MBC.LZ$orig.ident_clusters <- paste0(sample250529_UMAP.B.less_sub.MBC.LZ$seurat_clusters, "_", sample250529_UMAP.B.less_sub.MBC.LZ$orig.ident)
data1 <- DotPlot(sample250529_UMAP.B.less_sub.MBC.LZ, features = sample250529_markers_ref3.2, group.by = "orig.ident_clusters")$data
data1_wide <- dcast(data1, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data1_wide) <- data1_wide$features.plot
data1_wide <- data1_wide[,-1]
p1_new <- pheatmap(data1_wide,cluster_rows = FALSE,cluster_cols = FALSE)
#3.2.3 按组别绘制的热图
data2 <- DotPlot(sample250529_UMAP.B.less_sub.MBC.LZ, features = sample250529_markers_ref3.2, group.by = "labels_complete2")$data
data2_wide <- dcast(data2, features.plot ~ id, value.var = "avg.exp.scaled")
row.names(data2_wide) <- data2_wide$features.plot
data2_wide <- data2_wide[,-1]
p2_new <- pheatmap(data2_wide,cluster_rows = FALSE#,cluster_cols = FALSE
                   )
#3.2.4 将总体的两个热图绘图保存
pdf("sample250529_UMAP.B.less_sub.LZ_ref3.2_gene.expression.pheatmap[Ccr6-LZ.Vs.[Ccr6+LZ].pdf", width = 15, height = 10)
p1_new
grid.newpage()
p2_new
dev.off()




# 十一、拟时序分析 -------------------------------------------------------------------
## 实际上monocle3包能够提供完整的单细胞分析流程，也与前面一样需要质控、批次效应检查、粗降维、细降维的操作；
## 不过这一部分数据用Seurat包生成的结果导入，所以我们就忽略了这一包带来的运行问题。只关注我们最关注的内容。
## 拟时序基础设置 -------------------------------------------------------------------
##设置环境
rm(list=ls())
load(file="sample250529_UMAP.B.less2.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存
## 获取基础数据
sample250529_UMAP.B.less_df <- GetAssayData(sample250529_UMAP.B.less, assay="SCT", layer="counts")  #一个数值矩阵，其中行代表基因，列代表细胞
sample250529_cell_metadata <- sample250529_UMAP.B.less@meta.data #行是细胞，列是细胞属性（如细胞类型、培养条件、捕获日期等）
sample250529_gene_annotation <- data.frame(gene_short_name = rownames(sample250529_UMAP.B.less_df), row.names = rownames(sample250529_UMAP.B.less_df)) #一个数据框，其中行是特征（例如基因），列是基因属性，如生物类型、GC 含量等；应有gene_short_name这一列用于绘图
## 创建一个cell_data_set对象
sample250529_cds <- new_cell_data_set(sample250529_UMAP.B.less_df, cell_metadata=sample250529_cell_metadata, gene_metadata = sample250529_gene_annotation) #预处理
sample250529_cds <- preprocess_cds(sample250529_cds,method ="PCA") #粗降维PCA
sample250529_cds <- reduce_dimension(sample250529_cds, reduction_method="UMAP", preprocess_method="PCA") #降维
## 将seurat对象的UMAP导入
sample250529.B.less_int.embed <- Embeddings(sample250529_UMAP.B.less, reduction ="umap") #导入UMAP坐标点位置
sample250529.B.less_int.embed <- sample250529.B.less_int.embed[rownames(sample250529_cds@int_colData$reducedDims$UMAP),] #不知道这一步的意义何在？
sample250529_cds@int_colData$reducedDims$UMAP <- sample250529.B.less_int.embed #传入UMAP点的坐标数据

## 聚类、构建轨迹、选定起点 -------------------------------------------------------------------
## 聚类分区，不同分区细胞需要分出来，分别构建单独的轨迹分析
#聚类
sample250529_cds <- cluster_cells(sample250529_cds)
#绘图确认Partition唯一
sample250529.B.less_labels.complete2_monocle3_default <- plot_cells(sample250529_cds, reduction_method="UMAP", color_cells_by="labels_complete2", show_trajectory_graph=F) + 
  ggtitle("monole3_umap") #绘制默认图像
sample250529.B.less_labels.complete2_monocle3_cluster.Partition <- plot_cells(sample250529_cds, 
                                                                              color_cells_by = "partition", 
                                                                              show_trajectory_graph = FALSE) + 
  ggtitle("partitic") #确认Partition唯一，这保证细胞群体相对统一，有利于构建细胞轨迹
ggsave("sample250529_UMAP.B.less_labels.complete2_cluster.Partition.pdf", plot = sample250529.B.less_labels.complete2_monocle3_cluster.Partition, width = 10, height = 8)
## 构建细胞轨迹，此步较久
#构建细胞轨迹，此步较久
sample250529_cds <- learn_graph(sample250529_cds)
sample250529_cds1 <- sample250529_cds
sample250529_cds1 <- learn_graph(sample250529_cds1, learn_graph_control = list(euclidean_distance_ratio = 1.5, geodesic_distance_ratio = 0.2, nn.k = 35))
sample250529_cds2 <- sample250529_cds
sample250529_cds2 <- learn_graph(sample250529_cds2, learn_graph_control = list(euclidean_distance_ratio = 1.5, geodesic_distance_ratio = 0.2, nn.k = 35, orthogonal_proj_tip=TRUE, nn.metric="cosine"))
sample250529_cds3 <- sample250529_cds
sample250529_cds3 <- learn_graph(sample250529_cds3, learn_graph_control = list(euclidean_distance_ratio = 0.8, geodesic_distance_ratio = 0.8))
#画轨迹图
pdf("sample250529_UMAP.B.less_monocle3_learn.graph.pdf", width=8, height=10)
sample250529.B.less_labels.complete2_monocle3_learn.graph <- plot_cells(sample250529_cds, color_cells_by = "labels_complete2", label_groups_by_cluster = FALSE, 
                                                                        label_leaves = FALSE, label_branch_points=FALSE)
sample250529.B.less_labels.complete2_monocle3_learn.graph1 <- plot_cells(sample250529_cds1, color_cells_by = "labels_complete2", label_groups_by_cluster = FALSE, 
                                                                        label_leaves = FALSE, label_branch_points=FALSE)
sample250529.B.less_labels.complete2_monocle3_learn.graph2 <- plot_cells(sample250529_cds2, color_cells_by = "labels_complete2", label_groups_by_cluster = FALSE, 
                                                                         label_leaves = FALSE, label_branch_points=FALSE)
sample250529.B.less_labels.complete2_monocle3_learn.graph3 <- plot_cells(sample250529_cds3, color_cells_by = "labels_complete2", label_groups_by_cluster = FALSE, 
                                                                         label_leaves = FALSE, label_branch_points=FALSE)
print(sample250529.B.less_labels.complete2_monocle3_learn.graph)
print(sample250529.B.less_labels.complete2_monocle3_learn.graph1)
print(sample250529.B.less_labels.complete2_monocle3_learn.graph2)
print(sample250529.B.less_labels.complete2_monocle3_learn.graph3)
dev.off()
## 发育轨迹(拟时序)排列细胞。需要手动标注顺序。[选定起点]
sample250529_cds <- order_cells(sample250529_cds)
sample250529_cds1 <- order_cells(sample250529_cds1)
sample250529_cds2 <- order_cells(sample250529_cds2)
sample250529_cds3 <- order_cells(sample250529_cds3)
pdf("sample250529_UMAP.B.less_monocle3_pseudotime.pdf", width=8, height=10)
plot_cells(sample250529_cds, color_cells_by="pseudotime", label_cell_groups = FALSE,
             label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(sample250529_cds1, color_cells_by="pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(sample250529_cds2, color_cells_by="pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)
plot_cells(sample250529_cds3, color_cells_by="pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)
dev.off()
## 一种提取子集分析的可视化窗口操作方法：https://cole-trapnell-lab.github.io/monocle3/docs/trajectories/
sample250529_cds_sub <- choose_graph_segments(sample250529_cds)
## 三维绘图
#PCA降维绘图试试
sample250529_cds_3d <- reduce_dimension(sample250529_cds, max_components = 3) #降维到3，默认PCA（更合理），用reduction="UMAP"的问题是，这会与前面的绘图结果不太一致，得研究下怎么消除
sample250529_cds_3d <- cluster_cells(sample250529_cds_3d)
sample250529_cds_3d <- learn_graph(sample250529_cds_3d)
#sample250529_cds_3d <- order_cells(sample250529_cds_3d, root_pr_nodes=get_earliest_principal_node(sample250529_cds)) #官网说是指定起始点用的，不过似乎会报错
n_groups <- length(unique(colData(sample250529_cds_3d)$labels_complete2))
sample250529_cds_3d_plot_obj <- plot_cells_3d(sample250529_cds_3d, color_cells_by="labels_complete2", color_palette = viridis::viridis(n_groups))
print(sample250529_cds_3d_plot_obj)

##保存数据
save(sample250529_cds,sample250529_cds1,sample250529_cds2,sample250529_cds3, 
     file="sample250529_UMAP.B.less_cds.RData")
save(sample250529_cds_3d_plot_obj, file="sample250529_UMAP.B.less_cds.PCA.3d.plot_obj.RData")




## 寻找拟时轨迹差异基因 -------------------------------------------------------------------
# 此步是根据低维嵌入和主图，检测基因的差异表达。官网给出了两种方法去检测差异表达：
#   1.回归分析，即使用 fit_models()，更适用于临床数据，检测差异表达是否依赖于时间、治疗等变量；
#   2.图自相关分析[Graph-autocorrelation analysis]：使用 graph_test()，检测细胞簇的变化轨迹。

# 这里主要介绍后者，graph_test()实质上是借鉴了空问相关性分析中的一项强大技术--莫兰I检验[Moran's I test]。
#   莫兰I检验是一种多向、多维空间自相关的测量方法。
#   该统计量表明，轨迹上邻近位置的细胞对被测基因的表达水平是否相似或不同。
#   虽然莫兰I与皮尔逊相关性（pearson）的范围都在-1到1之间，但意义略有不同：
#     +1表示邻近细胞的表达完全相似；
#     0表示没有相关性；
#    -1表示邻近细胞没有相关性。

##找出轨迹差异基因
Track_genes <- graph_test(sample250529_cds, neighbor_graph="principal_graph", cores=8) #差异统计，此步非常久（单核运行，七万细胞需要2个小时）；服务器上核数多点运行更快；Windows系统只能单核运行
track.genes.strict <- row.names(subset(Track_genes, morans_I>0.25)) #阈值筛选，得到531个基因（似乎和上面有点对不上，差了100多个基因）
write.csv(Track_genes,"sample250529_UMAP.B.less_Track.genes.csv",row.names = F) #保存

##时序-基因聚类热图
#生成时序特征显著基因矩阵
num_clusters=3 #设置下面的K-means聚类的clusters数目
pt.matrix <- exprs(sample250529_cds)[match(track.genes.strict, rownames(rowData(sample250529_cds))),
                                     order(pseudotime(sample250529_cds))] #核心问题：前面的起点设置
pt.matrix <- t(apply(pt.matrix, 1, function(x){smooth.spline(x,df=3)$y})) #设置自由度参数使得平滑
pt.matrix <- t(apply(pt.matrix, 1, function(x){(x-mean(x))/sd(x)})) #对每个基因的表达值进行标准化
rownames(pt.matrix) <- track.genes.strict
#以时序为x轴，以时序特征显著基因为y轴，绘制基因聚类热图
color_heatmap <- rev(RColorBrewer::brewer.pal(11,"Spectral")) #画出来的色板
p <- ComplexHeatmap::Heatmap(
  pt.matrix, name = "z-score", show_row_names = T, show_column_names = F,
  col = circlize::colorRamp2(seq(from=-2,to=2,length=11),color_heatmap),
  row_names_gp = gpar(fontsize =6), row_title_rot= 0, km = num_clusters,
  cluster_rows = TRUE, #同时指定km和cluster_rows=TRUE，则使用k-means聚类；没有km，则使用层次聚类。
  cluster_row_slices = FALSE, #多个基因的横行布局类，即不重新排序，保持按簇1,2,3的顺序排列
  cluster_columns = FALSE, #列已经按伪时间排序了，所以不用对列进行聚类
  use_raster = TRUE #光栅化，F会让绘图变慢但更精细，没做挑选的话我觉得没必要
  )
pdf("sample250529_UMAP.B.less_Pseudotime.Genes.heatmap.pdf", width=8, height =30)
plot(p)
dev.off()

##在细胞簇注色情况下，单个基因的时序-基因表达图
#以时序为x轴，以基因表达情况为y轴，但散点横坐标可以为各簇细胞的seurat_clusters，相当于是指明在什么时间下什么细胞群体会高度表达该基因
#绘图需要的时间很长
morans.I_genes_top12 <- rownames(top_n(Track_genes, n=12, morans_I)) #基因展示，选取morans_I排名前12的基因，可自定义
morans.I_genes_top6 <- morans.I_genes_top12[1:6]
morans.I_genes_second6 <- morans.I_genes_top12[7:12]
six_of_important_genes <-  c("Il9r","Fas","Ptpn22","Tlr7","Zbtb18","Zbtb32")
pdf("sample250529_UMAP.B.less_Pseudotime.Genes.Jitterplot.pdf", width=8, height =8)
plot_genes_in_pseudotime(sample250529_cds[morans.I_genes_top6,], color_cells_by="seurat_clusters", min_expr=0.5, ncol=3)
plot_genes_in_pseudotime(sample250529_cds[morans.I_genes_second6,], color_cells_by="seurat_clusters", min_expr=0.5, ncol=3)
plot_genes_in_pseudotime(sample250529_cds[six_of_important_genes,], color_cells_by="seurat_clusters", min_expr=0.5, ncol=3)
dev.off()

##共表达基因模块图
#获取基因模块
morans.I.0.1_genelist <- row.names(subset(Track_genes, morans_I>0.1)) #共表达的阈值可以更宽松一点，相较于前面的时序绘图而言，这里找到1762个
gene.module_df <- find_gene_modules(sample250529_cds[morans.I.0.1_genelist,], resolution=1e-2, cores=6)
table(gene.module_df$module)
write.csv(gene.module_df,"sample250529_UMAP.B.less_Pseudotime.Genes.Module.csv", row.names = F)
#共表达基因模块热图
cell.group_df <- tibble(cell=row.names(colData(sample250529_cds)),
                     cell_group=colData(sample250529_cds)$labels_complete2)
heatmap_df <- aggregate_gene_expression(sample250529_cds, gene.module_df, cell.group_df)
row.names(heatmap_df) <- str_c("Module", row.names(heatmap_df)) #将Module列改为行名
module_pheatmap <- pheatmap::pheatmap(heatmap_df, scale="column", clustering_method="ward.D2")
ggsave("sample250529_UMAP.B.less_Pseudotime.GenesModule.pdf", plot = module_pheatmap, width =8, height = 6.5)
#散点图
gene.module_list <- lapply(1:length(unique(gene.module_df$module)), function(i){subset(gene.module_df, module==i)$id})
names(gene.module_list)<- paste0("module", 1:length(unique(gene.module_df$module))) #统一list子集的命名
load(file="sample250529_UMAP.B.less.RData")
sample250529_UMAP.B.less <- AddModuleScore(sample250529_UMAP.B.less, features=gene.module_list, name="module")
module_dotplot <- FeaturePlot(sample250529_UMAP.B.less, features = paste0("module", 1:length(unique(gene.module_df$module))), ncol=2) + 
  plot_layout() & scale_color_viridis_c(option ="C")
ggsave("sample250529_UMAP.B.less_PseudotimeGenes_Moduescore.pdf", plot = module_dotplot, width = 10, height = 80, limitsize = FALSE)
#针对各个模块的基因进行KEGG通路富集分析


# 针对WT单独做clusters比对、拟时序分析 -----------------------------------------------------------
## 拟时序基础设置 -------------------------------------------------------------------
##设置环境
rm(list=ls())
load(file="sample250529_UMAP.B.less_WT.and.KO.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存
## 获取基础数据
sample250529_UMAP.B.less.WT_df <- GetAssayData(sample250529_UMAP.B.less.WT, assay="SCT", layer="counts")  #一个数值矩阵，其中行代表基因，列代表细胞
sample250529_cell_metadata <- sample250529_UMAP.B.less.WT@meta.data #行是细胞，列是细胞属性（如细胞类型、培养条件、捕获日期等）
sample250529_gene_annotation <- data.frame(gene_short_name = rownames(sample250529_UMAP.B.less.WT_df), row.names = rownames(sample250529_UMAP.B.less.WT_df)) #一个数据框，其中行是特征（例如基因），列是基因属性，如生物类型、GC 含量等；应有gene_short_name这一列用于绘图
## 创建一个cell_data_set对象
sample250529_WTcds <- new_cell_data_set(sample250529_UMAP.B.less.WT_df, cell_metadata=sample250529_cell_metadata, gene_metadata = sample250529_gene_annotation) #预处理
sample250529_WTcds <- preprocess_cds(sample250529_WTcds,method ="PCA") #粗降维PCA
sample250529_WTcds <- reduce_dimension(sample250529_WTcds, reduction_method="UMAP", preprocess_method="PCA") #降维
## 将seurat对象的UMAP导入
sample250529.B.less.WT_int.embed <- Embeddings(sample250529_UMAP.B.less.WT, reduction ="umap") #导入UMAP坐标点位置
sample250529.B.less.WT_int.embed <- sample250529.B.less.WT_int.embed[rownames(sample250529_WTcds@int_colData$reducedDims$UMAP),] #不知道这一步的意义何在？
sample250529_WTcds@int_colData$reducedDims$UMAP <- sample250529.B.less.WT_int.embed #传入UMAP点的坐标数据



## 聚类、构建轨迹、选定起点 -------------------------------------------------------------------
## 聚类分区，不同分区细胞需要分出来，分别构建单独的轨迹分析
#聚类
sample250529_WTcds <- cluster_cells(sample250529_WTcds)
#绘图确认Partition唯一
sample250529.B.less.WT_labels.complete2_monocle3_default <- plot_cells(sample250529_WTcds, reduction_method="UMAP", color_cells_by="labels_complete2", show_trajectory_graph=F) + 
  ggtitle("monole3_umap") #绘制默认图像
sample250529.B.less.WT_labels.complete2_monocle3_cluster.Partition <- plot_cells(sample250529_WTcds, 
                                                                                 color_cells_by = "partition", 
                                                                                 show_trajectory_graph = FALSE) + 
  ggtitle("partitic") #确认Partition唯一，这保证细胞群体相对统一，有利于构建细胞轨迹
ggsave("sample250529_UMAP.B.less.WT_labels.complete2_cluster.Partition.pdf", plot = sample250529.B.less.WT_labels.complete2_monocle3_cluster.Partition, width = 10, height = 8)
## 构建细胞轨迹，此步较久
#构建细胞轨迹，此步较久
sample250529_WTcds <- learn_graph(sample250529_WTcds)
#画轨迹图
pdf("sample250529_UMAP.B.less.WT_monocle3_learn.graph.pdf", width=8, height=10)
sample250529.B.less.WT_labels.complete2_monocle3_learn.graph <- plot_cells(sample250529_WTcds, color_cells_by = "labels_complete2", label_groups_by_cluster = FALSE, 
                                                                           label_leaves = FALSE, label_branch_points=FALSE)
dev.off()
## 发育轨迹(拟时序)排列细胞。需要手动标注顺序。[选定起点]
sample250529_WTcds <- order_cells(sample250529_WTcds)
pdf("sample250529_UMAP.B.less.WT_monocle3_pseudotime.pdf", width=8, height=10)
plot_cells(sample250529_cds, color_cells_by="pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)
dev.off()
## 三维绘图
#PCA降维绘图试试
sample250529_WTcds_3d <- reduce_dimension(sample250529_WTcds, max_components = 3) #降维到3，默认PCA（更合理），用reduction="UMAP"的问题是，这会与前面的绘图结果不太一致，得研究下怎么消除
sample250529_WTcds_3d <- cluster_cells(sample250529_WTcds_3d)
sample250529_WTcds_3d <- learn_graph(sample250529_WTcds_3d)
n_groups <- length(unique(colData(sample250529_WTcds_3d)$labels_complete2))
sample250529_WTcds_3d_plot_obj <- plot_cells_3d(sample250529_WTcds_3d, color_cells_by="labels_complete2", color_palette = viridis::viridis(n_groups))
print(sample250529_WTcds_3d_plot_obj)

##保存数据
save(sample250529_WTcds,
     file="sample250529_UMAP.B.less.WT_cds.RData")
save(sample250529_WTcds_3d_plot_obj, file="sample250529_UMAP.B.less.WT_cds.PCA.3d.plot_obj.RData")







# 针对KO单独做clusters比对、拟时序分析 -----------------------------------------------------------
## 拟时序基础设置 -------------------------------------------------------------------
##设置环境
rm(list=ls())
load(file="sample250529_UMAP.B.less_KO.and.KO.RData")
plan(multisession, workers=8) #使用8核
options(future.globals.maxSize = 20 * 1024^3) #设置运行内存
## 获取基础数据
sample250529_UMAP.B.less.KO_df <- GetAssayData(sample250529_UMAP.B.less.KO, assay="SCT", layer="counts")  #一个数值矩阵，其中行代表基因，列代表细胞
sample250529_cell_metadata <- sample250529_UMAP.B.less.KO@meta.data #行是细胞，列是细胞属性（如细胞类型、培养条件、捕获日期等）
sample250529_gene_annotation <- data.frame(gene_short_name = rownames(sample250529_UMAP.B.less.KO_df), row.names = rownames(sample250529_UMAP.B.less.KO_df)) #一个数据框，其中行是特征（例如基因），列是基因属性，如生物类型、GC 含量等；应有gene_short_name这一列用于绘图
## 创建一个cell_data_set对象
sample250529_KOcds <- new_cell_data_set(sample250529_UMAP.B.less.KO_df, cell_metadata=sample250529_cell_metadata, gene_metadata = sample250529_gene_annotation) #预处理
sample250529_KOcds <- preprocess_cds(sample250529_KOcds,method ="PCA") #粗降维PCA
sample250529_KOcds <- reduce_dimension(sample250529_KOcds, reduction_method="UMAP", preprocess_method="PCA") #降维
## 将seurat对象的UMAP导入
sample250529.B.less.KO_int.embed <- Embeddings(sample250529_UMAP.B.less.KO, reduction ="umap") #导入UMAP坐标点位置
sample250529.B.less.KO_int.embed <- sample250529.B.less.KO_int.embed[rownames(sample250529_KOcds@int_colData$reducedDims$UMAP),] #不知道这一步的意义何在？
sample250529_KOcds@int_colData$reducedDims$UMAP <- sample250529.B.less.KO_int.embed #传入UMAP点的坐标数据



## 聚类、构建轨迹、选定起点 -------------------------------------------------------------------
## 聚类分区，不同分区细胞需要分出来，分别构建单独的轨迹分析
#聚类
sample250529_KOcds <- cluster_cells(sample250529_KOcds)
#绘图确认Partition唯一
sample250529.B.less.KO_labels.complete2_monocle3_default <- plot_cells(sample250529_KOcds, reduction_method="UMAP", color_cells_by="labels_complete2", show_trajectory_graph=F) + 
  ggtitle("monole3_umap") #绘制默认图像
sample250529.B.less.KO_labels.complete2_monocle3_cluster.Partition <- plot_cells(sample250529_KOcds, 
                                                                                 color_cells_by = "partition", 
                                                                                 show_trajectory_graph = FALSE) + 
  ggtitle("partitic") #确认Partition唯一，这保证细胞群体相对统一，有利于构建细胞轨迹
ggsave("sample250529_UMAP.B.less.KO_labels.complete2_cluster.Partition.pdf", plot = sample250529.B.less.KO_labels.complete2_monocle3_cluster.Partition, width = 10, height = 8)
## 构建细胞轨迹，此步较久
#构建细胞轨迹，此步较久
sample250529_KOcds <- learn_graph(sample250529_KOcds)
#画轨迹图
pdf("sample250529_UMAP.B.less.KO_monocle3_learn.graph.pdf", width=8, height=10)
sample250529.B.less.KO_labels.complete2_monocle3_learn.graph <- plot_cells(sample250529_KOcds, color_cells_by = "labels_complete2", label_groups_by_cluster = FALSE, 
                                                                           label_leaves = FALSE, label_branch_points=FALSE)
dev.off()
## 发育轨迹(拟时序)排列细胞。需要手动标注顺序。[选定起点]
sample250529_KOcds <- order_cells(sample250529_KOcds)
pdf("sample250529_UMAP.B.less.KO_monocle3_pseudotime.pdf", width=8, height=10)
plot_cells(sample250529_cds, color_cells_by="pseudotime", label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE)
dev.off()
## 三维绘图
#PCA降维绘图试试
sample250529_KOcds_3d <- reduce_dimension(sample250529_KOcds, max_components = 3) #降维到3，默认PCA（更合理），用reduction="UMAP"的问题是，这会与前面的绘图结果不太一致，得研究下怎么消除
sample250529_KOcds_3d <- cluster_cells(sample250529_KOcds_3d)
sample250529_KOcds_3d <- learn_graph(sample250529_KOcds_3d)
n_groups <- length(unique(colData(sample250529_KOcds_3d)$labels_complete2))
sample250529_KOcds_3d_plot_obj <- plot_cells_3d(sample250529_KOcds_3d, color_cells_by="labels_complete2", color_palette = viridis::viridis(n_groups))
print(sample250529_KOcds_3d_plot_obj)

##保存数据
save(sample250529_KOcds,
     file="sample250529_UMAP.B.less.KO_cds.RData")
save(sample250529_KOcds_3d_plot_obj, file="sample250529_UMAP.B.less.KO_cds.PCA.3d.plot_obj.RData")




# 提取GCMP进行进一步分析 -----------------------------------------------------------
sample250529_UMAP.B.GCMP = sample250529_UMAP.B.less[,sample250529_UMAP.B.less@meta.data$celltype.main %in% c("GC MP")]
save(sample250529_UMAP.B.GCMP, file="sample250529_UMAP.B.GCMP.RData")
##查看去除干不干净
pdf("sample250529_UMAP.B.GCMP_before.pdf", width=8, height=10)
DimPlot(sample250529_UMAP.B.GCMP, reduction = "umap", label=T)
DimPlot(sample250529_UMAP.B.GCMP, reduction = "umap", label=T, group.by ="orig.ident")
DimPlot(sample250529_UMAP.B.GCMP, reduction = "umap", label=T, group.by ="label.main")
dev.off()
#下游分析
#重新find markers一次
sample250529.GCMP_markers_wt.Vs.ko_lenient <- FindAllMarkers(sample250529_UMAP.B.GCMP, group.by="labels_complete2",only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) #宽松一点的阈值
sample250529_markers_cluster4.Vs.cluster1 <- FindMarkers(sample250529_UMAP.B.GCMP, ident.1=4, ident.2=1, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 
sample250529_markers_cluster4.Vs.cluster2 <- FindMarkers(sample250529_UMAP.B.GCMP, ident.1=4, ident.2=2, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 
sample250529_markers_cluster4.Vs.cluster11 <- FindMarkers(sample250529_UMAP.B.GCMP, ident.1=4, ident.2=11, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 
sample250529_markers_cluster4.Vs.cluster6 <- FindMarkers(sample250529_UMAP.B.GCMP, ident.1=4, ident.2=6, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) 

#看看wt和ko的细胞基因转录组情况
#然后两组的KEGG和GESA来一下
#后续再跑跑几个关键基因的t检验

