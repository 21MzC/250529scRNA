
# 下载R包 --------------------------------------------------------------------
install.packages("rjson")
BiocManager::install("limma", force = TRUE)
install.packages("glmnet")
BiocManager::install("hdf5r")
BiocManager::install('glmGamPoi')
BiocManager::install("SingleR")
BiocManager::install("celldex")
BiocManager::install("harmony")
BiocManager::install("GSVA")
BiocManager::install("AUCell")
install.packages("timeROC")
install.packages("survminer")
devtools::install_github('cole-trapnell-lab/monocle3') #有依赖包下载错误的话，试试BiocManager::install去单独下它的依赖包
remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', force = TRUE)
#####SCENIC更进一步的包，不过我也没装，太麻烦了
## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
devtools::install_github("aertslab/SCENIC") 
#####presto包也可以装，可以加速seurat的使用





# 导入R包 ----------------------------------------------------------------
###读取数据用，默认已安装
library(readr)
library(hdf5r)
library(patchwork)
##单细胞分析的包
library(glmGamPoi) #SeuratV5之后的SCT标准化依赖的包
library(sctransform) #SeuratV5之后的SCT标准化用的包
library(Seurat) #5.2.1，核心包，与4版本的代码有差别，官网https://satijalab.org/seurat/articles/get_started_v5_new
library(monocle3) #拟时序分析
library(SingleR) #自动注释包
library(celldex) #自动注释包的参考数据集单独成包
library(harmony) #批次效应
library(DoubletFinder) #找双细胞
##下游分析
library(rjson)
library(GEOquery)
library(limma) #微阵列和RNA-seq数据的差异表达分析
library(clusterProfiler) #功能富集分析
library(glmnet) #正则计算回归
library(GSVA) #基因集变异分析,单样本富集方法
##数据清洗，默认已安装
library(future) #并行计算
library(stringr) #字符处理
library(dplyr)
library(data.table)
library(Matrix)
library(tidyverse)
##绘图包，默认已安装
library(grid)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(patchwork) #图形排列对其用的，针对ggplot2
library(RColorBrewer)
library(tibble) #排版
##生存分析
library(survival)
library(timeROC)
library(survminer) #生存分析绘图





# 读取数据 --------------------------------------------------------------------
##基础准备
b_list <- list()
b_samples <- c("b2", "b10")
##导入，可改写为循环
b2.data <- Read10X(data.dir="C:/Users/Mingzhe Cai/Desktop/R_resource/250529scRNA/B2") #读取x10数据
b2 <- CreateSeuratObject(counts = b2.data, 
                              project = "b2", #样本名，后续要改循环的地方要用
                              min.cells = 3, 
                              min.features = 200) #创建Seurat对象，删除表达基因数小于200个的细胞，删除小于3个细胞表达的基因（质控）
b10.data <- Read10X(data.dir="C:/Users/Mingzhe Cai/Desktop/R_resource/250529scRNA/B10") #读取x10数据
b10 <- CreateSeuratObject(counts = b10.data, 
                         project = "b10", 
                         min.cells = 3, 
                         min.features = 200) #创建Seurat对象
b_list <- append(b_list, b2)
b_list <- append(b_list, b10)
##合并并保存RData
b_combined_2_10 <- merge(b_list[[1]], y=b_list[-1], add.cell.ids=b_samples ) #在meta.data的行名中添加样品id，
b_combined_new <- JoinLayers(b_combined_2_10) #b_combined_2_10@assays[["RNA"]]@layers有多个层次，不利于后续分析，本步合并
table(b_combined_new@meta.data$orig.ident) #查看样品名
save(b_combined_new, file="combined.RData")






# 质控1：按特征值质控 --------------------------------------------------------------
##清空环境，导入数据
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="combined.RData")
##获取线粒体基因表达情况，MT或者mt都是需要考虑的
b_combined_new[["percent.mt"]] <- PercentageFeatureSet(b_combined_new, pattern="^mt-") #有些会以大写MT开头，这里正则表达需要改一下
table(b_combined_new@meta.data$percent.mt) #检查线粒体情况的步骤；疑惑，质量这么好吗，还是cell ranger自带去除mt的效果？
##创建标准矩阵
b_combined_new <- NormalizeData(b_combined_new)
##质控前绘图 其实个人认为用ggsave储存更方便（方便后期改格式和进PS调整位置）
pdf("quality_control_before.pdf", width=8, height=10)
VlnPlot(b_combined_new, features =c("nFeature_RNA","nCount_RNA","percent.mt"), ncol =3)
FeatureScatter(b_combined_new, "nCount_RNA", "percent.mt", group.by="orig.ident" )
FeatureScatter(b_combined_new, "nCount_RNA", "nFeature_RNA", group.by="orig.ident" )
dev.off() #这次没有明显的双峰，不需要直接去
##数值查看，凭感觉去筛除极端值？这个质控范围怎么定的还需要读读文献
quantile(b_combined_new$nFeature_RNA, seq(0.01,0.1,0.01)) #显示前10%
quantile(b_combined_new$nFeature_RNA, seq(0.9,1,0.01))#显示90%-100% #这个样品nFeature设99%
quantile(b_combined_new$nCount_RNA, seq(0.01,0.1,0.01))
quantile(b_combined_new$nCount_RNA, seq(0.9,1,0.01)) #这里nCount_RNA设98%
quantile(b_combined_new$percent.mt, seq(0.01,0.1,0.01)) 
quantile(b_combined_new$percent.mt, seq(0.9,1,0.01)) #这个样品Mt设95%
##传递质控1后的矩阵
b_qc <- subset(b_combined_new, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 &
                 nCount_RNA < 23000 & percent.mt < 5)
##质控后绘图
pdf("quality_control_after.pdf", width=8, height=10)
VlnPlot(b_qc, features =c("nFeature_RNA","nCount_RNA","percent.mt"), ncol =3)
FeatureScatter(b_qc, "nCount_RNA", "percent.mt", group.by="orig.ident" )
FeatureScatter(b_qc, "nCount_RNA", "nFeature_RNA", group.by="orig.ident" )
dev.off()
##储存生成的文件
save(b_qc, file="qc.RData")





# 质控2：SCT标准化[替代标准化、去心化与高变检索三步] --------------------------------------------
#####质控2：SCT标准化[替代标准化、去心化与高变检索三步]
##清空环境，导入数据
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="qc.RData")
options(future.globals.maxSize = 8 * 1024^3) #设置运行内存
##分步运行
#标准化normalization
b_normalization <- NormalizeData(b_qc, 
                                 normalization.method = "LogNormalize") #标准化方法，包括LogNormalize、CLR、RC三种，默认为LogNormalize
#识别高变特征（特征选择）Identification of highly variable features (feature selection)
b_variable <- FindVariableFeatures(b_normalization, 
                                   selection.method = "vst", 
                                   nfeatures = 2000) #识别变化比较大的基因，默认值是2000（据说这样比较准），但SCT优势就是能找3000个
b_variable_top100 <- head(VariableFeatures(b_variable), 100)
write.table(top100, file="b_variable_top100.txt") #写出前100个突变基因并保存
b_variable_top20 <- head(VariableFeatures(b_variable), 20) #写出前10个突变基因
#绘图：平均表达量-SV变异关系图
pdf("average_VS_SV_top20.pdf", width=8, height=10)
VariableFeaturePlot(b_variable)
LabelPoints(plot = plot1, points = b_variable_top20, repel = TRUE)
dev.off()
#归一化scale
all.genes <- rownames(b_variable)
b_scale <- ScaleData(b_variable, features = all.genes)
save(b_scale, file="scale.RData")
##SCT标准化[实际上，如果不需要绘图的话，10万细胞以内直接使用SCT标准化会更快]
b_sct <- SCTransform(b_qc, 
                     vars.to.regress="percent.mt", #去除percent.mt的列（因为质控过可以不用了）
                     verbose=T) #显示进度条，可以设置成false
##储存生成的文件
save(b_sct, file="sct.RData")






# 线性降维分析（粗降） --------------------------------------------------------------
##清空环境，导入数据
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="sct.RData")
options(future.globals.maxSize = 8 * 1024^3) #设置运行内存
##PCA降维并绘图
b_PCA <- RunPCA(object = b_sct, 
               features = VariableFeatures(object=b_sct)) #设定特征对象，可以换用其他层次，会对结果产生相应的干扰
pdf("PCAPlot.pdf", width=8, height=10)
DimPlot(b_PCA, reduction = "pca")
dev.off()
#####确定维度数的三种方法
##1 指标法，这一指标没有后面两种重要常用
#1.1 确定每个PC 贡献所占百分比（一般应小于5%），用于过滤技术噪音/测序误差
stdev <- b_PCA[["pca"]]@stdev
pct <- stdev^2 / sum(stdev^2) * 100
pct
#1.2 主成分累积贡献（即计算每添加一个PC，累计贡献的百分比变为多少，一般应大于90%），保留主要生物学变异
cumu <-cumsum(pct)
cumu
#1.3 确定信息增益饱和点（一般应大于0.1%），确定信息增益饱和点（肘点）
diff <- c(0, abs(diff(pct)))
diff
##2.生物学信息法
##在特定的维度中如果发现熟悉的marker，则这个维度通常是应该被纳入考量的（有利于下游分析）
##具体而言有三种可视化方法VizDimReduction() 、 DimPlot() 和 DimHeatmap()，主要是展示形式不同，最推荐最后一个，信息量最大
#VizDimReduction()
pdf("PCA_VizDimReduction.pdf", width=16, height=20)
VizDimLoadings(b_PCA, 
               dims=1:6, #维度数是人为设置的，最大为50，这里为了避免重复下面没有弄相同范围的维度数了
               reduction="pca")
dev.off()
#DimHeatmap()
pdf("PCA_DimHeatmap.pdf", width=16, height=20)
DimHeatmap(b_PCA, dims = 7:20, cells = 500, balanced = TRUE)
dev.off()
##3.肘部图,本图是最重要的降维标准
##确定维度数时，可以选择拐点数，但一般不应小于6
pdf("ElbowPlot.pdf", width=8, height=10)
ElbowPlot(b_PCA, ndims=50) #感觉其实在10左右开始变平滑
dev.off()
##其实seuratV5官网还说到了一种监督下的方法
##储存生成的文件
save(b_PCA, file="PCA.RData")






# Harmony去除批次效应 -----------------------------------------------------------
#options(future.globals.maxSize = 8 * 1024^3) #设置运行内存
#b_PCA2 <- RunHarmony(b_PCA, group.by.vars="orig.ident", assay.use="SCT", max.iter.harmony = 20)
#table(b_PCA@meta.datasorig.ident)
#因为这里不确定批次效应的去处是否符合差异分析的需求，所以我没有去除；
#其实从UMAP的结果来看，批次差异并不明显





# UMAP降维聚类分析（细降） ----------------------------------------------------------
#####UMAP降维聚类分析（细降）
##清空环境，导入数据
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="PCA.RData")
options(future.globals.maxSize = 8 * 1024^3) #设置运行内存
##降维聚类
b_cluster <- FindNeighbors(b_PCA, dims = 1:15)
b_cluster <- FindClusters(b_cluster, resolution = 0.5)
#如果用了上面的去处批次效应的数据的话，需要添加参数
#b_cluster2 <- FindNeighbors(b_PCA2, reduction="harmony", dims = 1:15) %>% FindClusters(resolution = 0.5)
##UMAP降维并绘图
b_UMAP <- RunUMAP(b_cluster, dims = 1:15)
pdf("UMAPPlot.pdf", width=8, height=10)
DimPlot(b_UMAP, reduction = "umap", label=T) #调色参数cols=a，其中a是预设好的颜色vector
DimPlot(b_UMAP, reduction = "umap", label=F, group.by="orig.ident")#批次
FeaturePlot(b_UMAP, 
            features = c("Ptprc", "Sdc1", "Ighd", "Fas" #细胞marker，MBC
                         ))
FeaturePlot(b_UMAP, 
            features = c("Il9r","Cyb561a3","Cfp","Ptpn22","Mpeg1","Zfp385a")) #31个差异基因，第1-6个
FeaturePlot(b_UMAP, 
            features = c("Sdc4","Zmynd11","St8sia6","Tspan32","Tlr7","Naip5")) #31个差异基因，第7-12个
FeaturePlot(b_UMAP, 
            features = c("Pdia4","Ache",#"AhnaK",
                         "Cdkn1b","Smad3","9930111J21Rik1","Gpr174")) #31个差异基因，第13-19个，其中第15个基因"AhnaK"消失了
FeaturePlot(b_UMAP, 
            features = c("Cdh17","Fcgrt","Ptpre","Aff1","Lgals3bp","BC147527")) #31个差异基因，第20-25个
FeaturePlot(b_UMAP, 
            features = c("Serpina3f","Pglyrp1","Rnasel","Tnfrsf18","Ifit3","Pira2")) #31个差异基因，第26-31个
dev.off()
##储存生成的文件
saveRDS(b_UMAP, file = "UMAP.rds") #注意这里变格式了，后续读取格式也应改变





# 自动注释（粗定义群） --------------------------------------------------------------
#####自动注释，此步可跳过，但如果觉得样本不够纯净可以做做去除
##个人感觉主要是看样本纯不纯净，除非我们能找到亚群高度重合的论文数据集作为参考，否则意义不大
ref1_mus_immune <- celldex::ImmGenData() #传入参考基因集
b_SingleR <- GetAssayData(object=b_UMAP@assays[["SCT"]], layer="counts") #获取标准化矩阵层
b_SingleR.ref1 <- SingleR(test = b_SingleR, ref = ref1_mus_immune, 
                      labels = ref1_mus_immune$label.main) #这里ref和Labels可以用多个数据集综合注释，将传入参数变为list()即可
table(b_SingleR.ref1$labels) #查看细胞集情况
b_UMAP@meta.data$labels <- b_SingleR.ref1$labels #回传labels到UMAP中（不过懒得单独存为一个文件了，毕竟意义不大）
save(b_UMAP, file="b_UMAP_SingleR.RData")
#绘图
pdf("UMAP_SingleR.auto.annotation.pdf", width=8, height=10)
DimPlot(b_UMAP, group.by = "labels",reduction = "umap") #生成UMAP的自动注释绘图
dev.off()





# 双细胞筛除 -------------------------------------------------------------------
rm(list=ls())
Sys.setenv(R_MAX_NUM_DLLS=999)
load(file="b_UMAP_SingleR.RData")
options(future.globals.maxSize = 8 * 1024^3) #设置运行内存
#首先获得最佳的K值
#pK表示领域大小
sweep.res.list <- paramSweep(b_UMAP, PCs =1:15, sct=T)
sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
bcmvn <- find.pK(sweep.stats)
pk_best = bcmvn %>%
  dplyr::arrange(desc(BCmetric))%>%
  dplyr::pull(pK)%>%
  .[1]%>% as.character()%>% as.numeric()
#然后估算出双细胞群中，homotypic doublets的比例
annotations <- b_UMAP$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
print(homotypic.prop)
#双细胞占比为7%左右
nExp_poi <- round( 0.07 * nrow( b_UMAP@meta.data) )
nExp_poi.adj <- round( nExp_poi * (1-homotypic.prop) )
# 模拟出artificial doublet数量。不同取值对识别结果影响不大，默认为0.25
b_UMAP <- doubletFinder(b_UMAP, PCs=1:15,
                        pN=0.25, pK=pk_best, nExp=nExp_poi.adj, sct=T) #把默认值reuse.pANN=F加上反而会报错，不知道是为什么
colnames(b_UMAP@meta.data)
#将列名改为"Double_score"和"Is_Double"
colnames(b_UMAP@meta.data)
colnames(b_UMAP@meta.data)[length(colnames(b_UMAP@meta.data))-1] <- "Double_score"
colnames(b_UMAP@meta.data)[length(colnames(b_UMAP@meta.data))] <- "Is_Double"
#查看Doubletrinder分析结果
head(b_UMAP@meta.data[,c("Double_score","Is_Double")])
table(b_UMAP@meta.data["Is_Double"])
#绘制DoubletFinder分类的UMAP图并保存
pdf("UMAP_DoubletFinder.pdf", width=8, height=10)
DimPlot(b_UMAP, reduction="umap", group.by="Is_Double")
dev.off()





# 针对细胞周期进行了解 --------------------------------------------------------------
#获取G2M期相关基因
g2m_genes <- cc.genes$g2m.genes
g2m_genes <-CaseMatch(search=g2m_genes, match=rownames(b_UMAP))
#获取S期相关基因基因
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search=s_genes,match=rownames(b_UMAP)) #主要目的是匹配
#细胞周期阶段评分
b_UMAP <- CellCycleScoring(b_UMAP, g2m.features=g2m_genes, s.features=s_genes)
colnames(b_UMAP@meta.data)
table(b_UMAP$Phase)
#画图
pdf("UMAP_CellCycle.pdf", width=8, height=10)
DimPlot(b_UMAP, reduction= "umap", group.by="Phase")
dev.off()





# 手动注释（细定义群） --------------------------------------------------------------
##有监督注释有两种方法：看文献 + 查数据库，下面是两个推荐较多的数据库
# CellMarker, 网址是http://xteam.xbio.top/CellMarker/index.jsp
# PanglaoDB, 网址是https://panglaodb.se/index.html#google_vignette
##注释的时候要注意，很多基因是有别名的，比如什么CD45别名Ptprc之类的（打PTPRC都搜不出来）
# 推荐使用genecard来确认注释的基因的名称
# genecard，网址是https://www.genecards.org/

##设置环境
options(future.globals.maxSize = 8 * 1024^3) #设置运行内存

##有监督注释
#参考文献《Single-cell BCR and transcriptome analysis after influenza infection reveals spatiotemporal dynamics of antigen-specific B cells》f1d
b_markers_ref1 <- c("Cd1d1","Cd9", #边缘区（MZ）B 细胞【C6簇】——文献明确提及的经典标准
                   "Cr2","Plac8","Lars2", #边缘区（MZ）B 细胞【C6簇】——文献高表达的其他基因
                   "Ccr7","Ebf1","Cd74","Nr4a1", #Naive细胞【C1.2.4簇】——文献明确提及的经典标准
                   "Btg1","Cd79a","lghm","Fcmr","lghd","Sell","Fcer2a","Junb","Cd69",#Naive细胞【C1.2.4簇】——文献高表达的其他基因（3个簇，所以略有差异）
                   "Cd83", #似乎是一种比较交界状态的表达？
                   "Npm1","Eif5a","Ran","Mif","Eif4a1", #preGC【C3簇】——文献高表达的基因
                   "Ppia","Top1","Hmgb1","Hmgb2", #EarlyGC【C14簇】——包括上面的preGC的高表达基因，文献还有这些高表达的基因
                   "Hmces","Ptma", #GCDZ【G2M期C15簇、S期C10簇】——文献高表达的基因
                   "Tuba1b", "Top2a","Tubb5", #GCDZ【G2M期C16簇、S期C8簇】，文献除上一行外，还高表达的基因
                   "Aicda","Tnfrsf13c", "Actb", #GCLZ【S期C13簇、不明期C12簇】——文献高表达的基因
                   "Bcl6","Bach2","Zbtb20", "Mki67", "Foxo1","Ptprc","Apoe",#GCLZ/G2M【C11簇】和preMem【C9簇】——文献高表达的基因
                   "Vim", "Cd38", #Bmem【C5簇】——文献高表达的基因
                   "Slpi","Sdc1","Prdm1", #浆细胞【C7簇】——文献提及的经典标准
                   "Xbp1","Jchain" #浆细胞——文献高表达的其他基因
                   )
b_markers_ref1 <- c("Cr2","Cd1d1","Cd9","Lars2","Ccr7",
                    "Ebf1","Btg1","Cd79a","lghm","Fcmr",
                    "lghd","Cd74","Sell","Fcer2a","Junb",
                    "Nr4a1","Cd69","Cd83","Npm1","Eif5a",
                    "Ran","Mif","Eif4a1","Hmces","Ppia",
                    "Top1","Hmgb1","Hmgb2","Tuba1b", "Top2a",
                    "Ptma","Tubb5","Aicda","Tnfrsf13c", "Actb",
                    "Bcl6","Bach2","Zbtb20", "Mki67", "Foxo1",
                    "Ptprc","Apoe","Plac8","Vim", "Cd38",
                    "Slpi","Sdc1","Prdm1","Xbp1","Jchain") #按文献原有顺序重新排列
dotplot_b_UMAP_ref1 <- DotPlot(b_UMAP, features=b_markers_ref1, cols = c("lightblue","orange"),col.min=-1
                               ) + RotatedAxis()
ggsave("dotplot_b_UMAP_ref1.pdf", dotplot_b_UMAP_ref1, width=14, height =6)
vlnplot_b_UMAP_ref1<- VlnPlot(b_UMAP, features=b_markers_ref1, stack =T, flip =T) #+ NoLegend()
ggsave("vlnplot_b_UMAP_ref1.pdf", vlnplot_b_UMAP_ref1, width=14, height=6)
#定义细胞名
#b_UMAP$celltype.main <- recode(b_UMAP@meta.data[["seurat_clusters"]], "0"="T cells", "1"="B cells", "2"="T cells",)#后续依次类推，不同的簇可以根据需求定义为相同的细胞
#DimPlot(b_UMAP, reduction = "umap", label=T, group.by="celltype.main") #绘图，这里写celltype.main是因为还有可能要细化分型

##无监督注释
##找出各簇细胞的标志marker
#b_markers.cluster0 <- FindMarkers(b_UMAP, ident.1=0, ident.2=c(1:3)) #只看某一簇细胞（0）与其他指定簇（1到3）细胞的markers；若不指定ident.2，则为其余全部细胞
#b_markers <- FindAllMarkers(b_UMAP, only.pos=TRUE, logfc.threshold=1) #只看上调，一簇细胞与其他所有细胞簇的log2fc阈值为1
#严格
b_markers_all_strict <- FindAllMarkers(b_UMAP, only.pos=F, 
                                       logfc.threshold=1, 
                                       min.pct=0.5, #min_pct是两组细胞中至少应有一组的表达比例超过这一值
                                       min.diff.pct=0.3) #较为严格的阈值，min_diff_pct为两组表达比例差异
write.csv(b_markers_all_strict, file="SCT_markers_all_strict.csv")
#宽松
b_markers_all_lenient <- FindAllMarkers(b_UMAP, only.pos=F, logfc.threshold=0.5, min.pct=0.25, min.diff.pct=0.2) #宽松一点的阈值
b_markers_all_lenient$label <- ifelse(b_markers_all_lenient$avg_log2FC>0, "positive", "negative") #赋值
b_markers_all_lenient$label <- paste0(b_markers_all_lenient$cluster, "_", b_markers_all_lenient$label) #贴标签，方便找各簇上下调基因top
write.csv(b_markers_all_lenient, file="SCT_markers_all_lenient.csv")
#宽松p
b_markers_all_lenient <- arrange(b_markers_all_lenient, p_val_adj, by_group="label") #分label按p.adj排序
b_markers_all_lenient_ptop30 <- b_markers_all_lenient %>%
  group_by(label) %>%
  filter(row_number() <= 15 | row_number() > max(1, n()-15)) %>%
  ungroup() #按p值，提取每簇细胞上下调基因各top30
write.csv(b_markers_all_lenient_ptop30, file="SCT_markers_lenient_ptop30.csv") #保存
#宽松log2FC
b_markers_all_lenient <- arrange(b_markers_all_lenient, avg_log2FC, by_group="label")
b_markers_all_lenient_log2FCtop30 <- b_markers_all_lenient %>%
  group_by(label) %>%
  filter(row_number() <= 15 | row_number() > max(1, n()-15)) %>%
  ungroup() #按log2FC，提取每簇细胞上下调基因各top30
write.csv(b_markers_all_lenient_log2FCtop30, file="SCT_markers_lenient_log2FCtop30.csv")
##保存数据
save(b_markers_all_strict, b_markers_all_lenient, file="b_markers.RData")




# 拟时序分析 -------------------------------------------------------------------
##获取基础数据
b_UMAP_df <- GetAssayData(b_UMAP, assay="SCT", layer="counts") 
b_cell_metadata <- b_UMAP@meta.data
b_gene_annotation <- data.frame(gene_short_name = rownames(b_UMAP_df), row.names = rownames(b_UMAP_df)) #基因名作为细胞注释信息
##创建一个cell_data_set对象
b_cds <- new_cell_data_set(b_UMAP_df, cell_metadata=b_cell_metadata, gene_metadata = b_gene_annotation) #预处理
b_cds <- preprocess_cds(b_cds,method ="PCA") #cell_data_set对象生成
b_cds <- reduce_dimension(b_cds, reduction_method="UMAP", preprocess_method="PCA") #降维
##绘制默认图像
o1 <- plot_cells(b_cds, reduction_method="UMAP", color_cells_by="seurat_clusters",show_trajectory_graph=F) + 
  ggtitle("cds.umap") #绘制默认图像
print(o1)
##将seurat对象的UMAP导入
b_int.embed <- Embeddings(b_UMAP, reduction ="umap")
#排序
b_int.embed <- b_int.embed[rownames(b_cds@int_colData$reducedDims$UMAP),]
#导入
b_cds@int_colData$reducedDims$UMAP <- b_int.embed
#画图
o2 <- plot_cells(b_cds, reduction_method="UMAP", color_cells_by="seurat_clusters", show_trajectory_graph=F) + 
  ggtitle("monole3_umap") #绘制默认图像
print(o2)
p = o1 | o2
print(p)
ggsave("2.Reduction_compare.pdf",plot =p, width = 10, height = 5)
#聚类分区，不同分区购细胞会进行单独的轨迹分析
b_cds <- cluster_cells(b_cds)#画图
p1 <- plot_cells(b_cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + ggtitle("partitic")
ggsave("3.cluster_Partition.pdf", plot = p1, width =6, height = 5)
#构建细胞轨迹
b_cds <- learn_graph(b_cds, learn_graph_control = list(euclidean_distance_ratio = 0.8))
#画图
p3 = plot_cells(b_cds, color_cells_by = "partition", label_groups_by_cluster = FALSE, 
               label_leaves = FALSE, label_branch_points=FALSE)
#发育轨迹(拟时序)排列细胞。
b_cds <- order_cells(b_cds)
p=plot_cells(b_cds, color_cells_by="pseudotime", label_cell_groups = FALSE,
             label_leaves = FALSE, label_branch_points = FALSE)
ggsave("5.Trajectory_Pseudotime.pdf",plot = p, width =8, height = 6)#保存结果
saveRDS(b_cds,file="6.cds.rds")
#寻找拟时轨迹差异基因
#根据低维嵌入和主图测试基因的差异表达
#Monoc1e3引入了一种寻找此类基因的新方法，它借鉴了空问相关性分析中的一项强大技术--莫兰工检验莫兰工是一种多向、多维空间自相关的测量方法。
#该统计量可以告诉您，轨迹上邻近位置的细胞对被测基因的表达水平是否相似(或不同)。
#虽然皮尔逊相关性和莫兰工的范围都在-1到1之间，但对莫兰工的解释略有不同:4+1表示邻近细胞的表达完全相似:0表示没有相关性:-1表示邻近细胞没有相关性。
Track_genes <-graph_test(b_cds, neighbor_graph="principal_graph", cores=1)
#导出
write.csv(Track_genes,"7.Track_genes_a11.csv",row.names = F)
#保存
saveRDS(Track_genes,"8.Track_genes_a11.rds")
#拟时基因热图
genes <-row.names(subset(Track_genes, morans_I>0.25))
#选取合适的clusters
num_clusters=3
#画图
pt.matrix <-exprs(b_cds)[match(genes,rownames(rowData(b_cds))),order(pseudotime(b_cds))]
pt.matrix<-t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix<-t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix)<-genes
mycol=rev(RColorBrewer::brewer.pal(11,"Spectral"))
p=ComplexHeatmap::Heatmap(
  pt.matrix,name ="z-score",show_row_names =T, show_column_names = F,
  col=circlize::colorRamp2(seq(from=-2,to=2,length=11),mycol),
  row_names_gp = gpar(fontsize =6),row_title_rot= 0, km = num_clusters,
  cluster_rows=TRUE,cluster_row_slices =FALSE, cluster_columns=FALSE, use_raster=F)
pdf("9.PseudotimeGenes heatmap.pdf",width=8,height =8)
plot(p)
dev.off()
#基因展示，选取morans_I排名前10的基因，可自定义
genes =rownames(top_n(Track_genes,n=12,morans_I))#画图
p_new1 <- plot_genes_in_pseudotime(b_cds[genes,], color_cells_by="seurat_clusters", min_expr=0.5, ncol=3)
ggsave("10.Genes_Jitterplot.pdf",plot=p_new1,width =8, height = 6)
p_new2 <- plot_cells(b_cds, genes=genes, show_trajectory_graph=FALSE,
                     labe1_cell_groups=FALSE,labe1_leaves=FALSE)
p_new2$facet$params$ncol <- 3
ggsave("11.Genes_Featureplot.pdf", plot=p_new2,width = 12, height =9)
#寻找共表达基因模块
genelist <-row.names(subset(Track_genes,morans_I >0.1))
gene_module <- find_gene_modules(b_cds[genelist,],resolution=1e-2, cores = 6)
table(gene_module$module)
write.csv(gene_module,"12.PseudotimeGenes_Module.csv", row.names = F)
#热图
cell_group <-tibble(cell=row.names(colData(b_cds)),cell_group=colData(b_cds)$seurat_clusters)
agg_mat <-aggregate_gene_expression(b_cds, gene_module, cell_group)
row.names(agg_mat)<-str_c("Module",row.names(agg_mat))
p_pheatmap <- pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2")
ggsave("13.PseudotimeGenes_Module.pdf", plot = p_pheatmap, width =8, height = 6.5)
#散点图
ME.list <- lapply(1:length(unique(gene_module$module)), function(i){subset(gene_module, module==i)$id})
names(ME.list)<- paste0("module",1:length(unique(gene_module$module)))
b_UMAP2 <-AddModuleScore(b_UMAP2,features =ME.list, name = "module")
p_dot <- FeaturePlot(b_UMAP2, features = paste0("module", 1:length(unique(gene_module$module))), ncol=2)
p_dot <- p_dot+ plot_layout()&scale_color_viridis_c(option ="C")
ggsave("14.PseudotimeGenes_Moduescorel.pdf", plot = p_dot, width = 10, height = 80, limitsize = FALSE)
