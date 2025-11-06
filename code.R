library(Seurat)
library(tidyverse)
library(dplyr)
library(harmony)
library(DoubletFinder)
library(stats)
library(cowplot)
library(hdf5r)
library(rhdf5)
library(patchwork)
library(data.table)
library(limma)
library(devtools)
library(plyr)
library("scales")
library(ggsci)
library(ggplot2)
library(ggpubr)

setwd("/share/appspace_data/shared_groups/BGI/USERS/bgi_liuyang/TYM/data/H5/")
samples=list.files("./")
dir <- file.path('/share/appspace_data/shared_groups/BGI/USERS/bgi_liuyang/TYM/data/H5',samples)
names(dir)=str_split(dir,'/',simplify = T)[,11]
names(dir)=str_split(names(dir),'_',simplify = T)[,2]
scRNAlist <- list()
for(i in 1:length(dir)){
  A = Read10X_h5(filename = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(A, min.cells = 3, project = names(dir)[i],min.features =300)
  ##计算细胞中线粒体基因比例
  scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")
  colnames(scRNAlist[[i]]) <- paste(names(dir)[i], colnames(scRNAlist[[i]]), sep = "_")
}
scRNA <- Reduce(function(x, y) merge(x, y), scRNAlist)

col.num <- length(levels(scRNA@active.ident))
violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                  cols =rainbow(col.num), 
                  pt.size = 0.01, 
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

violin
scRNA <- subset(scRNA, subset = nFeature_RNA > 300 & nFeature_RNA < 7000 & percent.mt < 25 & nCount_RNA < 50000)

# SCTransform标准化
scRNA <- SCTransform(scRNA,assay = "RNA")
# 主成分分析（PCA）
scRNA <- RunPCA(scRNA, npcs = 50)
ElbowPlot(scRNA, ndims = 50)
# 运行Harmony进行批次校正
scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident", max.iter.harmony = 20)

pc.num=1:30
scRNA <- FindNeighbors(scRNA, dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 0.8)
scRNA <- RunUMAP(scRNA, dims = pc.num)
scRNA <- RunTSNE(scRNA, dims = pc.num)

options(repr.plot.width = 13, repr.plot.height = 6)
DimPlot(scRNA, reduction = "umap", group.by = c("orig.ident","seurat_clusters"),raster=FALSE)
DimPlot(scRNA, reduction = "tsne", group.by = c("orig.ident","seurat_clusters"),raster=FALSE)

DotPlot(scRNA, features = select_genes )+ coord_flip()+
  scale_y_discrete(limits = c(
    "18","22",
    "2","20","21",
    "8",
    "1","5","9","10","12","24",
    "0","3","4","6","7","13","14","15",
    "16","17","26","23","25",
    "19",
    "11"
  ))+
  scale_colour_gradientn(colours = c('#482878','#26828e','#6dcd59','#b4de2c','#fde725'))

celltype = read.table("/share/appspace_data/shared_groups/BGI/USERS/bgi_liuyang/TYM/plot/celltype.txt", header=T, sep="\t", check.names=F)
clusters <- scRNA@meta.data$seurat_clusters
scRNA@meta.data$CellType=celltype[match(clusters,celltype$ClusterID),'celltype']

p1 = DimPlot(scRNA, group.by="seurat_clusters",  label.size=5, reduction='umap',cols = allcolour,label = T)
p2 = DimPlot(scRNA, group.by="CellType",  label.size=5, reduction='umap',cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69"),label = T)
p = p1 +p2
p

A = scRNA@meta.data
a <- data.frame(table(A$seurat_clusters,A$orig.ident))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
#a$Var1 = factor(a$Var1,levels = Type$id)

a %>% 
  drop_na() %>%  #去掉所有空值，避免出错
  ggplot(aes(fill=Var2, y= percent, x = Var1)) + 
  geom_bar(position="fill", stat = "identity") +
  scale_fill_manual(values=allcolour)+ #设置填充的颜色 brewer.pal(9,"Set3")
  scale_y_continuous(labels = scales::percent) + #纵坐标变为百分比
  ylab("Percent(%)")+xlab("")+labs(fill="")+
  theme_gray()+
  rotate_x_text(90)+
  #coord_flip()+
  theme(legend.position = "right",
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank()
  ) -> p1

p1

celltype = read.table("/share/appspace_data/shared_groups/BGI/USERS/bgi_liuyang/TYM/plot/Type.txt", header=T, sep="\t", check.names=F)
clusters <- scRNA@meta.data$orig.ident
scRNA@meta.data$Group=celltype[match(clusters,celltype$ClusterID),'celltype']

a <- data.frame(table(A$seurat_clusters,A$Group))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
#a$Var1 = factor(a$Var1,levels = Type$id)

a %>% 
  drop_na() %>%  #去掉所有空值，避免出错
  ggplot(aes(fill=Var2, y= percent, x = Var1)) + 
  geom_bar(position="fill", stat = "identity") +
  scale_fill_manual(values=allcolour)+ #设置填充的颜色 brewer.pal(9,"Set3")
  scale_y_continuous(labels = scales::percent) + #纵坐标变为百分比
  ylab("Percent(%)")+xlab("")+labs(fill="")+
  theme_gray()+
  rotate_x_text(90)+
  #coord_flip()+
  theme(legend.position = "right",
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank()
  ) -> p1

p1

p=DimPlot(scRNA,group.by = "CellType",split.by = 'Group',cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69"))
p

a <- data.frame(table(A$orig.ident,A$CellType))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
#a$Var1 = factor(a$Var1,levels = Type$id)

a %>% 
  drop_na() %>%  #去掉所有空值，避免出错
  ggplot(aes(fill=Var2, y= percent, x = Var1)) + 
  geom_bar(position="fill", stat = "identity") +
  scale_fill_manual(values=allcolour)+ #设置填充的颜色 brewer.pal(9,"Set3")
  scale_y_continuous(labels = scales::percent) + #纵坐标变为百分比
  ylab("Percent(%)")+xlab("")+labs(fill="")+
  theme_gray()+
  rotate_x_text(90)+
  #coord_flip()+
  theme(legend.position = "right",
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank()
  ) -> p1

p1

a <- data.frame(table(A$Group,A$CellType))
a<- ddply(a,.(Var1),transform,percent=Freq/sum(Freq)*100) 
a$label = paste0(sprintf("%.1f", a$percent), "%")
#a$Var1 = factor(a$Var1,levels = Type$id)

a %>% 
  drop_na() %>%  #去掉所有空值，避免出错
  ggplot(aes(fill=Var2, y= percent, x = Var1)) + 
  geom_bar(position="fill", stat = "identity") +
  scale_fill_manual(values=allcolour)+ #设置填充的颜色 brewer.pal(9,"Set3")
  scale_y_continuous(labels = scales::percent) + #纵坐标变为百分比
  ylab("Percent(%)")+xlab("")+labs(fill="")+
  theme_gray()+
  rotate_x_text(90)+
  #coord_flip()+
  theme(legend.position = "right",
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),#去除背景
        panel.border = element_blank()
  ) -> p1

p1

allcolour <- c("#8DD3C7","#FFFFB3","#FB8072","#BEBADA","#80B1D3","#FDB462",
               "#05558b","#FCCDE5","#BC80BD","#D9D9D9","#FFFF00","#98FB98",
               "#EE82EE","#40E0D0","#B3DE69","#DC143C","#0000FF","#20B2AA",
               "#FFA500","#9370DB","#98FB98",
               "#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF",
               "#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C",
               "#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD",
               "#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4",
               "#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD",
               "#9932CC","#8B008B","#8B4513","#DEB887")

select_genes <- c("KRT18","TPPP3", #上皮
                  "CD68", #巨噬细胞
                  "FCGR3B", #中性粒细胞
                  "CD1C", #DC
                  "CLEC9A", #mDc
                  "LILRA4", #pDC
                  "TPSB2", #Mast
                  "CD3D", #T cell
                  "KLRD1", #NK
                  "MS4A1", #B 细胞
                  "IGHG4" #Plasma cell
                  ) 

options(repr.plot.width = 13, repr.plot.height = 6)
p=DotPlot(scRNA, features = select_genes )+ coord_flip()+
  scale_y_discrete(limits = c(
    "11","17","21",
    "0","1","2","3","5","6","10","12","13","14","15","16","18","22","28",
    "4",
    "19",
    "27",
    "29",
    "7","23","24","8",
    "9",
    "26",
    "20","25","30"
  ))+
  scale_colour_gradientn(colours = c('#482878','#26828e','#6dcd59','#b4de2c','#fde725'))
p

p=DotPlot(scRNA, features = select_genes )+ coord_flip()+
  scale_y_discrete(limits = c(
    "11","17","21",
    "0","1","2","3","5","6","10","12","13","14","15","16","18","22","28",
    "4",
    "19",
    "27",
    "29",
    "7","23","24","8",
    "9",
    "26",
    "20","25","30"
  ))+
  scale_colour_gradientn(colours = c('#482878','#26828e','#6dcd59','#b4de2c','#fde725'))

celltype = read.table("/share/appspace_data/shared_groups/BGI/USERS/bgi_liuyang/TYM/plot/celltype.txt", header=T, sep="\t", check.names=F)
clusters <- scRNA@meta.data$seurat_clusters
scRNA@meta.data$CellType=celltype[match(clusters,celltype$ClusterID),'celltype']

p1 = DimPlot(scRNA, group.by="seurat_clusters",  label.size=5, reduction='umap',cols = allcolour,label = T)
p2 = DimPlot(scRNA, group.by="CellType",  label.size=5, reduction='umap',cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69"),label = T)
p = p1 +p2
p

# 提取 UMAP 坐标和细胞类型信息
umap_data <- FetchData(scRNA, vars = c("UMAP_1", "UMAP_2", "CellType"))

# 计算每个 CellType 的坐标中位数（或均值）
library(dplyr)
centers <- umap_data %>%
  group_by(CellType) %>%
  summarise(
    UMAP_1 = median(UMAP_1),
    UMAP_2 = median(UMAP_2)
  )

p1 = DimPlot(scRNA, group.by="seurat_clusters",  label.size=5, reduction='umap',cols = allcolour,label = T)
p2=DimPlot(scRNA, group.by="CellType",  label.size=5, reduction='umap',cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69")) + geom_text_repel(
  data = centers,                          # 使用计算的中心坐标
  aes(x = UMAP_1, y = UMAP_2, label = CellType),
  size = 5,                                # 标签字体大小
  box.padding = 0.8,                       # 标签框周围的填充空间
  point.padding = 0.3,                     # 标签与对应点之间的最小间距
  force = 20,                              # 标签间的排斥力（值越大间距越大）
  max.overlaps = 100,                      # 允许的最大重叠次数（设为 Inf 显示全部标签）
  segment.color = "grey50",                # 连接线颜色
  segment.size = 0.5                       # 连接线粗细
)
p = p1 +p2
p

celltype = read.table("/share/appspace_data/shared_groups/BGI/USERS/bgi_liuyang/TYM/plot/Type.txt", header=T, sep="\t", check.names=F)
clusters <- scRNA@meta.data$orig.ident
scRNA@meta.data$Group=celltype[match(clusters,celltype$ClusterID),'celltype']

p=DimPlot(scRNA,group.by = "CellType",split.by = 'Group',cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69"),label = T)
p

# 假设原始Seurat对象为 seurat_obj，细胞类型存储在 metadata 列 "celltype"
macrophage <- subset(scRNA, subset = CellType == "Macrophages")

# SCTransform标准化
macrophage <- SCTransform(macrophage,assay = "RNA")
# 主成分分析（PCA）
macrophage <- RunPCA(macrophage, npcs = 50)
ElbowPlot(macrophage, ndims = 50)

pc.num=1:50
macrophage <- FindNeighbors(macrophage, dims = pc.num)
macrophage <- FindClusters(macrophage, resolution = 0.5)
macrophage <- RunUMAP(macrophage, dims = pc.num)
macrophage <- RunTSNE(macrophage, dims = pc.num)

options(repr.plot.width = 13, repr.plot.height = 6)
DimPlot(macrophage, reduction = "umap", group.by = c("orig.ident","seurat_clusters"),raster=FALSE)
DimPlot(macrophage, reduction = "tsne", group.by = c("orig.ident","seurat_clusters"),raster=FALSE)

select_genes <- c("CD68","CD11B", #所有巨噬细胞都表达的标记
                  "CD80","CD86","CD64","CD16","CD32","NOS2", #M1
                  "CD163","CD206", #M2
                  "S100A8","FCN1","CD14", #单核来源巨噬细胞
                  "FABP4","APOC1","MARC" #肺泡常驻巨噬细胞
                  ) 

p=DimPlot(macrophage, reduction = "umap", group.by = c("orig.ident","seurat_clusters"),raster=FALSE,cols = allcolour,label = T)
p

select_genes <- c(#"CD68","CD11B", #所有巨噬细胞都表达的标记
                  "SLC3A2",
                  "CD86",
                  "CCL2","CCL3","CXCL8",
                  "S100A8","FCN1","CD14", #单核来源巨噬细胞
                  "MRC1","TREM2","MSR1",
                  "FABP4","APOC1","MARCO" #肺泡常驻巨噬细胞
                  ) 

options(repr.plot.width = 13, repr.plot.height = 6)
p=DotPlot(macrophage, features = select_genes )+ coord_flip()+
  scale_y_discrete(limits = c(
    "1","3","4","10","12","16","20","21","22",
    "7","9","15",
    "0","2","5","6","8","11","13","14","17","18","19"
  ))+
  scale_colour_gradientn(colours = c('#482878','#26828e','#6dcd59','#b4de2c','#fde725'))
p

celltype = read.table("/share/appspace_data/shared_groups/BGI/USERS/bgi_liuyang/TYM/plot1/celltype.txt", header=T, sep="\t", check.names=F)
clusters <- macrophage@meta.data$seurat_clusters
macrophage@meta.data$CellType=celltype[match(clusters,celltype$ClusterID),'celltype']

options(repr.plot.width = 11, repr.plot.height = 6)
p = DimPlot(macrophage, group.by="CellType",  label.size=5, reduction='umap',cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69","#8B4513","#DEB887"),label = F,raster=FALSE)
p

macrophage$CellType = factor(macrophage$CellType, levels=c("Group1:M2 like anti-inflammatory (alveolar macrophage)","Group2:intermediate M1-M2 macrophage(alveolar macrophage)","Group3:pro-inflammatory M1 like (monocyte-derived macrophage)"))

options(repr.plot.width = 15.5, repr.plot.height = 13.5)
p=VlnPlot(macrophage, features = c("SLC3A2"),group.by = "CellType",cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69","#8B4513","#DEB887"), pt.size = 0)+
  theme(axis.text.x = element_text(angle = 75))+
geom_boxplot(alpha=0.1,outlier.size=0,size=0.3,width=0.3,fill="white")
p

# 提取SLC3A2基因的表达矩阵
slc3a2_expression <- GetAssayData(macrophage, slot = "counts", genes = "SLC3A2")
slc3a2 = as.data.frame(slc3a2_expression["SLC3A2",])
colnames(slc3a2) = "slc3a2"
meat = macrophage@meta.data
meat = meat[rownames(slc3a2),]
slc3a2$group = meat$CellType
slc3a2$class = meat$Group

options(repr.plot.width = 7.5, repr.plot.height = 18.5)
p=ggplot(slc3a2, aes(x = group, y = slc3a2, fill = group)) +
  geom_violin(trim = FALSE) +  # 绘制小提琴图
  stat_summary(
    fun = mean,
    geom = "crossbar",  # 使用 "crossbar" 或 "errorbar" 显示中位线
    width = 0.3,        # 调整宽度
    color = "black",    # 线条颜色
    size = 0.8          # 线条粗细
  ) +
theme_bw()+
  theme(panel.grid=element_blank())+
  #stat_compare_means(aes(group=`Pathologic Response`),method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  
  theme(legend.position = "right")+
  #ggtitle("") +
  ylab("SLC3A2")+
  #xlab("")+
  rotate_x_text(90)+
  ylim(0,15)+
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label.y = c(8,10,12),label = "p.signif")+
  scale_fill_manual(values= c("#8DD3C7" ,"#FFFFB3" , "#FB8072"))+
  theme(
    axis.title.x=element_blank(),
    axis.text.y = element_text(size = 15,face = "bold",color = "black"),
    axis.title.y = element_text(size = 15,face = "bold",color = "black"),
    axis.text.x = element_text(size = 15,face = "bold",color = "black"),
    legend.position = "none",
    
    plot.title = element_text(face = "bold",size=15,hjust = 0.5))
p

group=levels(factor(slc3a2$group))
slc3a2$group=factor(slc3a2$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

data <- slc3a2 %>% filter(slc3a2 <= 15)
options(repr.plot.width = 7.5, repr.plot.height = 18.5)
p=ggplot(data, aes(x = group, y = slc3a2, fill = group)) +
  geom_violin(trim = FALSE) +  # 绘制小提琴图
  stat_summary(
    fun = mean,
    geom = "crossbar",  # 使用 "crossbar" 或 "errorbar" 显示中位线
    width = 0.3,        # 调整宽度
    color = "red",    # 线条颜色
    size = 0.8          # 线条粗细
  ) +
theme_bw()+
  theme(panel.grid=element_blank())+
  #stat_compare_means(aes(group=`Pathologic Response`),method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  
  theme(legend.position = "right")+
  #ggtitle("") +
  ylab("SLC3A2")+
  #xlab("")+
  rotate_x_text(90)+
  ylim(0,20)+
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label.y = c(15,17,19),label = "p.signif")+
  scale_fill_manual(values= c("#8DD3C7" ,"#FFFFB3" , "#FB8072"))+
  theme(
    axis.title.x=element_blank(),
    axis.text.y = element_text(size = 15,face = "bold",color = "black"),
    axis.title.y = element_text(size = 15,face = "bold",color = "black"),
    axis.text.x = element_text(size = 15,face = "bold",color = "black"),
    legend.position = "none",
    
    plot.title = element_text(face = "bold",size=15,hjust = 0.5))
p

group=levels(factor(slc3a2$class))
slc3a2$class=factor(slc3a2$class, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

options(repr.plot.width = 5, repr.plot.height = 7.5)
p=ggplot(slc3a2, aes(x = class, y = slc3a2, fill = class)) +
  geom_violin(trim = FALSE) +  # 绘制小提琴图
  stat_summary(
    fun = mean,
    geom = "crossbar",  # 使用 "crossbar" 或 "errorbar" 显示中位线
    width = 0.3,        # 调整宽度
    color = "red",    # 线条颜色
    size = 0.8          # 线条粗细
  ) +
theme_bw()+
  theme(panel.grid=element_blank())+
  #stat_compare_means(aes(group=`Pathologic Response`),method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  
  theme(legend.position = "right")+
  #ggtitle("") +
  ylab("SLC3A2")+
  #xlab("")+
  rotate_x_text(90)+
  ylim(0,25)+
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label.y = c(19,21,23),label = "p.signif")+
  scale_fill_manual(values= c("#8DD3C7" ,"#FFFFB3" , "#FB8072"))+
  theme(
    axis.title.x=element_blank(),
    axis.text.y = element_text(size = 15,face = "bold",color = "black"),
    axis.title.y = element_text(size = 15,face = "bold",color = "black"),
    axis.text.x = element_text(size = 15,face = "bold",color = "black"),
    legend.position = "none",
    
    plot.title = element_text(face = "bold",size=15,hjust = 0.5))
p

p=DimPlot(macrophage,group.by = "CellType",split.by = 'Group',cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69"),label = F)
p

options(repr.plot.width = 16, repr.plot.height = 6)
p=DimPlot(macrophage,group.by = "CellType",split.by = 'Group',cols = c("#8DD3C7" ,"#FFFFB3" , "#FB8072" ,"#BEBADA","#80B1D3" ,"#FDB462" ,"#05558b",
  "#FCCDE5", "#BC80BD","#D9D9D9","#FFFF00","#98FB98", "#EE82EE", "#40E0D0","#B3DE69"),label = F)
p

data <- FetchData(macrophage, vars = c("SLC3A2", "orig.ident","Group"))
data %>%
  group_by(orig.ident) %>%
  filter(row_number() == 1)
data$Sample = data$orig.ident
data <- data %>%
  mutate(orig.ident = case_when(
    orig.ident == "C141" ~ "M1",
    orig.ident == "C142" ~ "M2",
    orig.ident == "C143" ~ "S1",
    orig.ident == "C144" ~ "M3",
    orig.ident == "C145" ~ "S2",
    orig.ident == "C146" ~ "S3",
    orig.ident == "C51" ~ "H1",
    orig.ident == "C52" ~ "H2",
    orig.ident == "C100" ~ "H3",
    orig.ident == "C148" ~ "S4",
    orig.ident == "C149" ~ "S5",
    orig.ident == "C152" ~ "S6",
    orig.ident == "SC249NORbal_fresh" ~ "H4",
    TRUE ~ orig.ident
  ))
p<-ggplot(data,aes(x=orig.ident,y=SLC3A2,fill=orig.ident))+
geom_violin(position=position_dodge(width=0.1),scale='width')+#小提琴
geom_boxplot(alpha=0.1,outlier.size=0,size=0.3,width=0.3,fill="white")+# 小提琴中的箱线图
#scale_fill_manual(values=col)+# 手动填充颜色
labs(x="",y="SLC3A2",color="Group")+
#scale_x_discrete(limits=c("BS","RS","RE","VE","SE","LE","P")) +
theme_classic()+
theme(axis.text.x=element_text(size=10,face="bold"),axis.text.y=element_text(size=10))+
theme(axis.title.y=element_text(size=12,face="bold"))+theme(axis.title.x=element_text(size=12))+
theme(legend.title=element_text(size=10),legend.text=element_text(size=8))
p

options(repr.plot.width = 14.5, repr.plot.height = 6.5)
p=ggplot(data, aes(x = orig.ident, y = SLC3A2, fill = orig.ident)) +
  
#geom_jitter(width = 0.2, height = 0.2, size = 1.5, alpha = 0.1) +  # 添加散点图  
geom_violin(trim = FALSE) +  # 绘制小提琴图
stat_summary(
    fun = mean,
    geom = "crossbar",  # 使用 "crossbar" 或 "errorbar" 显示中位线
    width = 0.3,        # 调整宽度
    color = "red",    # 线条颜色
    size = 0.8          # 线条粗细
  ) +
theme_bw()+
  theme(panel.grid=element_blank())+
  #stat_compare_means(aes(group=`Pathologic Response`),method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  
  theme(legend.position = "right")+
  #ggtitle("") +
  ylab("SLC3A2")+
  #xlab("")+
  rotate_x_text(90)+
  #ylim(0,15)+
  #stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label.y = c(8,10,12),label = "p.signif")+
  #scale_fill_manual(values= c("#8DD3C7" ,"#FFFFB3" , "#FB8072"))+
  theme(
    axis.title.x=element_blank(),
    axis.text.y = element_text(size = 15,face = "bold",color = "black"),
    axis.title.y = element_text(size = 15,face = "bold",color = "black"),
    axis.text.x = element_text(size = 15,face = "bold",color = "black"),
    legend.position = "none",
    
    plot.title = element_text(face = "bold",size=15,hjust = 0.5))
p

gene <- c('SLC3A2')
exp = macrophage[gene,]
exp = as.data.frame(exp@assays$RNA@counts)
exp = as.data.frame(t(exp))
A = macrophage@meta.data
B = cbind(A[,c(10)],exp)
colnames(B)[1] = "group"

data = melt(B,id.vars=c("group"))
colnames(data) = c("group","Gene","Expression")
group=levels(factor(data$group))
data$group=factor(data$group, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

data$group = factor(data$group, levels=c("Group1:M2 like anti-inflammatory (alveolar macrophage)","Group2:intermediate M1-M2 macrophage(alveolar macrophage)","Group3:pro-inflammatory M1 like (monocyte-derived macrophage)"))

options(repr.plot.width = 8.5, repr.plot.height = 17.5)
p=ggplot(data,aes(x = group , y=Expression, color=group))+
  stat_boxplot(geom="errorbar",width=0.5,size=1.3)+
  geom_boxplot(alpha=1,outlier.shape = NA,size=1.3,width=0.5,fatten=1)+
  theme_bw()+
  theme(panel.grid=element_blank())+
  #stat_compare_means(aes(group=`Pathologic Response`),method="wilcox.test",symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")+
  
  theme(legend.position = "right")+
  #ggtitle("") +
  ylab("SLC3A2")+
  #xlab("")+
  rotate_x_text(90)+
  ylim(0,15)+
  stat_compare_means(comparisons = my_comparisons,symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label.y = c(8,10,12),label = "p.signif")+
  scale_color_manual(values= c("#8DD3C7" ,"#FFFFB3" , "#FB8072"))+
  theme(
    axis.title.x=element_blank(),
    axis.text.y = element_text(size = 15,face = "bold",color = "black"),
    axis.title.y = element_text(size = 15,face = "bold",color = "black"),
    axis.text.x = element_text(size = 15,face = "bold",color = "black"),
    legend.position = "none",
    
    plot.title = element_text(face = "bold",size=15,hjust = 0.5))
p

options(repr.plot.width = 6.5, repr.plot.height = 6)
p=DimPlot(macrophage, group.by="seurat_clusters",  label.size=5, reduction='umap',cols = allcolour,label = T)
p

macrophage@meta.data <- macrophage@meta.data %>%
  mutate(orig.ident = case_when(
    orig.ident == "C141" ~ "M1",
    orig.ident == "C142" ~ "M2",
    orig.ident == "C143" ~ "S1",
    orig.ident == "C144" ~ "M3",
    orig.ident == "C145" ~ "S2",
    orig.ident == "C146" ~ "S3",
    orig.ident == "C51" ~ "H1",
    orig.ident == "C52" ~ "H2",
    orig.ident == "C100" ~ "H3",
    orig.ident == "C148" ~ "S4",
    orig.ident == "C149" ~ "S5",
    orig.ident == "C152" ~ "S6",
    orig.ident == "SC249NORbal_fresh" ~ "H4",
    TRUE ~ orig.ident
  ))

options(repr.plot.width = 6.5, repr.plot.height = 6)
p=DimPlot(macrophage, group.by="orig.ident",  label.size=5, reduction='umap',cols = allcolour,label = F)
p

scRNA@meta.data <- scRNA@meta.data %>%
  mutate(orig.ident = case_when(
    orig.ident == "C141" ~ "M1",
    orig.ident == "C142" ~ "M2",
    orig.ident == "C143" ~ "S1",
    orig.ident == "C144" ~ "M3",
    orig.ident == "C145" ~ "S2",
    orig.ident == "C146" ~ "S3",
    orig.ident == "C51" ~ "H1",
    orig.ident == "C52" ~ "H2",
    orig.ident == "C100" ~ "H3",
    orig.ident == "C148" ~ "S4",
    orig.ident == "C149" ~ "S5",
    orig.ident == "C152" ~ "S6",
    orig.ident == "SC249NORbal_fresh" ~ "H4",
    TRUE ~ orig.ident
  ))

options(repr.plot.width = 6.5, repr.plot.height = 6)
p=DimPlot(scRNA, group.by="orig.ident",  label.size=5, reduction='umap',cols = allcolour,label = F)
p

