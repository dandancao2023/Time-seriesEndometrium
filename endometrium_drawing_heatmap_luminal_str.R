rm(list=ls())

library('ComplexHeatmap')
library(gridExtra)
library(circlize)
library(tidyverse)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(patchwork)
library(purrr)
library(data.table)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(simplifyEnrichment)
library(locfit)
library(dendextend)

###############
# COLORS
if (T){
  group_col <- c(
    "#980043", #LH3
    "#dd1c77", #LH5
    "#df65b0", #LH7
    "#d7b5d8", #LH9
    "#fcc5c0", #LH11
    "#fbb4b9" #RIF
  )
  names(group_col) <- c("LH3","LH5","LH7","LH9","LH11","RIF")
  
  major_cell_type_col <- setNames(c(
    "#C06CAB", #Stromal
    "#016450", #Unciliate
    "#41ab5d", #Ciliated
    "#FDE400", #Endothelial
    "#899FD1", #NK/T
    "#be7669", #myeloid
    "#F47C2B", #B
    "#53bdda" #mast
  ),c("Stromal","Unciliated","Ciliated","Endothelial","NK/T","Myeloid","B","Mast"))
  epi_subtype_col <- c(
    "#005a32", #Glandular
    "#41b6c4", #Luminal
    "#41ab5d", #Ciliated
    "#225ea8", #Secretory
    "#4a1486", #EPCAM1+SOX9+
    "#a8ddb5", #LGR5+
    "#3690c0", #EMT
    "#807dba" #Proliferative
  )
  names(epi_subtype_col) <- c(
    "Glandular","Luminal","Ciliated","Secretory","EPCAM1+SOX9+","LGR5+","EMT","Proliferative"
  )
  str_subtype_col <- c(
    "#7a0177", #Str1
    "#ae017e", #Str2
    "#dd3497", #Str3
    "#ffaed7", #MET
    "#fcc5c0" #Proliferative
  )
  names(str_subtype_col) <- c(
    "Str1","Str2","Str3","MET","Proliferative"
  )
  
  NK_T_subtype_col <- c(
    "#0f5e9c", #NK1
    "#2389da", #NK2
    "#1ca3ec", #NK3
    "#5abcd8", #Cycling_NK
    "#090088", #NKT_like
    "#e54ed0", #Cycling_T
    "#9f45b0", #CD4_T
    "#ffe4f2", #CD4_Treg
    "#c79dd7", #CD8_T
    "#ef4f91", #CD8_CTL
    "#673888" #CD3-CD56+CD16+ 
  )
  names(NK_T_subtype_col) <- c(
    "NK1","NK2","NK3","Cycling_NK","NKT_like","Cycling_T","CD4_T","CD4_Treg","CD8_T","CD8_CTL","CD3-CD56+CD16+"
  )
  
  mye_subtype_col <- c(
    "#854442", #uM1
    "#a70000", #uM2
    "#ff7f00", #cDC2
    "#ff5252", #Neutrophil
    "#ffb515", #cDC_LAMP3
    "#ffe793", #cDC_CD24
    "#e08214", #cDC_cycling
    "#fee0b6", #pDC
    "#bf812d", #moDC
    "#ffbaba" #cDC1
  )
  names(mye_subtype_col) <- c(
    "uM1","uM2","cDC2","Neutrophil","cDC_LAMP3","cDC_CD24","cDC_cycling","pDC","moDC","cDC1"
  )
  cellcols <- c(major_cell_type_col,epi_subtype_col,str_subtype_col,NK_T_subtype_col,mye_subtype_col)
  GeneColor <- setNames(c("lightcoral","lightgoldenrod","black"),c("Transcription factor","Ligand-receptor gene","other"))
}
###################
# Methods
cluster_between_groups1 <- function(mat, factor) {
  
  if (!is.factor(factor)) {
    factor = factor(factor, levels = unique(factor))
  }
  
  dend_list = list()
  order_list = list()
  for(le in unique(levels(factor))) {
    m = mat[, factor == le, drop = FALSE]
    if (ncol(m) == 1) {
      order_list[[le]] = which(factor == le)
      dend_list[[le]] = structure(which(factor == le), class = "dendrogram", leaf = TRUE,
                                  height = 0, label = 1, members = 1)
    } else if(ncol(m) > 1) {
      hc1 = hclust(dist(match(colnames(m),names(factor))))
      dend_list[[le]] = reorder(as.dendrogram(hc1), wts = match(colnames(m),names(factor)), agglo.FUN = mean)
      order_list[[le]] = which(factor == le)[order.dendrogram(dend_list[[le]])]
      order.dendrogram(dend_list[[le]]) = order_list[[le]]
    }
    attr(dend_list[[le]], ".class_label") = le
  }
  
  parent = as.dendrogram(hclust(dist(t(sapply(order_list, function(x) rowMeans(mat[, x, drop = FALSE]))))))
  
  testparent <- reorder(parent,10^1:length(order_list),agglo.FUN = mean)
  
  dend_list = lapply(dend_list, function(dend) dendrapply(dend, function(node) {
    attr(node, "height") = 0
    node
  }))
  dend = merge_dendrogram(testparent, dend_list)
  order.dendrogram(dend) = unlist(order_list[order.dendrogram(testparent)])
  return(dend)
}


###################
# Heatmap
inputDir <- file.path(getwd(),"resources")
tfsDF <- read.csv(file.path(inputDir,"tfs_ABCDE220720.csv"))
lireDF <- read.csv(file.path(inputDir,"human_ligand_receptor_genelist.csv"))

celltypeFiles <- dir(file.path(inputDir,"TAgenes_subtype"))
celltypeFiles <- celltypeFiles[grepl("_mi.csv$",celltypeFiles)]
celltypeDFList <- lapply(celltypeFiles,function(fi) read.csv(file.path(inputDir,"TAgenes_subtype",fi))%>%mutate(celltype=gsub("_\\d+_mi.csv","",fi))%>%
                           arrange(desc(x)))
names(celltypeDFList) <- gsub("_\\d+_mi.csv","",celltypeFiles)
top50celltypeDF <- lapply(celltypeDFList,function(x) {x[1:50,] %>%`colnames<-`(c("Gene","value","celltype")) %>%
    mutate(tfs=ifelse(Gene %in% tfsDF$source,"Transcription factor",NA),lire=ifelse(Gene %in% lireDF$x,"Ligand-receptor gene",NA))%>%
    rowwise()%>%mutate(GeneInfo=paste(na.omit(c(tfs,lire)),collapse = ","))%>%
    mutate(GeneInfo=ifelse(GeneInfo=="","other",GeneInfo)) %>%mutate(GeneInfo=factor(GeneInfo,levels=names(GeneColor)))})

metaDir <- file.path(inputDir,"metadata")
CellClass <- setNames(c(rep("nkt",3),rep("epi",5),"str","nk",rep("nkt",4),"epi",rep("str",4),"t",rep("mye",2)),names(celltypeDFList))
oriDataRDS <- readRDS(file.path(inputDir,"major.harmony.all.rds"))
timeDFList <- lapply(dir(metaDir)[grepl("csv$",dir(metaDir))],function(fi) read.csv(file.path(metaDir,fi)))
names(timeDFList) <- gsub("json_","",gsub("_hvg.csv","",dir(metaDir)[grepl("csv$",dir(metaDir))]))
seudotimeCutoff <- c(-0.75,-0.25,0.25,0.75)

GeneCount <- 300

needplotsubcells <- c("Luminal","str") #"Glandular","Secretory",
highlightLRgenes <- lapply(needplotsubcells,function(x) top50celltypeDF[[x]]%>%pull(Gene,GeneInfo))
names(highlightLRgenes) <- needplotsubcells

pickGenes <- read.table(file.path(inputDir,"major_pickgenes.txt"),sep='\t',header = T)
pickGenes <- apply(pickGenes,2,function(x) x[x!=""])
#names(pickGenes) <- str_split(names(pickGenes),"_",simplify = T)[,1]
names(pickGenes) <- c("str","Luminal","Glandular","Secretory","epi")
needcluster <- T

plotGroupOrderHTfull <- lapply(needplotsubcells,function(subcell) {
  celltype <- CellClass[subcell]
  print(paste0("=========================running subcell: ",subcell,"============================"))
  if (!subcell%in% c("str","epi")){
    timeDF <- timeDFList[[celltype]]%>%mutate(group=factor(group,levels=c("LH3","LH5","LH7","LH9","LH11")),
                                              pseudo_group=factor(pseudo_group,levels=c("LH3","LH5","LH7","LH9","LH11")))%>%
      dplyr::filter(cell_type==subcell)
  }else{
    timeDF <- timeDFList[[celltype]]%>%mutate(group=factor(group,levels=c("LH3","LH5","LH7","LH9","LH11")),
                                              pseudo_group=factor(pseudo_group,levels=c("LH3","LH5","LH7","LH9","LH11")))
  }
  timeDF_oriOrder <- timeDF%>%group_by(group) %>%arrange(group,pseudotime) %>% as.data.frame() %>%
    `rownames<-`(.[,"cellid"]) %>%dplyr::select(group,pseudo_group)
  timeDF_psuOrder <- timeDF%>%group_by(pseudo_group) %>%arrange(pseudo_group,pseudotime) %>% as.data.frame() %>%
    `rownames<-`(.[,"cellid"]) %>%dplyr::select(pseudo_group,group)
  
  cellids_oriOrder <- rownames(timeDF_oriOrder)
  cellids_psuOrder <- rownames(timeDF_psuOrder)
  
  oriTimeColor <- group_col[names(group_col)!="RIF"]
  pseduTimeColor <- group_col[names(group_col)!="RIF"]
  if (needcluster){
    topGenes <- celltypeDFList[[subcell]] %>% dplyr::filter(x>=0.1) %>%pull(X)
  }else{
    topGenes <- pickGenes[[subcell]]
  }
  #topGenes <- celltypeDFList[[subcell]] %>% arrange(desc(x)) %>% .[1:GeneCount,"X"]
  exprDF <- as.data.frame(oriDataRDS@assays$RNA@data[topGenes,cellids_oriOrder])
  ###########
  # smoothing to denoise
  knn <- 10
  # calculate correlations between all cells
  x_t <- as.matrix(exprDF)
  cor_mat <- cor(x_t)
  # each column is ranked list of nearest neighbours of the column cell
  nhbr_mat <- apply(-cor_mat,1,rank,ties.method="random")
  idx_mat <- nhbr_mat<=knn
  avg_knn_mat <- sweep(idx_mat,2,colSums(idx_mat), "/")
  # calculate average over all kNNs
  imputed_mat <- x_t %*% avg_knn_mat
  exprDF_denoise <- imputed_mat
  ############
  
  
  #exprDF_scale <- t(scale(t(exprDF)))
  exprDF_denoise_scale <- t(scale(t(exprDF_denoise)))
  topn_order <- ifelse(dim(exprDF_denoise_scale)[2]<300,10,ceiling(dim(exprDF_denoise_scale)[2]/30) )
  
  for (cellorder in c("original","pseudo-time")){ #
    if (cellorder=="original"){
      cellids <- cellids_oriOrder
      column_split <- timeDF_oriOrder$group
      bottom_annotation <- HeatmapAnnotation(Group=anno_block(gp=gpar(fill=oriTimeColor,col="black"),
                                                              labels_gp=gpar(col="black",fontsize=11),show_name=T),
                                             gp=gpar(fontsize=11),
                                             show_annotation_name=T,show_legend=F,#annotation_height = unit(c(1.1),c("null")),
                                             annotation_name_side="left")
      column_title <- "Original collection date Order"
      
    }else{
      cellids <- cellids_psuOrder
      column_split <- timeDF_psuOrder$pseudo_group
      bottom_annotation <- HeatmapAnnotation(Group=timeDF_psuOrder$group,
                                             `Pseudo-Group`=anno_block(gp=gpar(fill=pseduTimeColor,col="black"),
                                                                       labels_gp = gpar(col="black",fontsize=11),show_name = T),
                                             col=list(Group=oriTimeColor),
                                             gp=gpar(fontsize=11),
                                             show_annotation_name=T,show_legend=F,annotation_height = unit(c(1,1.1),c("null","null")),
                                             annotation_name_side="left")
      column_title <- "Pseudo-Time Order"
    }
    
    topGenesOrder <- data.frame(gene=topGenes,max=colMeans(apply(exprDF_denoise_scale[topGenes,cellids],1,function(x) match(sort(x,decreasing = T)[1:topn_order],x)))) %>%
      arrange(max) %>%pull(gene)
    exprDF_scale_cellorder_sort <- exprDF_denoise_scale[topGenesOrder,cellids]
    
    if (needcluster){
      geneCluster <- cutree(hclust(dist(exprDF_scale_cellorder_sort)),k=5)
      
      
      integrateDF <- exprDF_scale_cellorder_sort %>% as.data.frame()%>%rownames_to_column("gene")%>%
        full_join(.,rownames_to_column(as.data.frame(geneCluster),"gene")) %>%
        pivot_longer(-c(gene,geneCluster)) %>%
        full_join(.,timeDF,by=c("name"="cellid")) %>%
        mutate(name=factor(name,levels = cellids))
      
      
      ###########
      # panel function with scatter plot and GO enrichment
      PanelFun <- function(index,nm) {
        tmpintegrateDF <- subset(integrateDF,geneCluster==nm)%>%
          group_by(name,geneCluster,pseudotime,group,LH_time,cell_type,pseudo_group) %>%dplyr::summarise(value=mean(value))%>%ungroup()
        
        gline <- tmpintegrateDF%>%arrange(parse_number(as.character(group))) %>% ggplot(.,aes(x=pseudotime,y=value))+ 
          geom_vline(aes(xintercept=seudotimeCutoff[1]),linetype="dashed",col="grey")+
          geom_vline(aes(xintercept=seudotimeCutoff[2]),linetype="dashed",col="grey")+
          geom_vline(aes(xintercept=seudotimeCutoff[3]),linetype="dashed",col="grey")+
          geom_vline(aes(xintercept=seudotimeCutoff[4]),linetype="dashed",col="grey")+
          geom_point(aes(color=group),size=0.3)+
          
          geom_smooth(col="black",span=0.2,method="loess")+
          annotate("segment",x=-.3,xend=.3,
                   y=min(tmpintegrateDF$value)-(max(tmpintegrateDF$value)-min(tmpintegrateDF$value))/3,
                   yend=min(tmpintegrateDF$value)-(max(tmpintegrateDF$value)-min(tmpintegrateDF$value))/3,
                   arrow=arrow(length = unit(0.03,"npc")))+
          annotate("text",x=0,y=min(tmpintegrateDF$value)-(max(tmpintegrateDF$value)-min(tmpintegrateDF$value))/3,label="Pseudo-time",hjust=0.5,vjust=-1,size=3)+
          coord_cartesian(ylim=c(min(tmpintegrateDF$value),
                                 max(tmpintegrateDF$value)),clip = "off")+
          scale_color_manual(
            name="Group",values=group_col,labels=names(group_col)
          )+
          labs(x="\n",title=paste0("Cluster ",nm," (",length(index)," Genes)"))+
          theme_classic()+theme(axis.text=element_blank(),axis.ticks=element_blank(),axis.title.y = element_blank(),axis.line.y = element_blank(),
                                plot.title = element_text(size=10),legend.position = "none",plot.margin = margin(0,0,0,0))
        
        gGO <- ggplot() + theme_void()
        
        g <- ggarrange(gline,gGO,nrow=1,common.legend = T,legend = "none",widths=c(1,2)) + theme(plot.margin = margin(0.1,0.1,0.1,0.1,"cm"))
        gd <- grid.grabExpr(print(g))
        pushViewport(viewport())
        grid.rect()
        grid.draw(gd)
        popViewport()
      }
      
      
      rowNameAnnotation <- rowAnnotation(
        Gene=anno_mark(at=match(highlightLRgenes[[subcell]],rownames(exprDF_scale_cellorder_sort)),
                       labels = highlightLRgenes[[subcell]],
                       labels_gp = gpar(col=GeneColor[names(highlightLRgenes[[subcell]])],
                                        fontface="bold",
                                        fontsize=9),
                       side="right",extend=unit(5,"mm")) ,
        
        GO=anno_link(align_to=geneCluster,which="row",panel_fun=PanelFun,size=unit(5,"cm"),width=unit(20,"cm"),
                     gap=unit(0.4,"cm"),
                     link_gp = gpar(fill=setNames(brewer.pal(9,"Purples")[unique(geneCluster)+1],1:5)))
      )
      reorderWeight <- data.frame(gene=names(geneCluster),cluster=geneCluster,
                                  col=max.col(exprDF_scale_cellorder_sort)) %>%group_by(cluster) %>%summarise(weight=median(col))
      
      
      
      
      dd <- cluster_between_groups1(t(exprDF_scale_cellorder_sort),factor(geneCluster,levels = pull(arrange(reorderWeight,desc(weight)),cluster)))
      leftAnnotation <- rowAnnotation(geneCluster=anno_block(gp=gpar(fill=brewer.pal(9,"Purples")[unique(geneCluster)+1])))
      rowsplit <- 5
      htheight <- unit(25,"cm")
    }else{
      rowNameAnnotation <- rowAnnotation(
        Gene=anno_mark(at=match(rownames(exprDF_scale_cellorder_sort),rownames(exprDF_scale_cellorder_sort)),
                       labels = rownames(exprDF_scale_cellorder_sort),
                       labels_gp = gpar(col=GeneColor[unlist(lapply(rownames(exprDF_scale_cellorder_sort),
                                                                    function(x) ifelse(x%in%tfsDF$X,"Transcription factor",
                                                                                       ifelse(x%in%lireDF$X,"Ligand-receptor gene","other"))))],
                                        fontface="bold",
                                        fontsize=9),
                       side="right",extend=unit(0,"mm"),link_width = unit(0, "mm"),link_gp = gpar(col=NA)) 
      )
      dd <- F
      leftAnnotation<-rowAnnotation(em=anno_empty(border=F,width = unit(1,"mm")))
      rowsplit <- NULL
      htheight <- unit(12,"cm")
    }
    rawQuantile <- c(0.05,0.95)
    
    setQuantile <- c(6/10,7/10,8/10,9/10)
    
    heatmap_col <- setNames(c(quantile(unlist(exprDF_scale_cellorder_sort),c(rawQuantile[1],setQuantile,rawQuantile[2]))),
                            c("navyblue","lightblue","lightgreen","yellow","coral","red"))
    
    col_fun <- colorRamp2(heatmap_col,names(heatmap_col) )
    
    
    
    ht <- Heatmap(exprDF_scale_cellorder_sort,
                  col = col_fun,
                  cluster_rows = dd,  #reorder(cluster_between_groups(t(exprDF_scale),geneCluster),reorderWeight$weight),
                  row_dend_reorder = F,
                  row_split = rowsplit,
                  show_heatmap_legend = F,
                  show_row_names = F,
                  left_annotation = leftAnnotation,
                  right_annotation = rowNameAnnotation,
                  row_dend_side = "left",
                  show_column_names = F,
                  column_title = column_title,
                  row_names_side = "right",
                  column_split = column_split,
                  column_gap=unit(0,"mm"),
                  bottom_annotation = bottom_annotation,
                  column_order = cellids,
                  row_title =NULL,
                  width=unit(20,"cm"),
                  height = htheight
    )
    
    mainLegend <- Legend(c(heatmap_col[1]-(heatmap_col[length(heatmap_col)]-heatmap_col[1])/10,
                           heatmap_col,
                           heatmap_col[length(heatmap_col)]+(heatmap_col[length(heatmap_col)]-heatmap_col[1])/10),
                         paste(c(paste0("<",sprintf("%.2f",heatmap_col[1])),sprintf("%.2f",heatmap_col),
                                 paste0(">",sprintf("%.2f",heatmap_col[length(heatmap_col)]))),
                               c("Quantile",sprintf("%.0f%%",c(rawQuantile[1],setQuantile,rawQuantile[2])*100),""),
                               sep="\n"),
                         title=paste0("Z-score of Normalized Expression of ",ifelse(needcluster,"MI score >0.1 ",""),"Genes (",dim(exprDF_scale_cellorder_sort)[1],")"),
                         col_fun=col_fun,
                         title_position = "topcenter",direction="horizontal",title_gp=gpar(fontface="bold"),legend_width = unit(1,"npc"))
    blocklgdDate <- Legend(title="Group based on date",at=names(oriTimeColor),legend_gp = gpar(fill=oriTimeColor),nrow=1,title_position="topcenter")
    blocklgdPseu <- Legend(title="Group based on pseudo-time",at=names(pseduTimeColor),legend_gp = gpar(fill=pseduTimeColor),nrow=1,title_position="topcenter")
    if (cellorder=="original"){
      blocklgd <- blocklgdDate
    }else{
      blocklgd <- packLegend(blocklgdDate,blocklgdPseu,direction = "horizontal",gap=unit(20,"mm"))
    }
    
    panelfunlgd1 <- Legend(title="p.adjust",labels=c("*","**","***"),graphics=list(
      function(x,y,w,h){grid.rect(x,y,w,h,gp=gpar(col="black",fill="grey",alpha=0.25))},
      function(x,y,w,h){grid.rect(x,y,w,h,gp=gpar(col="black",fill="grey",alpha=0.5))},
      function(x,y,w,h) {grid.rect(x,y,w,h,gp=gpar(col="black",fill="grey",alpha=1))}
    ), nrow=1,title_position="topcenter")
    panelfunlgd2 <- Legend(title="Ontology",at=c("BP","CC","MF"),legend_gp=gpar(fill=c("cornflowerblue","coral","lightseagreen")),nrow = 1,title_position = "topcenter")
    panelfunlgd <- packLegend(panelfunlgd1,panelfunlgd2,direction = "horizontal",gap=unit(20,"mm"))
    if(needcluster){
      genelgd <- Legend(title=ifelse(needcluster,"MI score top 50 genes","Genes"),labels=names(GeneColor),graphics = list(
        
        function(x,y,w,h) {
          grid.text("Gene",x,y,gp=gpar(col=GeneColor[1],fontface="bold"))
        },
        function(x,y,w,h) {
          grid.text("Gene",x,y,gp=gpar(col=GeneColor[2],fontface="bold"))
        },
        function(x,y,w,h) {
          grid.text("Gene",x,y,gp=gpar(col=GeneColor[3],fontface="bold"))
        }
        
      ), nrow=1,title_position="topcenter",grid_width = unit(2,"cm"))
    }else{
      genelgd <- Legend(title="",at=c(1),grid_height = unit(0,"mm"),grid_width = unit(0,"mm"),labels_gp = gpar(col=NA))
    }
    panelgenelgd <- packLegend(panelfunlgd,genelgd,direction = "horizontal",gap=unit(20,"mm"))
    
    grobs <- lapply(list(ht,mainLegend,blocklgd,genelgd),function(x) grid.grabExpr(draw(x)))
    lay <- rbind(c(1,1),
                 c(2,2),
                 c(3,3),
                 c(4,4))
    pdf(file.path(inputDir,"..","heatmap",paste0(subcell,"_",cellorder,"_complex",ifelse(needcluster,"","_pickuncluster"),".pdf")),height = 13,width=18)
    grid.arrange(grobs[[1]],grobs[[2]],grobs[[3]],grobs[[4]],
                 ncol=2,nrow=4,layout_matrix=lay,widths=c(25,25),heights=c(as.numeric(htheight),2,1,1),
                 top=textGrob(paste0(dim(timeDF)[1]," cells in group ",subcell),gp=gpar(fontsize=18,fontface="bold")))
    dev.off()
  }
})


