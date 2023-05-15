#https://www.datanovia.com/en/lessons/wilcoxon-test-in-r/#prerequisites
#Statistical tests
library(tidyverse)
library(rstatix)
library(ggpubr)
library(effsize)  
library(ggplot2)

setwd("~/new_frontal_official/hdWGCNA/Ast-subset")
seurat_obj <- readRDS("~/new_frontal_official/hdWGCNA/Ast-subset/Ast_hdWGCNA-seurat_obj-basic-end.rds")

modules <- c("Ast-M1", "Ast-M2", "Ast-M3", "Ast-M4", "Ast-M5", "Ast-M6", "Ast-M7", "Ast-M8", "Ast-M9", "Ast-M10")
modules <- data.frame(modules)

#########Kruskal_Wallis test#########
Kruskal_Wallis<- setNames(data.frame(matrix(ncol = 3, nrow = nrow(modules))), c("Module", "Kruskal-Wallis chi-squared", "p-value"))
for (i in 1:nrow(modules)){
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  #Kruskal-Wallis test
  Kruskal_Wallis_data<-kruskal.test(hME ~ diagnosis, data = hME)
  
  Kruskal_Wallis[i,1]<-modules[i,1]
  Kruskal_Wallis[i,2]<-Kruskal_Wallis_data$statistic
  Kruskal_Wallis[i,3]<-Kruskal_Wallis_data$p.value
}

write_csv(Kruskal_Wallis, "Kruskal_Wallis_data.csv")


#########Wilcoxon Rank Sum test#########
#E280A vs Sporadic  
Wilcoxon_rank_sum<- setNames(data.frame(matrix(ncol = 5, nrow = nrow(modules))), c("Module", "Wilcoxon p-value (BH adjusted)", "Delaney's A effect size", "Delaney's A magnitude", "Wilcoxon effect size"))

for (i in 1:nrow(modules)){
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  # Show a sample of the data by group
  set.seed(223)
  
  #signifiance test
  stat.test <- hME %>% 
    rstatix::wilcox_test(hME ~ diagnosis) %>%
    add_significance()
  stat.test
  W_effect_size<-hME %>% wilcox_effsize(hME ~ diagnosis)
  
  Wilcoxon_rank_sum[i,1] <- modules[i,1]
  Wilcoxon_rank_sum[i,2] <- stat.test[1,7]
  Wilcoxon_rank_sum[i,5] <- W_effect_size[1,4]
         
  #get AUC effect size
  Sporadic<-hME$hME[hME$diagnosis=='Sporadic AD']
  E280A<-hME$hME[hME$diagnosis=='E280A']
         
  Delaneys_A<-VD.A(new_frontal_official, E280A)
  Wilcoxon_rank_sum[i,3] <-Delaneys_A[["estimate"]]
  Wilcoxon_rank_sum[i,4] <-Delaneys_A[["magnitude"]]   
}

write_csv(Wilcoxon_rank_sum, "Wilcoxon_rank_sum_data.csv")

#########Graphs for Wilcoxon#########
for (i in 1:nrow(modules)){
  
  hME<-data.frame(diagnosis=seurat_obj@meta.data[["Diagnosis"]], hME=seurat_obj@meta.data[[paste(modules[i,1])]])
  rownames(hME)=rownames(seurat_obj@meta.data)
  
  #make boxplot
  bxp <-ggboxplot(
    hME, x = "diagnosis", y = "hME", 
    ylab = "hME", xlab = "Diganosis", add = "jitter"
  )
  
  #add to graph
  stat.test <- stat.test %>% add_xy_position(x = "diagnosis")
  bxp <-bxp + 
    stat_pvalue_manual(stat.test, tip.length = 0) +
    labs(title = paste(modules[i,1]))
  assign(paste('bxp_', print(modules[i,1]),sep = ""), bxp)
}


`bxp_Ast-M1`|`bxp_Ast-M2`|`bxp_Ast-M3`|`bxp_Ast-M4`|`bxp_Ast-M5`
`bxp_Ast-M6`|`bxp_Ast-M7`|`bxp_Ast-M8`|`bxp_Ast-M9`|`bxp_Ast-M10`



