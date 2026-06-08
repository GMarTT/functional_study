library(tidyverse)
library(tseries)
library(openxlsx)
library(FactoMineR)
library(factoextra)
library(ggpubr)
library(RColorBrewer)
library(pheatmap)
library(Hmisc)
library(corrplot)
library(DT)
library(ggrepel)
library(devtools)
library(FonctionsUtiles)
library(lmtest)
library(mixOmics)
library(bestNormalize)
library(kableExtra)
library(fda)
library(funData)
library(reshape2)
library(ggplotify)
library(patchwork)
library(RGCCA)
library(ClustVarLV)
library(ggnewscale)
library(kableExtra)
library(purrr)

# data
dat <- openxlsx::read.xlsx("C:/Users/g.martet/Documents/updatedata/dat.xlsx")
dim(dat)

# correspondance block
var_to_test <- colnames(dat %>%
                          dplyr::select(-annee, -region, -AMR_EHPAD_C3G_R, -AMR_EHPAD_FQ_R,
                                        -AMR_EHPAD_amoxi_clav_R, -AMR_ville_C3G_R,   
                                        -AMR_ville_amoxi_clav_R, -AMR_ville_FQ_R, -AMR_ES_FQ_R, 
                                        -AMR_ES_amoxi_clav_R, -AMR_ES_C3G_R, -log_AMR_ville_FQ_R))
correspondance_block <- data.frame(block = c( 
  rep("Animals", 6), rep("AMC", 20), rep("Humans",8),
  rep("Animals", 3), "AMC", rep("Humans", 4), "Animals", rep("Humans", 6),
  rep("Animals", 5), rep("Humans", 4), rep("Environment", 11),
  rep("Animals", 13)), 
  variable = var_to_test)

dim(correspondance_block)

# les fPC1 par Y
sum.var_AMR_EHPAD_FQ <- openxlsx::read.xlsx("C:/Users/g.martet/Documents/functional_study/sum.var_AMR_EHPAD_FQ_R.xlsx")
sum.var_AMR_ES_FQ <- openxlsx::read.xlsx("C:/Users/g.martet/Documents/functional_study/sum.var_AMR_ES_FQ_R.xlsx")
sum.var_AMR_ville_FQ <- openxlsx::read.xlsx("C:/Users/g.martet/Documents/functional_study/sum.var_AMR_ville_FQ_R.xlsx")
sum.var_AMR_EHPAD_C3G <- openxlsx::read.xlsx("C:/Users/g.martet/Documents/functional_study/sum.var_AMR_EHPAD_C3G_R.xlsx")
sum.var_AMR_ES_C3G <- openxlsx::read.xlsx("C:/Users/g.martet/Documents/functional_study/sum.var_AMR_ES_C3G_R.xlsx")
sum.var_AMR_ville_C3G <- openxlsx::read.xlsx("C:/Users/g.martet/Documents/functional_study/sum.var_AMR_ville_C3G_R.xlsx")

colnames(sum.var_AMR_EHPAD_FQ)[1] <- "region"
colnames(sum.var_AMR_ES_FQ)[1] <- "region"
colnames(sum.var_AMR_ville_FQ)[1] <- "region"
colnames(sum.var_AMR_EHPAD_C3G)[1] <- "region"
colnames(sum.var_AMR_ES_C3G)[1] <- "region"
colnames(sum.var_AMR_ville_C3G)[1] <- "region"

dfs <- list(sum.var_AMR_EHPAD_FQ, 
                                sum.var_AMR_ES_FQ, 
                                sum.var_AMR_ville_FQ, 
                                sum.var_AMR_EHPAD_C3G, 
                                sum.var_AMR_ES_C3G, 
                                sum.var_AMR_ville_C3G)

# Full join successif en détectant automatiquement les colonnes communes
unique_data <- merge(sum.var_AMR_EHPAD_FQ, sum.var_AMR_ES_FQ, 
      by = intersect(names(sum.var_AMR_EHPAD_FQ), names(sum.var_AMR_ES_FQ)))
unique_data <- merge(unique_data, sum.var_AMR_ville_FQ, 
                     by = intersect(names(unique_data), names(sum.var_AMR_ville_FQ)))
unique_data <- merge(unique_data, sum.var_AMR_EHPAD_C3G, 
                     by = intersect(names(unique_data), names(sum.var_AMR_EHPAD_C3G)))
unique_data <- merge(unique_data, sum.var_AMR_ES_C3G, 
                     by = intersect(names(unique_data), names(sum.var_AMR_ES_C3G)))
unique_data <- merge(unique_data, sum.var_AMR_ville_C3G, 
                     by = intersect(names(unique_data), names(sum.var_AMR_ville_C3G)))

sum(is.na(unique_data))


# unique_data <- Reduce(function(x, y) {
#   common_cols <- intersect(names(x), names(y))
#   dplyr::left_join(x, y, by = common_cols)
# }, dfs)

colnames(unique_data)
unique_data <- unique_data %>% 
  dplyr::select(-dplyr::ends_with(".1")) %>% 
  dplyr::select(-dplyr::ends_with(".2")) %>%
  dplyr::select(-dplyr::ends_with(".3"))

sum(is.na(unique_data))

setdiff(correspondance_block$variable, colnames(unique_data))
# variables qui ne sont jamais sélectionnées pour aucun Y

# la fonction
blocks <- levels(factor(correspondance_block$block))
clustvar_function <- function(clustering_Y = "unique",
                              method = "directional", 
                              K_choice, 
                              blocks){
  cls.var2 <- data.frame(variable = character(),
                         var2 = character())
  # clustering unique même pour chaque Y
  if (clustering_Y == "unique"){
    for (i in 1:length(blocks)) {
    var.clv <- unique_data[, which(colnames(unique_data) %in% intersect(var_to_test, colnames(unique_data)))]
    var.clv <- var.clv[, which(colnames(var.clv) %in% correspondance_block[correspondance_block$block == blocks[i], "variable"])]
    res.clv <- ClustVarLV::CLV(var.clv,
                               method = method,
                               sX = TRUE,
                               maxiter = 50,
                               graph = FALSE)
    if (is.na(K_choice[i])){
      d <- as.data.frame(res.clv$tabres)
      l <- c()
      for (i in 0:(nrow(d)-1)){
        l <- c(l, d$agg.crit.hac[nrow(d)-i] - d$agg.crit.hac[nrow(d)-(i+1)])}
      K = which.max(l)+1
      res.clv2 <- ClustVarLV::get_partition(res.clv, K = K)
      d <- data.frame(variable = names(res.clv2),
                      var2 = paste0("cluster_", res.clv2))
      cls.var2 <- rbind.data.frame(cls.var2, d)
    }
    else{
    res.clv2 <- ClustVarLV::get_partition(res.clv, K = K_choice[i])
    d <- data.frame(variable = names(res.clv2),
                    var2 = paste0("cluster_", res.clv2))
    cls.var2 <- rbind.data.frame(cls.var2, d)}
  }}
  else {
    for (k in 1:length(dfs)){
      clus_per_y <- list()
      cls.var2 <- data.frame(variable = character(),
                             var2 = character())
      for (i in 1:length(blocks)) {
      var.clv <- df[, which(colnames(dfs[[k]]) %in% intersect(var_to_test, colnames(df[[k]])))]
      var.clv <- var.clv[, which(colnames(var.clv) %in% correspondance_block[correspondance_block$block == blocks[i], "variable"])]
      res.clv <- ClustVarLV::CLV(var.clv,
                                 method = method,
                                 sX = TRUE,
                                 maxiter = 50,
                                 graph = FALSE)
      if (is.na(K_choice[i])){
        d <- as.data.frame(res.clv$tabres)
        l <- c()
        for (i in 0:(nrow(d)-1)){
          l <- c(l, d$agg.crit.hac[nrow(d)-i] - d$agg.crit.hac[nrow(d)-(i+1)])}
        K = which.max(l)+1
        res.clv2 <- ClustVarLV::get_partition(res.clv, K = K)
        d <- data.frame(variable = names(res.clv2),
                        var2 = paste0("cluster_", res.clv2))
        cls.var2 <- rbind.data.frame(cls.var2, d)
      }
      else{
        res.clv2 <- ClustVarLV::get_partition(res.clv, K = K_choice[i])
        d <- data.frame(variable = names(res.clv2),
                        var2 = paste0("cluster_", res.clv2))
        cls.var2 <- rbind.data.frame(cls.var2, d)}
      }
      clus_per_y <- c(clus_per_y, cls.var2)
      
    }
    return(clus_per_y)}
}

clustvar_function(clustering_Y = "by Y", 
                  method = "local", 
                  K_choice = c(4, 4, 4, 3),
                  blocks= blocks)
