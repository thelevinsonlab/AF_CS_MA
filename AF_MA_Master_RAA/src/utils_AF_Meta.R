library(tidyverse)

# Function for Multidimensional Scaling (MDS) plot
myMDS_AF <- function(FuncMatrix, Targets, main = "Title"){ 
  #So that samples are in the same order in the columns of the matrix
  library(limma)
  library(ggplot2)
  #and the id column of Targets
  
  # Run MDS and store results
  MDS_res = plotMDS(FuncMatrix, plot = FALSE)
  
  # Create dataframe for plotting
  myMDS_DF = data.frame(dim1 = MDS_res$x,
                        dim2 = MDS_res$y,
                        AF = Targets$condition)
  
  # Define color palette
  cbPalette = c("#E69F00", "#56B4E9")
  
  # Generate MDS plot
  MDS_plot = ggplot(myMDS_DF,aes(x = dim1, y = dim2)) + geom_point(aes(color = factor(AF))) +
    theme_bw() + scale_colour_manual(values=cbPalette) + theme(panel.grid.major = element_blank(),
                                                               panel.grid.minor = element_blank()) +
    ggtitle(main)
  
  print(MDS_plot)
}


get_tibble_union = function(meta_list, index_name){
  experiments = names(meta_list)
  meta_tible = tibble()
  
  for(e in experiments){
    meta_list[[e]][[index_name]] = mutate(meta_list[[e]][[index_name]],
                                          ExpID = e)
    
    meta_tible = dplyr::bind_rows(meta_list[[e]][[index_name]],
                                  meta_tible)
  }
  
  return(meta_tible)
}



# Function for running differential expression analysis using limma
run_AFlimma <- function(ExpMat, Targets){
  library(limma)
  
  # Fit a linear model to the expression data
  f = factor(Targets$condition, levels= c("AF","SR"))
  design = model.matrix(~0+f)
  colnames(design) = c("yes","no")
  fit = lmFit(ExpMat, design)
  
  #Define contrasts
  cont.matrix = makeContrasts(AF_sign = yes-no,
                              levels=design)
  
  #Empirical Bayes
  fit2 = contrasts.fit(fit, cont.matrix)
  fit2 = eBayes(fit2)
  
  # Obtain differentially expressed genes
  AF_results = as.data.frame(topTable(fit2,adjust.method = "BH",number = Inf)) %>% 
    tibble::rownames_to_column() %>%
    arrange(desc(abs(t))) %>% as_tibble()
  
  # Rename the first column
  colnames(AF_results)[1] = "ID"
  
  return(AF_results)
}

add_exp_id <- function(data, id) {
  data$TARGETS <- data$TARGETS %>% mutate(ExpID = id)
  return(data)
}

# We define now the get_all_limma_AF function here, which retrieves a matrix of certain statistics (p-values, t-values, logFC) from the differential expression results.
get_all_limma_AF <- function(meta_list, limma_column){
  
  sel_cols =  c("ID","ExpID", limma_column)
  
  limma_results = get_tibble_union(meta_list,"AF_limma") %>% 
    dplyr::select(all_of(sel_cols)) %>% 
    spread(sel_cols[2],sel_cols[3])
  
  limma_results_mat =  as.matrix(limma_results[,-1])
  rownames(limma_results_mat) = limma_results[[1]]
  
  return(limma_results_mat)
}

get_complete_gex_AF <- function(meta_list, gex_key = "GEX",
                                complete_targets){
  #METAheart object is a 
  
  all_genes = unique(unlist(lapply(meta_list,
                                   function(x) rownames(x[[gex_key]]))))
  
  gex_union = matrix(NA,nrow = length(all_genes),ncol = nrow(complete_targets))
  
  colnames(gex_union) = complete_targets$grl_id
  rownames(gex_union) = all_genes
  
  for(e in unique(complete_targets$ExpID)){
    fdf = filter(complete_targets,ExpID == e)
    gex = (meta_list[[e]])[[gex_key]]
    genes = rownames(gex)[rownames(gex) %in% all_genes]
    gex_union[genes,fdf$grl_id] = gex
  }
  
  return(gex_union)
  
}

run_anovastats_single = function(numeric_matrix, targets, 
                                 factor_a = "ExpID",
                                 pval = 0.05){
  library(sjstats)
  
  pca_anova_study = apply(numeric_matrix, 1, function(x, targets){
    pc_i = x
    factor_a_vect = factor(targets[[factor_a]])
    gene_aov = aov(pc_i ~ factor_a_vect)
    aov_stats = sjstats::anova_stats(gene_aov) 
    
  },targets = targets) %>% 
    bind_rows(.id = "PC") %>% 
    as_tibble() %>% 
    group_by(PC) %>% 
    gather(stats,value,-(PC:term)) %>% 
    spread(term,value) %>%
    ungroup()
  
  # here we identify PCs associated with the study label
  pcs_study = pca_anova_study %>% 
    filter(stats == "p.value" & factor_a_vect < pval)
  
  return(pcs_study)
}

run_anovastats = function(numeric_matrix, targets, 
                          factor_a = "HeartFailure", factor_b = "ExpID"){
  library(sjstats)
  library(fastmatch)
  
  anova_res = apply(numeric_matrix, 1, function(x, targets){
    
    gene_i = x
    factor_a_vect = factor(targets[[factor_a]])
    factor_b_vect = factor(targets[[factor_b]])
    gene_aov = aov(gene_i ~ factor_a_vect * factor_b_vect)
    aov_stats = sjstats::anova_stats(gene_aov) 
    
  },targets = targets)  %>% bind_rows(.id = "gene") %>% as_tibble() %>% 
    group_by(gene) %>% gather(stats,value,-(gene:term))  %>% spread(term,value) %>%
    ungroup()
  
  colnames(anova_res)[fmatch(c("factor_a_vect","factor_b_vect","factor_a_vect:factor_b_vect"),
                             colnames(anova_res))] = c(factor_a,factor_b,
                                                       paste(factor_a,factor_b,sep="_"))
  
  return(anova_res)
  
}


getRiskIndex_matrix = function(CoefMatrix, ExpressionMatrix){
  
  RiskMatrix = apply(CoefMatrix, 2, function(t_vector){
    
    t_red = na.omit(t_vector)
    t_red = t_red[names(t_red) %in% rownames(ExpressionMatrix)]
    ExpressionMatrix_red = ExpressionMatrix[names(t_red),]
    RiskE = t_red %*% ExpressionMatrix_red
    
  } )
  
  RiskMatrix = scale(RiskMatrix)
  rownames(RiskMatrix) = colnames(ExpressionMatrix)
  
  
  return(RiskMatrix)
}

getRisk_Stats_AF <- function(Experiment_List,genes,limma_t_mat){
  library(ROCR)
  Experiments = names(Experiment_List)
  names(Experiments) = Experiments
  Risk_Stats = lapply(Experiments, function(E){
    
    E_GEX = Experiment_List[[E]]$GEX
    RiskMatrix =  getRiskIndex_matrix(CoefMatrix = limma_t_mat[genes,],
                                      ExpressionMatrix = E_GEX) # Gets risk matrix
    
    SingleAUC = apply(RiskMatrix,2, function(Risk_vect){ #Each column is one estimation of the risk index
      Risk_vect = Risk_vect[Experiment_List[[E]]$TARGETS$sample]
      RiskDF = data.frame("RiskIX" = Risk_vect,
                          "AtrialFibrillation" = ifelse(Experiment_List[[E]]$TARGETS$condition == "AF",1,0),
                          stringsAsFactors = F)
      
      ROC_E = prediction(RiskDF$RiskIX, RiskDF$AtrialFibrillation) #Evaluate classification
      AUC_E = performance(ROC_E, measure = "auc")@y.values[[1]]
      return(AUC_E)
    })
    
    # Consensus score is calculated by taking the mean of each disease scores, without considering
    # the t values of the evaluated experiment
    ALL_Risk = rowMeans(RiskMatrix[Experiment_List[[E]]$TARGETS$sample,colnames(RiskMatrix)!=E])
    
    RiskDF = data.frame("RiskIX" = ALL_Risk,
                        "AtrialFibrillation" = ifelse(Experiment_List[[E]]$TARGETS$condition == "AF",
                                                      1,0),
                        stringsAsFactors = F)
    
    ROC_All = prediction(RiskDF$RiskIX, RiskDF$AtrialFibrillation)
    AUC_All = performance(ROC_All, measure = "auc")@y.values[[1]]
    
    RiskIX_results = list("RiskMatrix" = RiskMatrix,
                          "SingleAUC" = SingleAUC, #AUCs of predicting index experiment with experiments in t_mat
                          "AUC_All" = AUC_All) #AUCs of predicting index experiment with consensus of t_mat
    
    return(RiskIX_results)
    
  })
  return(Risk_Stats) #Index identical to Experiment_List
}

pairwise_ds_AF <- function(experiments, meta_list, ngenes =100, t_matrix){
  
  bind_rows(lapply(experiments, function(e){ #means that we get list of lists
    e_genes = (meta_list[[e]]$AF_limma[[1]])[1:ngenes] #Defining index-specific signature (predictor)
    e_genes = e_genes[e_genes %in% rownames(t_matrix)]
    e_ds_results = getRisk_Stats_AF(Experiment_List = meta_list,
                                    limma_t_mat = t_matrix,
                                    genes = e_genes) #Getting risk Stats
    
    AUC_n_louAUC = bind_rows(lapply(e_ds_results, function(Expment){
      #For each experiment in list (predicted) retrieve
      #performance of predictor Experiment
      return(list("single"= Expment$SingleAUC[e], #How predictor('s t values) predicted each experiment
                  "lou" = Expment$AUC_All)) #How predictor genes made the consensus better (all t's)
    }),.id = "PredictedExperiment")
    
    return(AUC_n_louAUC)
  }),.id = "PredictorExperiment") %>% arrange(PredictorExperiment,PredictedExperiment)
}

pairwise_ES_AF <- function(meta_list, ngenes = 200){
  library(fgsea)
  
  study_deg_list_up = lapply(meta_list, function(x){
    deg = dplyr::slice(x$AF_limma,
                       1:ngenes) %>%
      filter(t >0)
    return(deg[[1]])
  })
  
  study_deg_list_down = lapply(meta_list, function(x){
    deg = dplyr::slice(x$AF_limma,
                       1:ngenes) %>%
      filter(t <0)
    return(deg[[1]])
  })
  
  up_ES = lapply(experiments, function(x){
    
    stat_rank = meta_list[[x]][["AF_limma"]][["t"]]
    names(stat_rank) = meta_list[[x]][["AF_limma"]][["ID"]]
    stat_rank = sort(stat_rank)
    
    up_row = as_tibble(fgsea(pathways = study_deg_list_up,
                             stats = stat_rank,nperm = 1000)) %>%
      dplyr::select(pathway,ES)
  })
  
  up_ES = up_ES %>% 
    enframe("Reference") %>% unnest()
  
  colnames(up_ES) = c("Reference","DEG","ES")
  
  up_ES = mutate(up_ES, direction = "up")
  
  down_ES = lapply(experiments, function(x){
    
    stat_rank = meta_list[[x]][["AF_limma"]][["t"]]
    names(stat_rank) = meta_list[[x]][["AF_limma"]][["ID"]]
    stat_rank = sort(stat_rank)
    
    up_row = as_tibble(fgsea(pathways = study_deg_list_down,
                             stats = stat_rank,nperm = 1000)) %>%
      dplyr::select(pathway,ES)
  })
  
  down_ES = down_ES %>% 
    enframe("Reference") %>% unnest()
  
  colnames(down_ES) = c("Reference","DEG","ES")
  
  down_ES = mutate(down_ES, direction = "down")
  
  all_ES = bind_rows(down_ES, up_ES)
  
  all_ES = mutate(all_ES, number_genes = ngenes)
  
  return(all_ES)
}

run_fisher_meta_AF <- function(meta_list, n_missing = 6){
  library(survcomp)
  # Getting p-values from limma
  limma_pvals = get_all_limma_AF(meta_list = meta_list, "P.Value")
  
  # Use only genes that are present in all experiments (missing in n at most)
  limma_results_mat = limma_pvals[rowSums(is.na(limma_pvals))<=n_missing,]
  
  # Fisher combined test
  fisher_pvals = apply(limma_results_mat, 1, function(x){ 
    survcomp::combine.test(x, "fisher", na.rm = T)
  })
  
  fisher_pvals_adj = sort(p.adjust(fisher_pvals,"BH"))
  
  return(fisher_pvals_adj)
}

run_fisher_meta_AF <- function(meta_list, n_missing = 6){
  library(survcomp)
  # Getting p-values from limma
  limma_pvals = get_all_limma_AF(meta_list = meta_list, "P.Value")
  
  # Use only genes that are present in all experiments (missing in n at most)
  limma_results_mat = limma_pvals[rowSums(is.na(limma_pvals))<=n_missing,]
  
  # Fisher combined test
  fisher_pvals = apply(limma_results_mat, 1, function(x){ 
    survcomp::combine.test(x, "fisher", na.rm = T)
  })
  
  fisher_pvals_adj = sort(p.adjust(fisher_pvals,"BH"))
  
  return(fisher_pvals_adj)
}