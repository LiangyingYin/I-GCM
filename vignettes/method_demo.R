library(IGCM)

load("/mnt/md1/yujia/project/github_repository/IGCM/vignettes/Toydata_var200.RData") # Load toydata
all_feats <- colnames(X_assoc) # For demonstration purpose, feature dimension is low ~200
IGCM_res <- InvariantGCM(X_assoc,Y,envir) # Using the proposed method to identify direct causal variable set
assoc_feats_mat <- IGCM_res$candidate_features # Detected causal variable result matrix from the standard method: PC-simple
assoc_feats <- colnames(assoc_feats_mat) # Detected causal variables
IGCM_all_candidate_causal_set <- IGCM_res$all_candidate_causal_set # Extract all candidate variable set. Notably, the value indicates the end variable index, as the assoc_feats are ranked in descending order by variable importance

causal_set_index <- IGCM_all_candidate_causal_set[1] # Usually we take the 1st variable set as the direct final causal variable set
causal_feats <- assoc_feats[1:causal_set_index] # 

pa_variables # True causal variable set
intersect(pa_variables,assoc_feats) # Detected true causal from PCsimple
intersect(pa_variables,assoc_feats) # Detected true causal from I-GCMq()
