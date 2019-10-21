## ---------------------------
## Purpose of script: grouping contigs based on normalized coverage for svalbard permafrost metagenome
## Author: Yaxin Xue
## Email: xue.ethan@gmail.com
## ---------------------------
# change work_dir to your dir
work_dir = '/Users/yaxin/Workspace/Coverage-based-functional-analyisis-in-a-MAG-centric-view'
setwd(work_dir)

library(data.table)

# Grouping functions, KI and KD are defined as specical groups in PL correelated with depth
identify_normal_group <- function(cov_data, cov_cutoff){
  # Group UN means Abundannt in AL and some PL samples
  group = 'UN'
  al_cov = cov_data[1]
  pl_cov = cov_data[2:5]
  # Both
  if (all(cov_data >= cov_cutoff)){
    group = 'BO'
  }
  # Low abundance
  if (all(cov_data <= cov_cutoff)){
    group = 'LO'
  }
  # Active layer
  if (al_cov >= cov_cutoff && all(pl_cov <= cov_cutoff)){
    group = 'AL'
  }
  # At least one PL is higher than TH
  if (al_cov <= cov_cutoff && any(pl_cov >= cov_cutoff)){
    group = 'PL_SUB'
    if (all(pl_cov >= cov_cutoff)){
      group = 'PL_ALL'
    }
    else{
      # P1...P4
      pl_above_idx = which(as.numeric(pl_cov >= cov_cutoff)==1)
      if(length(pl_above_idx)==1){
        group = paste0('PL_P', pl_above_idx)
      }
    }
  }
  return(group)
}

identify_special_group <- function(pl_cov, sample_depth){
  group = 'UN'
  if (sd(pl_cov)!=0){
    pl_cor = cor(sample_depth, pl_cov)
    if (pl_cor >= 0.9){
      group = 'KI'
    }
    if (pl_cor <= -0.9){
      group = 'KD'
    }
  }
  return (group)
}


# define groups and sample depth
normal_groups = c('AL','BO','PL_SUB','PL_ALL','PL_P1','PL_P2','PL_P3','PL_P4','LO','UN')
special_groups = c('KI','KD')
sample_depth = c(1.10, 1.22, 1.35, 1.70)

#
ctg_cov_norm = fread('ctg_bin_cov.norm.csv')
# set median as coverage threshold
cov_cutoff = median(sapply(ctg_cov_norm[,3:7], median))
# identify normal groups
NORM_ctg_cov_wg = ctg_cov_norm[,c(1,2)]
NORM_ctg_cov_wg$group = apply(ctg_cov_norm[,3:7], 1, identify_normal_group, cov_cutoff = cov_cutoff)
NORM_ctg_cov_wg = cbind(NORM_ctg_cov_wg, ctg_cov_norm[,3:7])
# identify special groups
SPE_ctg_cov_wg = NORM_ctg_cov_wg[group %in% c('PL_SUB','PL_ALL')]
SPE_ctg_cov_wg$SPE_group = apply(SPE_ctg_cov_wg[,4:7],1, identify_special_group, sample_depth = sample_depth)
SPE_ctg_cov_wg = SPE_ctg_cov_wg[SPE_group %in% c('KI','KD')]

# save output with group column
# reorder columns with this order: Contig_ID,Bin_ID,group,S1,S2,....,SN

write.csv(NORM_ctg_cov_wg, 'NORM_group.csv',quote = F,row.names = F)
write.csv(SPE_ctg_cov_wg, 'SPE_group.csv',quote = F,row.names = F)
