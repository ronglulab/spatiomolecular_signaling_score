library(Seurat)
library(CellChat)

single_cell_idents <- function(seurat_obj, hspc_type){
  
  single_cells = subset(seurat_obj, ident = hspc_type)
  Idents(single_cells) <- colnames(single_cells@assays$RNA@counts)
  return(single_cells)
  
}

setup_cellchat <- function(dataset, varfeatpath){
  cellchat <- createCellChat(object = dataset)
  CellChatDB <- updateCellChatDB()
  cellchat@DB <- CellChatDB
  #CellChatDB.use = subsetDB(CellChatDB, search = c("NOTCH", "IL7"), key = "pathway_name")
  #cellchat@DB = CellChatDB.use
  cellchat <- subsetData(cellchat)
  
  
  load(varfeatpath)
  cellchat@var.features = varfeats
  rm(varfeats)
  #cellchat <- identifyOverExpressedGenes(cellchat)
  gc()
  
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  
  return(cellchat)
}

# LOAD IN ALL SINGLE-CELL RNA-SEQ DATA
load("P:/zthomas/Intercellular Interactions Project/Paper Files/code/scRNA-seq/hspc_seurat.rda")
load("P:/zthomas/Intercellular Interactions Project/Paper Files/code/scRNA-seq/bloodAndImmune_seurat.rda")
load("P:/zthomas/Intercellular Interactions Project/Paper Files/code/scRNA-seq/boneMarrowNiche_seurat.rda")

# SELECT HSPC TO DO CELLCHAT AT SINGLE-CELL LEVEL
hspc_for_single_cell = 'HSC'
single_cell_HSPC <- single_cell_idents(HSPC, hspc_for_single_cell)

# MERGE ALL DATA INTO ONE SEURAT OBJECT
bm_hspc <- merge(BoneMarrowNiche, HSPC)
bm_hspc_hm <- merge(bm_hspc, single_cell_HSPC)
All_bm_hspc <- merge(bm_hspc_hm, BloodAndImmune)
# NORMALIZE DATA
All_bm_hspc <- NormalizeData(All_bm_hspc)

# SET UP CELLCHAT OBJECT
cellchat <- setup_cellchat(All_bm_hspc, "P:/zthomas/Intercellular Interactions Project/scRNA-seq/varfeats/cellchat_varfeats.rda")
# REMOVE FILES YOU DON'T NEED, COMPUTING ACROSS THE CELLCHAT OBJECT WILL TAKE LOTS OF ROOM
rm(All_bm_hspc, bm_hspc, bm_hspc_hm, bm_niche, doi, hspc, true_bm)
gc()

# COMPUTE COMMUNICATION ON THE CELLCHAT OBJECT
cellchat <- computeCommunProb(cellchat)
cellchat <- computeCommunProbPathway(cellchat)

# SAVE OBJECT
save(cellchat, 'complete_cellchat_obj.rda')

# EXTRACT PATHWAYS BY SIGNALING TYPE (I.E. T CELLS RECEIVING FROM SINGLE CELL HSPCs)
#------------------------------------------------------------------------------------------------------------------
col_name = c( "Chondrocytes","EC-Arteriar","EC-Arteriolar","EC-Sinusoidal", "Fibroblasts","MSPC-Adipo","MSPC-Osteo","Myofibroblasts","Osteo", "Osteoblasts","Pericytes","Schwann-cells","Smooth-muscle",
              'HSC', 'MPP', 'GMP', 'CLP', 'MEP', 'lin-', 'CMP',
              'B cell','Dendritic cells','Eo/Baso prog.','Ery prog.','Ery/Mk prog.','Erythroblasts','Gran/Mono prog.','LMPPs','Mk prog.','Mono prog.','Monocytes','NK cells','Neutro prog.','Neutrophils','T cells','large pre-B.','pro-B','small pre-B.'
)

for(cluster_type in col_name){
  
  filename = str_replace_all(cluster_type, "[[:punct:]]", "")
  num = colnames(cellchat@netP$prob[,,1])
  
  layer = -1
  for(i in num){
    
    if(grepl(cluster_type, i, fixed = TRUE)){
      layer = i
    }
    
  }
  

  hsc_sending = cellchat@netP$prob[,which(num == layer),]
  hsc_receiving = cellchat@netP$prob[which(num == layer),,]
  
  hsc_s_filename = paste("SINGLECELL_to_", filename, ".csv", sep= "")
  hsc_r_filename = paste(filename, "_to_SINGLECELL.csv", sep="")
  
  write.csv(hsc_sending, file = hsc_s_filename)
  write.csv(hsc_receiving, file = hsc_r_filename)
  
}


# EXTRACT Different types of signaling (Cell Cell Contact, Secreted Signaling, ECM-Receptor)
#------------------------------------------------------------------------------------------------------------------

ssn = which(cellchat@LR$LRsig$annotation == "Secreted Signaling")
cccn = which(cellchat@LR$LRsig$annotation == "Cell-Cell Contact")
ecmr = which(cellchat@LR$LRsig$annotation == "ECM-Receptor")
# 
cc_ss <- cellchat
cc_ss@net$prob = cc_ss@net$prob[,,ssn]
cc_ss@net$pval = cc_ss@net$pval[,,ssn]
# 
cc_ccc <- cellchat
cc_ccc@net$prob <- cc_ccc@net$prob[,,cccn]
cc_ccc@net$pval <- cc_ccc@net$pval[,,cccn]

cc_ecmr <- cellchat
cc_ecmr@net$prob <- cc_ecmr@net$prob[,,ecmr ]
cc_ecmr@net$pval <- cc_ecmr@net$pval[,,ecmr ]
# 
# 
rm(cellchat)
gc()

cc_ss <- aggregateNet(cc_ss)
cc_ccc <- aggregateNet(cc_ccc)
cc_ecmr <- aggregateNet(cc_ecmr)
# 
a = cc_ccc@net$count
b = cc_ccc@net$weight
write.csv(a, file = paste(hspc_for_single_cell, "_ccc_counts.csv", sep=""))
write.csv(b, file = paste(hspc_for_single_cell, "_ccc_weights.csv", sep=""))
rm(a)
rm(b)
rm(cc_ccc)
gc()

a = cc_ss@net$count
b = cc_ss@net$weight
write.csv(a, file = paste(hspc_for_single_cell, "_ss_counts.csv", sep=""))
write.csv(b, file = paste(hspc_for_single_cell, "_ss_weights.csv", sep=""))
rm(a)
rm(b)
rm(cc_ss)
gc()

a = cc_ecmr@net$count
b = cc_ecmr@net$weight
write.csv(a, file = paste(hspc_for_single_cell, "_ecmr_counts.csv", sep=""))
write.csv(b, file = paste(hspc_for_single_cell, "_ecmr_weights.csv", sep=""))
rm(a)
rm(b)
rm(cc_ecmr)
gc()































