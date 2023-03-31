library(dplyr)

Sort_Unified_F <- read.delim("Sort_Unified_F.gff3",sep = "\t",header = FALSE)
names(Sort_Unified_F) <- c("chrloc","source","type","start","end","score","strand","phase","modified_attributes","attributes","gp_mod_attr","Overlap_Value")

rectify_bad_ensembl_models <- function(NCBI_Ens_XGC_Data){
  
  ###### for gp_mod_attr steps
  
  ### part_gene assigned
  
  if(unique(part_gene$modified_attributes) > 1 & unique(part_gene$source) == "NCBI" & 
     unique(part_gene$source) ==  "ensembl" & 
     unique(part_gene$source) == "Genbank" & 
     (length(part_gene[part_gene$type == "mRNA" & part_gene$source == "ensembl",]) > length(part_gene[part_gene$type == "mRNA" & part_gene$source == "NCBI",]) | length((part_gene[part_gene$type == "mRNA" & part_gene$source == "ensembl",])) > length(part_gene[part_gene$type == "mRNA" & part_gene$source == "Genbank",]))){
    
    other_ls[[length(other_ls)+1]] <- part_gene
    
    part_gene <- part_gene %>% filter(source == "NCBI" & source == "Genbank")
    
    ###if aaron wants to breakup all complicated models into 1:1 models,that is possible-
    
    #for(j in 1 : length(unique(part_gene$modified_attributes))){
    
    
    # filtered_attr_ls <- part_gene %>% filter(modified_attributes %in% unique(part_gene$modified_attributes)[[j]])
    
    ### and the rest to be taken from modgffsort 
    
    
    #} 
  }  

}
