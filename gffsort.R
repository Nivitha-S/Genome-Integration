
##########################################################

##### loading packages ############

library(dplyr)
library(tidyr)
library(stringr)

##########################################################


sort_NCBI_Ensembl_Data <- function(Processed_NCBI_Ens_Data,Mod_Data,rank_data) {

#Mod_Data <- rbind.data.frame(Processed_NCBI_Ens_Data[,c('chrloc','source','type','start','end','score','strand','phase','modified_attributes','attributes')],Unique_NCBI,Unique_Ens)  
  
filtered_df <- Mod_Data %>% left_join(rank_data,by = c("type")) 

## analyzing the architechure of the merged Data NCBI Ensembl 

gp_filtered_data <- filtered_df %>% group_by(chrloc,source,strand,modified_attributes) %>% arrange(rank,.by_group = TRUE) %>% ungroup() %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,rank,attributes) %>% distinct(.keep_all = TRUE) 

gp_filtered_NCBI_df <- gp_filtered_data %>% filter(source == "NCBI")

gp_filtered_NCBI_df$transcript <- sapply(gp_filtered_NCBI_df$attributes,function(x) 
  
  ifelse(str_detect(x,pattern = "transcript_id=.+"),
         
         gsub("transcript_id=","",str_extract(x,pattern = "transcript_id=.+")),
         
         gsub("Parent=.+-","",str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))
  
)


gp_filtered_NCBI_d <- gp_filtered_NCBI_df %>% mutate(
  
  
  trans =  ifelse( is.na(transcript),
                   gsub("Name=","",str_split(str_extract(attributes,pattern = "Name=.+;"),pattern = ";")[[1]][1]),
                   transcript
                   
  ) ) 


gp_filtered_NCBI_data <- gp_filtered_NCBI_d %>% group_by(chrloc,source,strand,modified_attributes) %>% arrange(trans,.by_group = TRUE) %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)

#########################################################################################


gp_filtered_Ensembl_df <- gp_filtered_data %>% filter(source == "ensembl")

gp_filtered_Ensembl_df$transcript <- sapply(gp_filtered_Ensembl_df$attributes,function(x)
  
  
  ifelse(str_detect(x,pattern = "transcript:.+"),
  gsub("transcript:","",str_split(str_extract(x,pattern = "transcript:.+;"),pattern = ";")[[1]][1]),
  gsub('Name=',"",str_split(str_extract(x,pattern = "Name=.+;"),pattern = ";")[[1]][1])       
   )

)

gp_filtered_Ensembl_d <- gp_filtered_Ensembl_df %>% group_by(chrloc,strand,modified_attributes) %>% arrange(transcript,.by_group = TRUE) %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)

#########################################################################################

COmbined_grped_NCBI_Ens_df <- rbind.data.frame(gp_filtered_NCBI_data,gp_filtered_Ensembl_d)

cbind_filtered_NCBI_data <- Processed_NCBI_Ens_Data %>% filter(source == "NCBI")

cbind_filtered_ensembl_data <- Processed_NCBI_Ens_Data %>% filter(source == "ensembl")

cbind_total_data <- cbind.data.frame(cbind_filtered_NCBI_data,cbind_filtered_ensembl_data)

names(cbind_total_data) <- c("chrloc_1","source1","type1","start1","end1","score1","strand1","phase1","modattr1","attr1","chrloc_2","source2","type2","start2","end2","score2","strand2","phase2","modattr2","attr2")

filtered_cbind_mod_attr_d <- cbind_total_data %>% select(modattr1,modattr2) %>% distinct(.keep_all = TRUE)

Combined_grped_NCBI_Ens_df <- COmbined_grped_NCBI_Ens_df %>% left_join(filtered_cbind_mod_attr_d,by = c("modified_attributes" = "modattr2"))

Combined_df <- Combined_grped_NCBI_Ens_df %>% mutate(
  
  
     gp_mod_attr = ifelse(is.na(modattr1),modified_attributes,modattr1)
  
  
)  

Final_grped_Data <- Combined_df %>% group_by(chrloc,gp_mod_attr) %>% arrange(gp_mod_attr,.by_group = TRUE) %>% ungroup() %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)

Final_grped_Data

}

######################################################################

#### Loading Data and callinf function sort_NCBI_Ensembl_Data #####

#########################################################################################

Processed_NCBI_Ens_75s_Raw_Data <- read.delim('Processed_NCBI_Ens_Fjoin_75s_Raw.gff3',header = FALSE,sep = '\t') 

names(Processed_NCBI_Ens_75s_Raw_Data) <- c('chrloc','source','type','start','end','score','strand','phase','modified_attributes','attributes','Overlap_Values')

rank_data <- data.frame(type = c("gene","ncRNA_gene","pseudogene","mRNA",
"transcript","primary_transcript","pseudogenic_transcript","exon","CDS",
"V_gene_segment","C_gene_segment","tRNA","rRNA","lnc_RNA","snoRNA","snRNA","scRNA",
"miRNA","Y_RNA","ncRNA","guide_RNA","match","origin_of_replication","D_loop"),rank = 1:24)

Unique_NCBI_Models <- read.delim('Unique_NCBI_Models_75s.gff3',header = FALSE,sep = '\t')

names(Unique_NCBI_Models) <- c('chrloc','source','type','start','end','score','strand','phase','modified_attributes','attributes')

Unique_Ens_Models <- read.delim('Unique_Ensembl_Models_75s.gff3',header = FALSE,sep = '\t')

names(Unique_Ens_Models) <- c('chrloc','source','type','start','end','score','strand','phase','modified_attributes','attributes')

#Unique NCBI and Unique Ensembl Data are retrieved through compute_NCBI_Ensembl_Fjoin_Model* Data

Processed_NCBI_Ens_75s_D <- rbind.data.frame(Processed_NCBI_Ens_75s_Raw_Data[,c('chrloc','source','type','start','end','score','strand','phase','modified_attributes','attributes')],Unique_NCBI_Models,Unique_Ens_Models)

Final_grped_Data <- sort_NCBI_Ensembl_Data(Processed_NCBI_Ens_75s_Raw_Data,Processed_NCBI_Ens_75s_D,rank_data)

Final_grped_D <- Final_grped_Data %>% filter(type == "exon") %>% select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% distinct(.keep_all = TRUE)

write.table(Final_grped_D,file = "NCBI_Ensembl_Mod_Final.gff3",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

write.table(Final_grped_Data,file = "Processed_NCBI_Ensembl_75s_Fjoin.gff3",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

##########################################################################################

########## Sorting out the exons/CDS sub-segments ##########

part_gene <- Final_grped_Data %>% filter(modified_attributes == "aasdh")

part_gene_1 <- Final_grped_Data %>% filter(modified_attributes == "znf185")

#######################################################################################

## sort NCBI Ensembl XGC Data 

sort_NCBI_Ensembl_XGC_Data <- function(Mod_Data,cbind_total_data,rank_data) {
  
## analyzing the architechure of the merged Data NCBI Ensembl 
  
gp_filtered_NCBI_df <- Mod_Data %>% filter(source == "NCBI") %>% distinct(.keep_all = TRUE)

filtered_NCBI_df <- gp_filtered_NCBI_df %>% left_join(rank_data,by = c("type")) 

## analyzing the architechure of the merged Data NCBI Ensembl 

gp_filtered_NCBI_data <- filtered_NCBI_df %>% group_by(chrloc,strand,modified_attributes) %>% 
arrange(rank,.by_group = TRUE) %>% ungroup() %>% 
select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,rank,attributes) %>% distinct(.keep_all = TRUE) 
  
gp_filtered_NCBI_data$transcript <- sapply(gp_filtered_NCBI_data$attributes,function(x) 
    
    ifelse(str_detect(x,pattern = "transcript_id=.+"),
           
           gsub("transcript_id=","",str_extract(x,pattern = "transcript_id=.+")),
           
           gsub("Parent=.+-","",str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))
    
  )
  
  
  gp_filtered_NCBI_d <- gp_filtered_NCBI_data %>% mutate(
    
    
    trans =  ifelse( is.na(transcript),
                     gsub("Name=","",str_split(str_extract(attributes,pattern = "Name=.+;"),pattern = ";")[[1]][1]),
                     transcript
                     
    ) ) 
  
  
  gp_filtered_NCBI <- gp_filtered_NCBI_d %>% group_by(chrloc,strand,modified_attributes) %>% arrange(trans,.by_group = TRUE) %>% select(chrloc,source,type,rank,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
  
  
  #########################################################################################
  
  gp_filtered_Ensembl_df <- Mod_Data %>% filter(source == "ensembl") %>% distinct(.keep_all = TRUE)
  
  filtered_Ens_df <- gp_filtered_Ensembl_df %>% left_join(rank_data,by = c("type")) 
  
  ## analyzing the architecture of the merged Data NCBI Ensembl 
  
  gp_filtered_Ens_data <- filtered_Ens_df %>% group_by(chrloc,strand,modified_attributes) %>% arrange(rank,.by_group = TRUE) %>% ungroup() %>% select(chrloc,source,type,rank,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE) 
  
  gp_filtered_Ens_data$transcript <- sapply(gp_filtered_Ens_data$attributes,function(x)
    
    
    ifelse(str_detect(x,pattern = "transcript:.+"),
           gsub("transcript:","",str_split(str_extract(x,pattern = "transcript:.+;"),pattern = ";")[[1]][1]),
           gsub('Name=',"",str_split(str_extract(x,pattern = "Name=.+;"),pattern = ";")[[1]][1])       
    )
    
  )
  
  gp_filtered_Ensembl_d <- gp_filtered_Ens_data %>% group_by(chrloc,strand,modified_attributes) %>% arrange(transcript,.by_group = TRUE) %>% select(chrloc,source,type,rank,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
  
  #########################################################################################
  
  gp_filtered_XGC_df <- Mod_Data %>% filter(source == "Genbank") %>% distinct(.keep_all = TRUE)
  
  filtered_XGC_df <- gp_filtered_XGC_df %>% left_join(rank_data,by = c("type")) 
  
  ## analyzing the architechure of the merged Data NCBI Ensembl 
  
  gp_filtered_XGC_data <- filtered_XGC_df %>% group_by(chrloc,strand,modified_attributes) %>% arrange(rank,.by_group = TRUE) %>% ungroup() %>% select(chrloc,source,type,rank,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE) 
  
  gp_filtered_XGC_data$transcript <- sapply(gp_filtered_XGC_data$attributes,function(x)
    
    ifelse(str_detect(x,pattern = "orig_transcript_id=.+"),
           gsub("orig_transcript_id=","",str_split(str_extract(x,pattern = "orig_transcript_id=.+;"),pattern = ";")[[1]][1]),
           gsub("Name=","",str_split(str_extract(x,pattern = "Name=.+;"),pattern = ";")[[1]][1])       
    )
    
  )
  
  gp_filtered_XGC_df <- gp_filtered_XGC_data %>% group_by(chrloc,strand,modified_attributes) %>% arrange(transcript,.by_group = TRUE) %>% select(chrloc,source,type,rank,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE) 
  
  #######################################################################################
  
  #Combined_grped_NCBI_Ens_XGC_df <- rbind.data.frame(gp_filtered_NCBI,gp_filtered_Ensembl_d,gp_filtered_XGC_df)
  
  rm(Mod_Data)
  
  #cbind_filtered_NE_ls <- list()
  #cbind_filtered_XGC_ls <- list()
  
   #for (i in 1 : nrow(Processed_NCBI_Ens_XGC_Data)) {
     
    # if(i %% 2 != 0){
       
     #  cbind_filtered_NE_ls[[length(cbind_filtered_NE_ls)+1]] <- Processed_NCBI_Ens_XGC_Data[i,c("chrloc","strand","modified_attributes")]
       
     #}
     
     #else if(i %% 2 == 0){
       
      # cbind_filtered_XGC_ls[[length(cbind_filtered_XGC_ls)+1]] <- Processed_NCBI_Ens_XGC_Data[i,c("chrloc","strand","modified_attributes")]
       
     #}
     
   #}
  
  #rm(Processed_NCBI_Ens_XGC_Data)
  
  #cbind_filtered_NE_df <- do.call("rbind.data.frame",cbind_filtered_NE_ls)
    
  #cbind_filtered_XGC_df <- do.call("rbind.data.frame",cbind_filtered_XGC_ls)
  
  #rm(cbind_filtered_NE_ls)
  
  #rm(cbind_filtered_XGC_ls)
  
  #cbind_total_data <- cbind.data.frame(cbind_filtered_XGC_df,cbind_filtered_NE_df)
  
  #names(cbind_total_data) <- c("chrloc1","strand1","modattr1","chrloc2","strand2","modattr2")
  
  #filtered_cbind_mod_attr_d <- cbind_total_data %>% filter(modattr1,chrloc2,strand2,modattr2) %>% distinct(.keep_all = TRUE)
  
  NCBI_d <- gp_filtered_NCBI %>% left_join(cbind_total_data,by = c("chrloc" = "chrloc2","strand" = "strand2","modified_attributes" = "modattr2"))
  
  Ens_d <- gp_filtered_Ensembl_d %>% left_join(cbind_total_data,by = c("chrloc" = "chrloc2","strand" = "strand2","modified_attributes" = "modattr2"))
  
  XGC_d <- gp_filtered_XGC_df %>% left_join(cbind_total_data,by = c("chrloc" = "chrloc2","strand" = "strand2","modified_attributes" = "modattr2"))
  
  NCBI_df <- NCBI_d %>% mutate(
    
   gp_mod_attr = ifelse(is.na(modattr1),modified_attributes,modattr1)
    
  )  
  
  Ens_df <- Ens_d %>% mutate(
    
    gp_mod_attr = ifelse(is.na(modattr1),modified_attributes,modattr1)
    
  )  
  
  XGC_df <- XGC_d %>% mutate(
    
    gp_mod_attr = ifelse(is.na(modattr1),modified_attributes,modattr1)
    
  )  
  
  Combined_df <- rbind.data.frame(NCBI_df,Ens_df,XGC_df)
  
  rm(NCBI_df)
  rm(Ens_df)
  rm(XGC_df)
  
  Final_grped_Data <- Combined_df %>% group_by(chrloc,strand,gp_mod_attr) %>% arrange(gp_mod_attr,.by_group = TRUE) %>% ungroup() %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes,gp_mod_attr) %>% distinct(.keep_all = TRUE)
  
  return(Final_grped_Data)
  
}

####################################################################################

#Complete_NCBI_Ens_XGC_df <- read.delim('Complete_NCBI_Ens_XGC_75s_50.gff3',header = FALSE,sep = '\t') 

NCBI_Ens_75s_XGC_50_Fjoin_Data <- read.delim('NCBI_Ens_75s_XGC_50.gff3',sep = '\t',header = FALSE)

cbind_total_data <- read.csv("cbind_total_data.csv")

#### arguments parsed from working directory #####

names(NCBI_Ens_75s_XGC_50_Fjoin_Data) <- c('OverlapVal','chrloc','source','type','start','end','score','strand','phase','modified_attributes','attributes')

names(Complete_NCBI_Ens_XGC_df) <- c('chrloc','source','type','start','end','score','strand','phase','modified_attributes','attributes')

######## obtaining the raw data from statistics file and working directory 

rank_data <- data.frame(type = c("gene","ncRNA_gene","pseudogene","mRNA","tRNA","rRNA","lnc_RNA","snoRNA","snRNA","scRNA",
                                 "miRNA","Y_RNA","ncRNA","guide_RNA","transcript","primary_transcript","pseudogenic_transcript","exon","CDS",
                                 "V_gene_segment","C_gene_segment","match","origin_of_replication","D_loop"),rank = 1:24)

Final_grped_Data <- sort_NCBI_Ensembl_XGC_Data(Complete_NCBI_Ens_XGC_df,cbind_total_data,rank_data)

Final_grped_Full_Data <- Final_grped_Data %>% select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% distinct(.keep_all = TRUE)

write.table(Final_grped_Data,file = "Final_grped_Data.gff3",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

write.table(Final_grped_Full_Data,file = "NCBI_Ensembl_Full_Data.gff3",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

##################################################################################

cbind_NCBI_Ens_XGC <- function(NCBI_Ens_XGC_Data){
  
  cbind_filtered_NE_ls <- list()
  cbind_filtered_XGC_ls <- list()
  
  for (i in 1 : nrow(NCBI_Ens_XGC_Data)) {
    
    if(i %% 2 != 0){
      
      cbind_filtered_NE_ls[[length(cbind_filtered_NE_ls)+1]] <- 
        NCBI_Ens_XGC_Data[i,]
      
      
    }
    
    else if(i %% 2 == 0){
      
      cbind_filtered_XGC_ls[[length(cbind_filtered_XGC_ls)+1]] <- 
        NCBI_Ens_XGC_Data[i,]
      
    }
    
  }
  
  cbind_filtered_NE_df <- do.call("rbind.data.frame",cbind_filtered_NE_ls)
  
  cbind_filtered_XGC_df <- do.call("rbind.data.frame",cbind_filtered_XGC_ls)
  
  names(cbind_filtered_XGC_df) <- c("OverlapVal1","chrloc1","source1","type1","start1","end1","score1","strand1","phase1","modattr1","attr1")
  
  names(cbind_filtered_NE_df) <- c("OverlapVal2","chrloc2","source2","type2","start2","end2","score2","strand2","phase2","modattr2","attr2")
  
  cbind_total_data <- cbind.data.frame(cbind_filtered_XGC_df,cbind_filtered_NE_df) %>% distinct(.keep_all = TRUE)  
  
  
}

#cbind_total_df <- cbind_NCBI_Ens_XGC(NCBI_Ens_75s_XGC_50_Fjoin_Data)

filtered_NCBI_Ens_XGC_ls <- list()
cbind_NCBI_Ens_XGC_ls <- list()

for (i in 1: length(unique(NCBI_Ens_75s_XGC_50_Fjoin_Data$chrloc))){
  
  filtered_NCBI_Ens_XGC_ls[[i]] <- NCBI_Ens_75s_XGC_50_Fjoin_Data %>% filter(chrloc %in% unique(NCBI_Ens_75s_XGC_50_Fjoin_Data$chrloc)[i])
  cbind_NCBI_Ens_XGC_ls[[i]] <- cbind_NCBI_Ens_XGC(filtered_NCBI_Ens_XGC_ls[[i]])
  
}

cbind_total_d <- do.call("rbind.data.frame",cbind_NCBI_Ens_XGC_ls)

cbind_total_df <- cbind_total_d[,c("chrloc2","strand2","modattr2","modattr1")] %>% distinct(.keep_all = FALSE) 

write.csv(cbind_total_df,file = "cbind_total_data.csv",row.names = FALSE,quote = FALSE)


count_NCBI_Ens_XGC <- function(cbind_NCBI_Ens_XGC_ls){
  
  cbind_total_data <- do.call("rbind.data.frame",cbind_NCBI_Ens_XGC_ls)
  
  names(cbind_total_data) <- c("OverlapVal1","chrloc_1","source1","type1","start1","end1","score1","strand1","phase1",
  "modattr1","attr1","OverlapVal2","chrloc_2","source2","type2","start2","end2","score2","strand2","phase2",
  "modattr2","attr2")
  
  NCBI_XGC_Data <- cbind_total_data %>% filter(source1 == "NCBI" & source2 == "Genbank" & 
  (type1 == "gene" | type1 == "pseudogene") & (type2 == "gene" | type2 == "pseudogene")) %>% distinct(.keep_all = TRUE) 
  
  Ens_XGC_Data <- cbind_total_data %>% filter(source1 == "ensembl" & source2 == "Genbank" & (type1 == "gene" | type1 == "pseudogene") & (type2 == "gene" | type2 == "pseudogene")) %>% distinct(.keep_all = TRUE) 
  
  ###### filtering out models that overlaps with NCBI ensembl and XGC ########

  filtered_NCBI_Ens_XGC <- NCBI_XGC_Data %>% full_join(Ens_XGC_Data,by = c("chrloc_2","source2","type2","start2","end2","score2","strand2","phase2","modattr2","attr2")) %>% filter(!is.na(modattr1.x) & !is.na(modattr1.y) & !is.na(modattr2)) %>% select(chrloc_1.x,start1.x,end1.x,strand1.x,modattr1.x,chrloc_1.y,start1.y,end1.y,strand1.y,modattr1.y,
  chrloc_2,start2,end2,strand2,modattr2) %>% distinct(.keep_all = TRUE)
    
  filtered_NCBI_Ens_XGC_df <- filtered_NCBI_Ens_XGC %>% mutate(
    
    Combined_attr_1 = paste0(start1.x,"_",end1.x),
    Combined_attr_2 = paste0(start1.y,"_",end1.y),
    Combined_attr_3 = paste0(start2,"_",end2)
    
  )
  
  Filtered_1_to_Many_NE_Model <- filtered_NCBI_Ens_XGC_df %>% 
    group_by(chrloc_1.x,start1.x,end1.x,strand1.x) %>% mutate(
    
    One_to_Many_Relationship_NE = paste0(unique(Combined_attr_2),collapse = ";"),
    One_to_Many_Relationship_NE_Count = length(unique(Combined_attr_2))
    
    
  ) 
  
  Filtered_Many_to_one_NE_Model <- Filtered_1_to_Many_NE_Model %>% 
    group_by(chrloc_1.y,start1.y,end1.y,strand1.y) %>% mutate(
      
      Many_to_One_Relationship_NE = paste0(unique(Combined_attr_1),collapse = ";"),
      Many_to_One_Relationship_NE_Count = length(unique(Combined_attr_1))
      
    ) 
  
  Filtered_1_to_Many_NX_Model <- Filtered_Many_to_one_NE_Model %>% 
    group_by(chrloc_1.x,start1.x,end1.x,strand1.x) %>% mutate(
      
      One_to_Many_Relationship_NX = paste0(unique(Combined_attr_3),collapse = ";"),
      One_to_Many_Relationship_NX_Count = length(unique(Combined_attr_3))
      
      
    ) 
  
  Filtered_Many_to_1_NX_Model <- Filtered_1_to_Many_NX_Model %>% 
    group_by(chrloc_2,start2,end2,strand2) %>% mutate(
      
      Many_to_One_Relationship_NX = paste0(unique(Combined_attr_1),collapse = ";"),
      Many_to_One_Relationship_NX_Count = length(unique(Combined_attr_1))
      
    ) 
  
  Filtered_1_to_Many_EX_Model <- Filtered_Many_to_1_NX_Model %>% 
    group_by(chrloc_1.y,start1.y,end1.y,strand1.y) %>% mutate(
      
      One_to_Many_Relationship_EX = paste0(unique(Combined_attr_3),collapse = ";"),
      One_to_Many_Relationship_EX_Count = length(unique(Combined_attr_3))
      
      
    ) 
  
  Filtered_Final_Model <- Filtered_1_to_Many_EX_Model %>% 
    group_by(chrloc_2,start2,end2,strand2) %>% mutate(
      
      Many_to_One_Relationship_EX = paste0(unique(Combined_attr_2),collapse = ";"),
      Many_to_One_Relationship_EX_Count = length(unique(Combined_attr_2))
      
    ) 

filtered_1_to_1_to_1 <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NE_Count == 1 & 
                                                Filtered_Final_Model$Many_to_One_Relationship_NE_Count == 1 &
                                                Filtered_Final_Model$One_to_Many_Relationship_NX_Count == 1 &
                                                Filtered_Final_Model$Many_to_One_Relationship_NX_Count == 1 &
                                                Filtered_Final_Model$One_to_Many_Relationship_EX_Count == 1 &
                                                Filtered_Final_Model$Many_to_One_Relationship_EX_Count == 1,]

filtered_complicated <-  Filtered_Final_Model[!(Filtered_Final_Model$One_to_Many_Relationship_NE_Count == 1 & 
                                                Filtered_Final_Model$Many_to_One_Relationship_NE_Count == 1 &
                                                Filtered_Final_Model$One_to_Many_Relationship_NX_Count == 1 &
                                                Filtered_Final_Model$Many_to_One_Relationship_NX_Count == 1 &
                                                Filtered_Final_Model$One_to_Many_Relationship_EX_Count == 1 &
                                                Filtered_Final_Model$Many_to_One_Relationship_EX_Count == 1),]

################################# NCBI - Ens ######################################

NCBI_Ens_One_to_Many <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NE_Count > 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_NE_Count == 1,] %>% ungroup() %>% select(chrloc_1.x,start1.x,end1.x,strand1.x,One_to_Many_Relationship_NE,
One_to_Many_Relationship_NE_Count,Many_to_One_Relationship_NE,Many_to_One_Relationship_NE_Count) %>% distinct(.keep_all = TRUE) 


NCBI_Ens_Many_to_One <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NE_Count == 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_NE_Count > 1,] %>% ungroup() %>% select(chrloc_1.y,start1.y,end1.y,strand1.y,One_to_Many_Relationship_NE,
           One_to_Many_Relationship_NE_Count,Many_to_One_Relationship_NE,Many_to_One_Relationship_NE_Count) %>% distinct(.keep_all = TRUE) 

NCBI_Ens_Many_to_Many <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NE_Count > 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_NE_Count > 1,] %>% ungroup() %>% select(chrloc_1.x,start1.x,end1.x,strand1.x,chrloc_1.y,start1.y,end1.y,strand1.y,One_to_Many_Relationship_NE,
    One_to_Many_Relationship_NE_Count,Many_to_One_Relationship_NE,Many_to_One_Relationship_NE_Count) %>% distinct(.keep_all = TRUE)

NCBI_Ens_One_to_One <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NE_Count == 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_NE_Count == 1,] %>% ungroup() %>% select(chrloc_1.x,start1.x,end1.x,strand1.x,chrloc_1.y,start1.y,end1.y,strand1.y,One_to_Many_Relationship_NE,
           One_to_Many_Relationship_NE_Count,Many_to_One_Relationship_NE,Many_to_One_Relationship_NE_Count) %>% distinct(.keep_all = TRUE)

############################## NCBI - XGC #########################

NCBI_XGC_One_to_Many <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NX_Count > 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_NX_Count == 1,] %>% ungroup() %>% select(chrloc_1.x,start1.x,end1.x,strand1.x,One_to_Many_Relationship_NX,
    One_to_Many_Relationship_NX_Count,Many_to_One_Relationship_NX,Many_to_One_Relationship_NX_Count) %>% distinct(.keep_all = TRUE) 

NCBI_XGC_Many_to_One <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NX_Count == 1 & 
Filtered_Final_Model$Many_to_One_Relationship_NX_Count > 1,] %>% ungroup() %>% select(chrloc_2,start2,end2,strand2,One_to_Many_Relationship_NX,
One_to_Many_Relationship_NX_Count,Many_to_One_Relationship_NX,Many_to_One_Relationship_NX_Count) %>% distinct(.keep_all = TRUE) 

NCBI_XGC_Many_to_Many <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NX_Count > 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_NX_Count > 1,] %>% ungroup() %>% select(chrloc_1.x,start1.x,end1.x,strand1.x,chrloc_2,start2,end2,strand2,One_to_Many_Relationship_NX,
     One_to_Many_Relationship_NX_Count,Many_to_One_Relationship_NX,Many_to_One_Relationship_NX_Count) %>% distinct(.keep_all = TRUE)

NCBI_XGC_One_to_One <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_NX_Count == 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_NX_Count == 1,] %>% ungroup() %>% select(chrloc_1.x,start1.x,end1.x,strand1.x,chrloc_2,start2,end2,strand2,One_to_Many_Relationship_NX,
           One_to_Many_Relationship_NX_Count,Many_to_One_Relationship_NX,Many_to_One_Relationship_NX_Count) %>% distinct(.keep_all = TRUE)

############################## Ens - XGC   #################################################

Ens_XGC_One_to_Many <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_EX_Count > 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_EX_Count == 1,] %>% ungroup() %>% select(chrloc_1.y,start1.y,end1.y,strand1.y,One_to_Many_Relationship_EX,
      One_to_Many_Relationship_EX_Count,Many_to_One_Relationship_EX,Many_to_One_Relationship_EX_Count) %>% distinct(.keep_all = TRUE) 


Ens_XGC_Many_to_One <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_EX_Count == 1 & 
  Filtered_Final_Model$Many_to_One_Relationship_EX_Count > 1,] %>% ungroup() %>% select(chrloc_2,start2,end2,strand2,One_to_Many_Relationship_EX,
      One_to_Many_Relationship_EX_Count,Many_to_One_Relationship_EX,Many_to_One_Relationship_EX_Count) %>% distinct(.keep_all = TRUE) 

Ens_XGC_Many_to_Many <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_EX_Count > 1 & 
    Filtered_Final_Model$Many_to_One_Relationship_EX_Count > 1,] %>% ungroup() %>% select(chrloc_1.y,start1.y,end1.y,strand1.y,chrloc_2,start2,end2,strand2,One_to_Many_Relationship_EX,
      One_to_Many_Relationship_EX_Count,Many_to_One_Relationship_EX,Many_to_One_Relationship_EX_Count) %>% distinct(.keep_all = TRUE)

Ens_XGC_One_to_One <-  Filtered_Final_Model[Filtered_Final_Model$One_to_Many_Relationship_EX_Count == 1 & 
Filtered_Final_Model$Many_to_One_Relationship_EX_Count == 1,] %>% ungroup() %>% select(chrloc_1.y,start1.y,end1.y,strand1.y,chrloc_2,start2,end2,strand2,One_to_Many_Relationship_EX,
One_to_Many_Relationship_EX_Count,Many_to_One_Relationship_EX,Many_to_One_Relationship_EX_Count) %>% distinct(.keep_all = TRUE)

########################################################################################


NCBI_Ens_XGC_Stats <- data.frame(type = c("Overlaped_NCBI_Ens_XGC","one_to_one_to_one","complicated_N_E_X",
"NCBI_Ens_One_to_One","NCBI_Ens_One_to_Many","NCBI_Ens_Many_to_One","NCBI_Ens_Many_to_Many",
"NCBI_XGC_One_to_One","NCBI_XGC_One_to_Many","NCBI_XGC_Many_to_One","NCBI_XGC_Many_to_Many",
"Ens_XGC_One_to_One","Ens_XGC_One_to_Many","Ens_XGC_Many_to_One","Ens_XGC_Many_to_Many"),
Count = c(

  nrow(filtered_1_to_1_to_1)+nrow(filtered_complicated),
  nrow(filtered_1_to_1_to_1),
  nrow(filtered_complicated),
  nrow(NCBI_Ens_One_to_One),
  nrow(NCBI_Ens_One_to_Many),
  nrow(NCBI_Ens_Many_to_One),
  nrow(NCBI_Ens_Many_to_Many),
  nrow(NCBI_XGC_One_to_One),
  nrow(NCBI_XGC_One_to_Many),
  nrow(NCBI_XGC_Many_to_One),
  nrow(NCBI_XGC_Many_to_Many),
  nrow(Ens_XGC_One_to_One),
  nrow(Ens_XGC_One_to_Many),
  nrow(Ens_XGC_Many_to_One),
  nrow(Ens_XGC_Many_to_Many)
  
))

return(list(NCBI_Ens_XGC_Stats,filtered_complicated,NCBI_Ens_One_to_Many,NCBI_Ens_Many_to_One,NCBI_Ens_Many_to_Many,NCBI_XGC_One_to_Many,NCBI_XGC_Many_to_One,NCBI_XGC_Many_to_Many,Ens_XGC_One_to_Many,Ens_XGC_Many_to_One,Ens_XGC_Many_to_Many))


}


NCBI_Ens_XGC_Stats <- count_NCBI_Ens_XGC(cbind_NCBI_Ens_XGC_ls)

write.csv(NCBI_Ens_XGC_Stats,file = "NCBI_Ens_XGC_Stats.csv",quote = FALSE,row.names = FALSE)


###########################################################################

############### Creation of Unified Gene Model ###########################

###########################################################################

###with the data analysis output of NCBI and Ensembl Output

#install.packages("data.table")
#library(data.table)

#part_gene_N <- part_gene_df %>% filter(type == "gene" | type == "mRNA")

#Unified_ls <- list()

for (j in 1 : length(unique(Final_grped_df$gp_mod_attr))) {
  
  part_gene <- Final_grped_df %>% filter(gp_mod_attr == unique(Final_grped_df$gp_mod_attr)[j])
  part_gene_id <- part_gene %>% filter(type == "gene")
  part_mRNA_id <- part_gene %>% filter(type == "mRNA")
  Non_gene_id <-  part_gene %>% filter(type != "gene" & type != "mRNA")
  
  #Mod_D <- rbind.data.frame(   part_gene_id[which.max(part_gene_id$Overlap_Value),],
   #                            part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),],
    #                              Non_gene_id)
  
    
  ## overlapping of gene id data with other exon and mrna features and eliminating the ones that do not overlap
  
  part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% select(chrloc,strand,start,end)
  Non_gene_id_df <- Non_gene_id[,c("chrloc","strand","start","end")]

  setDT(part_gene_id_df)
  setDT(Non_gene_id_df)
  setkey(Non_gene_id_df)
  
  overlapped_data <- foverlaps(part_gene_id_df,Non_gene_id_df,type = "any")
  
  selected_non_gene_ids <- Non_gene_id %>% 
    left_join(overlapped_data,by = c("chrloc","strand","start","end")) %>% 
    select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
    distinct(.keep_all = TRUE)

  #Unified_d <- Non_gene_id %>% 
  #left_join(overlapped_data,by = c("chrloc","strand","start","end")) %>% 
  #select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
  #distinct(.keep_all = TRUE)
  
  part_gene_id_f <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
    select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
    distinct(.keep_all = TRUE)
  
  part_mRNA_id_f <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>% 
    select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
    distinct(.keep_all = TRUE)
  
  Mod_df <- rbind.data.frame( part_gene_id_f,part_mRNA_id_f,selected_non_gene_ids)
 
 #Mod_ls[[j]] <- rbind.data.frame( 
  # part_gene_id[which.max(part_gene_id$Overlap_Value),],
   #part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),],
   #selected_non_gene_ids)
  
 #parent_id <- gsub("ID=","",str_split(str_extract(Unified_ls[[j]][Unified_ls[[j]]$type == "mRNA","attributes"],pattern = "ID=.+;"),pattern = ";")[[1]][1])

 parent_id <- gsub("ID=","",str_split(str_extract(Mod_df[Mod_df$type == "mRNA","attributes"],pattern = "ID=.+;"),pattern = ";")[[1]][1])
 
 #gene_id <- gsub("ID=","",str_split(str_extract(Unified_ls[[j]][Unified_ls[[j]]$type == "gene","attributes"],pattern = "ID=.+;"),pattern = ";")[[1]][1])

 #gene_id <- gsub("ID=","",str_split(str_extract(Mod_df[Mod_df$type == "gene","attributes"],pattern = "ID=.+;"),pattern = ";")[[1]][1])
 
 Non_Mod_df <-  Mod_df %>% 
 filter(type != "gene" & type != "mRNA")   
  
 Non_Mod_df <- Non_Mod_df %>% mutate(
   
   mod_attribute = gsub(str_split(str_extract(attributes,pattern = "Parent=.+;"),
   pattern = ";")[[1]][1],paste0("Parent=",parent_id),attributes),
   attr = if_else(source == "Genbank",
   gsub("Parent=",paste0("Parent=",parent_id,";Orgin="),attributes),
   mod_attribute)
   
 ) %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes = attr) %>% distinct(.keep_all = TRUE) 
  
 Unified_ls[[j]] <- rbind.data.frame(part_gene_id_f,part_mRNA_id_f,Non_Mod_df) %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
   
}

#Unified_Df <- do.call("rbind.data.frame",Unified_ls)

#Final_Unified_Df <- Unified_Df %>% 
#select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
#distinct(.keep_all = TRUE)

#Final_Unified_Df$source <- "Xenbase"

#write.table(Final_Unified_Df,file = "NCBI_Ensembl_XGC_Unified_Data_F.gff3",
#sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

######################################################################

##### An alternative strategy to create unified gene model catalogue:
####  aggregating exons across models by creating a customized gene model
####  with transcripts.Then,all the other features(Children) are further
####  aggregated under the respective customary gene and transcript.

######################################################################

Final_grped_Data <- read.delim("Final_grped_Data.gff3",sep = "\t",header = FALSE)

names(Final_grped_Data) <- c("chrloc","source","type","start","end","score","strand","phase","modified_attributes","attributes","gp_mod_attr")

Final_grped_df <- Final_grped_Data %>% mutate(
  
  
  Overlap_Value = end - start
  
  
)

Unified_M_ls <- list()
Unified_ls <- list()
Unified_F_ls <- list()

rank_data <- data.frame(type = c("gene","ncRNA_gene","pseudogene","mRNA","tRNA","rRNA","lnc_RNA","snoRNA","snRNA","scRNA",
"miRNA","Y_RNA","ncRNA","guide_RNA","transcript","primary_transcript","pseudogenic_transcript","exon","CDS",
"V_gene_segment","C_gene_segment","match","origin_of_replication","D_loop"),rank = 1:24)

########## correcting some entries while choosing the built models

Final_grped_df[(Final_grped_df$type == "CDS" & Final_grped_df$modified_attributes == "post11"),]$modified_attributes <- "LOC105945660"

#### 2nd entry ####

Correct_d <- Final_grped_df[Final_grped_df$source == "Genbank" & Final_grped_df$gp_mod_attr == "LOC101732650",]

## adding gene model

gene_d <- Correct_d[1,]

gene_d$type <- "gene"

gene_d$start <- min(Correct_d$start)

gene_d$end <- max(Correct_d$end)

gene_d$attributes <- paste0("ID=XENT-LOC101732650;locus_tag=XENTR_v10010934;gbkey=gene")
  
ID_p <- gsub("ID=","",str_split(str_extract(gene_d$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])

Correct_d$attributes <- sapply(Correct_d$attributes,function(x)
  
  gsub("Parent=",paste0("Parent=",ID_p,";Origin="),x)
  
  )

gene_d$Overlap_Value <- gene_d$end - gene_d$start

#gene_d$rank <- 1

Correct_df <- rbind.data.frame(gene_d,Correct_d)

Final_grped_df <- rbind.data.frame(Final_grped_df,gene_d)

abg4 <- Final_grped_df %>% filter(gp_mod_attr == "LOC101732650")

###########################

#### 3rd entry repeat!!

###########################

#for (i in 1 : length(unique(Final_grped_df$chrloc))){
  
#part_gene_chr <- Final_grped_df %>% filter(chrloc %in% Final_grped_df$chrloc[i])

#for (j in 1 : length(unique(Final_grped_df$gp_mod_attr))) {

part_gene <- part_gene_chr %>% group_by(chrloc,strand,gp_mod_attr) %>% filter(gp_mod_attr == unique(Final_grped_df$gp_mod_attr)[j])

if(((nrow(part_gene[part_gene$type == "gene",]) >= 2) |
   (nrow(part_gene[part_gene$type == "pseudogene",]) >= 2) |
   (nrow(part_gene[part_gene$type == "ncRNA_gene",]) >= 2)) & 
   (length(unique(part_gene[part_gene$type == "gene" | 
    part_gene$type == "pseudogene" | 
    part_gene$type == "ncRNA_gene",]$type) %in% 
    c("gene","pseudogene","ncRNA_gene")) == 1 )){

  part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | 
            type == "ncRNA_gene")
  
  part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
  type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
  type == "snRNA" | type == "scRNA" | type == "miRNA" |  
  type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
  type == "transcript" | type == "pseudogenic_transcript" | 
  type == "primary_transcript")
  
  Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
    select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
    distinct(.keep_all = TRUE)
  
  #### finding the start and end positions for customary gene/transcript model
  
  part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
   select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
    distinct(.keep_all = TRUE)
  
  part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
    select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
    distinct(.keep_all = TRUE)
  
  part_gene_id_df$start <- min(part_gene$start)
  
  part_gene_id_df$end <- max(part_gene$end)
  
  part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x) 
    
    gsub("ID=","gID=",x)
    
  )
  
  if(part_gene_id_df$type == "pseudogene"){
    
    part_gene_id_df$attributes <- paste0("ID=pseudogene-",paste0(unique(part_gene$modified_attributes),
     collapse = "_"),";",paste0(unique(part_gene_id$attributes),collapse = ";"))
    
  }
  
  if(part_gene_id_df$type == "ncRNA_gene"){
    
    part_gene_id_df$attributes <- paste0("ID=ncrnagene-",paste0(unique(part_gene$modified_attributes),
    collapse = "_"),";",paste0(unique(part_gene_id$attributes),collapse = ";")) 
    
    
  }
  
  if(part_gene_id_df$type == "gene"){
    
  part_gene_id_df$attributes <- paste0("ID=gene-",paste0(unique(part_gene$modified_attributes),
  collapse = "_"),";",paste0(unique(part_gene_id$attributes),collapse = ";"))    
  
  }
  
  part_gene_id_df$modified_attributes <- paste0(unique(part_gene_id$modified_attributes),collapse = ";")    
 
  
  if(nrow(part_mRNA_id_df) == 1){
    
    part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
      
      gsub("ID=","transID=",x)
      
    )
    
    part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
      
      gsub("Parent=","Origin=",x)
      
    )
    
    if(part_gene_id_df$type == "gene"){
    
    part_mRNA_id_df$attributes <- paste0("ID=rna-",paste0(unique(part_gene$modified_attributes),
    collapse = "_"),";","Parent=gene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";",
    paste0(unique(part_mRNA_id$attributes),collapse = ";"))
    
    }
    
    if(part_gene_id_df$type == "pseudogene"){
      
      part_mRNA_id_df$attributes <- paste0("ID=rna-p-",paste0(unique(part_gene$modified_attributes),
       collapse = "_"),";","Parent=pseudogene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";",
        paste0(unique(part_mRNA_id$attributes),collapse = ";"))
      
    }
    
    if(part_gene_id_df$type == "ncRNA_gene"){
      
      part_mRNA_id_df$attributes <- paste0("ID=rna-nc-",paste0(unique(part_gene$modified_attributes),
        collapse = "_"),";","Parent=ncrnagene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";",
          paste0(unique(part_mRNA_id$attributes),collapse = ";"))
      
    }
    
    part_mRNA_id_df$start <- min(part_gene$start)
    
    part_mRNA_id_df$end <- max(part_gene$end)
    
    part_mRNA_id_df$modified_attributes <- paste0(unique(part_mRNA_id$modified_attributes),collapse = ";")    
    
    parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
    
  }
  
  else if(nrow(part_mRNA_id_df) < 1){
  
    parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
    
  }
  

  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x) gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
  
  Unified_M_ls[[j]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
    select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
    distinct(.keep_all = TRUE) 

  }
  
 else if((nrow(part_gene[part_gene$type == "gene" | 
  part_gene$type == "pseudogene" | 
  part_gene$type == "ncRNA_gene",]) < 2) | 
  (nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
   nrow(part_gene[(part_gene$type == "pseudogene" | 
   part_gene$type == "ncRNA_gene"),]) >= 1)){
     
  if(nrow(part_gene[part_gene$type == "gene" | 
          part_gene$type == "pseudogene" | 
          part_gene$type == "ncRNA_gene",]) < 2){
    
    part_gene_df <- part_gene %>% filter(type == "gene"| type == "pseudogene" | type == "ncRNA_gene") %>%
      select(chrloc,source,type,start,end,score,phase,strand,modified_attributes,attributes) 
    
    if(nrow(part_gene_df) == 1){
    
    #part_gene_df <- part_gene %>% filter(type == "gene"| type == "pseudogene" | type == "ncRNA_gene") %>%
     # select(chrloc,source,type,start,end,score,phase,strand,modified_attributes,attributes) 
  
    part_trans_M <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
    type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
    type == "snRNA" | type == "scRNA" | type == "miRNA" |  
    type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
    type == "transcript" | type == "pseudogenic_transcript" | 
    type == "primary_transcript")
    
    part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
      
      chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
      
    ) %>% distinct(.keep_all = TRUE)
    
    if(length((part_gene_df$type == "pseudogene" | part_gene_df$type == "ncRNA_gene")) == 1){
      
     if(nrow(part_trans_M_df) >= 1){
        
       if(part_gene_df$type == "pseudogene"){
        
        #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 

        part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";gID="),x))
        
        part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=pseudogene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Origin="),x))
        
        part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";transID="),x))
        
          
        }
        
        else if(part_gene_df$type == "ncRNA_gene"){
          
          #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
          
          part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";gID="),x))
          
          part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=ncrnagene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Origin="),x))
          
          part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-nc-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";transID="),x))
          
        }
       
        
        part_gene_df$start <- min(part_gene$start)
        
        part_gene_df$end <- max(part_gene$end)
         
        part_trans_M_df$start <- part_gene_df$start
       
        part_trans_M_df$end <- part_gene_df$end 
        
        Non_gene_id <- part_gene  %>% filter(!(type %in% part_gene_df$type | 
                       type %in% part_trans_M$type ))  %>% select(
                                                     
        chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                     
        ) %>% distinct(.keep_all = TRUE)
        
        parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
        
        Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
          gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
        
        part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
          select(chrloc,source,type,start,end,score,                                                                           
          strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)
        
      }
      
      else if(nrow(part_trans_M) < 1){
        
        if(part_gene_df$type == "pseudogene"){
          
          #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
          
          part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";gId="),x))
          
          #parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";"))
        
          }
        
        else if(part_gene_df$type == "ncRNA_gene"){
          
          #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
          
          part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";gId="),x))
          
          #parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";"))
          
          
        }
        
        Non_gene_id <- part_gene  %>% filter(!(type == part_gene_df$type))  %>% select(
          
          chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
          
        ) %>% distinct(.keep_all = TRUE)
        
        part_gene_df$start <- min(part_gene$start)
        
        part_gene_df$end <- max(part_gene$end)
        
        parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
        
        Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
          gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
        
        part_gene <- rbind.data.frame(part_gene_df,Non_gene_id) %>% 
          select(chrloc,source,type,start,end,score,                                                                           
          strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)
        
      }
      
    }
    
    if(length(part_gene_df$type == "gene") == 1){
      
       if(nrow(part_trans_M_df) < 1){
       
      #ab <- part_gene_df[part_gene_df$type == "gene",]$attributes 
      
      #orginal_gene_id <- str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]
      
      part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) 
        
        gsub("ID=",paste0("ID=gene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";gID="),
        x)
      )
        
      parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
      
      Non_gene_id <- part_gene  %>% filter(!(type == part_gene_df$type))  %>% select(
                                                 
       chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                 
       ) %>% distinct(.keep_all = TRUE)
      
      part_gene_df$start <- min(part_gene$start)
      
      part_gene_df$end <- max(part_gene$end)
      
      Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
        gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
      
      part_gene <- rbind.data.frame(part_gene_df,Non_gene_id) %>% 
       select(chrloc,source,type,start,end,score,                                                                           
       strand,phase,modified_attributes,attributes) %>% 
        distinct(.keep_all = TRUE)
      
       }
       
       else if(nrow(part_trans_M_df) >= 1){
       
         part_gene_df$attributes <- sapply(part_gene_df$attributes, function(x)
           
           gsub("ID=",paste0("ID=gene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";gId="),x)
           
         )
         
         
         part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
           
           gsub("Parent=",paste0("Parent=gene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Origin="),x)
           
         )
         
         part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
           
           gsub("ID=",paste0("ID=rna-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Origin="),x)
           
         )
         
           
        parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,
        pattern = "ID=.+;"),pattern = ";")[[1]][1])
        
        part_gene_df$start <- min(part_gene$start)
        
        part_gene_df$end <- max(part_gene$end)
        
        part_trans_M_df$start <- part_gene_df$start
        
        part_trans_M_df$end <- part_gene_df$end
        
        Non_gene_id <- part_gene  %>% filter(!(type %in% part_gene_df$type | type %in% part_trans_M$type))  %>% select(
          
          chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
          
        ) %>% distinct(.keep_all = TRUE)
        
        Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
          gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
        
        part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
          select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)
        
       }
    }
    
    }
    
    else if(nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) < 1 ){
      
      if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) < 1 )){
      
      org_parent_id <- gsub("Parent=","",str_split(str_extract(part_gene$attributes[1],pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
      
      gene_id <- paste0("ID=gene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";","OrigingID:",org_parent_id,";","gbkey=gene")
      
      parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
      
      mRNA_id <- paste0("ID=rna-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Parent=",parent_id,";")
      
      pid <- gsub("ID=","",str_split(str_extract(mRNA_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
        
      part_gene$attributes <- sapply(part_gene$attributes,function(x) 
        
        gsub("Parent=",paste0("Parent=",pid,";Origin="),x)
        
      ) 
      
      mod_gene_id <- part_gene[1,]
      
      mod_mRNA_id <- part_gene[1,]
      
      mod_gene_id$type <- "gene"
      mod_mRNA_id$type <- "mRNA"
      
      mod_gene_id$start <- min(part_gene$start)
      mod_gene_id$end <- max(part_gene$end)
      
      mod_mRNA_id$start <- min(part_gene$start)
      mod_mRNA_id$end <- max(part_gene$end)
      
      mod_gene_id$attributes <- gene_id
      mod_mRNA_id$attributes <- mRNA_id
      
      part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,part_gene) %>% select(chrloc,source,type,start,end,score,
      strand,phase,modified_attributes,attributes) %>% 
        distinct(.keep_all = TRUE)
      
      }
      
      else if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) >= 1 )){
        
        mod_mRNA_id <- part_gene[which.max(part_gene$type == "mRNA"),]
        
        org_parent_id <- gsub("Parent=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
        
        gene_id <- paste0("ID=gene-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";","OrigingID:",org_parent_id,";","gbkey=gene")
        
        parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
        
        #mRNA_id <- paste0("ID=rna-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Parent=",parent_id,";")
        
        mod_gene_id <- mod_mRNA_id       
         
        mod_gene_id$type <- "gene"
        mod_mRNA_id$type <- "mRNA"
        
        mod_gene_id$start <- min(part_gene$start)
        mod_gene_id$end <- max(part_gene$end)
        
        mod_mRNA_id$start <- min(part_gene$start)
        mod_mRNA_id$end <- max(part_gene$end)
        
        mod_gene_id$attributes <- gene_id
        
        mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
        mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("ID=",paste0("ID=rna-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Origin="),x))
          
        p_id <- gsub("ID=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "ID.+;"),pattern = ";")[[1]][1])
          
        Non_gene_id <- part_gene %>% filter(!(type %in% mod_gene_id$type | type %in% mod_mRNA_id$type))
        
        Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
          
          gsub("Parent=",paste0("Parent=",p_id,";Origin="),x)
          
          )
        
        part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,Non_gene_id) %>% select(chrloc,source,type,start,end,score,
         strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)        
        
      }
      
      if(unique(part_gene$source) == "NCBI"){
        
        part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,
         strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)        
        
        
      }
    
    }
    
    Unified_M_ls[[j]] <- part_gene %>% select(chrloc,source,type,start,end,score,
    strand,phase,modified_attributes,attributes) %>% 
    distinct(.keep_all = TRUE)
    
}
   
 else if(nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
         nrow(part_gene[(part_gene$type == "pseudogene" | 
          part_gene$type == "ncRNA_gene"),]) >= 1){
    
  part_gene_No <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene") %>% 
   group_by(type) %>% summarise( Count = n(), source_t = paste0(source,collapse = ";"),mod_attr = paste0(unique(modified_attributes),collapse = ";") )
     
     for(i in 1 : length(part_gene_No$Count)){
       
      if(part_gene_No$Count[i] == 1){
       
      source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
         
      part_gene_Mod_attr <- unique(str_split(part_gene_No$mod_attr[i],pattern = ";")[[1]])
        
      part_gene_type <- part_gene_No[i,]$type    
         
      #if(length(unique(part_gene_No$mod_attr)) == 1){
       
      part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
       
       #}
       
       #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene_No$mod_attr)) == nrow(part_gene_No))){
         
       #part_gene_M <- part_gene %>% filter( (source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
       
       #}
       
       #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene_No$mod_attr)) != nrow(part_gene_No))){
         
        #part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
         
       #}
       
     part_trans_M <- part_gene_M %>% filter(type == "mRNA" | type == "tRNA" | 
     type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | 
     type == "snRNA" | type == "scRNA" | type == "miRNA" |  
     type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
     type == "transcript" | type == "pseudogenic_transcript" | 
     type == "primary_transcript")
       
     part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
         
      chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
         
       ) %>% distinct(.keep_all = TRUE)
       
       part_gene_M_df <- part_gene_M %>% filter(type == part_gene_type) %>% select(
         
         chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
         
       ) %>% distinct(.keep_all = TRUE)
       
       
    if(part_gene_type == "pseudogene" | part_gene_type == "ncRNA_gene"){
         
         if(nrow(part_trans_M_df) == 1){
           
           if(part_gene_type == "pseudogene"){
             
             #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
             
             part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";gID="),x))
             
             part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=pseudogene-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";Origin="),x))
             
             part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";transID="),x))
             
             parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
             
           }
           
           else if(part_gene_type == "ncRNA_gene"){
             
             #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
             
             part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-",paste0(unique(part_gene_M$modified_attributes)),";gID="),x))
             
             part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=ncrnagene-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";Origin="),x))
             
             part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-ncg-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";transID="),x))
             
             parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
             
             
           }
           
           part_gene_M_df$start <- min(part_gene_M$start)
           part_gene_M_df$end <- max(part_gene_M$end)
           part_trans_M_df$start <- min(part_gene_M$start)
           part_trans_M_df$end <- max(part_gene_M$end)
           
           Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | 
           type %in% part_trans_M$type ))  %>% select(
                                                        
            chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                        
            ) %>% distinct(.keep_all = TRUE)
           
           
           Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
             gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
           
           part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
            select(chrloc,source,type,start,end,score,                                                                           
            strand,phase,modified_attributes,attributes) %>% 
            distinct(.keep_all = TRUE)
           
         }
         
         else if(nrow(part_trans_M_df) != 1){
           
           if(part_gene_type == "pseudogene"){
             
             #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
             
             part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";gID="),x))
             
           }
           
           else if(part_gene_type == "ncRNA_gene"){
             
             #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
             
             part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";gID="),x))             
             
           }
           
           
           part_gene_M_df$start <- min(part_gene_M$start)
           part_gene_M_df$end <- max(part_gene_M$end)

           
           Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type))  %>% select(
             
             chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
             
           ) %>% distinct(.keep_all = TRUE)
           
           parent_id_M <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
           
           Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
             gsub("Parent=",paste0("Parent=",parent_id_M,";Origin="),x))
           
           part_gene_M <- rbind.data.frame(part_gene_M_df,Non_gene_id) %>% 
             select(chrloc,source,type,start,end,score,                                                                           
             strand,phase,modified_attributes,attributes) %>% 
             distinct(.keep_all = TRUE)
           
         }
         
       }
       
       if(part_gene_type == "gene"){
        
         if(nrow(part_trans_M_df) < 1){
           
           part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=gene-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";gID="),x))
             
           parent_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,
            pattern = "ID=.+;"),pattern = ";")[[1]][1])
           
           part_gene_M_df$start <- min(part_gene_M$start)
           part_gene_M_df$end <- max(part_gene_M$end)
           
           Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type))  %>% select(
             
             chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
             
           ) %>% distinct(.keep_all = TRUE)
           
           Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
             gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
           
           part_gene_M <- rbind.data.frame(part_gene_M_df,Non_gene_id) %>% 
           select(chrloc,source,type,start,end,score,                                                                           
            strand,phase,modified_attributes,attributes) %>% 
             distinct(.keep_all = TRUE)
           
         }
         
         else if(nrow(part_trans_M_df) >= 1){
           
           part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) 
            gsub("ID=",paste0("ID=gene-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";gID="),x))
      
           part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
             
             gsub("Parent=",paste0("Parent=gene-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";Origin="),x)
             
           )
            
           part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
             
             gsub("ID=",paste0("ID=rna-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),";transID="),x)
             
           )    
      
           
           parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,
                        pattern = "ID=.+;"),pattern = ";")[[1]][1])
           
           part_gene_M_df$start <- min(part_gene_M$start)
           part_gene_M_df$end <- max(part_gene_M$end)
           part_trans_M_df$start <- min(part_gene_M$start)
           part_trans_M_df$end <- max(part_gene_M$end)
           
           Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | type %in% part_trans_M_df$type))  %>% select(
             
             chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
             
           ) %>% distinct(.keep_all = TRUE)
           
           Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
             gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
           
           part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
             select(chrloc,source,type,start,end,score,                                                                           
             strand,phase,modified_attributes,attributes) %>% 
             distinct(.keep_all = TRUE)
           
         }
       }  
       
        Unified_ls[[i]] <- part_gene_M %>% select(chrloc,source,type,start,end,score,
         strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)
      
       }   
  
      else if(part_gene_No$Count[i] > 1){
         
       source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
        
       part_gene_type <- part_gene_No[i,]$type    
        
       part_gene_Mod_attr <- unique(part_gene[
       part_gene$type == part_gene_type,]$modified_attributes)
        
       #if(length(unique(part_gene_No$mod_attr)) == 1){
         
         part_gene_Mod <- part_gene %>% filter((source %in% source_tot ) & (modified_attributes %in% part_gene_Mod_attr))
         
       #}
       
       #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene[(part_gene$source %in% source_tot) & (part_gene$modified_attributes %in% part_gene_Mod_attr),]$source)) == length(unique(source_tot)))){
         
         #part_gene_Mod <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
         
       #}
       
       #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene[(part_gene$source %in% source_tot) & (part_gene$modified_attributes %in% part_gene_Mod_attr) ,]$source)) == length(unique(source_tot)))){
         
         #part_gene_Mod <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
         
       #}
       
       part_gene_id <- part_gene_Mod %>% filter(type == part_gene_type)
      
       part_mRNA_id <- part_gene_Mod %>% filter(type == "mRNA" | type == "tRNA"| 
       type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
       type == "snRNA" | type == "scRNA" | type == "miRNA" |  
       type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
       type == "transcript" | type == "pseudogenic_transcript" | 
       type == "primary_transcript")
      
       Non_gene_id <-  part_gene_Mod %>% filter(!(type == part_gene_type | type %in% part_mRNA_id$type)) %>% 
         select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
       
      #### finding the start and end positions for customary gene/transcript model
      
     part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
     select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
      
     part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
     select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
      
     part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x) 
        
        gsub("ID=","gID=",x)
        
      )
      
     
      part_gene_id_df$start <- min(part_gene_Mod$start)
      
      part_gene_id_df$end <- max(part_gene_Mod$end)
      
      if(part_gene_type == "pseudogene"){
      
      part_gene_id_df$attributes <- 
      paste0("ID=pseudogene-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),";",
      paste0(unique(part_gene_id$attributes),collapse = ";"))    
      
      
      }
      
      if(part_gene_type == "ncRNA_gene"){
        
        part_gene_id_df$attributes <- 
          paste0("ID=ncrnagene-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),";",
                 paste0(unique(part_gene_id$attributes),collapse = ";"))    
        
        
      }
      
      if(part_gene_type == "gene"){
        
        part_gene_id_df$attributes <- 
          paste0("ID=gene-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),";",
                 paste0(unique(part_gene_id$attributes),collapse = ";"))    
        
        
      }
      
      
    part_gene_id_df$modified_attributes <- paste0(unique(part_gene_id$modified_attributes),collapse = ";")    
      
    if(nrow(part_mRNA_id_df) >= 1){
      
      part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
          
        gsub("ID=","transID=",x)
          
      )
        
      part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
          
        gsub("Parent=","Origin=",x)
          
      )
        
        part_mRNA_id_type <- part_mRNA_id_df$type
        
        part_mRNA_id_df$start <- min(part_gene_Mod$start)
        
        part_mRNA_id_df$end <- max(part_gene_Mod$end)
        
        if(part_gene_type == "pseudogene"){
        
        part_mRNA_id_df$attributes <- paste0("ID=rna-p-",paste0(unique(part_gene_Mod$modified_attributes),
         collapse = "_"),";","Parent=pseudogene-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),
          ";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
        
        }
        
        if(part_gene_type == "ncRNA_gene"){
          
          part_mRNA_id_df$attributes <- paste0("ID=rna-nc-",paste0(unique(part_gene_Mod$modified_attributes),
          collapse = "_"),";","Parent=ncrnagene-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),
          ";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
          
        }
        
        if(part_gene_type == "gene"){
          
          part_mRNA_id_df$attributes <- paste0("ID=rna-",paste0(unique(part_gene_Mod$modified_attributes),
           collapse = "_"),";","Parent=gene-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),
           ";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
          
        }
        
        part_mRNA_id_df$modified_attributes <- paste0(unique(part_mRNA_id$modified_attributes),collapse = ";")    
        
        parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
      
      }
      
      else if(nrow(part_mRNA_id_df) != 1){
        
        parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),
                     pattern = ";")[[1]][1])
        
      }
      
      Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
        
        
        gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x)
        
        ) 
      
    
      
      Unified_ls[[i]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>%
      select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
       distinct(.keep_all = TRUE)
      
         }
       
      }
     
     Unified_M_ls[[j]] <- do.call("rbind.data.frame",Unified_ls) 
     
    } 
  
  }

 else {
   
   Unified_M_ls[[j]] <- part_gene %>% select(chrloc,source,type,start,end,score,
   strand,phase,modified_attributes,attributes) %>% 
   distinct(.keep_all = TRUE)
   
   }
 
}

#}  
  
#Unified_M_Df <- do.call("rbind.data.frame",Unified_M_ls)


########### passing each chromosome gene models and further grouping them according to
### their grouping variable (XGC gp variable)/declaring a function and re-routing the 
### filtered models according to chromosome location accordingly

unified_NCBI_Ens_XGC <- function(NCBI_Ens_XGC_Data){
  
  rank_data <- data.frame(type = c("gene","ncRNA_gene","pseudogene","mRNA","tRNA","rRNA","lnc_RNA","snoRNA","snRNA","scRNA",
   "miRNA","Y_RNA","ncRNA","guide_RNA","transcript","primary_transcript","pseudogenic_transcript","exon","CDS",
    "V_gene_segment","C_gene_segment","match","origin_of_replication","D_loop"),rank = 1:24)
  
  for (j in 1 : length(unique(NCBI_Ens_XGC_Data$gp_mod_attr))) {
    
    part_gene_d <- NCBI_Ens_XGC_Data %>% filter(gp_mod_attr == unique(NCBI_Ens_XGC_Data$gp_mod_attr)[j])
    
    chr_d <- unique(part_gene_d$chrloc)  
    
    for (k in 1 : length(unique(part_gene_d$strand))) {
      
      part_gene <- part_gene_d %>% filter(
        
        strand == unique(part_gene_d$strand)[k]
        
      )
    
  if(unique(part_gene$strand) == "+"){
    
    if((unique(part_gene$gp_mod_attr) == "LOC101734477" | unique(part_gene$gp_mod_attr) == "LOC116406437" | unique(part_gene$gp_mod_attr) == "smyd5")){
      
      Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
       strand,phase,modified_attributes,attributes) %>% 
        distinct(.keep_all = TRUE)
      
    }
    
    else if(unique(part_gene$gp_mod_attr) == "U4" | 
    unique(part_gene$gp_mod_attr) == "xtr" | 
    unique(part_gene$gp_mod_attr) == "RNaseP_nuc" | 
    unique(part_gene$gp_mod_attr) == "U2" | 
    unique(part_gene$gp_mod_attr) == "XENTR_v10003163" |
    unique(part_gene$gp_mod_attr) == "trnag-ccc" |
    unique(part_gene$gp_mod_attr) == "trnar-ucu" |
    unique(part_gene$gp_mod_attr) == "5S_rRNA" |
    unique(part_gene$gp_mod_attr) == "trnap-agg" |
    unique(part_gene$gp_mod_attr) == "trnad-guc" |
    unique(part_gene$gp_mod_attr) == "trnav-aac" |
    unique(part_gene$gp_mod_attr) == "trnas-uga" |
    unique(part_gene$gp_mod_attr) == "trnat-ugu" |
    unique(part_gene$gp_mod_attr) == "trnaa-cgc"|
    unique(part_gene$gp_mod_attr) == "trnae-uuc" |
    unique(part_gene$gp_mod_attr) == "U7" |
    unique(part_gene$gp_mod_attr) == "trnaq-cug" |
    unique(part_gene$gp_mod_attr) == "trnas-aga" |
    unique(part_gene$gp_mod_attr) == "trnaq-uug" 
    ){
      
      part_df <- part_gene %>% left_join(rank_data,by = c("type"))
      
      if(length(unique(part_df[part_df$type == "tRNA",]$type)) >= 1){
        
        part_df <- part_df %>% arrange(rank)
        
        part_df$transcript <- sapply(part_df$attributes,function(x) 
          
          ifelse(str_detect(x,pattern = "Parent=.+"),
                 
          gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1]),
                 
          gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
          
        )
        
        part_df$attributes <- as.vector(part_df$attributes)
        
        part_df$transcript <- as.vector(part_df$transcript)
        
        part_df <- part_df %>% group_by(chrloc,strand,modified_attributes) %>% 
        arrange(transcript) %>%
        ungroup() %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
          distinct(.keep_all = TRUE)
        
        Unified_M_ls[[k]] <- part_df %>% select(chrloc,source,type,start,end,score,
         strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)
        
      }
      
      if(length(unique(part_df[part_df$type == "ncRNA_gene" | part_df$type == "miRNA",]$type)) >= 1){
        
        #part_df <- part_df %>% left_join(rank_data,by = c("type"))
        
        part_df <- part_df %>% arrange(rank)
      
        part_df$transcript <- sapply(part_df$attributes,function(x) 
        
          ifelse(str_detect(x,pattern = "Parent=.+"),
               
           gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1]),
               
           gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
        
        )
        
        part_df <- part_df %>% arrange(transcript) %>% 
          select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
        
        Unified_M_ls[[k]] <- part_df %>% select(chrloc,source,type,start,end,score,
         strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)
        
      } 
      
      else{
      
      Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
       strand,phase,modified_attributes,attributes) %>% 
       distinct(.keep_all = TRUE)
    
      }
    }
    
    else if(((nrow(part_gene[part_gene$type == "gene",]) >= 2) |
      (nrow(part_gene[part_gene$type == "pseudogene",]) >= 2) |
      (nrow(part_gene[part_gene$type == "ncRNA_gene",]) >= 2)) & 
       (length(unique(part_gene[part_gene$type == "gene" | 
        part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]$type) %in% 
        c("gene","pseudogene","ncRNA_gene")) == 1 )){
      
      part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
      
      part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
                                             type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                             type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                             type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                             type == "transcript" | type == "pseudogenic_transcript" | 
                                             type == "primary_transcript")
      
      Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
        select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
        distinct(.keep_all = TRUE)
      
      #### finding the start and end positions for customary gene/transcript model
      
      part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
        select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
        distinct(.keep_all = TRUE)
      
      part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
        select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
        distinct(.keep_all = TRUE)
      
      part_gene_id_df$start <- min(part_gene$start)
      
      part_gene_id_df$end <- max(part_gene$end)
      
      part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x) 
        
        gsub("ID=","gID=",x)
        
      )
      
      if(part_gene_id_df$type == "pseudogene"){
        
        part_gene_id_df$attributes <- paste0("ID=pseudogene-p-",paste0(unique(part_gene$modified_attributes),
                                                                     collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_gene_id$attributes),collapse = ";"))
        
      }
      
      if(part_gene_id_df$type == "ncRNA_gene"){
        
        part_gene_id_df$attributes <- paste0("ID=ncrnagene-p-",paste0(unique(part_gene$modified_attributes),
        collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_gene_id$attributes),collapse = ";")) 
        
        
      }
      
      if(part_gene_id_df$type == "gene"){
        
        part_gene_id_df$attributes <- paste0("ID=gene-p-",paste0(unique(part_gene$modified_attributes),
                                                               collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_gene_id$attributes),collapse = ";"))    
        
      }
      
      part_gene_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
      
      
      if(nrow(part_mRNA_id_df) == 1){
        
        part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
          
          gsub("ID=","transID=",x)
          
        )
        
        part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
          
          gsub("Parent=","Origin=",x)
          
        )
        
        if(part_gene_id_df$type == "gene"){
          
          part_mRNA_id_df$attributes <- paste0("ID=rna-gp-",paste0(unique(part_gene$modified_attributes),
            collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=gene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
             paste0(unique(part_mRNA_id$attributes),collapse = ";"))
          
        }
        
        if(part_gene_id_df$type == "pseudogene"){
          
          part_mRNA_id_df$attributes <- paste0("ID=rna-pseudo-p-",paste0(unique(part_gene$modified_attributes),
          collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=pseudogene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
          paste0(unique(part_mRNA_id$attributes),collapse = ";"))
          
        }
        
        if(part_gene_id_df$type == "ncRNA_gene"){
          
         part_mRNA_id_df$attributes <- paste0("ID=rna-nc-p-",paste0(unique(part_gene$modified_attributes),
         collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=ncrnagene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
         paste0(unique(part_mRNA_id$attributes),collapse = ";"))
          
        }
        
        part_mRNA_id_df$start <- min(part_gene$start)
        
        part_mRNA_id_df$end <- max(part_gene$end)
        
        part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
        
        parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
        
      }
      
      else if(nrow(part_mRNA_id_df) < 1){
        
        parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
        
      }
      
      
      Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x) gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
      
      if(nrow(Non_gene_id) >= 1){
        
        
        Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
        
        
      }
      
      Unified_M_ls[[k]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
        select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
        distinct(.keep_all = TRUE) 
      
    }
    
    else if((nrow(part_gene[part_gene$type == "gene" | 
             part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) < 2) | 
            (nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
             (nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1 |
             nrow(part_gene[(part_gene$type == "ncRNA_gene"),]) >= 1 ))){
      
      if(nrow(part_gene[part_gene$type == "gene" | 
                        part_gene$type == "pseudogene" |
                        part_gene$type == "ncRNA_gene",]) < 2){
        
        part_gene_df <- part_gene %>% filter(type == "gene"| type == "pseudogene"| type == "ncRNA_gene") %>%
          select(chrloc,source,type,start,end,score,phase,strand,modified_attributes,attributes) 
        
        if(nrow(part_gene_df) == 1){
          
          #part_gene_df <- part_gene %>% filter(type == "gene"| type == "pseudogene" | type == "ncRNA_gene") %>%
          # select(chrloc,source,type,start,end,score,phase,strand,modified_attributes,attributes) 
          
          part_trans_M <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
                                                 type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                 type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                 type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                 type == "transcript" | type == "pseudogenic_transcript" | 
                                                 type == "primary_transcript")
          
          part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
            
            chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
            
          ) %>% distinct(.keep_all = TRUE)
          
          if(length((part_gene_df$type == "pseudogene")) == 1 | length((part_gene_df$type == "ncRNA_gene")) == 1){
            
            if(nrow(part_trans_M_df) >= 1){
              
              if(part_gene_df$type == "pseudogene"){
                
                #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                
                part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gID="),x))
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=pseudogene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";Origin="),x))
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-pseudo-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";transID="),x))
                
                
              }
              
              else if(part_gene_df$type == "ncRNA_gene"){
                
                gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                
                part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gID="),x))
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=ncrnagene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";Origin="),x))
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-nc-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";transID="),x))
                
              }
              
              
              part_gene_df$start <- min(part_gene$start)
              
              part_gene_df$end <- max(part_gene$end)
              
              part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
              
              part_trans_M_df$start <- part_gene_df$start
              
              part_trans_M_df$end <- part_gene_df$end 
              
              part_trans_M_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
              
              Non_gene_id <- part_gene  %>% filter(!(type %in% part_gene_df$type | 
              type %in% part_trans_M$type ))  %>% select(
                                                         
              chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                         
              ) %>% distinct(.keep_all = TRUE)
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
              
              if(nrow(Non_gene_id) >= 1){
                
                
                Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
                
                
              }
              
              
              part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                strand,phase,modified_attributes,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            else if(nrow(part_trans_M) < 1){
              
              if(part_gene_df$type == "pseudogene"){
                
                #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                
                part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gId="),x))
                
                #parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";"))
                
              }
              
              else if(part_gene_df$type == "ncRNA_gene"){
                
                #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                
               part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gId="),x))
                
                #parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";"))
                
                
              }
              
              Non_gene_id <- part_gene  %>% filter(!(type == part_gene_df$type))  %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              ) %>% distinct(.keep_all = TRUE)
              
              part_gene_df$start <- min(part_gene$start)
              
              part_gene_df$end <- max(part_gene$end)
              
              part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
              
              if(nrow(Non_gene_id) >= 1){
                
                
                Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
                
                
              }
              
              part_gene <- rbind.data.frame(part_gene_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                strand,phase,modified_attributes,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
          }
          
          if(length(part_gene_df$type == "gene") == 1){
            
            if(nrow(part_trans_M_df) < 1){
              
              #ab <- part_gene_df[part_gene_df$type == "gene",]$attributes 
              
              #orginal_gene_id <- str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]
              
              part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) 
                
                gsub("ID=",paste0("ID=gene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gID="),
                     x)
              )
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
              
              Non_gene_id <- part_gene  %>% filter(!(type == part_gene_df$type))  %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              ) %>% distinct(.keep_all = TRUE)
              
              if(nrow(Non_gene_id) >= 1){
                
                
                Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
                
                
              }
              
              part_gene_df$start <- min(part_gene$start)
              
              part_gene_df$end <- max(part_gene$end)
              
              part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
              
              part_gene <- rbind.data.frame(part_gene_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                strand,phase,modified_attributes,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            else if(nrow(part_trans_M_df) >= 1){
              
              part_gene_df$attributes <- sapply(part_gene_df$attributes, function(x)
                
                gsub("ID=",paste0("ID=gene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gId="),x)
                
              )
              
              
              part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
                
                gsub("Parent=",paste0("Parent=gene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";Origin="),x)
                
              )
              
              part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
                
                gsub("ID=",paste0("ID=rna-gp-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";Origin="),x)
                
              )
              
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,
                                                               pattern = "ID=.+;"),pattern = ";")[[1]][1])
              
              part_gene_df$start <- min(part_gene$start)
              
              part_gene_df$end <- max(part_gene$end)
              
              part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
              
              part_trans_M_df$start <- part_gene_df$start
              
              part_trans_M_df$end <- part_gene_df$end
              
              part_trans_M_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
              
              Non_gene_id <- part_gene  %>% filter(!(type %in% part_gene_df$type | type %in% part_trans_M$type))  %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              ) %>% distinct(.keep_all = TRUE)
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
              
              if(nrow(Non_gene_id) >= 1){
                
                
                Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
                
                
              }
              
              part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
          }
          
        }
        
        else if(nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene",]) < 1 ){
          
          if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) < 1 )){
            
            org_parent_id <- gsub("Parent=","",str_split(str_extract(part_gene$attributes[1],pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
            
            mod_gene_id <- part_gene[1,]
            
            mod_mRNA_id <- part_gene[1,]
            
            mod_gene_id$type <- "gene"
            mod_mRNA_id$type <- "mRNA"
            
            mod_gene_id$start <- min(part_gene$start)
            mod_gene_id$end <- max(part_gene$end)
            
            mod_mRNA_id$start <- min(part_gene$start)
            mod_mRNA_id$end <- max(part_gene$end)
            
            gene_id <- paste0("ID=gene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",mod_gene_id$start,"_",mod_gene_id$end,";","OrigingID:",org_parent_id,";","gbkey=gene")
            
            parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
            
            mRNA_id <- paste0("ID=rna-gp-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",mod_gene_id$start,"_",mod_gene_id$end,";Parent=",parent_id,";")
            
            pid <- gsub("ID=","",str_split(str_extract(mRNA_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
            
            part_gene$attributes <- sapply(part_gene$attributes,function(x) 
              
              gsub("Parent=",paste0("Parent=",pid,";Origin="),x)
              
            )
            
            mod_gene_id$attributes <- gene_id
            mod_mRNA_id$attributes <- mRNA_id
            
            mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
            mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
            
            part_gene$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
            
            part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,part_gene) %>% select(chrloc,source,type,start,end,score,
                                                                                        strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          else if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) >= 1 )){
            
            mod_mRNA_id <- part_gene[which.max(part_gene$type == "mRNA"),]
            
            org_parent_id <- gsub("Parent=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
            
            #mRNA_id <- paste0("ID=rna-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Parent=",parent_id,";")
            
            mod_gene_id <- mod_mRNA_id       
            
            mod_gene_id$type <- "gene"
            mod_mRNA_id$type <- "mRNA"
            
            mod_gene_id$start <- min(part_gene$start)
            mod_gene_id$end <- max(part_gene$end)
            
            mod_mRNA_id$start <- min(part_gene$start)
            mod_mRNA_id$end <- max(part_gene$end)
            
            gene_id <- paste0("ID=gene-p-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",mod_gene_id$start,"_",mod_gene_id$end,";","OrigingID:",org_parent_id,";","gbkey=gene")
            
            parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
            
            mod_gene_id$attributes <- gene_id
            
            mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
            mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("ID=",paste0("ID=rna-gp-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",mod_gene_id$start,"_",mod_gene_id$end,";Origin="),x))
            
            p_id <- gsub("ID=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "ID.+;"),pattern = ";")[[1]][1])
            
            Non_gene_id <- part_gene %>% filter(!(type %in% mod_gene_id$type | type %in% mod_mRNA_id$type))
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
              
              gsub("Parent=",paste0("Parent=",p_id,";Origin="),x)
              
            )
            
            mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
            mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
            
            part_gene$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
            
            part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,Non_gene_id) %>% select(chrloc,source,type,start,end,score,
                                                                                          strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)        
            
          }
          
          if(unique(part_gene$source) == "NCBI"){
            
            part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,
             strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)        
            
            
          }
          
        }
        
        Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
          strand,phase,modified_attributes,attributes) %>% 
          distinct(.keep_all = TRUE)
        
      }
      
      else if(nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
             (nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1 |
              nrow(part_gene[(part_gene$type == "ncRNA_gene"),]) >= 1)){
        
        part_gene_No <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene") %>% 
          group_by(type) %>% summarise( Count = n(), source_t = paste0(source,collapse = ";"),mod_attr = paste0(unique(modified_attributes),collapse = ";") )
        
        for(i in 1 : length(part_gene_No$Count)){
          
          if(part_gene_No$Count[i] == 1){
            
            source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
            
            part_gene_Mod_attr <- unique(str_split(part_gene_No$mod_attr[i],pattern = ";")[[1]])
            
            part_gene_type <- part_gene_No[i,]$type    
            
            #if(length(unique(part_gene_No$mod_attr)) == 1){
            
            part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
            
            #}
            
            #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene_No$mod_attr)) == nrow(part_gene_No))){
            
            #part_gene_M <- part_gene %>% filter( (source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
            
            #}
            
            #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene_No$mod_attr)) != nrow(part_gene_No))){
            
            #part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
            
            #}
            
            part_trans_M <- part_gene_M %>% filter(type == "mRNA" | type == "tRNA" | 
                                                     type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | 
                                                     type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                     type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                     type == "transcript" | type == "pseudogenic_transcript" | 
                                                     type == "primary_transcript")
            
            part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df <- part_gene_M %>% filter(type == part_gene_type) %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes,Overlap_Value
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df <- part_gene_M_df[which.max(part_gene_M_df$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df$type <- part_gene_type
            
            if(part_gene_type == "pseudogene" | part_gene_type == "ncRNA_gene"){
              
              if(nrow(part_trans_M_df) == 1){
                
                part_gene_M_df$start <- min(part_gene_M$start)
                part_gene_M_df$end <- max(part_gene_M$end)
                part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                part_trans_M_df$start <- min(part_gene_M$start)
                part_trans_M_df$end <- max(part_gene_M$end)
                part_trans_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                
                if(part_gene_type == "pseudogene"){
                  
                  #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                  
                  part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=pseudogene-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";Origin="),x))
                  
                  part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-pseudo-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";transID="),x))
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                }
                
                else if(part_gene_type == "ncRNA_gene"){
                  
                  #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-p-",paste0(unique(part_gene_M$modified_attributes)),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                  
                  part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=ncrnagene-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";Origin="),x))
                  
                  part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-ncg-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";transID="),x))
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  
                }
                
                Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | 
                                                           type %in% part_trans_M$type ))  %>% select(
                                                             
                                                             chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                             
                                                           ) %>% distinct(.keep_all = TRUE)
                
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                
                if(nrow(Non_gene_id) >= 1){
                  
                  
                  Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                  
                  
                }
                
                part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                  strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
              else if(nrow(part_trans_M_df) != 1){
                
                part_gene_M_df$start <- min(part_gene_M$start)
                part_gene_M_df$end <- max(part_gene_M$end)
                part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                
                
                if(part_gene_type == "pseudogene"){
                  
                  #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                  
                }
                
                else if(part_gene_type == "ncRNA_gene"){
                  
                  #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))             
                  
                }
                
                Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type))  %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                  
                ) %>% distinct(.keep_all = TRUE)
                
                parent_id_M <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub("Parent=",paste0("Parent=",parent_id_M,";Origin="),x))
                
                if(nrow(Non_gene_id) >= 1){
                  
                  
                Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")   
                  
                }
                
                part_gene_M <- rbind.data.frame(part_gene_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                  strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
            }
            
            if(part_gene_type == "gene"){
              
              if(nrow(part_trans_M_df) < 1){
                
                part_gene_M_df$start <- min(part_gene_M$start)
                part_gene_M_df$end <- max(part_gene_M$end)
                
                part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=gene-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
             
                parent_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,
                 pattern = "ID=.+;"),pattern = ";")[[1]][1])
                
                part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                
                Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type))  %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                  
                ) %>% distinct(.keep_all = TRUE)
                
                if(nrow(Non_gene_id) >= 1){
                  
                  
                  Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                  
                  
                }
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                
                part_gene_M <- rbind.data.frame(part_gene_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                  strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
              else if(nrow(part_trans_M_df) >= 1){
                
                part_gene_M_df$start <- min(part_gene_M$start)
                part_gene_M_df$end <- max(part_gene_M$end)
                part_trans_M_df$start <- min(part_gene_M$start)
                part_trans_M_df$end <- max(part_gene_M$end)
                
                part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) 
                  gsub("ID=",paste0("ID=gene-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
                  
                  gsub("Parent=",paste0("Parent=gene-p-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";Origin="),x)
                  
                )
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
                  
                  gsub("ID=",paste0("ID=rna-gp-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";transID="),x)
                  
                )
                
                parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,
                pattern = "ID=.+;"),pattern = ";")[[1]][1])
                
                part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                part_trans_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                
                Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | type %in% part_trans_M_df$type))  %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                  
                ) %>% distinct(.keep_all = TRUE)
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                  
                  
                }
                
                
                part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                  strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
            }  
            
            Unified_ls[[i]] <- part_gene_M %>% select(chrloc,source,type,start,end,score,
              strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }   
          
          else if(part_gene_No$Count[i] > 1){
            
            source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
            
            part_gene_type <- part_gene_No[i,]$type    
            
            part_gene_Mod_attr <- unique(part_gene[
              part_gene$type == part_gene_type,]$modified_attributes)
            
            #if(length(unique(part_gene_No$mod_attr)) == 1){
            
            part_gene_Mod <- part_gene %>% filter((source %in% source_tot ) & (modified_attributes %in% part_gene_Mod_attr))
            
            #}
            
            #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene[(part_gene$source %in% source_tot) & (part_gene$modified_attributes %in% part_gene_Mod_attr),]$source)) == length(unique(source_tot)))){
            
            #part_gene_Mod <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
            
            #}
            
            #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene[(part_gene$source %in% source_tot) & (part_gene$modified_attributes %in% part_gene_Mod_attr) ,]$source)) == length(unique(source_tot)))){
            
            #part_gene_Mod <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
            
            #}
            
            part_gene_id <- part_gene_Mod %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
            
            part_mRNA_id <- part_gene_Mod %>% filter(type == "mRNA" | type == "tRNA"| 
                                                       type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                       type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                       type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                       type == "transcript" | type == "pseudogenic_transcript" | 
                                                       type == "primary_transcript")
            
            Non_gene_id <-  part_gene_Mod %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
              select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
            
            #### finding the start and end positions for customary gene/transcript model
            
            part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
              select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
            
            part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
              select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
            
            part_gene_id_df$type <- part_gene_type
            
            part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x) 
              
              gsub("ID=","gID=",x)
              
            )
            
            part_gene_id_df$start <- min(part_gene_Mod$start)
            
            part_gene_id_df$end <- max(part_gene_Mod$end)
            
            part_gene_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
            
            if(part_gene_type == "pseudogene"){
              
              part_gene_id_df$attributes <- 
                paste0("ID=pseudogene-p-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
                 paste0(unique(part_gene_id$attributes),collapse = ";"))    
              
            }
            
            if(part_gene_type == "ncRNA_gene"){
              
              part_gene_id_df$attributes <- 
               paste0("ID=ncrnagene-p-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
                paste0(unique(part_gene_id$attributes),collapse = ";"))    
              
            }
            
            if(part_gene_type == "gene"){
              
              part_gene_id_df$attributes <- 
                paste0("ID=gene-p-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
                       paste0(unique(part_gene_id$attributes),collapse = ";"))    
              
              
            }
            
            
            part_gene_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")    
            
            if(nrow(part_mRNA_id_df) >= 1){
              
              part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
                
                gsub("ID=","transID=",x)
                
              )
              
              part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
                
                gsub("Parent=","Origin=",x)
                
              )
              
              part_mRNA_id_type <- part_mRNA_id_df$type
              
              part_mRNA_id_df$start <- min(part_gene_Mod$start)
              
              part_mRNA_id_df$end <- max(part_gene_Mod$end)
              
              part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
              
              if(part_gene_type == "pseudogene"){
                
                part_mRNA_id_df$attributes <- paste0("ID=rna-pseudo-p-",paste0(unique(part_gene_Mod$modified_attributes),
                 collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=pseudogene-p-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,
                ";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
                
              }
              
              if(part_gene_type == "ncRNA_gene"){
                
                part_mRNA_id_df$attributes <- paste0("ID=rna-nc-p-",paste0(unique(part_gene_Mod$modified_attributes),
                collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=ncrnagene-p-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,
                ";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
                
              }
              
              if(part_gene_type == "gene"){
                
                part_mRNA_id_df$attributes <- paste0("ID=rna-gp-",paste0(unique(part_gene_Mod$modified_attributes),
                 collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=gene-p-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,
                  ";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
                
              }
              
              #part_mRNA_id_df$modified_attributes <- paste0(unique(part_mRNA_id$modified_attributes),collapse = ";")    
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
              
            }
            
            else if(nrow(part_mRNA_id_df) != 1){
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),
                                                   pattern = ";")[[1]][1])
              
            }
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
              
              
              gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x)
              
            ) 
            
            if(nrow(Non_gene_id) >= 1){
              
              
              Non_gene_id$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
              
              
            }
            
            Unified_ls[[i]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>%
              select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
        }
        
        Unified_M_ls[[k]] <- do.call("rbind.data.frame",Unified_ls) 
        
      } 
      
    }
    
    else {
      
      Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
       strand,phase,modified_attributes,attributes) %>% 
        distinct(.keep_all = TRUE)
      
    }
    
    
  }
      
   if(unique(part_gene$strand) == "-"){
     
            if(unique(part_gene$gp_mod_attr) == "U4" | 
             unique(part_gene$gp_mod_attr) == "xtr" | 
             unique(part_gene$gp_mod_attr) == "RNaseP_nuc" | 
             unique(part_gene$gp_mod_attr) == "U2" | 
             unique(part_gene$gp_mod_attr) == "XENTR_v10003163" |
             unique(part_gene$gp_mod_attr) == "trnag-ccc" |
             unique(part_gene$gp_mod_attr) == "trnar-ucu" |
             unique(part_gene$gp_mod_attr) == "5S_rRNA" |
             unique(part_gene$gp_mod_attr) == "trnap-agg" |
             unique(part_gene$gp_mod_attr) == "trnad-guc" |
             unique(part_gene$gp_mod_attr) == "trnav-aac" |
             unique(part_gene$gp_mod_attr) == "trnas-uga" |
             unique(part_gene$gp_mod_attr) == "trnat-ugu" |
             unique(part_gene$gp_mod_attr) == "trnaa-cgc"|
             unique(part_gene$gp_mod_attr) == "trnae-uuc" |
             unique(part_gene$gp_mod_attr) == "U7" |
             unique(part_gene$gp_mod_attr) == "trnaq-cug" |
             unique(part_gene$gp_mod_attr) == "trnas-aga" |
             unique(part_gene$gp_mod_attr) == "trnaq-uug" 
     ){
       
              part_df <- part_gene %>% left_join(rank_data,by = c("type"))
              
              if(length(unique(part_df[part_df$type == "tRNA",]$type)) >= 1){
                
               part_df <- part_df %>% arrange(rank)
                
               part_df$transcript <- sapply(part_df$attributes,function(x) 
                  
                ifelse(str_detect(x,pattern = "Parent=.+"),
                         
                  gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1]),
                         
                  gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
                  
                )
                
                part_df$attributes <- as.vector(part_df$attributes)
                
                part_df$transcript <- as.vector(part_df$transcript)
                
                part_df <- part_df %>% group_by(chrloc,strand,modified_attributes) %>% 
                  arrange(transcript) %>%
                  ungroup() %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
                  distinct(.keep_all = TRUE)
                
                Unified_M_ls[[k]] <- part_df %>% select(chrloc,source,type,start,end,score,
                 strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
              if(length(unique(part_df[part_df$type == "ncRNA_gene" | part_df$type == "miRNA",]$type)) >= 1){
                
                #part_df <- part_df %>% left_join(rank_data,by = c("type"))
                
                part_df <- part_df %>% arrange(rank)
                
                part_df$transcript <- sapply(part_df$attributes,function(x) 
                  
                  ifelse(str_detect(x,pattern = "Parent=.+"),
                         
                         gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1]),
                         
                         gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
                  
                )
                
                part_df <- part_df %>% arrange(transcript) %>% 
                  select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                
                Unified_M_ls[[k]] <- part_df %>% select(chrloc,source,type,start,end,score,
                  strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              } 
              
              else{
                
                Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
                 strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
            
      }
    
     
    else if(((nrow(part_gene[part_gene$type == "gene",]) >= 2) |
         (nrow(part_gene[part_gene$type == "pseudogene",]) >= 2) |
         (nrow(part_gene[part_gene$type == "ncRNA_gene",]) >= 2)) & 
        (length(unique(part_gene[part_gene$type == "gene" | 
          part_gene$type == "pseudogene" | 
          part_gene$type == "ncRNA_gene" ,]$type) %in% 
          c("gene","pseudogene","ncRNA_gene")) == 1 )){
       
       part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
       
       part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
        type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
        type == "snRNA" | type == "scRNA" | type == "miRNA" |  
        type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
        type == "transcript" | type == "pseudogenic_transcript" | 
        type == "primary_transcript")
       
       Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
         select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
         distinct(.keep_all = TRUE)
       
       #### finding the start and end positions for customary gene/transcript model
       
       part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
         select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
         distinct(.keep_all = TRUE)
       
       part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
         select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>%
         distinct(.keep_all = TRUE)
       
       part_gene_id_df$start <- min(part_gene$start)
       
       part_gene_id_df$end <- max(part_gene$end)
       
       part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x) 
         
         gsub("ID=","gID=",x)
         
       )
       
       if(part_gene_id_df$type == "pseudogene"){
         
         part_gene_id_df$attributes <- paste0("ID=pseudogene-n-",paste0(unique(part_gene$modified_attributes),
         collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_gene_id$attributes),collapse = ";"))
         
       }
       
       if(part_gene_id_df$type == "ncRNA_gene"){
         
        part_gene_id_df$attributes <- paste0("ID=ncrnagene-n-",paste0(unique(part_gene$modified_attributes),
        collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_gene_id$attributes),collapse = ";")) 
         
         
       }
       
       if(part_gene_id_df$type == "gene"){
         
         part_gene_id_df$attributes <- paste0("ID=gene-n-",paste0(unique(part_gene$modified_attributes),
         collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_gene_id$attributes),collapse = ";"))    
         
       }
       
       part_gene_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
       
       
       if(nrow(part_mRNA_id_df) == 1){
         
         part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
           
           gsub("ID=","transID=",x)
           
         )
         
         part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
           
           gsub("Parent=","Origin=",x)
           
         )
         
         if(part_gene_id_df$type == "gene"){
           
           part_mRNA_id_df$attributes <- paste0("ID=rna-gn-",paste0(unique(part_gene$modified_attributes),
            collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=gene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
             paste0(unique(part_mRNA_id$attributes),collapse = ";"))
           
         }
         
         if(part_gene_id_df$type == "pseudogene"){
           
           part_mRNA_id_df$attributes <- paste0("ID=rna-pseudo-n-",paste0(unique(part_gene$modified_attributes),
              collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=pseudogene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
                paste0(unique(part_mRNA_id$attributes),collapse = ";"))
           
         }
         
         if(part_gene_id_df$type == "ncRNA_gene"){
           
           part_mRNA_id_df$attributes <- paste0("ID=rna-nc-n-",paste0(unique(part_gene$modified_attributes),
           collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=ncrnagene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
           paste0(unique(part_mRNA_id$attributes),collapse = ";"))
           
         }
         
         part_mRNA_id_df$start <- min(part_gene$start)
         
         part_mRNA_id_df$end <- max(part_gene$end)
         
         part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
         
         parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
         
       }
       
       else if(nrow(part_mRNA_id_df) < 1){
         
         parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
         
       }
       
       
       Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x) gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
       
       if(nrow(Non_gene_id) >= 1){
         
         
         Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
         
         
       }
       
       
       Unified_M_ls[[k]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
         select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
         distinct(.keep_all = TRUE) 
       
     }
     
     else if((nrow(part_gene[part_gene$type == "gene" | 
              part_gene$type == "pseudogene" | 
              part_gene$type == "ncRNA_gene",]) < 2) | 
             (nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
              (nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1 |
             nrow(part_gene[(part_gene$type == "ncRNA_gene"),]) >= 1))
             ){
       
       if(nrow(part_gene[part_gene$type == "gene" | 
                         part_gene$type == "pseudogene" | 
                         part_gene$type == "ncRNA_gene",]) < 2){
         
         part_gene_df <- part_gene %>% filter(type == "gene"| type == "pseudogene" | type == "ncRNA_gene") %>%
           select(chrloc,source,type,start,end,score,phase,strand,modified_attributes,attributes) 
         
         if(nrow(part_gene_df) == 1){
           
           #part_gene_df <- part_gene %>% filter(type == "gene"| type == "pseudogene" | type == "ncRNA_gene") %>%
           # select(chrloc,source,type,start,end,score,phase,strand,modified_attributes,attributes) 
           
           part_trans_M <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
            type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
             type == "snRNA" | type == "scRNA" | type == "miRNA" |  
              type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
               type == "transcript" | type == "pseudogenic_transcript" | 
                type == "primary_transcript")
           
           part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
             
             chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
             
           ) %>% distinct(.keep_all = TRUE)
           
           if(length((part_gene_df$type == "pseudogene")) == 1 |length((part_gene_df$type == "ncRNA_gene")) == 1 ){
             
             if(nrow(part_trans_M_df) >= 1){
               
               if(part_gene_df$type == "pseudogene"){
                 
                 #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                 
                 part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gID="),x))
                 
                 part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=pseudogene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";Origin="),x))
                 
                 part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-pseudo-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";transID="),x))
                 
                 
               }
               
               else if(part_gene_df$type == "ncRNA_gene"){
                 
                 #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                 
                 part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gID="),x))
                 
                 part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=ncrnagene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";Origin="),x))
                 
                 part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-nc-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";transID="),x))
                 
               }
               
               part_gene_df$start <- min(part_gene$start)
               
               part_gene_df$end <- max(part_gene$end)
               
               part_trans_M_df$start <- part_gene_df$start
               
               part_trans_M_df$end <- part_gene_df$end 
               
               part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
               
               part_trans_M_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")

               Non_gene_id <- part_gene  %>% filter(!(type %in% part_gene_df$type | 
               type %in% part_trans_M$type ))  %>% select(
                                                          
                 chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                          
                ) %>% distinct(.keep_all = TRUE)
               
               parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
               
               Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                 gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
               
               if(nrow(Non_gene_id) >= 1){
                 
                 
                 Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                 
                 
               }
               
               part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
                 select(chrloc,source,type,start,end,score,                                                                           
                 strand,phase,modified_attributes,attributes) %>% 
                 distinct(.keep_all = TRUE)
               
             }
             
             else if(nrow(part_trans_M) < 1){
               
               if(part_gene_df$type == "pseudogene"){
                 
                 #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                 
                 part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gId="),x))
                 
                 #parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";"))
                 
               }
               
               else if(part_gene_df$type == "ncRNA_gene"){
                 
                 #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                 
                part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gId="),x))
                 
                 #parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";"))
                 
                 
               }
               
               Non_gene_id <- part_gene  %>% filter(!(type == part_gene_df$type))  %>% select(
                 
                 chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                 
               ) %>% distinct(.keep_all = TRUE)
               
               part_gene_df$start <- min(part_gene$start)
               
               part_gene_df$end <- max(part_gene$end)
               
               part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
               
               parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
               
               Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                 gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
               
               if(nrow(Non_gene_id) >= 1){
                 
                 
                 Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                 
                 
               }
               
               part_gene <- rbind.data.frame(part_gene_df,Non_gene_id) %>% 
                 select(chrloc,source,type,start,end,score,                                                                           
                  strand,phase,modified_attributes,attributes) %>% 
                 distinct(.keep_all = TRUE)
               
             }
             
           }
           
           if(length(part_gene_df$type == "gene") == 1){
             
             if(nrow(part_trans_M_df) < 1){
               
               #ab <- part_gene_df[part_gene_df$type == "gene",]$attributes 
               
               #orginal_gene_id <- str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]
               
               part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) 
                 
                 gsub("ID=",paste0("ID=gene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gID="),
                      x)
               )
               
               part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
               
               parent_id <- gsub("ID=","",str_split(str_extract(part_gene_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
               
               Non_gene_id <- part_gene  %>% filter(!(type == part_gene_df$type))  %>% select(
                 
                 chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                 
               ) %>% distinct(.keep_all = TRUE)
               
               part_gene_df$start <- min(part_gene$start)
               
               part_gene_df$end <- max(part_gene$end)
               
               if(nrow(Non_gene_id) >= 1){
                 
                 
                 Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                 
                 
               }
               
               Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                 gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
               
               part_gene <- rbind.data.frame(part_gene_df,Non_gene_id) %>% 
                 select(chrloc,source,type,start,end,score,                                                                           
                        strand,phase,modified_attributes,attributes) %>% 
                 distinct(.keep_all = TRUE)
               
             }
             
             else if(nrow(part_trans_M_df) >= 1){
               
               part_gene_df$attributes <- sapply(part_gene_df$attributes, function(x)
                 
                 gsub("ID=",paste0("ID=gene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";gId="),x)
                 
               )
               
               
               part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
                 
                 gsub("Parent=",paste0("Parent=gene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";Origin="),x)
                 
               )
               
               part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
                 
                 gsub("ID=",paste0("ID=rna-gn-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_df$start,"_",part_gene_df$end,";Origin="),x)
                 
               )
               
               
               parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,
                pattern = "ID=.+;"),pattern = ";")[[1]][1])
               
               part_gene_df$start <- min(part_gene$start)
               
               part_gene_df$end <- max(part_gene$end)
               
               part_trans_M_df$start <- part_gene_df$start
               
               part_trans_M_df$end <- part_gene_df$end
               
               part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
               
               part_trans_M_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
               
               Non_gene_id <- part_gene  %>% filter(!(type %in% part_gene_df$type | type %in% part_trans_M$type))  %>% select(
                 
               chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                 
               ) %>% distinct(.keep_all = TRUE)
               
               Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                 gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
               
               if(nrow(Non_gene_id) >= 1){
                 
                 
                 Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                 
                 
               }
               
               
               part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
                 select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
                 distinct(.keep_all = TRUE)
               
             }
           }
           
         }
         
         else if(nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) < 1 ){
           
           if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) < 1 )){
             
             org_parent_id <- gsub("Parent=","",str_split(str_extract(part_gene$attributes[1],pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
             
             mod_gene_id <- part_gene[1,]
             
             mod_mRNA_id <- part_gene[1,]
             
             mod_gene_id$type <- "gene"
             mod_mRNA_id$type <- "mRNA"
             
             mod_gene_id$start <- min(part_gene$start)
             mod_gene_id$end <- max(part_gene$end)
             
             mod_mRNA_id$start <- min(part_gene$start)
             mod_mRNA_id$end <- max(part_gene$end)
             
             mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             
             gene_id <- paste0("ID=gene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",mod_gene_id$start,"_",mod_gene_id$end,";","OrigingID:",org_parent_id,";","gbkey=gene")
             
             parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
             
             mRNA_id <- paste0("ID=rna-gn-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",mod_gene_id$start,"_",mod_gene_id$end,";Parent=",parent_id,";")
             
             pid <- gsub("ID=","",str_split(str_extract(mRNA_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
             
             part_gene$attributes <- sapply(part_gene$attributes,function(x) 
               
               gsub("Parent=",paste0("Parent=",pid,";Origin="),x)
               
             )
             
             mod_gene_id$attributes <- gene_id
             mod_mRNA_id$attributes <- mRNA_id
             
             mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             
             part_gene$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             
             part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,part_gene) %>% select(chrloc,source,type,start,end,score,
              strand,phase,modified_attributes,attributes) %>% 
               distinct(.keep_all = TRUE)
             
           }
           
           else if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) >= 1 )){
             
             mod_mRNA_id <- part_gene[which.max(part_gene$type == "mRNA"),]
             
             org_parent_id <- gsub("Parent=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
             
             #mRNA_id <- paste0("ID=rna-",paste0(unique(part_gene$modified_attributes),collapse = "_"),";Parent=",parent_id,";")
             
             mod_gene_id <- mod_mRNA_id       
             
             mod_gene_id$type <- "gene"
             mod_mRNA_id$type <- "mRNA"
             
             mod_gene_id$start <- min(part_gene$start)
             mod_gene_id$end <- max(part_gene$end)
             
             mod_mRNA_id$start <- min(part_gene$start)
             mod_mRNA_id$end <- max(part_gene$end)
             
             mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             
             gene_id <- paste0("ID=gene-n-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",mod_gene_id$start,"_",mod_gene_id$end,";","OrigingID:",org_parent_id,";","gbkey=gene")
             
             parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
             
             mod_gene_id$attributes <- gene_id
             
             mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
             mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("ID=",paste0("ID=rna-gn-",paste0(unique(part_gene$modified_attributes),collapse = "_"),"_",chr_d,"_",mod_gene_id$start,"_",mod_gene_id$end,";Origin="),x))
             
             p_id <- gsub("ID=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "ID.+;"),pattern = ";")[[1]][1])
             
             Non_gene_id <- part_gene %>% filter(!(type %in% mod_gene_id$type | type %in% mod_mRNA_id$type))
             
             Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
               
               gsub("Parent=",paste0("Parent=",p_id,";Origin="),x)
               
             )
             
             mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             
             if(nrow(Non_gene_id) >= 1){
             
             Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
             
             
             }
             
             part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,Non_gene_id) %>% select(chrloc,source,type,start,end,score,
              strand,phase,modified_attributes,attributes) %>% 
               distinct(.keep_all = TRUE)        
             
           }
           
           if(unique(part_gene$source) == "NCBI"){
             
             part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,
              strand,phase,modified_attributes,attributes) %>% 
               distinct(.keep_all = TRUE)        
             
             
           }
           
         }
         
         Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
          strand,phase,modified_attributes,attributes) %>% 
           distinct(.keep_all = TRUE)
         
       }
       
       else if(nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
        (nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1 |
        nrow(part_gene[(part_gene$type == "ncRNA_gene"),]) >= 1)){
         
         part_gene_No <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene") %>% 
           group_by(type) %>% summarise( Count = n(), source_t = paste0(source,collapse = ";"),mod_attr = paste0(unique(modified_attributes),collapse = ";") )
         
         for(i in 1 : length(part_gene_No$Count)){
           
           if(part_gene_No$Count[i] == 1){
             
             source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
             
             part_gene_Mod_attr <- unique(str_split(part_gene_No$mod_attr[i],pattern = ";")[[1]])
             
             part_gene_type <- part_gene_No[i,]$type    
             
             #if(length(unique(part_gene_No$mod_attr)) == 1){
             
             part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
             
             #}
             
             #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene_No$mod_attr)) == nrow(part_gene_No))){
             
             #part_gene_M <- part_gene %>% filter( (source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
             
             #}
             
             #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene_No$mod_attr)) != nrow(part_gene_No))){
             
             #part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
             
             #}
             
             part_trans_M <- part_gene_M %>% filter(type == "mRNA" | type == "tRNA" | 
                                                      type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | 
                                                      type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                      type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                      type == "transcript" | type == "pseudogenic_transcript" | 
                                                      type == "primary_transcript")
             
             part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
               
               chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
               
             ) %>% distinct(.keep_all = TRUE)
             

             part_gene_M_df <- part_gene_M %>% filter(type == part_gene_type) %>% select(
               
               chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes,Overlap_Value
               
             ) %>% distinct(.keep_all = TRUE)
             
             part_gene_M_df <- part_gene_M_df[which.max(part_gene_M_df$Overlap_Value),] %>% select(
               
               chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
               
             ) %>% distinct(.keep_all = TRUE)
             
             part_gene_M_df$type <- part_gene_type
             
             if(part_gene_type == "pseudogene" | part_gene_type == "ncRNA_gene"){
               
               if(nrow(part_trans_M_df) == 1){
                 
                 part_gene_M_df$start <- min(part_gene_M$start)
                 part_gene_M_df$end <- max(part_gene_M$end)
                 part_trans_M_df$start <- min(part_gene_M$start)
                 part_trans_M_df$end <- max(part_gene_M$end)
                 part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                 
                 
                 if(part_gene_type == "pseudogene"){
                   
                   #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                   
                   part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                   
                   part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=pseudogene-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";Origin="),x))
                   
                   part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-pseudo-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";transID="),x))
                   
                   parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                   
                 }
                 
                 else if(part_gene_type == "ncRNA_gene"){
                   
                   #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                   
                   part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-n-",paste0(unique(part_gene_M$modified_attributes)),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                   
                   part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",paste0("Parent=ncrnagene-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";Origin="),x))
                   
                   part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",paste0("ID=rna-ncg-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";transID="),x))
                   
                   parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                   
                   
                 }
                 
                 
                 part_trans_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")

                 Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | 
                  type %in% part_trans_M$type ))  %>% select(
                                                              
                  chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                              
                  ) %>% distinct(.keep_all = TRUE)
                 
                 
                 Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                   gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                 
                 if(nrow(Non_gene_id) >= 1){
                   
                   
                   Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                   
                   
                 }
                 
                 part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                   select(chrloc,source,type,start,end,score,                                                                           
                   strand,phase,modified_attributes,attributes) %>% 
                   distinct(.keep_all = TRUE)
                 
               }
               
               else if(nrow(part_trans_M_df) != 1){
                 
                 part_gene_M_df$start <- min(part_gene_M$start)
                 part_gene_M_df$end <- max(part_gene_M$end)
                 part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                 
                 if(part_gene_type == "pseudogene"){
                   
                   #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                   
                   part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=pseudogene-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                   
                 }
                 
                 else if(part_gene_type == "ncRNA_gene"){
                   
                   #gene_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                   
                   part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=ncrnagene-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))             
                   
                 }
                 
                 
                 Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type))  %>% select(
                   
                   chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                   
                 ) %>% distinct(.keep_all = TRUE)
                 
                 parent_id_M <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
                 
                 Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                   gsub("Parent=",paste0("Parent=",parent_id_M,";Origin="),x))
                 
                 if(nrow(Non_gene_id) >= 1){
                   
                   
                   Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                   
                   
                 }
                 
                 part_gene_M <- rbind.data.frame(part_gene_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                  strand,phase,modified_attributes,attributes) %>% 
                   distinct(.keep_all = TRUE)
                 
               }
               
             }
             
             if(part_gene_type == "gene"){
               
               if(nrow(part_trans_M_df) < 1){
                 
                 part_gene_M_df$start <- min(part_gene_M$start)
                 part_gene_M_df$end <- max(part_gene_M$end)
                 part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                 
                 part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",paste0("ID=gene-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                 
                 parent_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,
                 pattern = "ID=.+;"),pattern = ";")[[1]][1])
                 
                 Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type))  %>% select(
                   
                   chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                   
                 ) %>% distinct(.keep_all = TRUE)
                 
                 Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                   gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                 
                 if(nrow(Non_gene_id) >= 1){
                   
                   Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                   
                   
                 }
                 
                 part_gene_M <- rbind.data.frame(part_gene_M_df,Non_gene_id) %>% 
                   select(chrloc,source,type,start,end,score,                                                                           
                   strand,phase,modified_attributes,attributes) %>% 
                   distinct(.keep_all = TRUE)
                 
               }
               
               else if(nrow(part_trans_M_df) >= 1){
                 
                 part_gene_M_df$start <- min(part_gene_M$start)
                 part_gene_M_df$end <- max(part_gene_M$end)
                 part_trans_M_df$start <- min(part_gene_M$start)
                 part_trans_M_df$end <- max(part_gene_M$end)
                 part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ';')
                 
                 part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) 
                   gsub("ID=",paste0("ID=gene-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";gID="),x))
                 
                 part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
                   
                   gsub("Parent=",paste0("Parent=gene-n-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";Origin="),x)
                   
                 )
                 
                 part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes, function(x)
                   
                   gsub("ID=",paste0("ID=rna-gn-",paste0(unique(part_gene_M$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_M_df$start,"_",part_gene_M_df$end,";transID="),x)
                   
                 )    
                 
                 
                 parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,
                  pattern = "ID=.+;"),pattern = ";")[[1]][1])
                 
                 
                 part_trans_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ';')
                 
                 Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | type %in% part_trans_M_df$type))  %>% select(
                   
                   chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                   
                 ) %>% distinct(.keep_all = TRUE)
                 
                 Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                   gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                 
                 if(nrow(Non_gene_id) >= 1){
                   
                   
                   Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ';')
                   
                   
                 }
                 
                 part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                   select(chrloc,source,type,start,end,score,                                                                           
                    strand,phase,modified_attributes,attributes) %>% 
                   distinct(.keep_all = TRUE)
                 
               }
             }  
             
             Unified_ls[[i]] <- part_gene_M %>% select(chrloc,source,type,start,end,score,
              strand,phase,modified_attributes,attributes) %>% 
               distinct(.keep_all = TRUE)
             
           }   
           
           else if(part_gene_No$Count[i] > 1){
             
             source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
             
             part_gene_type <- part_gene_No[i,]$type    
             
             part_gene_Mod_attr <- unique(part_gene[
               part_gene$type == part_gene_type,]$modified_attributes)
             
             #if(length(unique(part_gene_No$mod_attr)) == 1){
             
             part_gene_Mod <- part_gene %>% filter((source %in% source_tot ) & (modified_attributes %in% part_gene_Mod_attr))
             
             #}
             
             #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene[(part_gene$source %in% source_tot) & (part_gene$modified_attributes %in% part_gene_Mod_attr),]$source)) == length(unique(source_tot)))){
             
             #part_gene_Mod <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
             
             #}
             
             #else if((length(unique(part_gene_No$mod_attr)) > 1) & (length(unique(part_gene[(part_gene$source %in% source_tot) & (part_gene$modified_attributes %in% part_gene_Mod_attr) ,]$source)) == length(unique(source_tot)))){
             
             #part_gene_Mod <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
             
             #}
             
             part_gene_id <- part_gene_Mod %>% filter(type == "gene" | type == "pseudogene")
             
             part_mRNA_id <- part_gene_Mod %>% filter(type == "mRNA" | type == "tRNA"| 
                                                        type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                        type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                        type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                        type == "transcript" | type == "pseudogenic_transcript" | 
                                                        type == "primary_transcript")
             
             Non_gene_id <-  part_gene_Mod %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
               select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
             
             #### finding the start and end positions for customary gene/transcript model
             
             part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
               select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
             
             part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
               select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
             
             part_gene_id_df$type <- part_gene_type
             
             part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x) 
               
               gsub("ID=","gID=",x)
               
             )
             
             part_gene_id_df$start <- min(part_gene_Mod$start)
             
             part_gene_id_df$end <- max(part_gene_Mod$end)
             
             part_gene_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
             
             if(part_gene_type == "pseudogene"){
               
               part_gene_id_df$attributes <- 
                 paste0("ID=pseudogene-n-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
                  paste0(unique(part_gene_id$attributes),collapse = ";"))    
               
               
             }
             
             if(part_gene_type == "ncRNA_gene"){
               
               part_gene_id_df$attributes <- 
                 paste0("ID=ncrnagene-n-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
                        paste0(unique(part_gene_id$attributes),collapse = ";"))    
               
               
             }
             
             
             if(part_gene_type == "gene"){
               
               part_gene_id_df$attributes <- 
                 paste0("ID=gene-n-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",
                        paste0(unique(part_gene_id$attributes),collapse = ";"))    
               
               
             }
             
             
             #part_gene_id_df$modified_attributes <- paste0(unique(part_gene_id$modified_attributes),collapse = ";")    
             
             if(nrow(part_mRNA_id_df) >= 1){
               
               part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
                 
                 gsub("ID=","transID=",x)
                 
               )
               
               
               part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
                 
                 gsub("Parent=","Origin=",x)
                 
               )
               
               part_mRNA_id_type <- part_mRNA_id_df$type
               
               part_mRNA_id_df$start <- min(part_gene_Mod$start)
               
               part_mRNA_id_df$end <- max(part_gene_Mod$end)
               
               if(part_gene_type == "pseudogene"){
                 
                 part_mRNA_id_df$attributes <- paste0("ID=rna-pseudo-n-",paste0(unique(part_gene_Mod$modified_attributes),
                  collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=pseudogene-n-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),
                  "_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
                 
               }
               
               if(part_gene_type == "ncRNA_gene"){
                 
                 part_mRNA_id_df$attributes <- paste0("ID=rna-nc-n-",paste0(unique(part_gene_Mod$modified_attributes),
                  collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=ncrnagene-n-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),
                  "_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
                 
               }
               
               if(part_gene_type == "gene"){
                 
                 part_mRNA_id_df$attributes <- paste0("ID=rna-gn-",paste0(unique(part_gene_Mod$modified_attributes),
                 collapse = "_"),"_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";","Parent=gene-n-",paste0(unique(part_gene_Mod$modified_attributes),collapse = "_"),
                 "_",chr_d,"_",part_gene_id_df$start,"_",part_gene_id_df$end,";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
                 
               }
               
               part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")    
               
               
               parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
               
             }
             
             else if(nrow(part_mRNA_id_df) != 1){
               
               parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),
               pattern = ";")[[1]][1])
               
             }
             
             Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
               
               
               gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x)
               
             ) 
             
             if(nrow(Non_gene_id) >= 1){
               
               
               Non_gene_id$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")    
               
               
             }
             
             Unified_ls[[i]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>%
               select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
               distinct(.keep_all = TRUE)
             
           }
           
         }
         
         Unified_M_ls[[k]] <- do.call("rbind.data.frame",Unified_ls) 
         
       } 
       
     }
     
     else {
       
       Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
        strand,phase,modified_attributes,attributes) %>% 
         distinct(.keep_all = TRUE)
       
      }
     
     
     }
  
    }
    
  Unified_F_ls[[j]] <- do.call("rbind.data.frame",Unified_M_ls)
  
  }

  Unified_d <- do.call("rbind.data.frame",Unified_F_ls)
  return(Unified_d)

}


filtered_Unified_NCBI_Ens_XGC_ls <- list()
Unified_NCBI_Ens_XGC_ls <- list()

for (i in 1: length(unique(Final_grped_df$chrloc))){
  
  filtered_Unified_NCBI_Ens_XGC_ls[[i]] <- Final_grped_df %>% filter(chrloc == unique(Final_grped_df$chrloc)[i])
  Unified_NCBI_Ens_XGC_ls[[i]] <- unified_NCBI_Ens_XGC(filtered_Unified_NCBI_Ens_XGC_ls[[i]])
  
}

Final_Unified_df <- do.call("rbind.data.frame",Unified_NCBI_Ens_XGC_ls)

Final_Unified_df$attributes <- as.vector(Final_Unified_df$attributes)

Final_Unified_d <- Final_Unified_df %>% distinct(.keep_all = TRUE)

Test_Chr1_Unified_File <- Filtered_df %>% filter(chrloc == "Chr1") %>% 
 select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
 distinct(.keep_all = TRUE)

Test_Chr1_Unified_File$source <- "Xenbase"

Final_Unified_M_Df <- Final_Unified_d %>% 
  select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
  distinct(.keep_all = TRUE)

Final_Unified_M_Df$source <- "Xenbase"

write.table(Final_Unified_M_Df,file = "NCBI_Ensembl_XGC_Unified_M_Final.gff3",
sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

write.table(Test_Chr1_Unified_File,file = "Test_Chr1_NCBI_Ens_XGC_Unified_File.gff3",
sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file


########## removing duplicates ###################

filtered_rep_Unified_d <- Final_Unified_d %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene") %>%
 group_by(modified_attributes) %>%
  summarise(
  
    mod_attr_c = n()
  
)

filtered_df <- filtered_rep_Unified_d[filtered_rep_Unified_d$mod_attr_c > 1,]

filtered_d_2 <- filtered_df[filtered_df$mod_attr_c == 2,]

filtered_d_2 <- filtered_df[filtered_df$mod_attr_c > 2,]


#### checking source ####

#Filtered_Data <- Final_Unified_d[!(Final_Unified_d$modified_attributes %in% filtered_d_2$modified_attributes) & (Final_Unified_d$type %in% c("gene","mRNA","exon")),]

#gp_mod_att <- unique(Final_grped_df[Final_grped_df$modified_attributes %in% filtered_d_2$modified_attributes,]$gp_mod_attr)

########################## 

g_trnaq_cug <- Filtered_df %>% filter(modified_attributes == "trnaq-cug")

g_trnaq_cug$ID <- sapply(g_trnaq_cug$attributes,function(x) gsub("ID=","",
                                                                 
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_trnaq_cug$Parent <- sapply(g_trnaq_cug$attributes,function(x) gsub("Parent=","",
                                                                     
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))  
  
#########################  

g_trnat_ugu <- Final_Unified_df %>% filter(modified_attributes == "trnat-ugu")

g_trnat_ugu$ID <- sapply(g_trnat_ugu$attributes,function(x) gsub("ID=","",
                                                                 
  str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_trnat_ugu$Parent <- sapply(g_trnat_ugu$attributes,function(x) gsub("Parent=","",
                                                                     
  str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

############## checking for double repeats ##########

g_abhd17b_o <- Final_grped_df %>% filter(modified_attributes == "abhd17b" | gp_mod_attr == "abhd17b")

g_abhd17b_o$ID <- sapply(g_abhd17b_o$attributes,function(x) gsub("ID=","",
                                                               
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_abhd17b_o$Parent <- sapply(g_abhd17b_o$attributes,function(x) gsub("Parent=","",
                                                                       
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

g_XENTR_v10001772_abhd <- Final_Unified_d %>% filter(modified_attributes == "abhd17b;XENTR_v10001772")

g_abhd17b <- Final_Unified_d %>% filter(modified_attributes == "abhd17b")

###checking for duplicate values after removal of dup

g_mod_abhd17b <- Filtered_df %>% filter(modified_attributes == "abhd17b")

g_mod_abhd17b_XENTR_v10001772 <- Filtered_df %>% filter(modified_attributes == "abhd17b;XENTR_v10001772")

g_mod_abhd17b_XENTR_v10001772$ID <- sapply(g_mod_abhd17b_XENTR_v10001772$attributes,function(x) gsub("ID=","",
                                                                 
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_mod_abhd17b_XENTR_v10001772$Parent <- sapply(g_mod_abhd17b_XENTR_v10001772$attributes,function(x) gsub("Parent=","",
                                                                     
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

#################################################################

g_trnac_gca <- Filtered_df %>% filter(modified_attributes == "trnac-gca")

g_trnac_gca$ID <- sapply(g_trnac_gca$attributes,function(x) gsub("ID=","",
                                                                                                     
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_trnac_gca$Parent <- sapply(g_trnac_gca$attributes,function(x) gsub("Parent=","",
                                                                                                         
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

##############################################################

g_U4 <- Filtered_df %>% filter(modified_attributes == "U4")

g_U4$ID <- sapply(g_U4$attributes,function(x) gsub("ID=","",
                                                                 
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_U4$Parent <- sapply(g_U4$attributes,function(x) gsub("Parent=","",
                                                                     
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

##################################################################

##### removing XENTR ones #####

Final_Unified_M <- Final_Unified_d %>% mutate(
  
  mod_attr_split = str_split(modified_attributes,pattern = ";"),
  
  #rep_ids = ifelse(rep_val == "TRUE",str_split(modified_attributes,pattern = ";"),modified_attributes)
  
)

Final_Unified_M$mod_attr_c <- lapply(Final_Unified_M$mod_attr_split, function(x)
  
  length(x)
  
  )

#rep_val_values <- unique(Final_Unified_M$rep_ids)

Final_df <- Final_Unified_M %>% mutate(
  
  #mod_count = length(mod_attr_split),
  rep_val = ifelse((str_detect(modified_attributes,pattern = "XENTR") & mod_attr_c >= 2),"TRUE","FALSE"),
  #rep_final = ifelse(rep_val == "TRUE",rep_ids[!(str_detect(rep_ids,"XENTR"))],modified_attributes)
  
  
)

Final_Unified_M_Values <- Final_df[Final_df$rep_val == "TRUE",]

rep_val_values <- lapply(Final_Unified_M_Values$mod_attr_split, function(x)
  
  x[!(str_detect(x,"XENTR"))]
  
  )

rep_v <- unique(unlist(rep_val_values))

g_jade1 <- Final_Unified_d %>% filter(modified_attributes == "jade1")

Filtered_df <- Final_Unified_d %>% filter(!(modified_attributes %in% rep_v))

#Unified_Final_D <- Filtered_df %>% select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% distinct(.keep_all = TRUE)

#Unified_D <- read.delim("NCBI_Ens_XGC_Unified_Final.gff3",sep = "\t",header = FALSE)

Test_Chr1_Unified_File <- Filtered_df %>% filter(chrloc == "Chr1") %>% 
  select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
  distinct(.keep_all = TRUE)

Test_Chr1_Unified_File$source <- "Xenbase"

#Unified_Final_D$source <- "Xenbase"

#write.table(Unified_Final_D,file = "NCBI_Ens_XGC_Unified_Final.gff3",
#sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

write.table(Test_Chr1_Unified_File,file = "Test_Chr1_NCBI_Ens_XGC_Unified_File.gff3",
sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

######### sorting data ########

#Unified_D <- read.delim("NCBI_Ens_XGC_Unified_Final.gff3",sep = "\t",header = FALSE)

#names(Unified_D) <- c("chrloc","source","type","start","end","score","strand","phase","attributes")

#Sort_Unified_F <- Sorted_Unified_Data %>% arrange(rank) %>% select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% distinct(.keep_all = TRUE)

UD_s <- str_sort(unique(Filtered_df$chrloc),numeric = TRUE)

rank_d <- data.frame(chr = UD_s, rank = 1:77)

Sorted_Unified_Data <- Filtered_df %>% left_join(rank_d,by = c("chrloc" = "chr"))

Sort_Unified_F <- Sorted_Unified_Data %>% arrange(rank) %>% 
  select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
  distinct(.keep_all = TRUE)

Sort_Unified_F$source <- "Xenbase"

write.table(Sort_Unified_F,file = "Final_NCBI_Ens_XGC_Unified_Model.gff3",
sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing into final file

#####################################################

#Final_Unified_df <- Final_Unified_df %>% filter(!(modified_attributes %in% filtered_df$modified_attributes))

######################################################################

#####################################################################

#Unified_M_Df$ID <- sapply(Unified_M_Df$attributes,function(x) gsub("ID=","",
                                                                   
 #   str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

#Unified_M_Df$Parent <- sapply(Unified_M_Df$attributes,function(x) gsub("Parent=","",
                                                                       
 #   str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

#error_file <- read.table("GFF3_error_uniq.txt",header = TRUE,sep = "\t")

#names(error_file) <- c("description")  

#error_file$description <- sapply(error_file$description,function(x) str_trim(x,c("both")))

#Missing_P <- Unified_M_Df %>% filter(ID %in% error_file$description)  

#####################################################################

Test_Chr1_Unified_File <- Sorted_Unified_Data %>% filter(chrloc == "Chr1") %>% 
  select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
  distinct(.keep_all = TRUE)

Test_Chr1_Unified_File$source <- "Xenbase"

Test_Chr1_Unified_File$ID <- sapply(Test_Chr1_Unified_File$attributes,function(x) gsub("ID=","",
                                                                   
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

Test_Chr1_Unified_File$Parent <- sapply(Test_Chr1_Unified_File$attributes,function(x) gsub("Parent=","",
                                                                       
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

error_file <- read.table("GFF3_error_uniq.txt",header = TRUE,sep = "\t")

names(error_file) <- c("description")  

error_file$description  <- sapply(error_file$description, function(x)
  
  gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1])
  
  )

#Missing_P <- Test_Chr1_Unified_File %>% filter(Parent %in% error_file$description)  

#Test_Unified_Gene_Models <- Unified_M_Df %>% filter(modified_attributes == "aacs") %>%
 # select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
  #distinct(.keep_all = TRUE)

#Test_Unified_Gene_Models$source <- "Xenbase"

#write.table(Test_Unified_Gene_Models,file = "Test_NCBI_Ens_XGC_Unified_Models.gff3",
#sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

#### checking the custom models of some complicated genes overall and chromosome in particular

abv1 <- Final_grped_df %>% filter(modified_attributes == "b3gnt4" | modified_attributes == "pycr3")

g_pycr3 <- Final_Unified_df %>% filter(modified_attributes == "pycr3")

abv <- Final_Unified_df %>% filter(modified_attributes == "b3gnt4;pycr3" | modified_attributes == "b3gnt4;b3gnt4" | modified_attributes == "pycr3;pycr3") %>% 
distinct(.keep_all = TRUE)

abv$ID <- sapply(abv$attributes,function(x) gsub("ID=","",
                                                   
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

abv$Parent <- sapply(abv$attributes,function(x) gsub("Parent=","",
                                                       
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

abg1 <- Final_Unified_df %>% filter(modified_attributes == "5S_rRNA") %>% distinct(.keep_all = TRUE)

######################################################################

abg3 <- Test_Chr1_Unified_File %>% filter(modified_attributes == "trnat-ugu")

abg3_o <- Final_Unified_df %>% filter(modified_attributes == "trnat-ugu")

abg3_o$ID <- sapply(abg3_o$attributes,function(x) gsub("ID=","",
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

abg3_o$Parent <- sapply(abg3_o$attributes,function(x) gsub("Parent=","",
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

####### upon observation, the complicated models are properly assigned according to 
### their respective strands,chromosomes and start and end.

#################### long ncRNA genes

ncRNA_D <- Final_Unified_d %>% filter(modified_attributes == "U4" | modified_attributes == "xtr-mir-101a-1" | modified_attributes == "xtr-mir-146b" | modified_attributes == "smyd5" | modified_attributes == "RNAaseP_nuc" | modified_attributes == "LOC116406437" | modified_attributes == "LOC101734477")

ncRNA_D_o <- Final_grped_df %>% filter(modified_attributes == "U4" | modified_attributes == "xtr-mir-101a-1" | modified_attributes == "xtr-mir-146b" | modified_attributes == "smyd5" | modified_attributes == "RNAaseP_nuc" | modified_attributes == "LOC116406437" | modified_attributes == "LOC101734477")

