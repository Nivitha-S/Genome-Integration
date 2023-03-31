########## NCBI Ens XGC Unified statistics ############

Sort_Unified_F <- read.delim("Sort_Unified_F.gff3",sep = "\t",header = FALSE)

names(Sort_Unified_F) <- c("chrloc","source","type","start","end","score","strand","phase","modified_attributes","attributes","gp_mod_attr","Overlap_Value")


Sort_Unified_F_grped <- Sort_Unified_F %>% group_by(chrloc,strand,gp_mod_attr) %>% mutate(
  
  
  f_mod_attr = paste0(unique(modified_attributes),collapse = ";")
  
)

Final_Unified_M <- Sort_Unified_F_grped %>% mutate(
  
  mod_attr_split = str_split(f_mod_attr,pattern = ";")
  
)

Final_Unified_M$mod_attr_c <- lapply(Final_Unified_M$mod_attr_split, function(x)
  
  length(x)
  
)

Final_df <- Final_Unified_M %>% mutate(
  
  rep_val = ifelse((str_detect(f_mod_attr,pattern = "XENTR") & mod_attr_c >= 2),"TRUE","FALSE")
  
)

Final_Unified_M_Values <- Final_df[Final_df$rep_val == "TRUE",]

rep_val_values <- lapply(Final_Unified_M_Values$mod_attr_split, function(x)
  
  x[!(str_detect(x,"XENTR"))]
  
)

rep_v <- unique(unlist(rep_val_values))

Filtered_df <- Final_Unified_M %>% filter(!(f_mod_attr %in% rep_v))

#Filtered_Gene_Data <- Filtered_df %>% filter(type == "gene")

#trna_mod_atr <- data.frame(ids = Filtered_df[Filtered_df$type == "tRNA",]$f_mod_attr)

#Filtered_Gene_df <- Filtered_Gene_Data %>% filter(!(f_mod_attr %in% unique(trna_mod_atr$ids)))

#Com_Filtered_df <- Filtered_df %>% left_join(
  
 # Combined_gp_var_df,by = c("chrloc","strand","f_mod_attr" = "modified_attributes"))

#Com_Filtered_d <- Com_Filtered_df %>% mutate(
  
 # Final_Mod_Attr = ifelse(is.na(Combined_Mod_Attr),f_mod_attr,Combined_Mod_Attr)
  
#) %>% select(chrloc,source,type,start = start.x,end = end.x,score,strand,phase,f_mod_attr,attributes,Final_Mod_Attr) %>% distinct(.keep_all = TRUE)

######## Computing NCBI Ens XGC Statistics ################

Unified_gene_final <- Filtered_df %>% filter(type == "gene") %>% group_by(chrloc,strand,f_mod_attr) %>% mutate(
  
  Total_Source = paste0(source,collapse = ";")
  
) 

Unified_gene_final$Total_Source_Split <- sapply(Unified_gene_final$Total_Source, function(x)
  
  str_split(x,pattern = ";")
  
)

Unified_gene_final$NCBI_Count <- sapply(Unified_gene_final$Total_Source_Split, function(x)
  
  as.vector(unlist(length(x[x == "NCBI"])))  
  
)

Unified_gene_final$Genbank_Count <- sapply(Unified_gene_final$Total_Source_Split, function(x)
  
  as.vector(unlist(length(x[x == "Genbank"])))  
  
)

Unified_gene_final$ensembl_Count <- sapply(Unified_gene_final$Total_Source_Split, function(x)
  
  as.vector(unlist(length(x[x == "ensembl"])))  
  
)

Unified_gene_final$NCBI_Count <- as.vector(Unified_gene_final$NCBI_Count)
Unified_gene_final$Genbank_Count <- as.vector(Unified_gene_final$Genbank_Count)
Unified_gene_final$ensembl_Count <- as.vector(Unified_gene_final$ensembl_Count)

UF <- Unified_gene_final %>% filter(Final_Mod_Attr == unique(Unified_gene_final$Final_Mod_Attr)[1])

######## NCBI & Ens & XGC #########

one_one_one_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 1 & 
                     Unified_gene_final$Genbank_Count == 1 & Unified_gene_final$ensembl_Count == 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
###16058 

complicated_count <- Unified_gene_final[(Unified_gene_final$NCBI_Count >= 1 & 
                     Unified_gene_final$Genbank_Count >= 1 & Unified_gene_final$ensembl_Count >= 1) & 
                     !(Unified_gene_final$NCBI_Count == 1 & Unified_gene_final$Genbank_Count == 1 & 
                        Unified_gene_final$ensembl_Count == 1),] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
###1998 

######## NCBI & XGC #########

NCBI_XGC_one_to_one_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 1 & 
                              Unified_gene_final$Genbank_Count == 1 & 
                               Unified_gene_final$ensembl_Count == 0,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##1867 

NCBI_XGC_Complicated <- Unified_gene_final[(Unified_gene_final$NCBI_Count >= 1 & 
                         Unified_gene_final$Genbank_Count >= 1 & 
                           Unified_gene_final$ensembl_Count == 0) & !(Unified_gene_final$NCBI_Count == 1 & 
                            Unified_gene_final$Genbank_Count == 1 & Unified_gene_final$ensembl_Count == 0),] %>% 
                             distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##157

NCBI_XGC_one_many_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 1 & 
                            Unified_gene_final$Genbank_Count > 1 & 
                              Unified_gene_final$ensembl_Count == 0,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##105

NCBI_XGC_many_one_count <- Unified_gene_final[Unified_gene_final$NCBI_Count > 1 & 
                            Unified_gene_final$Genbank_Count == 1 & Unified_gene_final$ensembl_Count == 0,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##28

NCBI_XGC_many_many_count <- Unified_gene_final[Unified_gene_final$NCBI_Count > 1 & 
Unified_gene_final$Genbank_Count > 1 & Unified_gene_final$ensembl_Count == 0,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")

##24

######## NCBI & Ens #########

NCBI_Ens_one_to_one_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 1 & 
                              Unified_gene_final$Genbank_Count == 0 & 
                                Unified_gene_final$ensembl_Count == 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")

##282

NCBI_Ens_Complicated <- Unified_gene_final[(Unified_gene_final$NCBI_Count >= 1 & Unified_gene_final$Genbank_Count == 0 & 
                         Unified_gene_final$ensembl_Count >= 1) & 
                         !(Unified_gene_final$NCBI_Count == 1 & 
                            Unified_gene_final$Genbank_Count == 0 & 
                             Unified_gene_final$ensembl_Count == 1),] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##25

NCBI_Ens_one_to_many_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 1 & 
                              Unified_gene_final$Genbank_Count == 0 & 
                              Unified_gene_final$ensembl_Count > 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")

##20

NCBI_Ens_many_to_one_count <- Unified_gene_final[Unified_gene_final$NCBI_Count > 1 & 
                               Unified_gene_final$Genbank_Count == 0 & 
                                Unified_gene_final$ensembl_Count == 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##3

NCBI_Ens_many_to_many_count <- Unified_gene_final[Unified_gene_final$NCBI_Count > 1 & 
                                Unified_gene_final$Genbank_Count == 0 & 
                                  Unified_gene_final$ensembl_Count > 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##2

######## Ens & XGC #########

Ens_XGC_one_to_one_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 0 & 
                             Unified_gene_final$Genbank_Count == 1 & 
                              Unified_gene_final$ensembl_Count == 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##442

Ens_XGC_Complicated <- Unified_gene_final[(Unified_gene_final$NCBI_Count == 0 & 
                                             Unified_gene_final$Genbank_Count >= 1 & 
                                             Unified_gene_final$ensembl_Count >= 1) & 
                                            !(Unified_gene_final$NCBI_Count == 0 & 
                                                Unified_gene_final$Genbank_Count == 1 & 
                                                Unified_gene_final$ensembl_Count == 1),] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##30

Ens_XGC_one_to_many_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 0 & 
                              Unified_gene_final$Genbank_Count == 1 & Unified_gene_final$ensembl_Count > 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##24

Ens_XGC_many_to_one_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 0 & 
                              Unified_gene_final$Genbank_Count > 1 & Unified_gene_final$ensembl_Count == 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##3

Ens_XGC_many_to_many_count <- Unified_gene_final[Unified_gene_final$NCBI_Count == 0 & 
                              Unified_gene_final$Genbank_Count > 1 & Unified_gene_final$ensembl_Count > 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
##3

######## Unique Models ######

Unique_NCBI_Models <- Unified_gene_final[Unified_gene_final$NCBI_Count >= 1 & 
                       Unified_gene_final$Genbank_Count == 0 & 
                        Unified_gene_final$ensembl_Count == 0,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
###5552

Unique_UCB_Models <- Unified_gene_final[Unified_gene_final$NCBI_Count == 0 & 
                      Unified_gene_final$Genbank_Count >= 1 & 
                       Unified_gene_final$ensembl_Count == 0,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
###3168

Unique_Ens_Models <- Unified_gene_final[Unified_gene_final$NCBI_Count == 0 & 
                      Unified_gene_final$Genbank_Count == 0 & 
                        Unified_gene_final$ensembl_Count >= 1,] %>% distinct(.keep_all = TRUE) %>% summarise(count = n(),.groups = "keep")
###4267


