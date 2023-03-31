
library(dplyr)
library(stringr)

Sort_Unified_F <- read.delim("Final_Modified_NCBI_Ens_XGC_Unified_Model.gff3",sep = "\t",header = FALSE)

names(Sort_Unified_F) <- c("chrloc","source","type","start","end","score","strand","phase","attributes")

Final_grped_Data <- read.delim("Final_grped_Data.gff3",sep = "\t",header = FALSE)

names(Final_grped_Data) <- c("chrloc","source","type","start","end","score","strand","phase","modified_attributes","attributes","gp_mod_attr")

gene_Final_grped_Data <- Final_grped_Data %>% filter(type == "gene")

####### XB-ID >Entrez-ID>gene symbol>Ensembl-ID ###########

gene_Sort_Unified_F <- Sort_Unified_F %>% filter(type == "gene")

gene_Sort_Unified_F$geneid <- sapply(gene_Sort_Unified_F$attributes, function(x)
  
  ifelse(str_detect(x,pattern = "geneid="),
  gsub("geneid=","",str_extract(x,pattern = "geneid=.+")),
  gsub("gene=","",str_split(str_extract(x,pattern = "gene=.+;"),pattern = ";")[[1]][1])

))

gene_Sort_Unified_F$XB_ID <- sapply(gene_Sort_Unified_F$attributes, function(x)
  
  gsub("ID=","",
  str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1])
  
)

gene_Sort_Unified_F$Ensembl_ID <- sapply(gene_Sort_Unified_F$attributes, function(x)
  
  ifelse(str_detect(x,pattern = "geneid="),
  gsub("geneid=","",str_extract(x,pattern = "geneid=.+")),
  gsub("gene_id=","",
  str_split(str_extract(x,pattern = "gene_id=.+;"),pattern = ";")[[1]][1]))
  
)

gene_Sort_Unified_F$Dbxref <- sapply(gene_Sort_Unified_F$attributes, function(x)
  
  gsub("Dbxref=","",
  str_split(str_extract(x,pattern = "Dbxref=.+;"),pattern = ";")[[1]][1])

)

gene_Sort_Unified_F$Entrez_ID <- sapply(gene_Sort_Unified_F$Dbxref, function(x)
  
  ifelse(str_detect(x,pattern = "GeneID:"),gsub("GeneID:","",
  str_split(x,pattern = ",")[[1]][1]),NA)
  
)

gene_Sort_Unified_F$XB_Page_ID <- sapply(gene_Sort_Unified_F$Dbxref, function(x)
  
  ifelse(str_detect(x,pattern = "Xenbase:"),gsub("Xenbase:","",
  str_split(x,pattern = ",")[[1]][2]),NA)
  
)


xref_mappings_data <- xref_mappings_d %>% filter(!is.na(Entrez_ID) | !is.na(geneid) | !is.na(Ensembl_ID) | !is.na(XB_Page_ID))

xref_mappings_data <- xref_mappings_data[,c("XB_ID","Entrez_ID","geneid","Ensembl_ID","XB_Page_ID")] %>% distinct(.keep_all = TRUE)

xref_mappings_data$gene_id <- sapply(xref_mappings_data$geneid, function(x)
  
  paste0(unique(str_split(x,pattern = ";")[[1]]),collapse = ";")

)

xref_mappings_data$Entrez_ID <- sapply(xref_mappings_data$Entrez_ID, function(x)
  
  paste0(unique(str_split(x,pattern = ";")[[1]]),collapse = ";")
  
  
)

xref_mappings_data$Ensembl_ID <- sapply(xref_mappings_data$Ensembl_ID, function(x)
  
  paste0(unique(str_split(x,pattern = ";")[[1]]),collapse = ";")
  
)

xref_mappings_data <- xref_mappings_data[,c("XB_ID","Entrez_ID","gene_id","Ensembl_ID","XB_Page_ID")] %>% distinct(.keep_all = TRUE)

write.csv(xref_mappings_data,file = "xref_mappings.csv",quote = FALSE,row.names = FALSE)

#########################  read in grped data ############################


