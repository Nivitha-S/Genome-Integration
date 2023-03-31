library(dplyr)
library(stringr)

gtf_raw <- read.delim("Final_NCBI_Ens_XGC_Unified_Model.gtf",sep = "\t",header = FALSE)

gtf_raw$V10 <- sapply(gtf_raw$V9,function(x) str_extract(x,pattern = "[[:space:]]*\\[provisional.*\\]"))

#gtf_raw$V9 <- sapply(gtf_raw$V9,function(x) gsub("[[:space:]]*\\[provisional.*\\]","_provisional",x))

gtf_raw$V9 <- sapply(gtf_raw$V9, function(x) gsub("[[:space:]]*\\[provisional*","_provisional",x))

gtf_raw$V9 <- sapply(gtf_raw$V9, function(x) gsub("\\]","",x))

gtf_raw$V9 <- sapply(gtf_raw$V9, function(x) gsub("[[:space:]]*provisional*","_provisional",x))

filtered_ids <- gtf_raw %>% filter(!is.na(V10))

gtf_raw$V9 <- sapply(gtf_raw$V9,function(x) 
  ifelse(str_detect(x,pattern = "gene_name[[:space:]]+;"),
    gsub("[[:space:]]*gene_name[[:space:]]*;","",x),x))

gtf_raw$V9 <- sapply(gtf_raw$V9,function(x) 
  ifelse(str_detect(x,pattern = "transcript_name[[:space:]]+;"),
         gsub("[[:space:]]*transcript_name[[:space:]]*;","",x),x))

#####################################################################

gtf_raw$V9 <- sapply(gtf_raw$V9,function(x) 
ifelse(str_detect(x,pattern = "[[:space:]]*gene_id[[:space:]]*(.)*;"),
gsub(str_split(str_extract(x,
pattern = "[[:space:]]*gene_id[[:space:]]*(.)*;"),pattern = ";")[[1]][1],
paste0(" gene_id ",'"',gsub("[[:space:]]*gene_id[[:space:]]*","",
str_split(str_extract(x,
pattern = "[[:space:]]*gene_id[[:space:]]*(.)*;"),pattern = ";")[[1]][1]),'"'),x),x))

gtf_raw$V9 <- sapply(gtf_raw$V9,function(x) 
ifelse(str_detect(x,pattern = "[[:space:]]*gene_name[[:space:]]*(.)*;"),
gsub(str_split(str_extract(x,
pattern = "[[:space:]]*gene_name[[:space:]]*(.)*;"),pattern = ";")[[1]][1],
paste0(" gene_name ",'"',gsub("[[:space:]]*gene_name[[:space:]]*","",
str_split(str_extract(x,
pattern = "[[:space:]]*gene_name[[:space:]]*(.)*;"),pattern = ";")[[1]][1]),'"'),x),x))

gtf_raw$V9 <- sapply(gtf_raw$V9,function(x) 
ifelse(str_detect(x,pattern = "[[:space:]]*transcript_id[[:space:]]*(.)*;"),
gsub(str_split(str_extract(x,
pattern = "[[:space:]]*transcript_id[[:space:]]*(.)*;"),pattern = ";")[[1]][1],
paste0(" transcript_id ",'"',gsub("[[:space:]]*transcript_id[[:space:]]*","",
str_split(str_extract(x,
pattern = "[[:space:]]*transcript_id[[:space:]]*(.)*;"),pattern = ";")[[1]][1]),'"'),x),x))

gtf_raw$V9 <- sapply(gtf_raw$V9,function(x) 
ifelse(str_detect(x,pattern = "[[:space:]]*transcript_name[[:space:]]*(.)*;"),
gsub(str_split(str_extract(x,
pattern = "[[:space:]]*transcript_name[[:space:]]*(.)*;"),pattern = ";")[[1]][1],
paste0(" transcript_name ",'"',gsub("[[:space:]]*transcript_name[[:space:]]*","",
str_split(str_extract(x,
pattern = "[[:space:]]*transcript_name[[:space:]]*(.)*;"),pattern = ";")[[1]][1]),'"'),x),x))

gtf_raw <- gtf_raw %>% select(V1,V2,V3,V4,V5,V6,V7,V8,V9)


gtf_raw$V10 <- sapply(gtf_raw$V9,function(x)
  
  str_extract(x,pattern = "\\|")
  
)

gtf_raw_d <- gtf_raw %>% filter(!is.na(V10))

gtf_raw_df <- gtf_raw %>% filter(is.na(V10))

gtf_raw_df <- gtf_raw_df %>% select(V1,V2,V3,V4,V5,V6,V7,V8,V9)

write.table(gtf_raw_df,file = "Final_Mod_Unified_F.gtf",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

