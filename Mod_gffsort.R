
Final_grped_Data <- read.delim("Final_grped_Data.gff3",sep = "\t",header = FALSE)

names(Final_grped_Data) <- c("chrloc","source","type","start","end","score","strand","phase","modified_attributes","attributes","gp_mod_attr")

Final_grped_df <- Final_grped_Data %>% mutate(
  
  Overlap_Value = end - start
  
)

#Unified_ls <- list()
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

filtered_gene_grped_Data <- Final_grped_df %>% filter(type == "gene" | 
type == "pseudogene" | type == "ncRNA_gene") %>% 
group_by(chrloc) %>% 
arrange(start,end,.by_group = TRUE) %>% 
select(chrloc,strand,gp_mod_attr) %>%  
distinct(.keep_all = TRUE)  

filtered_gene_grped_Data$rank <- 1: nrow(filtered_gene_grped_Data)

rank_d <- filtered_gene_grped_Data %>% select(chrloc,strand,gp_mod_attr,rank)

sorted_final_grped_data <- Final_grped_df %>% left_join(rank_d,by = c("chrloc","strand","gp_mod_attr"))

sorted_grped_data_f <- sorted_final_grped_data %>% arrange(rank) %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes,gp_mod_attr,Overlap_Value) 

UD_s <- str_sort(unique(sorted_grped_data_f$chrloc),numeric = TRUE)

rank_d <- data.frame(chr = UD_s, rank = 1:77)

Sorted_Unified_Data <- sorted_grped_data_f %>% left_join(rank_d,by = c("chrloc" = "chr"))

Sort_Unified_F <- Sorted_Unified_Data %>% arrange(rank) %>% 
select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes,gp_mod_attr,Overlap_Value) %>%
distinct(.keep_all = TRUE)

############# correcting nodal5 entries and sox17a entries ##############

g_sox17a_NCBI_data <- Sort_Unified_F %>% filter(modified_attributes == "sox17a")

g_sox17b_Ens_data <- Sort_Unified_F %>% filter(modified_attributes == "sox17b.1")

g_sox17a_XGC_data <- Sort_Unified_F %>% filter(start == 115156640 & end == 115158916)

filtered_data <- Sort_Unified_F %>% filter(gp_mod_attr %in% c("LOC100490072","LOC101734157","LOC101734231") | (modified_attributes %in% "nodal5.2" & gp_mod_attr %in% "nodal5"))

ab <- Sort_Unified_F %>% filter(gp_mod_attr %in% c("LOC100490072","LOC101734157","LOC101734231") | (modified_attributes %in% "nodal5.2" & gp_mod_attr %in% "nodal5"))

write.table(Sort_Unified_F,file = "Sort_Unified_F.gff3",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

###################################################################

Sort_Unified_F <- read.delim("Sort_Unified_F.gff3",sep = "\t",header = FALSE)

Sort_Unified_F <- Sort_Unified_F %>% distinct(.keep_all = TRUE)

names(Sort_Unified_F) <- c("chrloc","source","type","start","end","score","strand","phase","modified_attributes","attributes","gp_mod_attr","Overlap_Value")

########### passing each chromosome gene models and further grouping them according to
### their grouping variable (XGC gp variable)/declaring a function and re-routing the 
### filtered models according to chromosome location accordingly

gene_mrna_count <- 0
exon_count <- 0
cds_count <- 0
V_gene_segment_count <- 0
C_gene_segment_count <- 0
match_count <- 0
or_count <- 0
D_loop_count <- 0
Non_trans_count <- 0
other_ls <- list()
#Unified_F_ls <- list()
#Unified_ls <- list()

unified_NCBI_Ens_XGC <- function(NCBI_Ens_XGC_Data){
  
  Unified_F_ls <- list()
  
  rank_data <- data.frame(type = c("gene","ncRNA_gene","pseudogene","mRNA","tRNA","rRNA","lnc_RNA","snoRNA","snRNA","scRNA",
                                   "miRNA","Y_RNA","ncRNA","guide_RNA","transcript","primary_transcript","pseudogenic_transcript","exon","CDS",
                                   "V_gene_segment","C_gene_segment","match","origin_of_replication","D_loop"),rank = 1:24)
  
  for (j in 1 : length(unique(NCBI_Ens_XGC_Data$gp_mod_attr))) {
    
    print(j)
    
    part_gene_d <- NCBI_Ens_XGC_Data %>% filter(gp_mod_attr == unique(NCBI_Ens_XGC_Data$gp_mod_attr)[j])
    
    chr_d <- unique(part_gene_d$chrloc)  
    
    Unified_M_ls <- list()
    
    for (k in 1 : length(unique(part_gene_d$strand))) {
      
      part_gene <- part_gene_d %>% filter(
        
        strand == unique(part_gene_d$strand)[k]
        
      )
      
      if(unique(part_gene$strand) == "+"){
        
        if((unique(part_gene$gp_mod_attr) == "LOC101734477" | unique(part_gene$gp_mod_attr) == "LOC116406437" | unique(part_gene$gp_mod_attr) == "smyd5")){
          
          if(unique(part_gene$gp_mod_attr) == "LOC101734477"){
            
            g_LOC101734477_ls <- list()
            
            part_gene_LOC101734477 <- part_gene %>% filter(modified_attributes == "LOC101734477") %>% distinct(.keep_all = TRUE)
            
            part_gene_f_c <- part_gene_LOC101734477 %>% filter(!(type == "gene" | type == "mRNA" | type == "transcript")) %>% 
              group_by(type) %>% summarise(count = n())
            
            gene_mrna_count <- gene_mrna_count + 1
            exon_i_count <- exon_count + 1
            exon_count <- exon_count + part_gene_f_c[part_gene_f_c$type == "exon",]$count
            cds_i_count <- cds_count + 1
            cds_count <- cds_count + part_gene_f_c[part_gene_f_c$type == "CDS",]$count   
            
            g_LOC101734477_ls[[length(g_LOC101734477_ls)+1]] <- gene_mrna_count
            g_LOC101734477_ls[[length(g_LOC101734477_ls)+1]] <- paste0(exon_i_count,";",exon_count)
            g_LOC101734477_ls[[length(g_LOC101734477_ls)+1]] <- paste0(cds_i_count,";",cds_count)
            
            Unified_M_ls[[k]] <- part_gene_LOC101734477 %>% select(chrloc,source,type,start,end,score,
                                                                   strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          else if(unique(part_gene$gp_mod_attr) == "LOC116406437"){
            
            g_LOC116406437_ls <- list()
            
            part_gene_LOC116406437 <- part_gene %>% filter(modified_attributes == "LOC116406437") %>% distinct(.keep_all = TRUE)
            
            part_gene_f_c <- part_gene_LOC116406437 %>% filter(!(type == "gene" | type == "mRNA" | type == "transcript")) %>% 
              group_by(type) %>% summarise(count = n())
            
            gene_mrna_count <- gene_mrna_count + 1
            
            exon_i_count <- exon_count + 1
            exon_count <- exon_count + part_gene_f_c[part_gene_f_c$type == "exon",]$count
            cds_i_count <- cds_count + 1
            cds_count <- cds_count + part_gene_f_c[part_gene_f_c$type == "CDS",]$count   
            
            g_LOC116406437_ls[[length(g_LOC116406437_ls)+1]] <- gene_mrna_count
            g_LOC116406437_ls[[length(g_LOC116406437_ls)+1]] <- paste0(exon_i_count,";",exon_count)
            g_LOC116406437_ls[[length(g_LOC116406437_ls)+1]] <- paste0(cds_i_count,";",cds_count)
            
            Unified_M_ls[[k]] <- part_gene_LOC116406437 %>% select(chrloc,source,type,start,end,score,
              strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          else if(unique(part_gene$gp_mod_attr) == "smyd5"){
            
            smyd5_count_ls <- list()
            
            part_gene_smyd5 <- part_gene %>% filter(modified_attributes == "smyd5") %>% distinct(.keep_all = TRUE)
            
            part_gene_f_c <- part_gene_smyd5 %>% filter(!(type == "gene" | type == "mRNA" | type == "transcript")) %>% 
              group_by(type) %>% summarise(count = n())
            
            gene_mrna_count <- gene_mrna_count + 1
            exon_i_count <- exon_count + 1
            exon_count <- exon_count + part_gene_f_c[part_gene_f_c$type == "exon",]$count
            cds_i_count <- cds_count + 1
            cds_count <- cds_count + part_gene_f_c[part_gene_f_c$type == "CDS",]$count   
            
            smyd5_count_ls[[length(smyd5_count_ls)+1]] <- gene_mrna_count
            smyd5_count_ls[[length(smyd5_count_ls)+1]] <- paste0(exon_i_count,";",exon_count)
            smyd5_count_ls[[length(smyd5_count_ls)+1]] <- paste0(cds_i_count,";",cds_count)
            
            Unified_M_ls[[k]] <- part_gene_smyd5 %>% select(chrloc,source,type,start,end,score,
                                                            strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)           
            
          }
          
        }
        
        else if(((nrow(part_gene[part_gene$type == "gene",]) >= 2) &
                 (length(unique(part_gene[part_gene$type == "gene",]$type) %in% 
                         c("gene")) == 1 ) &
                 nrow(part_gene[part_gene$type == "tRNA",]) == 0) |
                ((nrow(part_gene[part_gene$type == "pseudogene",]) >= 2) & 
                 (length(unique(part_gene[part_gene$type == "pseudogene",]$type) %in%  
                         c("pseudogene")) == 1))){
          
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
          
          Non_gene_f_c <- Non_gene_id %>% 
            group_by(type) %>% summarise(count = n())
          
          gene_mrna_count <- gene_mrna_count + 1
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
            
            exon_i_count <- exon_count + 1  
            exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
            
            cds_i_count <- cds_count + 1   
            cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
            
            match_i_count <- match_count + 1
            match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
            
            V_gene_segment_i_count <- V_gene_segment_count + 1
            V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
            
            C_gene_segment_i_count <- C_gene_segment_count + 1
            C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
            
            or_i_count <- or_count + 1
            or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
            
            D_loop_i_count <- D_loop_count + 1
            D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
          }
          
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
          
          part_gene_id_df$attributes <- paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,
          ";",paste0(unique(part_gene_id$attributes),collapse = ";"))
          
          part_gene_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
          
          if(nrow(part_mRNA_id_df) == 1){
            
            part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
              
              gsub("ID=","transID=",x)
              
            )
            
            part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
              
              gsub("Parent=","Origin=",x)
              
            )
            
            part_mRNA_id_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                 gene_mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                 gene_mrna_count,";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))
            
            part_mRNA_id_df$start <- min(part_gene$start)
            
            part_mRNA_id_df$end <- max(part_gene$end)
            
            part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
            
            parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
            
            if(nrow(Non_gene_id) >= 1){
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
              
              Non_gene_id$type_count <- NA
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                
              }
              
              Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                
                paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                
              )
              
              Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
              
              Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                
                paste0("ID=XB-",x)
                
              )
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                
                gsub("ID=","OriginID=",x)
                
                
              )
              
              Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
              
              Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
              
            }
            
            Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
            
            Unified_M_ls[[k]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
              select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE) 
            
          }
          
          else if(nrow(part_mRNA_id_df) < 1){
            
            part_mRNA_id_df <- part_gene_id_df 
            
            part_mRNA_id_df$type <- "transcript"
            
            part_mRNA_id_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                 gene_mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                 gene_mrna_count,";")
            
            part_mRNA_id_df$start <- part_gene_id_df$start
            
            part_mRNA_id_df$end <- part_gene_id_df$end 
            
            part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
            
            parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
            
            if(nrow(Non_gene_id) >= 1){
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
              
              Non_gene_id$type_count <- NA
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                
              }
              
              Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                
                paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                
              )
              
              Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
              
              Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                
                paste0("ID=XB-",x)
                
              )
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                
                gsub("ID=","OriginID=",x)
                
                
              )
              
              Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
              
              Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
              
            }
            
            Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
            
            Unified_M_ls[[k]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
              select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE) 
            
          }
          
        }
        
        else if((nrow(part_gene[part_gene$type == "gene" | 
                                part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) < 2 &
                 nrow(part_gene[part_gene$type == "tRNA",]) == 0) | 
                (nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
                 nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1 & 
                 nrow(part_gene[part_gene$type == "ncRNA_gene" | 
                                part_gene$type == "tRNA",]) == 0)){
          
          if(nrow(part_gene[part_gene$type == "gene" | 
                            part_gene$type == "pseudogene" |
                            part_gene$type == "ncRNA_gene",]) < 2){
            
            part_gene_df <- part_gene %>% filter(type == "gene"| type == "pseudogene"| type == "ncRNA_gene") %>%
              select(chrloc,source,type,start,end,score,phase,strand,modified_attributes,attributes) 
            
            if(nrow(part_gene_df) == 1){
              
              part_trans_M <- part_gene %>% filter(type == "mRNA" | 
                                                     type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                     type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                     type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                     type == "transcript" | type == "pseudogenic_transcript" | 
                                                     type == "primary_transcript")
              
              part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              ) %>% distinct(.keep_all = TRUE)
              
              part_gene_df$start <- min(part_gene$start)
              
              part_gene_df$end <- max(part_gene$end)
              
              part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
              
              Non_gene_id <- part_gene  %>% filter(!(type %in% part_gene_df$type | 
                                                       type %in% part_trans_M$type ))  %>% select(
                                                         chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                       ) %>% distinct(.keep_all = TRUE)
              
              Non_gene_f_c <- Non_gene_id %>% 
                group_by(type) %>% summarise(count = n())
              
              gene_mrna_count <- gene_mrna_count + 1
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                exon_i_count <- exon_count + 1  
                exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                cds_i_count <- cds_count + 1    
                cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                match_i_count <- match_count + 1
                match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                V_gene_segment_i_count <- V_gene_segment_count + 1
                V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                C_gene_segment_i_count <- C_gene_segment_count + 1
                C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                or_i_count <- or_count + 1
                or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                D_loop_i_count <- D_loop_count + 1
                D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                
              }
              
              part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",
                                                                                         paste0("ID=XBXT10g",
                                                                                                paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";gID="),x))
              
              if(nrow(part_trans_M_df) >= 1){
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) 
                  gsub("Parent=",paste0("Parent=XBXT10g",
                                        paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                        gene_mrna_count,";Origin="),x))
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",
                                                                                                 paste0("ID=XB-mRNA",
                                                                                                        paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";transID="),x))
                
                part_trans_M_df$start <- part_gene_df$start
                
                part_trans_M_df$end <- part_gene_df$end 
                
                part_trans_M_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                    gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    
                    gsub("ID=","OriginID=",x)
                    
                    
                  )
                  
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                         strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
              else if(nrow(part_trans_M_df) < 1){
                
                part_trans_M_df <- part_gene_df 
                
                part_trans_M_df$type <- "transcript"
                
                part_trans_M_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                     gene_mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                     gene_mrna_count,";")
                
                part_trans_M_df$start <- part_gene_df$start
                
                part_trans_M_df$end <- part_gene_df$end 
                
                part_trans_M_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                    gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    
                    gsub("ID=","OriginID=",x)
                    
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                         strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
                
              }
              
              Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                        strand,phase,modified_attributes,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            else if(nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) < 1 ){
              
              if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) < 1 )){
                
                Non_gene_f_c <- part_gene %>% 
                  group_by(type) %>% summarise(count = n())
                
                gene_mrna_count <- gene_mrna_count + 1
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                org_parent_id <- gsub("Parent=","",str_split(str_extract(part_gene$attributes[1],pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
                
                mod_gene_id <- part_gene[1,]
                
                mod_mRNA_id <- part_gene[1,]
                
                mod_gene_id$type <- "gene"
                mod_mRNA_id$type <- "mRNA"
                
                mod_gene_id$start <- min(part_gene$start)
                mod_gene_id$end <- max(part_gene$end)
                
                mod_mRNA_id$start <- min(part_gene$start)
                mod_mRNA_id$end <- max(part_gene$end)
                
                gene_id <- paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                  gene_mrna_count,";gbkey=gene")
                
                parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
                
                mRNA_id <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                  gene_mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                  gene_mrna_count,";")
                
                pid <- gsub("ID=","",str_split(str_extract(mRNA_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
                
                part_gene$attributes <- sapply(part_gene$attributes,function(x) 
                  
                  gsub("Parent=",paste0("Parent=",pid,";Origin="),x)
                  
                )
                
                part_gene$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  part_gene[part_gene$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  part_gene[part_gene$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  part_gene[part_gene$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  part_gene[part_gene$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  part_gene[part_gene$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  part_gene[part_gene$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  part_gene[part_gene$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                part_gene$type_mod <- sapply(part_gene$type_count, function(x)
                  
                  paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                  
                )
                
                part_gene$type_f <- paste0(casefold(part_gene$type),part_gene$type_mod)
                
                part_gene$type_final <- sapply(part_gene$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                part_gene$attributes <- sapply(part_gene$attributes, function(x)
                  
                  
                  gsub("ID=","OriginID=",x)
                  
                  
                )
                
                part_gene$attributes <- paste0(part_gene$type_final,part_gene$attributes)
                
                part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                mod_gene_id$attributes <- gene_id
                mod_mRNA_id$attributes <- mRNA_id
                
                mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                part_gene$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                mod_gene_id <- mod_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                mod_mRNA_id <- mod_mRNA_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,part_gene) %>% select(chrloc,source,type,start,end,score,
                                                                                            strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
              else if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) >= 1 )){
                
                mod_mRNA <- part_gene %>% filter(type == "mRNA")
                
                mod_mRNA_id <- part_gene[which.max(part_gene$type == "mRNA"),]
                
                org_parent_id <- gsub("Parent=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
                
                mod_gene_id <- mod_mRNA_id       
                
                mod_gene_id$type <- "gene"
                mod_mRNA_id$type <- "mRNA"
                
                mod_gene_id$start <- min(part_gene$start)
                mod_gene_id$end <- max(part_gene$end)
                
                mod_mRNA_id$start <- min(part_gene$start)
                mod_mRNA_id$end <- max(part_gene$end)
                
                Non_gene_id <- part_gene %>% filter(!(type %in% mod_gene_id$type | type %in% mod_mRNA$type))
                
                Non_gene_f_c <- Non_gene_id %>% 
                  group_by(type) %>% summarise(count = n())
                
                gene_mrna_count <- gene_mrna_count + 1
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                gene_id <- paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                  gene_mrna_count,";gbkey=gene")
                
                parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
                
                mod_gene_id$attributes <- gene_id
                
                mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("ID=",paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                                                                      gene_mrna_count,";gID="),x))
                
                p_id <- gsub("ID=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "ID.+;"),pattern = ";")[[1]][1])
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                  
                  gsub("Parent=",paste0("Parent=",p_id,";Origin="),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  gsub("ID=","OriginID=",x)
                  
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                mod_gene_id <- mod_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                mod_mRNA_id <- mod_mRNA_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,Non_gene_id) %>% select(chrloc,source,type,start,end,score,
                                                                                              strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)        
                
              }
              
              if(unique(part_gene$source) == "NCBI"){
                
                if(nrow(part_gene[part_gene$type == "mRNA" | part_gene$type == "tRNA" |
                                  part_gene$type == "rRNA" | part_gene$type == "lnc_RNA" | part_gene$type == "snoRNA" |
                                  part_gene$type == "snRNA" | part_gene$type == "scRNA" | part_gene$type == "miRNA" |
                                  part_gene$type == "Y_RNA" | part_gene$type == "ncRNA" | part_gene$type == "guide_RNA" |
                                  part_gene$type == "transcript" | part_gene$type == "primary_transcript" | part_gene$type == "pseudogenic_transcript",]) >= 1){
                  
                  Non_trans_count <- Non_trans_count + 1
                  
                  part_gene$attributes <- sapply(part_gene$attributes, function(x)
                    
                    gsub("ID=",paste0("ID=XB-trans",paste0(rep(0,6-nchar(Non_trans_count)),collapse = ""),Non_trans_count,";OriginID="),x)
                    
                  )
                  
                  part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                    strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)   
                  
                } 
                
                else if(nrow(part_gene[part_gene$type == "exon" | part_gene$type == "CDS" |
                                       part_gene$type == "match" | part_gene$type == "V_gene_segment" |
                                       part_gene$type == "C_gene_segment"|part_gene$type == "origin_of_replication" |
                                       part_gene$type == "D_loop",]) >= 1){
                  
                  Non_gene_f_c <- part_gene %>%
                    group_by(type) %>% summarise(count = n())
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    exon_i_count <- exon_count + 1
                    exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    cds_i_count <- cds_count + 1    
                    cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    match_i_count <- match_count + 1
                    match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    V_gene_segment_i_count <- V_gene_segment_count + 1
                    V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    C_gene_segment_i_count <- C_gene_segment_count + 1
                    C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    or_i_count <- or_count + 1
                    or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    D_loop_i_count <- D_loop_count + 1
                    D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                    
                  }
                  
                  part_gene$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    part_gene[part_gene$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    part_gene[part_gene$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    part_gene[part_gene$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    part_gene[part_gene$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    part_gene[part_gene$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    part_gene[part_gene$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    part_gene[part_gene$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  part_gene$type_mod <- sapply(part_gene$type_count, function(x)
                    
                    paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                    
                  )
                  
                  part_gene$type_f <- paste0(casefold(part_gene$type),part_gene$type_mod)
                  
                  part_gene$type_final <- sapply(part_gene$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  part_gene$attributes <- sapply(part_gene$attributes, function(x)
                    
                    
                    gsub("ID=","OriginID=",x)
                    
                    
                  )
                  
                  part_gene$attributes <- paste0(part_gene$type_final,part_gene$attributes)
                  
                  part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                    strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)        
                }
                
              }
              
            }
            
            Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                      strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          else if(nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
                  nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1){
            
            Unified_ls <- list()
            
            part_gene_No <- part_gene %>% filter(type == "gene" | type == "pseudogene") %>% 
              group_by(type) %>% summarise( Count = n(), source_t = paste0(source,collapse = ";"),
                                            mod_attr = paste0(unique(modified_attributes),collapse = ";") )
            
            for(i in 1 : length(part_gene_No$Count)){
              
              if(part_gene_No$Count[i] == 1){
                
                gene_mrna_count <- gene_mrna_count + 1
                
                source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
                
                part_gene_Mod_attr <- unique(str_split(part_gene_No$mod_attr[i],pattern = ";")[[1]])
                
                part_gene_type <- part_gene_No[i,]$type    
                
                part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
                
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
                
                part_gene_M_df$start <- min(part_gene_M$start)
                part_gene_M_df$end <- max(part_gene_M$end)
                part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                
                part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",
                                                                                               paste0("ID=XBXT10g",
                                                                                                      paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";gID="),x))
                
                Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | 
                                                           type %in% part_trans_M$type ))  %>% select(
                                                             chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                           ) %>% distinct(.keep_all = TRUE)
                
                Non_gene_f_c <- Non_gene_id %>% 
                  group_by(type) %>% summarise(count = n())
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1  
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                if(nrow(part_trans_M_df) == 1){
                  
                  part_trans_M_df$start <- min(part_gene_M$start)
                  part_trans_M_df$end <- max(part_gene_M$end)
                  part_trans_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                  
                  part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",
                                                                                                   paste0("Parent=XBXT10g",
                                                                                                          paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";Origin="),x))
                  
                  part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",
                                                                                                   paste0("ID=XB-mRNA",
                                                                                                          paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";transID="),x))
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  if(nrow(Non_gene_id) >= 1){  
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                    select(chrloc,source,type,start,end,score,                                                                           
                           strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)
                  
                }
                
                else if(nrow(part_trans_M_df) != 1){
                  
                  part_trans_M_df <- part_gene_M_df
                  
                  part_trans_M_df$type <- "transcript"
                  
                  part_trans_M_df$start <- min(part_gene_M$start)
                  part_trans_M_df$end <- max(part_gene_M$end)
                  part_trans_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                  
                  part_trans_M_df$attributes <- 
                    paste0("ID=XB-mRNA",
                           paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";","Parent=XBXT10g",
                           paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";")
                  
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                    select(chrloc,source,type,start,end,score,                                                                           
                           strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)
                  
                }
                
                Unified_ls[[i]] <- part_gene_M %>% select(chrloc,source,type,start,end,score,
                                                          strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }   
              
              else if(part_gene_No$Count[i] > 1){
                
                gene_mrna_count <- gene_mrna_count + 1
                
                source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
                
                part_gene_type <- part_gene_No[i,]$type    
                
                part_gene_Mod_attr <- unique(part_gene[
                  part_gene$type == part_gene_type,]$modified_attributes)
                
                part_gene_Mod <- part_gene %>% filter((source %in% source_tot ) & (modified_attributes %in% part_gene_Mod_attr))
                
                part_gene_id <- part_gene_Mod %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
                
                part_mRNA_id <- part_gene_Mod %>% filter(type == "mRNA" | type == "tRNA"| 
                                                           type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                           type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                           type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                           type == "transcript" | type == "pseudogenic_transcript" | 
                                                           type == "primary_transcript")
                
                Non_gene_id <-  part_gene_Mod %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
                  select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
                
                Non_gene_f_c <- Non_gene_id %>% 
                  group_by(type) %>% summarise(count = n())
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                
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
                
                part_gene_id_df$attributes <- 
                  paste0("ID=XBXT10g",
                         paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",paste0(unique(part_gene_id$attributes),collapse = ";"))    
                
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
                  
                  part_mRNA_id_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";","Parent=XBXT10g",
                                                       paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",
                                                       paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene_M <- rbind.data.frame(part_gene_M_df,part_mRNA_id_df,Non_gene_id) %>% 
                    select(chrloc,source,type,start,end,score,                                                                           
                           strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)
                  
                  
                }
                
                else if(nrow(part_mRNA_id_df) != 1){
                  
                  part_mRNA_id_df <- part_gene_id_df
                  
                  part_mRNA_id_df$type <- "transcript"
                  
                  part_mRNA_id_df$start <- min(part_gene_Mod$start)
                  part_mRNA_id_df$end <- max(part_gene_Mod$end)
                  part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
                  
                  part_mRNA_id_df$attributes <- 
                    paste0("ID=XB-mRNA",
                           paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";","Parent=XBXT10g",
                           paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";")
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene_M <- rbind.data.frame(part_gene_M_df,part_mRNA_id_df,Non_gene_id) %>% 
                    select(chrloc,source,type,start,end,score,                                                                           
                           strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)
                  
                }
                
                Unified_ls[[i]] <- part_gene_M %>%
                  select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
            }
            
            Unified_M_ls[[k]] <- do.call("rbind.data.frame",Unified_ls) 
            
          } 
          
        }
        
        else if(((nrow(part_gene[part_gene$type == "ncRNA_gene",]) > 1) &
                 (length(unique(part_gene[part_gene$type == "ncRNA_gene",]$type) %in% 
                         c("ncRNA_gene")) == 1 )) | (nrow(part_gene[part_gene$type == "gene" | 
                                                                    part_gene$type == "pseudogene",]) >= 1 & 
                                                     nrow(part_gene[part_gene$type == "ncRNA_gene",]) >= 1) |
                (nrow(part_gene[part_gene$type == "tRNA",]) >= 1)) {
          
          other_ls[[length(other_ls)+ 1]] <- part_gene
          
          part_df <- part_gene %>% left_join(rank_data,by = c("type"))
          
          if(length(unique(part_df[part_df$type == "tRNA",]$type)) >= 1 | ((length(unique(part_df[part_df$type == "tRNA",]$type)) >= 1) & (length(unique(part_df[part_df$type == "ncRNA_gene" | part_df$type == "miRNA",]$type)) >= 1))){
            
            #part_df <- part_df %>% arrange(rank)
            
            part_df$ID <- sapply(part_df$attributes,function(x) gsub("ID=","",
            str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
            
            part_df$Parent <- sapply(part_df$attributes,function(x) gsub("Parent=","",
            str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))
            
            filtered_rna_gene_ids <- part_df %>% filter(type == "mRNA" |
            type == "tRNA"| type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" |
            type == "snRNA" | type == "scRNA" | type == "miRNA" | type == "Y_RNA" |
            type == "ncRNA" | type == "guide_RNA" | type == "transcript" | 
            type == "primary_transcript" |type == "pseudogenic_transcript") %>% select(ID,Parent)
            
            part_d <- part_df %>% left_join(filtered_rna_gene_ids,by = c("Parent" = "ID"))
            
            part_d <- part_d %>% mutate(
              
              transcript = ifelse(is.na(Parent.y),Parent,Parent.y),
              final_transcript = ifelse(is.na(transcript),ID,transcript)
              
              
            ) %>% arrange(final_transcript) %>% select(chrloc,source,type,start,end,
            score,strand,phase,modified_attributes,attributes,final_transcript,Overlap_Value)
            
            filtered_mod_attr_ls <- list()
            final_mod_attr_ls <- list()
            
            for (l in 1 : length(unique(part_d$final_transcript))) {
              
              filtered_mod_attr_ls[[l]] <- part_d %>% filter(final_transcript %in% unique(part_d$final_transcript)[l]) 
              
              gene_filtered_mod_attr <- filtered_mod_attr_ls[[l]] %>% filter(type == "gene" | 
              type == "pseudogene" | type == "ncRNA_gene")
              
              transcript_filtered_mod_attr <- filtered_mod_attr_ls[[l]] %>% filter(type == "mRNA" |
              type == "tRNA"| type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" |
              type == "snRNA" | type == "scRNA" | type == "miRNA" | type == "Y_RNA" |
              type == "ncRNA" | type == "guide_RNA" | type == "transcript" | 
              type == "primary_transcript" |type == "pseudogenic_transcript")
              
              
              gene_filtered_mod_attr_df <- gene_filtered_mod_attr[which.max(gene_filtered_mod_attr$Overlap_Value),] %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              )
              
              transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr[which.max(transcript_filtered_mod_attr$Overlap_Value),] %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              )
              
              Non_gene_id <- filtered_mod_attr_ls[[l]] %>% filter(!(type %in% gene_filtered_mod_attr$type | type %in% transcript_filtered_mod_attr$type)) %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              )
              
              Non_gene_f_c <- Non_gene_id %>% 
                group_by(type) %>% summarise(count = n())
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                exon_i_count <- exon_count + 1 
                exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                cds_i_count <- cds_count + 1   
                cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                match_i_count <- match_count + 1
                match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                V_gene_segment_i_count <- V_gene_segment_count + 1
                V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                C_gene_segment_i_count <- C_gene_segment_count + 1
                C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                or_i_count <- or_count + 1
                or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                D_loop_i_count <- D_loop_count + 1
                D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                
              }
              
              if(nrow(gene_filtered_mod_attr_df) == 1){
                
                gene_mrna_count <- gene_mrna_count + 1
                
                gene_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                gene_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                
                gene_filtered_mod_attr_df$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                
                gene_filtered_mod_attr$attributes <- sapply(gene_filtered_mod_attr$attributes, function(x)
                  
                  gsub("ID=","gID=",x)
                  
                )
                
                gene_filtered_mod_attr_df$attributes <- paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",paste0(unique(gene_filtered_mod_attr$attributes),collapse = ";"))
                
                if(nrow(transcript_filtered_mod_attr_df) == 1 ){
                  
                  transcript_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                  transcript_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                  
                  transcript_filtered_mod_attr_df$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                  
                  transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                    
                    gsub("ID=","transID=",x)
                    
                  )
                  
                  transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                    
                    gsub("Parent=","Origin=",x)
                    
                  )
                  
                  transcript_filtered_mod_attr_df$attributes <- paste0("ID=XB-rna",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",
                  "Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",
                   paste0(unique(transcript_filtered_mod_attr$attributes),collapse = ";"))
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(
                    transcript_filtered_mod_attr_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])   
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                  
                  final_mod_attr_ls[[l]] <- rbind.data.frame(gene_filtered_mod_attr_df,transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                    
                  chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                    
                  )
                  
                }
                
                else if(nrow(transcript_filtered_mod_attr_df) < 1){
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(gene_filtered_mod_attr_df$attributes,pattern = "ID=.+;"))[[1]][1])   
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                    
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                  
                  final_mod_attr_ls[[l]] <- rbind.data.frame(gene_filtered_mod_attr_df,transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                    
                    chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                    
                  )
                  
                }
                
              }
              
              if(nrow(gene_filtered_mod_attr_df) < 1){
                
                if( nrow(transcript_filtered_mod_attr_df) == 1 ){
                  
                  Non_trans_count <- Non_trans_count + 1
                  
                  transcript_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                  transcript_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                  
                  transcript_filtered_mod_attr_df$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                  
                  transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                    
                    gsub("ID=","transID=",x)
                    
                  )
                  
                  transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                    
                    gsub("Parent=","Origin=",x)
                    
                  )
                  
                  transcript_filtered_mod_attr_df$attributes <- 
                    
                    sapply(transcript_filtered_mod_attr_df$attributes, function(x)
                      
                      gsub("ID=",paste0("ID=XB-trans",paste0(rep(0,6-nchar(Non_trans_count)),collapse = ""),Non_trans_count,";transID="),x)    
                      
                    )
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(
                    transcript_filtered_mod_attr_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])   
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                  
                  final_mod_attr_ls[[l]] <- rbind.data.frame(transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                    
                    chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                    
                  )
                  
                }
                
                else if(nrow(transcript_filtered_mod_attr_df) < 1){
                  
                  #parent_id <- gsub("ID=","",str_split(str_extract(gene_filtered_mod_attr_df$attributes,pattern = "ID=.+;"))[[1]][1])   
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    #Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    #gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                    
                    #Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                    
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                  
                  final_mod_attr_ls[[l]] <- rbind.data.frame(Non_gene_id) %>% select(
                    
                    chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                    
                  )
                  
                }
                
              }  
              
            }
            
            final_mod_attr_df <- do.call("rbind.data.frame",final_mod_attr_ls)
            
            Unified_M_ls[[k]] <- final_mod_attr_df %>% select(chrloc,source,type,start,end,score,
             strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          if((length(unique(part_df[part_df$type == "ncRNA_gene" | part_df$type == "miRNA",]$type)) >= 1 & nrow(part_df[part_df$type == 'tRNA',]) == 0) | 
             (nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene",]) >= 1 & 
              nrow(part_gene[part_gene$type == "ncRNA_gene",]) > 1)){
            
            part_df <- part_df %>% arrange(rank)
            
            final_gb_ls <- list()
            
            if (nrow(part_df[part_df$source == "Genbank",]) >= 1){
              
              filtered_genbank_ls <- part_df %>% filter(source == "Genbank")
              
              filtered_genbank_ls <- filtered_genbank_ls %>% arrange(gp_mod_attr)
              
              filtered_gb_ls <- list()
              
              for (h in 1 : length(unique(filtered_genbank_ls$gp_mod_attr))) {
                
                gene_mrna_count <- gene_mrna_count + 1
                
                filtered_gb_ls[[h]] <- filtered_genbank_ls %>% filter(gp_mod_attr == unique(filtered_genbank_ls$gp_mod_attr)[h]) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                  
                )
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes,function(x) 
                  gsub("ID=",paste0("ID=XBXT10g",
                                    paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";gID="),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  gsub("ID=",paste0("ID=XB-mRNA",
                                    paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";transID="),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  gsub("Parent=",paste0("Parent=XBXT10g",
                                        paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";Origin="),x))
                
                if(nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]) >= 1){
                  
                  parent_id <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))  
                  
                }
                
                else if (nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]) < 1){
                  
                  
                  parent_id <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))  
                  
                }
                
                Non_gene_id <- filtered_gb_ls[[h]] %>% filter(!(type == "gene" | type == "mRNA"))
                
                Non_gene_f_c <- Non_gene_id %>% 
                  group_by(type) %>% summarise(count = n())
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1  
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x)
                    
                  )
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    
                    gsub("ID=","OriginID=",x)
                    
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                final_gb_ls[[h]] <- rbind.data.frame(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",],
                                                     filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",],
                                                     Non_gene_id) %>% select(
                                                       
                                                       chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes                                       
                                                       
                                                     )
                
              }
              
            }
            
            final_gb_df <- do.call("rbind.data.frame",final_gb_ls)
            
            part_df <- part_df %>% filter(!(source %in% "Genbank"))
            
            part_df$transcript <- sapply(part_df$attributes,function(x) 
              
              ifelse(str_detect(x,pattern = "Parent=.+"),
                     
                     gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1]),
                     
                     gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
              
            )
            
            part_df <- part_df %>% arrange(transcript)
            
            filtered_gene_ls <- list()
            filtered_M_ls <- list()
            
            for(k in 1: length(unique(part_df$transcript))){
              
              gene_mrna_count <- gene_mrna_count + 1
              
              filtered_gene_ls[[k]] <- part_df %>% filter(transcript %in% unique(part_df$transcript)[k]) %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              )
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",]$attributes,function(x) 
                gsub("ID=",paste0("ID=XBXT10g",
                                  paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";gID="),x))
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                      filtered_gene_ls[[k]]$type == "miRNA" | 
                                      filtered_gene_ls[[k]]$type == "rRNA" | 
                                      filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                      filtered_gene_ls[[k]]$type == "snoRNA" | 
                                      filtered_gene_ls[[k]]$type == "scRNA" | 
                                      filtered_gene_ls[[k]]$type == "Y_RNA" | 
                                      filtered_gene_ls[[k]]$type == "ncRNA" | 
                                      filtered_gene_ls[[k]]$type == "guide_RNA" |
                                      filtered_gene_ls[[k]]$type == "transcript" | 
                                      filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "miRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "rRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "snoRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "Y_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "guide_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "transcript" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes,function(x) 
                                                                                                                                                                                 gsub("Parent=",paste0("Parent=XBXT10g",
                                                                                                                                                                                                       paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";Origin="),x))
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | filtered_gene_ls[[k]]$type == "miRNA" | 
                                      filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                      filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                      filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" |
                                      filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                      filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes,function(x) 
                                                                                                                                                                                 gsub("ID=",paste0("ID=XB-rna",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";transID="),x))
              if(nrow(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                            filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                            filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                            filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                            filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                            filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]) >= 1){                                                                                                                                              
                
                parent_id <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                            filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                            filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                            filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                            filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                            filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
                
              }
              
              else if (nrow(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                  filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                  filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                  filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                  filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                  filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]) < 1){
                
                parent_id <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
                
                
              }
              
              Non_gene_id <- filtered_gene_ls[[k]] %>% filter(!(type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | type == "snRNA" | 
                                                                  type == "scRNA" | type == "miRNA" | type == "transcript" | type == "primary_transcript" |
                                                                  type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "ncRNA_gene" | type == "mRNA" | type == "gene" | type == "pseudogene"
                                                                
              ))
              
              Non_gene_f_c <- Non_gene_id %>% 
                group_by(type) %>% summarise(count = n())
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                exon_i_count <- exon_count + 1  
                exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                cds_i_count <- cds_count + 1    
                cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                match_i_count <- match_count + 1
                match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                V_gene_segment_i_count <- V_gene_segment_count + 1
                V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                C_gene_segment_i_count <- C_gene_segment_count + 1
                C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                or_i_count <- or_count + 1
                or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                D_loop_i_count <- D_loop_count + 1
                D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                
              }
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  gsub("ID=","OriginID=",x)
                  
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
              
              filtered_M_ls[[k]] <- rbind.data.frame(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",],
                                                     filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | filtered_gene_ls[[k]]$type == "miRNA" | 
                                                                             filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                             filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                             filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" |
                                                                             filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                                             filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",],Non_gene_id) %>% select(
                                                                               
                                                                               chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                                               
                                                                             )
              
            }
            
            Final_mod_d <- do.call("rbind.data.frame",filtered_M_ls)
            
            if(nrow(final_gb_df) >= 1){
              
              Final_mod_d <- rbind.data.frame(final_gb_df,Final_mod_d)
              
            }
            
            Unified_M_ls[[k]] <- Final_mod_d %>% select(chrloc,source,type,start,end,score,
                                                        strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
        } 
        
        else{
          
          Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                    strand,phase,modified_attributes,attributes) %>% 
            distinct(.keep_all = TRUE)
          
        }
        
      }
      
      if(unique(part_gene$strand) == "-"){
        
        if((unique(part_gene$gp_mod_attr) == "LOC101734477" | unique(part_gene$gp_mod_attr) == "LOC116406437" | unique(part_gene$gp_mod_attr) == "smyd5")){
          
          if(unique(part_gene$gp_mod_attr) == "LOC101734477"){
            
            g_LOC101734477_ls <- list()
            
            part_gene_LOC101734477 <- part_gene %>% filter(modified_attributes == "LOC101734477") %>% distinct(.keep_all = TRUE)
            
            part_gene_f_c <- part_gene_LOC101734477 %>% filter(!(type == "gene" | type == "mRNA" | type == "transcript")) %>% 
              group_by(type) %>% summarise(count = n())
            
            gene_mrna_count <- gene_mrna_count + 1
            exon_i_count <- exon_count + 1
            exon_count <- exon_count + part_gene_f_c[part_gene_f_c$type == "exon",]$count
            cds_i_count <- cds_count + 1
            cds_count <- cds_count + part_gene_f_c[part_gene_f_c$type == "CDS",]$count   
            
            g_LOC101734477_ls[[length(g_LOC101734477_ls)+1]] <- gene_mrna_count
            g_LOC101734477_ls[[length(g_LOC101734477_ls)+1]] <- paste0(exon_i_count,";",exon_count)
            g_LOC101734477_ls[[length(g_LOC101734477_ls)+1]] <- paste0(cds_i_count,";",cds_count)
            
            Unified_M_ls[[k]] <- part_gene_LOC101734477 %>% select(chrloc,source,type,start,end,score,
                                                                   strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          else if(unique(part_gene$gp_mod_attr) == "LOC116406437"){
            
            g_LOC116406437_ls <- list()
            
            part_gene_LOC116406437 <- part_gene %>% filter(modified_attributes == "LOC116406437") %>% distinct(.keep_all = TRUE)
            
            part_gene_f_c <- part_gene_LOC116406437 %>% filter(!(type == "gene" | type == "mRNA" | type == "transcript")) %>% 
              group_by(type) %>% summarise(count = n())
            
            gene_mrna_count <- gene_mrna_count + 1
            
            exon_i_count <- exon_count + 1
            exon_count <- exon_count + part_gene_f_c[part_gene_f_c$type == "exon",]$count
            cds_i_count <- cds_count + 1
            cds_count <- cds_count + part_gene_f_c[part_gene_f_c$type == "CDS",]$count   
            
            g_LOC116406437_ls[[length(g_LOC116406437_ls)+1]] <- gene_mrna_count
            g_LOC116406437_ls[[length(g_LOC116406437_ls)+1]] <- paste0(exon_i_count,";",exon_count)
            g_LOC116406437_ls[[length(g_LOC116406437_ls)+1]] <- paste0(cds_i_count,";",cds_count)
            
            Unified_M_ls[[k]] <- part_gene_LOC116406437 %>% select(chrloc,source,type,start,end,score,
                                                                   strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          else if(unique(part_gene$gp_mod_attr) == "smyd5"){
            
            smyd5_count_ls <- list()
            
            part_gene_smyd5 <- part_gene %>% filter(modified_attributes == "smyd5") %>% distinct(.keep_all = TRUE)
            
            part_gene_f_c <- part_gene_smyd5 %>% filter(!(type == "gene" | type == "mRNA" | type == "transcript")) %>% 
              group_by(type) %>% summarise(count = n())
            
            gene_mrna_count <- gene_mrna_count + 1
            exon_i_count <- exon_count + 1
            exon_count <- exon_count + part_gene_f_c[part_gene_f_c$type == "exon",]$count
            cds_i_count <- cds_count + 1
            cds_count <- cds_count + part_gene_f_c[part_gene_f_c$type == "CDS",]$count   
            
            smyd5_count_ls[[length(smyd5_count_ls)+1]] <- gene_mrna_count
            smyd5_count_ls[[length(smyd5_count_ls)+1]] <- paste0(exon_i_count,";",exon_count)
            smyd5_count_ls[[length(smyd5_count_ls)+1]] <- paste0(cds_i_count,";",cds_count)
            
            Unified_M_ls[[k]] <- part_gene_smyd5 %>% select(chrloc,source,type,start,end,score,
                                                            strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)           
            
          }
          
        }
        
        else if(((nrow(part_gene[part_gene$type == "gene",]) >= 2) &
                 (length(unique(part_gene[part_gene$type == "gene",]$type) %in% 
                         c("gene")) == 1 ) &
                 nrow(part_gene[part_gene$type == "tRNA",]) == 0) |
                ((nrow(part_gene[part_gene$type == "pseudogene",]) >= 2) & 
                 (length(unique(part_gene[part_gene$type == "pseudogene",]$type) %in%  
                         c("pseudogene")) == 1))){
          
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
          
          Non_gene_f_c <- Non_gene_id %>% 
            group_by(type) %>% summarise(count = n())
          
          gene_mrna_count <- gene_mrna_count + 1
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
            
            exon_i_count <- exon_count + 1  
            exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
            
            cds_i_count <- cds_count + 1   
            cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
            
            match_i_count <- match_count + 1
            match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
            
            V_gene_segment_i_count <- V_gene_segment_count + 1
            V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
            
            C_gene_segment_i_count <- C_gene_segment_count + 1
            C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
            
            or_i_count <- or_count + 1
            or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
            
            D_loop_i_count <- D_loop_count + 1
            D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
          }
          
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
          
          part_gene_id_df$attributes <- paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,
                                               ";",paste0(unique(part_gene_id$attributes),collapse = ";"))
          
          part_gene_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
          
          if(nrow(part_mRNA_id_df) == 1){
            
            part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
              
              gsub("ID=","transID=",x)
              
            )
            
            part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
              
              gsub("Parent=","Origin=",x)
              
            )
            
            part_mRNA_id_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                 gene_mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                 gene_mrna_count,";",paste0(unique(part_mRNA_id$attributes),collapse = ";"))
            
            part_mRNA_id_df$start <- min(part_gene$start)
            
            part_mRNA_id_df$end <- max(part_gene$end)
            
            part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
            
            parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
            
            if(nrow(Non_gene_id) >= 1){
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
              
              Non_gene_id$type_count <- NA
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                
              }
              
              Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                
                paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                
              )
              
              Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
              
              Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                
                paste0("ID=XB-",x)
                
              )
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                
                gsub("ID=","OriginID=",x)
                
                
              )
              
              Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
              
              Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
              
            }
            
            Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
            
            Unified_M_ls[[k]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
              select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE) 
            
          }
          
          else if(nrow(part_mRNA_id_df) < 1){
            
            part_mRNA_id_df <- part_gene_id_df 
            
            part_mRNA_id_df$type <- "transcript"
            
            part_mRNA_id_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                 gene_mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                 gene_mrna_count,";")
            
            part_mRNA_id_df$start <- part_gene_id_df$start
            
            part_mRNA_id_df$end <- part_gene_id_df$end 
            
            part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
            
            parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
            
            if(nrow(Non_gene_id) >= 1){
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
              
              Non_gene_id$type_count <- NA
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                
              }
              
              Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                
                paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                
              )
              
              Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
              
              Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                
                paste0("ID=XB-",x)
                
              )
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                
                gsub("ID=","OriginID=",x)
                
                
              )
              
              Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
              
              Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
              
            }
            
            Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
            
            Unified_M_ls[[k]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
              select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE) 
            
          }
          
        }
        
        else if((nrow(part_gene[part_gene$type == "gene" | 
                                part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) < 2 &
                 nrow(part_gene[part_gene$type == "tRNA",]) == 0) | 
                (nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
                 nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1 & 
                 nrow(part_gene[part_gene$type == "ncRNA_gene" | 
                                part_gene$type == "tRNA",]) == 0)){
          
          if(nrow(part_gene[part_gene$type == "gene" | 
                            part_gene$type == "pseudogene" |
                            part_gene$type == "ncRNA_gene",]) < 2){
            
            part_gene_df <- part_gene %>% filter(type == "gene"| type == "pseudogene"| type == "ncRNA_gene") %>%
              select(chrloc,source,type,start,end,score,phase,strand,modified_attributes,attributes) 
            
            if(nrow(part_gene_df) == 1){
              
              part_trans_M <- part_gene %>% filter(type == "mRNA" | 
                                                     type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                     type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                     type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                     type == "transcript" | type == "pseudogenic_transcript" | 
                                                     type == "primary_transcript")
              
              part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              ) %>% distinct(.keep_all = TRUE)
              
              part_gene_df$start <- min(part_gene$start)
              
              part_gene_df$end <- max(part_gene$end)
              
              part_gene_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
              
              Non_gene_id <- part_gene  %>% filter(!(type %in% part_gene_df$type | 
                                                       type %in% part_trans_M$type ))  %>% select(
                                                         chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                       ) %>% distinct(.keep_all = TRUE)
              
              Non_gene_f_c <- Non_gene_id %>% 
                group_by(type) %>% summarise(count = n())
              
              gene_mrna_count <- gene_mrna_count + 1
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                exon_i_count <- exon_count + 1  
                exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                cds_i_count <- cds_count + 1    
                cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                match_i_count <- match_count + 1
                match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                V_gene_segment_i_count <- V_gene_segment_count + 1
                V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                C_gene_segment_i_count <- C_gene_segment_count + 1
                C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                or_i_count <- or_count + 1
                or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                D_loop_i_count <- D_loop_count + 1
                D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                
              }
              
              part_gene_df$attributes <- sapply(part_gene_df$attributes,function(x) gsub("ID=",
                                                                                         paste0("ID=XBXT10g",
                                                                                                paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";gID="),x))
              
              if(nrow(part_trans_M_df) >= 1){
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) 
                  gsub("Parent=",paste0("Parent=XBXT10g",
                                        paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                        gene_mrna_count,";Origin="),x))
                
                part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",
                                                                                                 paste0("ID=XB-mRNA",
                                                                                                        paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";transID="),x))
                
                part_trans_M_df$start <- part_gene_df$start
                
                part_trans_M_df$end <- part_gene_df$end 
                
                part_trans_M_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                    gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    
                    gsub("ID=","OriginID=",x)
                    
                    
                  )
                  
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                         strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
              else if(nrow(part_trans_M_df) < 1){
                
                part_trans_M_df <- part_gene_df 
                
                part_trans_M_df$type <- "transcript"
                
                part_trans_M_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                     gene_mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                     gene_mrna_count,";")
                
                part_trans_M_df$start <- part_gene_df$start
                
                part_trans_M_df$end <- part_gene_df$end 
                
                part_trans_M_df$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                    gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    
                    gsub("ID=","OriginID=",x)
                    
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                part_gene <- rbind.data.frame(part_gene_df,part_trans_M_df,Non_gene_id) %>% 
                  select(chrloc,source,type,start,end,score,                                                                           
                         strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
                
              }
              
              Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                        strand,phase,modified_attributes,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            else if(nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) < 1 ){
              
              if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) < 1 )){
                
                Non_gene_f_c <- part_gene %>% 
                  group_by(type) %>% summarise(count = n())
                
                gene_mrna_count <- gene_mrna_count + 1
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                org_parent_id <- gsub("Parent=","",str_split(str_extract(part_gene$attributes[1],pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
                
                mod_gene_id <- part_gene[1,]
                
                mod_mRNA_id <- part_gene[1,]
                
                mod_gene_id$type <- "gene"
                mod_mRNA_id$type <- "mRNA"
                
                mod_gene_id$start <- min(part_gene$start)
                mod_gene_id$end <- max(part_gene$end)
                
                mod_mRNA_id$start <- min(part_gene$start)
                mod_mRNA_id$end <- max(part_gene$end)
                
                gene_id <- paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                  gene_mrna_count,";gbkey=gene")
                
                parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
                
                mRNA_id <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                  gene_mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                  gene_mrna_count,";")
                
                pid <- gsub("ID=","",str_split(str_extract(mRNA_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
                
                part_gene$attributes <- sapply(part_gene$attributes,function(x) 
                  
                  gsub("Parent=",paste0("Parent=",pid,";Origin="),x)
                  
                )
                
                part_gene$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  part_gene[part_gene$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  part_gene[part_gene$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  part_gene[part_gene$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  part_gene[part_gene$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  part_gene[part_gene$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  part_gene[part_gene$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  part_gene[part_gene$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                part_gene$type_mod <- sapply(part_gene$type_count, function(x)
                  
                  paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                  
                )
                
                part_gene$type_f <- paste0(casefold(part_gene$type),part_gene$type_mod)
                
                part_gene$type_final <- sapply(part_gene$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                part_gene$attributes <- sapply(part_gene$attributes, function(x)
                  
                  
                  gsub("ID=","OriginID=",x)
                  
                  
                )
                
                part_gene$attributes <- paste0(part_gene$type_final,part_gene$attributes)
                
                part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                mod_gene_id$attributes <- gene_id
                mod_mRNA_id$attributes <- mRNA_id
                
                mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                part_gene$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                mod_gene_id <- mod_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                mod_mRNA_id <- mod_mRNA_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,part_gene) %>% select(chrloc,source,type,start,end,score,
                                                                                            strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
              else if((unique(part_gene$source) == "Genbank") & (nrow(unique(part_gene[part_gene$type == "mRNA",])) >= 1 )){
                
                mod_mRNA <- part_gene %>% filter(type == "mRNA")
                
                mod_mRNA_id <- part_gene[which.max(part_gene$type == "mRNA"),]
                
                org_parent_id <- gsub("Parent=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "Parent=.+;"),pattern = ";")[[1]][1])  
                
                mod_gene_id <- mod_mRNA_id       
                
                mod_gene_id$type <- "gene"
                mod_mRNA_id$type <- "mRNA"
                
                mod_gene_id$start <- min(part_gene$start)
                mod_gene_id$end <- max(part_gene$end)
                
                mod_mRNA_id$start <- min(part_gene$start)
                mod_mRNA_id$end <- max(part_gene$end)
                
                Non_gene_id <- part_gene %>% filter(!(type %in% mod_gene_id$type | type %in% mod_mRNA$type))
                
                Non_gene_f_c <- Non_gene_id %>% 
                  group_by(type) %>% summarise(count = n())
                
                gene_mrna_count <- gene_mrna_count + 1
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                gene_id <- paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                  gene_mrna_count,";gbkey=gene")
                
                parent_id <- gsub("ID=","",str_split(str_extract(gene_id,pattern = "ID=.+;"),pattern = ";")[[1]][1])  
                
                mod_gene_id$attributes <- gene_id
                
                mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                mod_mRNA_id$attributes <- sapply(mod_mRNA_id$attributes,function(x) gsub("ID=",paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),
                                                                                                      gene_mrna_count,";gID="),x))
                
                p_id <- gsub("ID=","",str_split(str_extract(mod_mRNA_id$attributes,pattern = "ID.+;"),pattern = ";")[[1]][1])
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                  
                  gsub("Parent=",paste0("Parent=",p_id,";Origin="),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  gsub("ID=","OriginID=",x)
                  
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                mod_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                mod_mRNA_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")
                
                mod_gene_id <- mod_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                mod_mRNA_id <- mod_mRNA_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                part_gene <- rbind.data.frame(mod_gene_id,mod_mRNA_id,Non_gene_id) %>% select(chrloc,source,type,start,end,score,
                                                                                              strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)        
                
              }
              
              if(unique(part_gene$source) == "NCBI"){
                
                if(nrow(part_gene[part_gene$type == "mRNA" | part_gene$type == "tRNA" |
                                  part_gene$type == "rRNA" | part_gene$type == "lnc_RNA" | part_gene$type == "snoRNA" |
                                  part_gene$type == "snRNA" | part_gene$type == "scRNA" | part_gene$type == "miRNA" |
                                  part_gene$type == "Y_RNA" | part_gene$type == "ncRNA" | part_gene$type == "guide_RNA" |
                                  part_gene$type == "transcript" | part_gene$type == "primary_transcript" | part_gene$type == "pseudogenic_transcript",]) >= 1){
                  
                  Non_trans_count <- Non_trans_count + 1
                  
                  part_gene$attributes <- sapply(part_gene$attributes, function(x)
                    
                    gsub("ID=",paste0("ID=XB-trans",paste0(rep(0,6-nchar(Non_trans_count)),collapse = ""),Non_trans_count,";OriginID="),x)
                    
                  )
                  
                  part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                    strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)   
                  
                } 
                
                else if(nrow(part_gene[part_gene$type == "exon" | part_gene$type == "CDS" |
                                       part_gene$type == "match" | part_gene$type == "V_gene_segment" |
                                       part_gene$type == "C_gene_segment"|part_gene$type == "origin_of_replication" |
                                       part_gene$type == "D_loop",]) >= 1){
                  
                  Non_gene_f_c <- part_gene %>%
                    group_by(type) %>% summarise(count = n())
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    exon_i_count <- exon_count + 1
                    exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    cds_i_count <- cds_count + 1    
                    cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    match_i_count <- match_count + 1
                    match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    V_gene_segment_i_count <- V_gene_segment_count + 1
                    V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    C_gene_segment_i_count <- C_gene_segment_count + 1
                    C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    or_i_count <- or_count + 1
                    or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    D_loop_i_count <- D_loop_count + 1
                    D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                    
                  }
                  
                  part_gene$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    part_gene[part_gene$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    part_gene[part_gene$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    part_gene[part_gene$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    part_gene[part_gene$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    part_gene[part_gene$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    part_gene[part_gene$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    part_gene[part_gene$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  part_gene$type_mod <- sapply(part_gene$type_count, function(x)
                    
                    paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                    
                  )
                  
                  part_gene$type_f <- paste0(casefold(part_gene$type),part_gene$type_mod)
                  
                  part_gene$type_final <- sapply(part_gene$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  part_gene$attributes <- sapply(part_gene$attributes, function(x)
                    
                    
                    gsub("ID=","OriginID=",x)
                    
                    
                  )
                  
                  part_gene$attributes <- paste0(part_gene$type_final,part_gene$attributes)
                  
                  part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                    strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)        
                }
                
              }
              
            }
            
            Unified_M_ls[[k]] <- part_gene %>% select(chrloc,source,type,start,end,score,
                                                      strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          else if(nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
                  nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1){
            
            Unified_ls <- list()
            
            part_gene_No <- part_gene %>% filter(type == "gene" | type == "pseudogene") %>% 
              group_by(type) %>% summarise( Count = n(), source_t = paste0(source,collapse = ";"),
                                            mod_attr = paste0(unique(modified_attributes),collapse = ";") )
            
            for(i in 1 : length(part_gene_No$Count)){
              
              if(part_gene_No$Count[i] == 1){
                
                gene_mrna_count <- gene_mrna_count + 1
                
                source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
                
                part_gene_Mod_attr <- unique(str_split(part_gene_No$mod_attr[i],pattern = ";")[[1]])
                
                part_gene_type <- part_gene_No[i,]$type    
                
                part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
                
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
                
                part_gene_M_df$start <- min(part_gene_M$start)
                part_gene_M_df$end <- max(part_gene_M$end)
                part_gene_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                
                part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) gsub("ID=",
                                                                                               paste0("ID=XBXT10g",
                                                                                                      paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";gID="),x))
                
                Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | 
                                                           type %in% part_trans_M$type ))  %>% select(
                                                             chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                           ) %>% distinct(.keep_all = TRUE)
                
                Non_gene_f_c <- Non_gene_id %>% 
                  group_by(type) %>% summarise(count = n())
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1  
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                if(nrow(part_trans_M_df) == 1){
                  
                  part_trans_M_df$start <- min(part_gene_M$start)
                  part_trans_M_df$end <- max(part_gene_M$end)
                  part_trans_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                  
                  part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("Parent=",
                                                                                                   paste0("Parent=XBXT10g",
                                                                                                          paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";Origin="),x))
                  
                  part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) gsub("ID=",
                                                                                                   paste0("ID=XB-mRNA",
                                                                                                          paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";transID="),x))
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  if(nrow(Non_gene_id) >= 1){  
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                    select(chrloc,source,type,start,end,score,                                                                           
                           strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)
                  
                }
                
                else if(nrow(part_trans_M_df) != 1){
                  
                  part_trans_M_df <- part_gene_M_df
                  
                  part_trans_M_df$type <- "transcript"
                  
                  part_trans_M_df$start <- min(part_gene_M$start)
                  part_trans_M_df$end <- max(part_gene_M$end)
                  part_trans_M_df$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                  
                  part_trans_M_df$attributes <- 
                    paste0("ID=XB-mRNA",
                           paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";","Parent=XBXT10g",
                           paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";")
                  
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(part_gene_M$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                    select(chrloc,source,type,start,end,score,                                                                           
                           strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)
                  
                }
                
                Unified_ls[[i]] <- part_gene_M %>% select(chrloc,source,type,start,end,score,
                                                          strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }   
              
              else if(part_gene_No$Count[i] > 1){
                
                gene_mrna_count <- gene_mrna_count + 1
                
                source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
                
                part_gene_type <- part_gene_No[i,]$type    
                
                part_gene_Mod_attr <- unique(part_gene[
                  part_gene$type == part_gene_type,]$modified_attributes)
                
                part_gene_Mod <- part_gene %>% filter((source %in% source_tot ) & (modified_attributes %in% part_gene_Mod_attr))
                
                part_gene_id <- part_gene_Mod %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
                
                part_mRNA_id <- part_gene_Mod %>% filter(type == "mRNA" | type == "tRNA"| 
                                                           type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                           type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                           type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                           type == "transcript" | type == "pseudogenic_transcript" | 
                                                           type == "primary_transcript")
                
                Non_gene_id <-  part_gene_Mod %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
                  select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) 
                
                Non_gene_f_c <- Non_gene_id %>% 
                  group_by(type) %>% summarise(count = n())
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                
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
                
                part_gene_id_df$attributes <- 
                  paste0("ID=XBXT10g",
                         paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",paste0(unique(part_gene_id$attributes),collapse = ";"))    
                
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
                  
                  part_mRNA_id_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";","Parent=XBXT10g",
                                                       paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",
                                                       paste0(unique(part_mRNA_id$attributes),collapse = ";"))   
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene_M <- rbind.data.frame(part_gene_M_df,part_mRNA_id_df,Non_gene_id) %>% 
                    select(chrloc,source,type,start,end,score,                                                                           
                           strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)
                  
                  
                }
                
                else if(nrow(part_mRNA_id_df) != 1){
                  
                  part_mRNA_id_df <- part_gene_id_df
                  
                  part_mRNA_id_df$type <- "transcript"
                  
                  part_mRNA_id_df$start <- min(part_gene_Mod$start)
                  part_mRNA_id_df$end <- max(part_gene_Mod$end)
                  part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
                  
                  part_mRNA_id_df$attributes <- 
                    paste0("ID=XB-mRNA",
                           paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";","Parent=XBXT10g",
                           paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";")
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  part_gene_M <- rbind.data.frame(part_gene_M_df,part_mRNA_id_df,Non_gene_id) %>% 
                    select(chrloc,source,type,start,end,score,                                                                           
                           strand,phase,modified_attributes,attributes) %>% 
                    distinct(.keep_all = TRUE)
                  
                }
                
                Unified_ls[[i]] <- part_gene_M %>%
                  select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% 
                  distinct(.keep_all = TRUE)
                
              }
              
            }
            
            Unified_M_ls[[k]] <- do.call("rbind.data.frame",Unified_ls) 
            
          } 
          
        }
        
        else if(((nrow(part_gene[part_gene$type == "ncRNA_gene",]) > 1) &
                 (length(unique(part_gene[part_gene$type == "ncRNA_gene",]$type) %in% 
                         c("ncRNA_gene")) == 1 )) | (nrow(part_gene[part_gene$type == "gene" | 
                                                                    part_gene$type == "pseudogene",]) >= 1 & 
                                                     nrow(part_gene[part_gene$type == "ncRNA_gene",]) >= 1) |
                (nrow(part_gene[part_gene$type == "tRNA",]) >= 1)) {
          
          other_ls[[length(other_ls)+ 1]] <- part_gene
          
          part_df <- part_gene %>% left_join(rank_data,by = c("type"))
          
          if(length(unique(part_df[part_df$type == "tRNA",]$type)) >= 1 | ((length(unique(part_df[part_df$type == "tRNA",]$type)) >= 1) & (length(unique(part_df[part_df$type == "ncRNA_gene" | part_df$type == "miRNA",]$type)) >= 1))){
            
            #part_df <- part_df %>% arrange(rank)
            
            part_df$ID <- sapply(part_df$attributes,function(x) gsub("ID=","",
                                                                     str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
            
            part_df$Parent <- sapply(part_df$attributes,function(x) gsub("Parent=","",
                                                                         str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))
            
            filtered_rna_gene_ids <- part_df %>% filter(type == "mRNA" |
                                                          type == "tRNA"| type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" |
                                                          type == "snRNA" | type == "scRNA" | type == "miRNA" | type == "Y_RNA" |
                                                          type == "ncRNA" | type == "guide_RNA" | type == "transcript" | 
                                                          type == "primary_transcript" |type == "pseudogenic_transcript") %>% select(ID,Parent)
            
            part_d <- part_df %>% left_join(filtered_rna_gene_ids,by = c("Parent" = "ID"))
            
            part_d <- part_d %>% mutate(
              
              transcript = ifelse(is.na(Parent.y),Parent,Parent.y),
              final_transcript = ifelse(is.na(transcript),ID,transcript)
              
              
            ) %>% arrange(final_transcript) %>% select(chrloc,source,type,start,end,
                                                       score,strand,phase,modified_attributes,attributes,final_transcript,Overlap_Value)
            
            filtered_mod_attr_ls <- list()
            final_mod_attr_ls <- list()
            
            for (l in 1 : length(unique(part_d$final_transcript))) {
              
              filtered_mod_attr_ls[[l]] <- part_d %>% filter(final_transcript %in% unique(part_d$final_transcript)[l]) 
              
              gene_filtered_mod_attr <- filtered_mod_attr_ls[[l]] %>% filter(type == "gene" | 
                                                                               type == "pseudogene" | type == "ncRNA_gene")
              
              transcript_filtered_mod_attr <- filtered_mod_attr_ls[[l]] %>% filter(type == "mRNA" |
                                                                                     type == "tRNA"| type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" |
                                                                                     type == "snRNA" | type == "scRNA" | type == "miRNA" | type == "Y_RNA" |
                                                                                     type == "ncRNA" | type == "guide_RNA" | type == "transcript" | 
                                                                                     type == "primary_transcript" |type == "pseudogenic_transcript")
              
              
              gene_filtered_mod_attr_df <- gene_filtered_mod_attr[which.max(gene_filtered_mod_attr$Overlap_Value),] %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              )
              
              transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr[which.max(transcript_filtered_mod_attr$Overlap_Value),] %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              )
              
              Non_gene_id <- filtered_mod_attr_ls[[l]] %>% filter(!(type %in% gene_filtered_mod_attr$type | type %in% transcript_filtered_mod_attr$type)) %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              )
              
              Non_gene_f_c <- Non_gene_id %>% 
                group_by(type) %>% summarise(count = n())
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                exon_i_count <- exon_count + 1 
                exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                cds_i_count <- cds_count + 1   
                cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                match_i_count <- match_count + 1
                match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                V_gene_segment_i_count <- V_gene_segment_count + 1
                V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                C_gene_segment_i_count <- C_gene_segment_count + 1
                C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                or_i_count <- or_count + 1
                or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                D_loop_i_count <- D_loop_count + 1
                D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                
              }
              
              if(nrow(gene_filtered_mod_attr_df) == 1){
                
                gene_mrna_count <- gene_mrna_count + 1
                
                gene_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                gene_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                
                gene_filtered_mod_attr_df$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                
                gene_filtered_mod_attr$attributes <- sapply(gene_filtered_mod_attr$attributes, function(x)
                  
                  gsub("ID=","gID=",x)
                  
                )
                
                gene_filtered_mod_attr_df$attributes <- paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",paste0(unique(gene_filtered_mod_attr$attributes),collapse = ";"))
                
                if(nrow(transcript_filtered_mod_attr_df) == 1 ){
                  
                  transcript_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                  transcript_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                  
                  transcript_filtered_mod_attr_df$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                  
                  transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                    
                    gsub("ID=","transID=",x)
                    
                  )
                  
                  transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                    
                    gsub("Parent=","Origin=",x)
                    
                  )
                  
                  transcript_filtered_mod_attr_df$attributes <- paste0("ID=XB-rna",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",
                                                                       "Parent=XBXT10g",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";",
                                                                       paste0(unique(transcript_filtered_mod_attr$attributes),collapse = ";"))
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(
                    transcript_filtered_mod_attr_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])   
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                  
                  final_mod_attr_ls[[l]] <- rbind.data.frame(gene_filtered_mod_attr_df,transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                    
                    chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                    
                  )
                  
                }
                
                else if(nrow(transcript_filtered_mod_attr_df) < 1){
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(gene_filtered_mod_attr_df$attributes,pattern = "ID=.+;"))[[1]][1])   
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                    
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                  
                  final_mod_attr_ls[[l]] <- rbind.data.frame(gene_filtered_mod_attr_df,transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                    
                    chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                    
                  )
                  
                }
                
              }
              
              if(nrow(gene_filtered_mod_attr_df) < 1){
                
                if( nrow(transcript_filtered_mod_attr_df) == 1 ){
                  
                  Non_trans_count <- Non_trans_count + 1
                  
                  transcript_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                  transcript_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                  
                  transcript_filtered_mod_attr_df$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                  
                  transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                    
                    gsub("ID=","transID=",x)
                    
                  )
                  
                  transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                    
                    gsub("Parent=","Origin=",x)
                    
                  )
                  
                  transcript_filtered_mod_attr_df$attributes <- 
                    
                    sapply(transcript_filtered_mod_attr_df$attributes, function(x)
                      
                      gsub("ID=",paste0("ID=XB-trans",paste0(rep(0,6-nchar(Non_trans_count)),collapse = ""),Non_trans_count,";transID="),x)    
                      
                    )
                  
                  parent_id <- gsub("ID=","",str_split(str_extract(
                    transcript_filtered_mod_attr_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])   
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                  
                  final_mod_attr_ls[[l]] <- rbind.data.frame(transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                    
                    chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                    
                  )
                  
                }
                
                else if(nrow(transcript_filtered_mod_attr_df) < 1){
                  
                  #parent_id <- gsub("ID=","",str_split(str_extract(gene_filtered_mod_attr_df$attributes,pattern = "ID=.+;"))[[1]][1])   
                  
                  if(nrow(Non_gene_id) >= 1){
                    
                    #Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    #gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                    
                    Non_gene_id$type_count <- NA
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                      
                    }
                    
                    if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                      
                      Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                      
                    }
                    
                    Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                      
                      paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                      
                    )
                    
                    Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                    
                    Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                      
                      paste0("ID=XB-",x)
                      
                    )
                    
                    Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                      
                      
                      gsub("ID=","OriginID=",x)
                      
                      
                    )
                    
                    Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                    
                    Non_gene_id$modified_attributes <- paste0(unique(filtered_mod_attr_ls[[l]]$modified_attributes),collapse = ";")
                    
                    #Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                    
                    
                  }
                  
                  Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes)
                  
                  final_mod_attr_ls[[l]] <- rbind.data.frame(Non_gene_id) %>% select(
                    
                    chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                    
                  )
                  
                }
                
              }  
              
            }
            
            final_mod_attr_df <- do.call("rbind.data.frame",final_mod_attr_ls)
            
            Unified_M_ls[[k]] <- final_mod_attr_df %>% select(chrloc,source,type,start,end,score,
                                                              strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
          if((length(unique(part_df[part_df$type == "ncRNA_gene" | part_df$type == "miRNA",]$type)) >= 1 & nrow(part_df[part_df$type == 'tRNA',]) == 0) | 
             (nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene",]) >= 1 & 
              nrow(part_gene[part_gene$type == "ncRNA_gene",]) > 1)){
            
            part_df <- part_df %>% arrange(rank)
            
            final_gb_ls <- list()
            
            if (nrow(part_df[part_df$source == "Genbank",]) >= 1){
              
              filtered_genbank_ls <- part_df %>% filter(source == "Genbank")
              
              filtered_genbank_ls <- filtered_genbank_ls %>% arrange(gp_mod_attr)
              
              filtered_gb_ls <- list()
              
              for (h in 1 : length(unique(filtered_genbank_ls$gp_mod_attr))) {
                
                gene_mrna_count <- gene_mrna_count + 1
                
                filtered_gb_ls[[h]] <- filtered_genbank_ls %>% filter(gp_mod_attr == unique(filtered_genbank_ls$gp_mod_attr)[h]) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                  
                )
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes,function(x) 
                  gsub("ID=",paste0("ID=XBXT10g",
                                    paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";gID="),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  gsub("ID=",paste0("ID=XB-mRNA",
                                    paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";transID="),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  gsub("Parent=",paste0("Parent=XBXT10g",
                                        paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";Origin="),x))
                
                if(nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]) >= 1){
                  
                  parent_id <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))  
                  
                }
                
                else if (nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]) < 1){
                  
                  
                  parent_id <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))  
                  
                }
                
                Non_gene_id <- filtered_gb_ls[[h]] %>% filter(!(type == "gene" | type == "mRNA"))
                
                Non_gene_f_c <- Non_gene_id %>% 
                  group_by(type) %>% summarise(count = n())
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  exon_i_count <- exon_count + 1  
                  exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  cds_i_count <- cds_count + 1    
                  cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  match_i_count <- match_count + 1
                  match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  V_gene_segment_i_count <- V_gene_segment_count + 1
                  V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  C_gene_segment_i_count <- C_gene_segment_count + 1
                  C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  or_i_count <- or_count + 1
                  or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  D_loop_i_count <- D_loop_count + 1
                  D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                  
                }
                
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x)
                    
                  )
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    
                    gsub("ID=","OriginID=",x)
                    
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                
                final_gb_ls[[h]] <- rbind.data.frame(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",],
                                                     filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",],
                                                     Non_gene_id) %>% select(
                                                       
                                                       chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes                                       
                                                       
                                                     )
                
              }
              
            }
            
            final_gb_df <- do.call("rbind.data.frame",final_gb_ls)
            
            part_df <- part_df %>% filter(!(source %in% "Genbank"))
            
            part_df$transcript <- sapply(part_df$attributes,function(x) 
              
              ifelse(str_detect(x,pattern = "Parent=.+"),
                     
                     gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1]),
                     
                     gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
              
            )
            
            part_df <- part_df %>% arrange(transcript)
            
            filtered_gene_ls <- list()
            filtered_M_ls <- list()
            
            for(k in 1: length(unique(part_df$transcript))){
              
              gene_mrna_count <- gene_mrna_count + 1
              
              filtered_gene_ls[[k]] <- part_df %>% filter(transcript %in% unique(part_df$transcript)[k]) %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                
              )
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",]$attributes,function(x) 
                gsub("ID=",paste0("ID=XBXT10g",
                                  paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";gID="),x))
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                      filtered_gene_ls[[k]]$type == "miRNA" | 
                                      filtered_gene_ls[[k]]$type == "rRNA" | 
                                      filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                      filtered_gene_ls[[k]]$type == "snoRNA" | 
                                      filtered_gene_ls[[k]]$type == "scRNA" | 
                                      filtered_gene_ls[[k]]$type == "Y_RNA" | 
                                      filtered_gene_ls[[k]]$type == "ncRNA" | 
                                      filtered_gene_ls[[k]]$type == "guide_RNA" |
                                      filtered_gene_ls[[k]]$type == "transcript" | 
                                      filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "miRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "rRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "snoRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "Y_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "guide_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "transcript" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes,function(x) 
                                                                                                                                                                                 gsub("Parent=",paste0("Parent=XBXT10g",
                                                                                                                                                                                                       paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";Origin="),x))
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | filtered_gene_ls[[k]]$type == "miRNA" | 
                                      filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                      filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                      filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" |
                                      filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                      filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes,function(x) 
                                                                                                                                                                                 gsub("ID=",paste0("ID=XB-rna",paste0(rep(0,6-nchar(gene_mrna_count)),collapse = ""),gene_mrna_count,";transID="),x))
              if(nrow(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                            filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                            filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                            filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                            filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                            filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]) >= 1){                                                                                                                                              
                
                parent_id <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                            filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                            filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                            filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                            filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                            filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
                
              }
              
              else if (nrow(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                  filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                  filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                  filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                  filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                  filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]) < 1){
                
                parent_id <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
                
                
              }
              
              Non_gene_id <- filtered_gene_ls[[k]] %>% filter(!(type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | type == "snRNA" | 
                                                                  type == "scRNA" | type == "miRNA" | type == "transcript" | type == "primary_transcript" |
                                                                  type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "ncRNA_gene" | type == "mRNA" | type == "gene" | type == "pseudogene"
                                                                
              ))
              
              Non_gene_f_c <- Non_gene_id %>% 
                group_by(type) %>% summarise(count = n())
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                exon_i_count <- exon_count + 1  
                exon_count <- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                cds_i_count <- cds_count + 1    
                cds_count <- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                match_i_count <- match_count + 1
                match_count <- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                V_gene_segment_i_count <- V_gene_segment_count + 1
                V_gene_segment_count <- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                C_gene_segment_i_count <- C_gene_segment_count + 1
                C_gene_segment_count <- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                or_i_count <- or_count + 1
                or_count <- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                D_loop_i_count <- D_loop_count + 1
                D_loop_count <- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                
              }
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x,";")
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  gsub("ID=","OriginID=",x)
                  
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
              
              filtered_M_ls[[k]] <- rbind.data.frame(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",],
                                                     filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | filtered_gene_ls[[k]]$type == "miRNA" | 
                                                                             filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                             filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                             filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" |
                                                                             filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                                             filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",],Non_gene_id) %>% select(
                                                                               
                                                                               chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes
                                                                               
                                                                             )
              
            }
            
            Final_mod_d <- do.call("rbind.data.frame",filtered_M_ls)
            
            if(nrow(final_gb_df) >= 1){
              
              Final_mod_d <- rbind.data.frame(final_gb_df,Final_mod_d)
              
            }
            
            Unified_M_ls[[k]] <- Final_mod_d %>% select(chrloc,source,type,start,end,score,
                                                        strand,phase,modified_attributes,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
        } 
        
        else{
          
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

for (i in 1: length(unique(Sort_Unified_F$chrloc))){
  
  filtered_Unified_NCBI_Ens_XGC_ls[[i]] <- Sort_Unified_F %>% filter(chrloc == unique(Sort_Unified_F$chrloc)[i])
  Unified_NCBI_Ens_XGC_ls[[i]] <- unified_NCBI_Ens_XGC(filtered_Unified_NCBI_Ens_XGC_ls[[i]])
  
}

Final_Unified_df <- do.call("rbind.data.frame",Unified_NCBI_Ens_XGC_ls)

Final_Unified_df$attributes <- as.vector(Final_Unified_df$attributes)

Final_Unified_d <- Final_Unified_df %>% distinct(.keep_all = TRUE)

write.table(Final_Unified_d,file = "Final_Unified_d.gff3",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

###############  changing attributes for LOC and smy #############

Final_Unified_d[Final_Unified_d$type == "gene" & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes <- sapply(Final_Unified_d[Final_Unified_d$type == "gene" & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes,function(x)

  gsub("ID=",paste0("ID=XBXT10g",paste0(rep(0,6-nchar(g_LOC116406437_ls[[1]])),collapse = ""),g_LOC116406437_ls[[1]]),x)                                                        

 )

Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes <- sapply(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes,function(x)

 gsub("Parent=",paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(g_LOC116406437_ls[[1]])),collapse = ""),g_LOC116406437_ls[[1]]),x)                                                        

)

Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes <- sapply(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes,function(x)

 gsub("ID=",paste0("ID=XB-rna",paste0(rep(0,6-nchar(g_LOC116406437_ls[[1]])),collapse = ""),g_LOC116406437_ls[[1]]),x)                                                        

)

parent_ids <- gsub("ID=","",str_split(str_extract(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes,pattern = "ID=.+;"))[[1]][1])

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "tRNA") & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes <- sapply(Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes,function(x)

 gsub("Parent=",paste0("Parent=",parent_ids,";Origin=",x))                                

)

range_exon_LOC116406437  <- str_split(g_LOC116406437_ls[[2]],pattern = ";")[[1]][1] : str_split(g_LOC116406437_ls[[2]],pattern = ";")[[1]][2]

Final_Unified_d$type_count <- NA

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437" & Final_Unified_d$type == "exon",]$type_count <- range_exon_LOC116406437

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437" & Final_Unified_d$type == "exon",]$attributes <- sapply(
  
  Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437" & Final_Unified_d$type == "exon",]$attributes,function(x)

  gsub("ID=",paste0("ID=XB-","exon",paste0(rep(0,6-nchar(Final_Unified_d$type_count)),collapse = ""),Final_Unified_d$type_count,";OriginID="),x)                                
  
)

range_CDS_LOC116406437  <- str_split(g_LOC116406437_ls[[3]],pattern = ";")[[1]][1] : str_split(g_LOC116406437_ls[[3]],pattern = ";")[[1]][2]

Final_Unified_d[!(Final_Unified_d$type == "gene" | 
 Final_Unified_d$type == "mRNA" | 
   Final_Unified_d$type == "transcript") & 
   Final_Unified_d$modified_attributes == "LOC116406437" & 
   Final_Unified_d$type == "CDS",]$type_count <- range_CDS_LOC116406437

  
Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437" & Final_Unified_d$type == "CDS",]$attributes <- sapply(
  
  Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437" & Final_Unified_d$type == "CDS",]$attributes,function(x)
    
    gsub("ID=",paste0("ID=XB-","CDS",paste0(rep(0,6-nchar(Final_Unified_d$type_count)),collapse = ""),Final_Unified_d$type_count,";OriginID="),x)                                
  
)

############### changing LOC116406437 ##############

Final_Unified_d[Final_Unified_d$type == "gene" & Final_Unified_d$modified_attributes == "LOC101734477",]$attributes <- sapply(Final_Unified_d[Final_Unified_d$type == "gene" & Final_Unified_d$modified_attributes == "LOC101734477",]$attributes,function(x)
  
  gsub("ID=",paste0("ID=XBXT10g",paste0(rep(0,6-nchar(g_LOC101734477_ls[[1]])),collapse = ""),g_LOC101734477_ls[[1]]),x)                                                        
  
)

Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477",]$attributes <- sapply(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477",]$attributes,function(x)
  
  gsub("Parent=",paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(g_LOC101734477_ls[[1]])),collapse = ""),g_LOC101734477_ls[[1]]),x)                                                        
  
)

Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477",]$attributes <- sapply(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC116406437",]$attributes,function(x)
  
  gsub("ID=",paste0("ID=XB-rna",paste0(rep(0,6-nchar(g_LOC101734477_ls[[1]])),collapse = ""),g_LOC101734477_ls[[1]]),x)                                                        
  
)

parent_ids <- gsub("ID=","",str_split(str_extract(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477",]$attributes,pattern = "ID=.+;"))[[1]][1])

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477",]$attributes <- sapply(Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477",]$attributes,function(x)
  
  gsub("Parent=",paste0("Parent=",parent_ids,";Origin="),x)                                
  
)

range_exon_LOC101734477  <- str_split(g_LOC101734477_ls[[2]],pattern = ";")[[1]][1] : str_split(g_LOC101734477_ls[[2]],pattern = ";")[[1]][2]

#Final_Unified_d$type_count <- NA

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477" & Final_Unified_d$type == "exon",]$type_count <- range_exon_LOC101734477

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477" & Final_Unified_d$type == "exon",]$attributes <- sapply(
  
  Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477" & Final_Unified_d$type == "exon",]$attributes,function(x)
    
    gsub("ID=",paste0("ID=XB-","exon",paste0(rep(0,6-nchar(Final_Unified_d$type_count)),collapse = ""),Final_Unified_d$type_count,";OriginID="),x)                                
  
)

range_CDS_LOC101734477  <- str_split(g_LOC101734477_ls[[3]],pattern = ";")[[1]][1] : str_split(g_LOC101734477_ls[[3]],pattern = ";")[[1]][2]

Final_Unified_d[!(Final_Unified_d$type == "gene" | 
                    Final_Unified_d$type == "mRNA" | 
                    Final_Unified_d$type == "transcript") & 
                  Final_Unified_d$modified_attributes == "LOC101734477" & 
                  Final_Unified_d$type == "CDS",]$type_count <- range_CDS_LOC101734477


Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477" & Final_Unified_d$type == "CDS",]$attributes <- sapply(
  
  Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "LOC101734477" & Final_Unified_d$type == "CDS",]$attributes,function(x)
    
    gsub("ID=",paste0("ID=XB-","CDS",paste0(rep(0,6-nchar(Final_Unified_d$type_count)),collapse = ""),Final_Unified_d$type_count,";OriginID="),x)                                
  
)

################################# smyd5 change ###############

Final_Unified_d[Final_Unified_d$type == "gene" & Final_Unified_d$modified_attributes == "smyd5",]$attributes <- sapply(Final_Unified_d[Final_Unified_d$type == "gene" & Final_Unified_d$modified_attributes == "smyd5",]$attributes,function(x)
  
  gsub("ID=",paste0("ID=XBXT10g",paste0(rep(0,6-nchar(g_smyd5_ls[[1]])),collapse = ""),g_smyd5_ls[[1]]),x)                                                        
  
)

Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5",]$attributes <- sapply(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5",]$attributes,function(x)
  
  gsub("Parent=",paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(g_smyd5_ls[[1]])),collapse = ""),g_smyd5_ls[[1]]),x)                                                        
  
)

Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5",]$attributes <- sapply(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5",]$attributes,function(x)
  
  gsub("ID=",paste0("ID=XB-rna",paste0(rep(0,6-nchar(g_smyd5_ls[[1]])),collapse = ""),g_smyd5_ls[[1]]),x)                                                        
  
)

parent_ids <- gsub("ID=","",str_split(str_extract(Final_Unified_d[(Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5",]$attributes,pattern = "ID=.+;"))[[1]][1])

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5",]$attributes <- sapply(Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5",]$attributes,function(x)
  
  gsub("Parent=",paste0("Parent=",parent_ids,";Origin="),x)                                
  
)

range_exon_smyd5  <- str_split(g_smyd5_ls[[2]],pattern = ";")[[1]][1] : str_split(g_smyd5_ls[[2]],pattern = ";")[[1]][2]

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5" & Final_Unified_d$type == "exon",]$type_count <- range_exon_smyd5

Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5" & Final_Unified_d$type == "exon",]$attributes <- sapply(
  
  Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5" & Final_Unified_d$type == "exon",]$attributes,function(x)
    
    gsub("ID=",paste0("ID=XB-","exon",paste0(rep(0,6-nchar(Final_Unified_d$type_count)),collapse = ""),Final_Unified_d$type_count,";OriginID="),x)                                
  
)

range_CDS_smyd5  <- str_split(g_smyd5_ls[[3]],pattern = ";")[[1]][1] : str_split(g_smyd5_ls[[3]],pattern = ";")[[1]][2]

Final_Unified_d[!(Final_Unified_d$type == "gene" | 
                    Final_Unified_d$type == "mRNA" | 
                    Final_Unified_d$type == "transcript") & 
                  Final_Unified_d$modified_attributes == "smyd5" & 
                  Final_Unified_d$type == "CDS",]$type_count <- range_CDS_smyd5


Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5" & Final_Unified_d$type == "CDS",]$attributes <- sapply(
  
  Final_Unified_d[!(Final_Unified_d$type == "gene" | Final_Unified_d$type == "mRNA" | Final_Unified_d$type == "transcript") & Final_Unified_d$modified_attributes == "smyd5" & Final_Unified_d$type == "CDS",]$attributes,function(x)
    
    gsub("ID=",paste0("ID=XB-","CDS",paste0(rep(0,6-nchar(Final_Unified_d$type_count)),collapse = ""),Final_Unified_d$type_count,";OriginID="),x)                                
  
)

               
          ########### removing XENTR ones ##########

Final_Unified_d <- read.delim("Final_Unified_d.gff3",sep = "\t",header = FALSE)

names(Final_Unified_d) <- c("chrloc","source","type","start","end","score","strand","phase","modified_attributes","attributes")

Final_Unified_M <- Final_Unified_d %>% mutate(
  
  mod_attr_split = str_split(modified_attributes,pattern = ";")
  
)

Final_Unified_M$mod_attr_c <- lapply(Final_Unified_M$mod_attr_split, function(x)
  
  length(x)
  
)

Final_df <- Final_Unified_M %>% mutate(
  
  rep_val = ifelse((str_detect(modified_attributes,pattern = "XENTR") & mod_attr_c >= 2),"TRUE","FALSE")

)

Final_Unified_M_Values <- Final_df[Final_df$rep_val == "TRUE",]

rep_val_values <- lapply(Final_Unified_M_Values$mod_attr_split, function(x)
  
  x[!(str_detect(x,"XENTR"))]
  
)

rep_v <- unique(unlist(rep_val_values))

Filtered_df <- Final_Unified_M %>% filter(!(modified_attributes %in% rep_v))

#write.table(Filtered_df,file = "Filtered_df.gff3",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

################# removing repetitive complicated models ##############

Filtered_Gene_Data <- Filtered_df %>% filter(type == "gene")

trna_mod_atr <- data.frame(ids = Filtered_df[Filtered_df$type == "tRNA",]$modified_attributes)

Filtered_Gene_df <- Filtered_Gene_Data %>% filter(!(modified_attributes %in% unique(trna_mod_atr$ids)))

#Filtered_Gene_d <- Filtered_Gene_df %>% filter(!(modified_attributes == "smyd5" | modified_attributes == "LOC101734477" | modified_attributes == "LOC116406437"))  

filtered_id_list <- unique(unlist(Filtered_Gene_df$mod_attr_split))

##########################################

##### finding overlaps 

library(data.table)

g_1 <- Filtered_Gene_df %>% select(chrloc,strand,start,end)
g_2 <- Filtered_Gene_df %>% select(chrloc,strand,start,end)

setDT(g_1)
setDT(g_2)
setkey(g_1)

overlap_g_d <- foverlaps(g_2,g_1,type = c("within"))

overlap_g_data <- overlap_g_d %>% left_join(Filtered_Gene_df,by = c("chrloc","strand","start","end"))

added_Overlap_g_data <- overlap_g_data %>% left_join(Filtered_Gene_df,by = c("chrloc","strand","i.start" = "start","i.end" = "end"))

filtered_added_Overlap_g_d <- added_Overlap_g_data %>% filter(!(modified_attributes.x == modified_attributes.y))

############ grouping variables to map them to add in the grp variable #########

gp_filtered_added_Overlap_g_d <- filtered_added_Overlap_g_d %>% group_by(chrloc,strand,start,end) %>% mutate(
  
 Combined_Mod_Attr = paste0(paste0(unique(unlist(mod_attr_split.x)),collapse = ";"),";",paste0(unique(unlist(mod_attr_split.y)),collapse = ";"))
  
) %>% select(chrloc,strand,start,end,modified_attributes.x,i.start,i.end,modified_attributes.y,Combined_Mod_Attr) %>% distinct(.keep_all = TRUE) 

gp_mod_attr_2 <- gp_filtered_added_Overlap_g_d %>% select(chrloc,strand,start = i.start,end = i.end,modified_attributes = modified_attributes.y,Combined_Mod_Attr) %>% distinct(.keep_all = TRUE)

gp_mod_attr_1 <- gp_filtered_added_Overlap_g_d %>% select(chrloc,strand,start,end,modified_attributes = modified_attributes.x,Combined_Mod_Attr) %>% distinct(.keep_all = TRUE)

Combined_gp_var_df <- rbind(gp_mod_attr_1,gp_mod_attr_2) %>% distinct(.keep_all = TRUE)

#######################################################

Com_Filtered_df <- Filtered_df %>% left_join(

Combined_gp_var_df,by = c("chrloc","strand","modified_attributes"))

Com_Filtered_d <- Com_Filtered_df %>% mutate(
  
  Final_Mod_Attr = ifelse(is.na(Combined_Mod_Attr),modified_attributes,Combined_Mod_Attr)
  
) %>% select(chrloc,source,type,start = start.x,end = end.x,score,strand,phase,modified_attributes,attributes,Final_Mod_Attr) %>% distinct(.keep_all = TRUE)

Com_Filtered_d <- Com_Filtered_d %>% mutate(
  
  Overlap_Value = end - start
  
)

#################################

### 2 out of 3 strategy #######


##incase if NCBI and ensembl has only one or lesser transcripts compared to XGC,
## to ensembl,go for only ncbi and xgc merging

###rectifying gene model error###

Com_Filtered_d[Com_Filtered_d$modified_attributes == "smyd5" & Com_Filtered_d$start >= 190777503 & Com_Filtered_d$end <= 190784498,]$Final_Mod_Attr <- "smyd5"

Com_Filtered_d[Com_Filtered_d$modified_attributes == "LOC101734477" & Com_Filtered_d$start >= 190776610 & Com_Filtered_d$end <= 190784506,]$Final_Mod_Attr <- "smyd5"

Com_Filtered_d[Com_Filtered_d$modified_attributes == "LOC101734477" & Com_Filtered_d$start >= 190777332 & Com_Filtered_d$end <= 190781859,]$Final_Mod_Attr <- "smyd5"

Com_Filtered_d[Com_Filtered_d$Final_Mod_Attr == "LOC101734477;smyd5",]$Final_Mod_Attr <- "smyd5;LOC116406512;LOC116406437;LOC116408425;LOC105946267" 

Com_Filtered_d[Com_Filtered_d$Final_Mod_Attr == "smyd5;LOC116406512;LOC116406437;LOC116408425;LOC105946267",]$Final_Mod_Attr <- "LOC100145764;ENSXETG00000033287;smyd5;LOC116406512" 

Com_Filtered_d <- Com_Filtered_d %>% distinct(.keep_all = TRUE)

Any_NA_com <- Com_Filtered_d %>% filter(is.na(Final_Mod_Attr))

write.table(Com_Filtered_d,file = "Com_Filtered_d.gff3",quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

Com_Filtered_d <- read.delim("Com_Filtered_d.gff3",sep = "\t",header = FALSE)

names(Com_Filtered_d) <- c("chrloc","source","type","start","end","score","strand","phase","modified_attributes","attributes","Final_Mod_Attr","Overlap_Value")

############################################################

############# re-naming of models ################

###########################################################

gene_count <- 0
mrna_count <- 0
exon_count <- 0
cds_count <- 0
V_gene_segment_count <- 0
C_gene_segment_count <- 0
match_count <- 0
or_count <- 0
D_loop_count <- 0

Sort_Unified_Data <- function(Unified_Models){
  
Final_gene_ls <- list()
  
for (u in 1 : length(unique(Unified_Models$Final_Mod_Attr))){

print(u)  
  
part_gene_d <- Unified_Models %>% filter(Final_Mod_Attr %in% unique(Unified_Models$Final_Mod_Attr)[u])
  
filtered_strand_ls <- list()  

for (s in 1 : length(unique(part_gene_d$strand))){
  
  part_gene <- part_gene_d %>% filter(strand %in% unique(part_gene_d$strand)[s])
  
  if(unique(part_gene$strand) == "+"){
    
    if(nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | 
                      part_gene$type == "ncRNA_gene",]) < 2 &
       nrow(part_gene[part_gene$type == "mRNA" | part_gene$type == "tRNA" | 
                      part_gene$type == "rRNA" | part_gene$type == "lnc_RNA" | part_gene$type == "snoRNA" | 
                      part_gene$type == "snRNA" | part_gene$type == "scRNA" | part_gene$type == "miRNA" |  
                      part_gene$type == "Y_RNA" | part_gene$type == "ncRNA" | part_gene$type == "guide_RNA" | 
                      part_gene$type == "primary_transcript" |
                      part_gene$type == "transcript" | part_gene$type == "pseudogenic_transcript",]) < 2){
      
      #if(nrow(part_gene_id) <= 1 & nrow(part_mRNA_id) <= 1){
      
      part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
      
      part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
                                             type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                             type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                             type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "primary_transcript" |
                                             type == "transcript" | type == "pseudogenic_transcript" ) 
      
      Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
        select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
        distinct(.keep_all = TRUE)
      
      Non_gene_f_c <- Non_gene_id %>% 
        group_by(type) %>% summarise(count = n())
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
        
        exon_i_count <- exon_count + 1  
        exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
        
        cds_i_count <- cds_count + 1   
        cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
        
        match_i_count <- match_count + 1
        match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
        
        V_gene_segment_i_count <- V_gene_segment_count + 1
        V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
        
        C_gene_segment_i_count <- C_gene_segment_count + 1
        C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
        
        or_i_count <- or_count + 1
        or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
        
        D_loop_i_count <- D_loop_count + 1
        D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
      }
      
      if(nrow(part_gene_id) == 1 & nrow(part_mRNA_id) <= 1){  
        
        if(nrow(part_gene_id) == 1){
          
          gene_count <<- gene_count + 1
          
          part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x)
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                 paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x  
                 
            ))
          
        }
        
        if(nrow(part_mRNA_id) == 1){
          
          mrna_count <<- mrna_count + 1
          
          part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x)
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                 paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x 
                 
            ))
          
          part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x)
            
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                 paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x  
                 
            ))
          
          #part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
          
          # paste0(x,";geneid=",paste0(unique(part_gene$Final_mod_attr),collapse = ";")) 
          
          #)
          
        }
        
        #### finding the start and end positions for customary gene/transcript model
        
        if(nrow(part_mRNA_id) == 1){
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
        }
        
        else if(nrow(part_mRNA_id) < 1 & nrow(part_gene_id) == 1){
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
        }
        
        if(nrow(Non_gene_id) >= 1){
          
          #if(length(parent_id) == 1){
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
          
          #}
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                   gsub("Parent=(.)*\\|","",x),x)  
            
          )
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            ifelse(str_detect(x,pattern = "\\|Parent="),
            gsub("\\|Parent=",";Parent=",x),x)
            
          )
          
          Non_gene_id$type_count <- NA
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
            
          }
          
          Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
            
            ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
            
          )
          
          Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
          
          Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
            
            paste0("ID=XB-",x)
            
          )
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            
            ifelse(str_detect(x,pattern = "ID=.+;"),
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),
            pattern = ";")[[1]][1],"",x),paste0(";",x))
            
          )
          
          Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
          
        }
        
        Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_gene_id <- part_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_mRNA_id <- part_mRNA_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id,part_mRNA_id,Non_gene_id)
        
      }
      
      if(nrow(part_gene_id) == 0 & nrow(part_mRNA_id) <= 1){
        
        if(nrow(part_mRNA_id) == 1){
          
          mrna_count <<- mrna_count + 1
          
          part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x)
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                 paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x 
                 
            ))
          
          #part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x)
          
          # gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
          #  paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x  
          
          # ))
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
          if(nrow(Non_gene_id) >= 1){
            
            #if(length(parent_id) == 1){
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
            
            #}
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
              gsub("Parent=(.)*\\|","",x),x)  
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "\\|Parent="),
              gsub("\\|Parent=",";Parent=",x),x)
              
            )
            
            Non_gene_id$type_count <- NA
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
              
            }
            
            Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
              
              ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
              
            )
            
            Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
            
            Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
              
              paste0("ID=XB-",x)
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              
            ifelse(str_detect(x,pattern = "ID=.+;"),
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),
            pattern = ";")[[1]][1],"",x),paste0(";",x))
              
            )
            
            Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
            
          }
          
        }
        
        if(nrow(part_mRNA_id) == 0){
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            gsub(paste0(str_split(str_extract(x,pattern = "Parent=.+;"),
                                  pattern = ";")[[1]][1],";"),"",x)
            
          )
          
          if(nrow(Non_gene_id) >= 1){
            
            Non_gene_id$type_count <- NA
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
              
            }
            
            Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
              
              ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
              
            )
            
            Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
            
            Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
              
              paste0("ID=XB-",x)
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              
              ifelse(str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),paste0(";",x))
              
            )
            
            Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
            
          }
          
        }
        
        part_gene_id <- part_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_mRNA_id <- part_mRNA_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id,part_mRNA_id,Non_gene_id)
        
      }
      
    }  
    
    else if ((nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) > 1) | (
      nrow(part_gene[part_gene$type == "mRNA" | part_gene$type == "tRNA" | 
                     part_gene$type == "rRNA" | part_gene$type == "lnc_RNA" | part_gene$type == "snoRNA" | 
                     part_gene$type == "snRNA" | part_gene$type == "scRNA" | part_gene$type == "miRNA" |  
                     part_gene$type == "Y_RNA" | part_gene$type == "ncRNA" | part_gene$type == "guide_RNA" | 
                     part_gene$type == "primary_transcript" |
                     part_gene$type == "transcript" | part_gene$type == "pseudogenic_transcript",]) > 1)){
      
      if(nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
         nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1 & 
         nrow(part_gene[part_gene$type == "ncRNA_gene" | 
                        part_gene$type == "tRNA",]) == 0){
        
        Unified_ls <- list()
        
        part_gene_No <- part_gene %>% filter(type == "gene" | type == "pseudogene") %>% 
          group_by(type) %>% summarise( Count = n(), source_t = paste0(source,collapse = ";"),
                                        mod_attr = paste0(unique(modified_attributes),collapse = ";") )
        
        for(i in 1 : length(part_gene_No$Count)){
          
          if(part_gene_No$Count[i] == 1){
            
            gene_count <<- gene_count + 1
            
            source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
            
            part_gene_Mod_attr <- unique(str_split(part_gene_No$mod_attr[i],pattern = ";")[[1]])
            
            part_gene_type <- part_gene_No[i,]$type    
            
            part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
            
            part_trans_M <- part_gene_M %>% filter(type == "mRNA" | type == "tRNA" | 
                                                     type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | 
                                                     type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                     type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                     type == "transcript" | type == "pseudogenic_transcript" | 
                                                     type == "primary_transcript")
            
            part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df <- part_gene_M %>% filter(type == part_gene_type) %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes,Overlap_Value
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df <- part_gene_M_df[which.max(part_gene_M_df$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df$type <- part_gene_type
            
            part_gene_M_df$start <- min(part_gene_M$start)
            part_gene_M_df$end <- max(part_gene_M$end)
            part_gene_M_df$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
            
            part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) 
              gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                   paste0("ID=XBXT10g",
                          paste0(rep(0,6-nchar(gene_count)),collapse = ""),
                          gene_count),x))
            
            part_gene_M_df$attributes <- sapply(
              
              part_gene_M_df$attributes,function(x)
                paste0(x,";geneid=",paste0(unique(part_gene_M$Final_Mod_Attr),
                                           collapse = ";"))  
              
            )
            
            Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | 
                                                       type %in% part_trans_M$type)) %>% select(
                                                         chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                                                       ) %>% distinct(.keep_all = TRUE)
            
            Non_gene_f_c <- Non_gene_id %>% 
              group_by(type) %>% summarise(count = n())
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              exon_i_count <- exon_count + 1  
              exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              cds_i_count <- cds_count + 1    
              cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              match_i_count <- match_count + 1
              match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              V_gene_segment_i_count <- V_gene_segment_count + 1
              V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              C_gene_segment_i_count <- C_gene_segment_count + 1
              C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              or_i_count <- or_count + 1
              or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              D_loop_i_count <- D_loop_count + 1
              D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
              
            }
            
            if(nrow(part_trans_M_df) == 1){
              
              mrna_count <<- mrna_count + 1
              
              part_trans_M_df$start <- min(part_gene_M$start)
              part_trans_M_df$end <- max(part_gene_M$end)
              part_trans_M_df$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
              
              part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                     paste0("Parent=XBXT10g",
                            paste0(rep(0,6-nchar(gene_count)),collapse = ""),
                            gene_count),x))
              
              part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                     paste0("ID=XB-mRNA",
                            paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
                            mrna_count),x))
              
              part_trans_M_df$attributes <- sapply(
                
                part_trans_M_df$attributes,function(x)
                  paste0(x,";geneid=",paste0(unique(part_gene_M$Final_Mod_Attr),
                                             collapse = ";"))  
                
              )
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
              
              if(nrow(Non_gene_id) >= 1){  
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                       paste0("Parent=",parent_id),x))
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x))
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(unique(part_gene_M$Final_Mod_Attr),
                                               collapse = ";"))  
                  
                )
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M_df <- part_gene_M_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_trans_M_df <- part_trans_M_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                       strand,phase,Final_Mod_Attr,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            else if(nrow(part_trans_M_df) != 1){
              
              #mrna_count <- mrna_count + 1
              
              #part_trans_M_df <- part_gene_M_df
              
              #part_trans_M_df$type <- "transcript"
              
              #part_trans_M_df$start <- min(part_gene_M$start)
              #part_trans_M_df$end <- max(part_gene_M$end)
              #part_trans_M_df$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
              
              #part_trans_M_df$attributes <- 
              # paste0("ID=XB-mRNA",
              #paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
              #mrna_count,";",
              #"Parent=XBXT10g",
              #paste0(rep(0,6-nchar(gene_count)),collapse = ""),
              #gene_count,";")
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x))
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                  paste0(x,";geneid=",paste0(unique(part_gene_M$Final_Mod_Attr),
                  collapse = ";"))  
                  
                )
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M_df <- part_gene_M_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_trans_M_df <- part_trans_M_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                       strand,phase,Final_Mod_Attr,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            Unified_ls[[i]] <- part_gene_M %>% select(chrloc,source,type,start,end,score,
                                                      strand,phase,Final_Mod_Attr,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }   
          
          else if(part_gene_No$Count[i] > 1){
            
            gene_count <<- gene_count + 1
            
            source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
            
            part_gene_type <- part_gene_No[i,]$type    
            
            part_gene_Mod_attr <- unique(part_gene[
              part_gene$type == part_gene_type,]$modified_attributes)
            
            part_gene_Mod <- part_gene %>% filter((source %in% source_tot ) & (modified_attributes %in% part_gene_Mod_attr))
            
            part_gene_id <- part_gene_Mod %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
            
            part_mRNA_id <- part_gene_Mod %>% filter(type == "mRNA" | type == "tRNA"| 
                                                       type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                       type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                       type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                       type == "transcript" | type == "pseudogenic_transcript" | 
                                                       type == "primary_transcript")
            
            Non_gene_id <-  part_gene_Mod %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
              select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) 
            
            Non_gene_f_c <- Non_gene_id %>% 
              group_by(type) %>% summarise(count = n())
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              exon_i_count <- exon_count + 1
              exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              cds_i_count <- cds_count + 1    
              cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              match_i_count <- match_count + 1
              match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              V_gene_segment_i_count <- V_gene_segment_count + 1
              V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              C_gene_segment_i_count <- C_gene_segment_count + 1
              C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              or_i_count <- or_count + 1
              or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              D_loop_i_count <- D_loop_count + 1
              D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
              
            }
            
            
            #### finding the start and end positions for customary gene/transcript model
            
            part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
              select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) 
            
            part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
              select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) 
            
            part_gene_id_df$type <- part_gene_type
            
            #part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x) 
            
            # gsub("ID=","gID=",x)
            
            #)
            
            part_gene_id_df$start <- min(part_gene_Mod$start)
            
            part_gene_id_df$end <- max(part_gene_Mod$end)
            
            part_gene_id_df$Final_Mod_Attr <- paste0(unique(part_gene_Mod$Final_Mod_Attr),collapse = ";")
            
            part_gene_id_df$attributes <- 
              sapply(part_gene_id_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                     paste0("ID=XBXT10g",
                            paste0(rep(0,6-nchar(gene_count)),collapse = ""),
                            gene_count),x))
            
            part_gene_id_df$attributes <- sapply(
              
              part_gene_id_df$attributes,function(x)
                paste0(x,";geneid=",paste0(unique(part_gene_Mod$Final_Mod_Attr),
                                           collapse = ";"))  
              
            )
            
            if(nrow(part_mRNA_id_df) >= 1){
              
              #part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
              
              # gsub("ID=","transID=",x)
              
              #)
              
              #part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
              
              # gsub("Parent=","Origin=",x)
              
              #)
              
              mrna_count <<- mrna_count + 1
              
              part_mRNA_id_type <- part_mRNA_id_df$type
              
              part_mRNA_id_df$start <- min(part_gene_Mod$start)
              
              part_mRNA_id_df$end <- max(part_gene_Mod$end)
              
              #part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
              
              part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                     paste0("ID=XB-mRNA",
                            paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
                            mrna_count),x))  
              
              part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                     paste0("Parent=XBXT10g",
                            paste0(rep(0,6-nchar(gene_count)),collapse = ""),
                            gene_count),x))
              
              part_gene_id_df$attributes <- sapply(
                
                part_gene_id_df$attributes,function(x)
                  paste0(x,";geneid=",paste0(unique(part_gene_Mod$Final_Mod_Attr),
                                             collapse = ";"))  
                
              )
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                       paste0("Parent=",parent_id),x))
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                  gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                  pattern = ";")[[1]][1],"",x),paste0(";",x))
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id$Final_Mod_Attr <- paste0(unique(part_gene_Mod$Final_Mod_Attr),collapse = ";")
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(unique(part_gene_Mod$Final_Mod_Attr),
                                               collapse = ";"))  
                  
                )
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                       strand,phase,Final_Mod_Attr,attributes) %>% 
                distinct(.keep_all = TRUE)
              
              
            }
            
            else if(nrow(part_mRNA_id_df) != 1){
              
              #mrna_count <- mrna_count + 1
              
              #part_mRNA_id_df <- part_gene_id_df
              
              #part_mRNA_id_df$type <- "transcript"
              
              #part_mRNA_id_df$start <- min(part_gene_Mod$start)
              #part_mRNA_id_df$end <- max(part_gene_Mod$end)
              #part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
              
              #part_mRNA_id_df$attributes <- 
              #paste0("ID=XB-mRNA",
              #paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count,";","Parent=XBXT10g",
              #paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count,";")
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                       paste0("Parent=",parent_id),x))
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x))
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id$Final_Mod_Attr <- paste0(unique(part_gene_Mod$Final_Mod_Attr),collapse = ";")
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(unique(part_gene_Mod$Final_Mod_Attr),
                                               collapse = ";"))  
                  
                )
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                       strand,phase,Final_Mod_Attr,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            Unified_ls[[i]] <- part_gene_M %>%
              select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
        }
        
        filtered_strand_ls[[s]] <- do.call("rbind.data.frame",Unified_ls) 
        
      }
      
      else if(((nrow(part_gene[part_gene$type == "ncRNA_gene",]) > 1) &
               (length(unique(part_gene[part_gene$type == "ncRNA_gene",]$type) %in% 
                       c("ncRNA_gene")) == 1 )) | (nrow(part_gene[part_gene$type == "gene" | 
                                                                  part_gene$type == "pseudogene",]) >= 1 & 
                                                   nrow(part_gene[part_gene$type == "ncRNA_gene",]) >= 1) |
              (nrow(part_gene[part_gene$type == "tRNA",]) >= 1)) {
        
        #other_ls[[length(other_ls)+ 1]] <- part_gene
        
        #part_df <- part_gene %>% left_join(rank_data,by = c("type"))
        
        if(length(unique(part_gene[part_gene$type == "tRNA",]$type)) >= 1 | 
           ((length(unique(part_gene[part_gene$type == "tRNA",]$type)) >= 1) & 
            (length(unique(part_gene[part_gene$type == "ncRNA_gene" | 
                                     part_gene$type == "miRNA",]$type)) >= 1))){
          
          part_gene$ID <- sapply(part_gene$attributes,function(x) gsub("ID=","",
                                                                       str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
          
          part_gene$Parent <- sapply(part_gene$attributes,function(x) gsub("Parent=","",
                                                                           str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))
          
          filtered_rna_gene_ids <- part_gene %>% filter(type == "mRNA" |
                                                          type == "tRNA"| type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" |
                                                          type == "snRNA" | type == "scRNA" | type == "miRNA" | type == "Y_RNA" |
                                                          type == "ncRNA" | type == "guide_RNA" | type == "transcript" | 
                                                          type == "primary_transcript" |type == "pseudogenic_transcript") %>% select(ID,Parent)
          
          part_d <- part_gene %>% left_join(filtered_rna_gene_ids,by = c("Parent" = "ID"))
          
          part_d <- part_d %>% mutate(
            
            transcript = ifelse(is.na(Parent.y),Parent,Parent.y),
            final_transcript = ifelse(is.na(transcript),ID,transcript)
            
            
          ) %>% arrange(final_transcript) %>% select(chrloc,source,type,start,end,
           score,strand,phase,modified_attributes,attributes,Final_Mod_Attr,final_transcript,Overlap_Value)
          
          filtered_mod_attr_ls <- list()
          final_mod_attr_ls <- list()
          
          for (l in 1 : length(unique(part_d$final_transcript))) {
            
            filtered_mod_attr_ls[[l]] <- part_d %>% filter(final_transcript %in% unique(part_d$final_transcript)[l]) 
            
            gene_filtered_mod_attr <- filtered_mod_attr_ls[[l]] %>% filter(type == "gene" | 
                                                                             type == "pseudogene" | type == "ncRNA_gene")
            
            transcript_filtered_mod_attr <- filtered_mod_attr_ls[[l]] %>% filter(type == "mRNA" |
                                                                                   type == "tRNA"| type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" |
                                                                                   type == "snRNA" | type == "scRNA" | type == "miRNA" | type == "Y_RNA" |
                                                                                   type == "ncRNA" | type == "guide_RNA" | type == "transcript" | 
                                                                                   type == "primary_transcript" |type == "pseudogenic_transcript")
            
            gene_filtered_mod_attr_df <- gene_filtered_mod_attr[which.max(gene_filtered_mod_attr$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            )
            
            transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr[which.max(transcript_filtered_mod_attr$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            )
            
            Non_gene_id <- filtered_mod_attr_ls[[l]] %>% filter(!(type %in% gene_filtered_mod_attr$type | type %in% transcript_filtered_mod_attr$type)) %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            )
            
            Non_gene_f_c <- Non_gene_id %>% 
              group_by(type) %>% summarise(count = n())
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              exon_i_count <- exon_count + 1 
              exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              cds_i_count <- cds_count + 1   
              cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              match_i_count <- match_count + 1
              match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              V_gene_segment_i_count <- V_gene_segment_count + 1
              V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              C_gene_segment_i_count <- C_gene_segment_count + 1
              C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              or_i_count <- or_count + 1
              or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              D_loop_i_count <- D_loop_count + 1
              D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
              
            }
            
            if(nrow(gene_filtered_mod_attr_df) == 1){
              
              gene_count <<- gene_count + 1
              
              gene_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
              gene_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
              
              gene_filtered_mod_attr_df$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
              
              #gene_filtered_mod_attr$attributes <- sapply(gene_filtered_mod_attr$attributes, function(x)
              
              # gsub("ID=","gID=",x)
              
              #)
              
              gene_filtered_mod_attr_df$attributes <- sapply(
                
                gene_filtered_mod_attr_df$attributes,function(x)
                  gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                       paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),
                       x))
              
              gene_filtered_mod_attr_df$attributes <- sapply(
                
                gene_filtered_mod_attr_df$attributes,function(x)
                  paste0(x,";geneid=",paste0(filtered_mod_attr_ls[[l]]$Final_Mod_Attr,
                                             collapse = ";"))  
                
              )
              
              if(nrow(transcript_filtered_mod_attr_df) == 1 ){
                
                mrna_count <<- mrna_count + 1
                
                transcript_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                transcript_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                
                transcript_filtered_mod_attr_df$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                
                transcript_filtered_mod_attr_df$attributes <-  sapply(
                  
                  transcript_filtered_mod_attr_df$attributes,function(x)
                    gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                         paste0("ID=XB-RNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),
                         x))
                
                transcript_filtered_mod_attr_df$attributes <-  sapply(
                  
                  transcript_filtered_mod_attr_df$attributes,function(x)
                    gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                         paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),
                         x))
                
                transcript_filtered_mod_attr_df$attributes <- sapply(
                  
                  transcript_filtered_mod_attr_df$attributes,function(x)
                    paste0(x,";geneid=",paste0(filtered_mod_attr_ls[[l]]$Final_Mod_Attr,
                                               collapse = ";"))  
                  
                )
                
                parent_id <- gsub("ID=","",str_split(str_extract(
                  transcript_filtered_mod_attr_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])   
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                           gsub("Parent=(.)*\\|","",x),x)  
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "\\|Parent="),
                           gsub("\\|Parent=",";Parent=",x),x)
                    
                  )
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                    
                    
                    ifelse(str_detect(x,pattern = "ID=.+;"),
                           gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                          pattern = ";")[[1]][1],"",x),paste0(";",x)    
                           
                    )
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                  
                  Non_gene_id$attributes <- sapply(
                    
                    Non_gene_id$attributes,function(x)
                      paste0(x,";geneid=",paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),
                                                 collapse = ";"))  
                    
                  )
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                gene_filtered_mod_attr_df <- gene_filtered_mod_attr_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                final_mod_attr_ls[[l]] <- rbind.data.frame(gene_filtered_mod_attr_df,transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                  
                )
                
              }
              
              else if(nrow(transcript_filtered_mod_attr_df) < 1){
                
                parent_id <- gsub("ID=","",str_split(str_extract(gene_filtered_mod_attr_df$attributes,pattern = "ID=.+;"))[[1]][1])   
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                           gsub("Parent=(.)*\\|","",x),x)  
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "\\|Parent="),
                           gsub("\\|Parent=",";Parent=",x),x)
                    
                  )
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                    
                    
                    ifelse(str_detect(x,pattern = "ID=.+;"),
                           gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                          pattern = ";")[[1]][1],"",x),paste0(";",x)   
                           
                    )
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                  
                  Non_gene_id$attributes <- sapply(
                    
                    Non_gene_id$attributes,function(x)
                      paste0(x,";geneid=",paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),
                                                 collapse = ";"))  
                    
                  )
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                gene_filtered_mod_attr_df <- gene_filtered_mod_attr_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                final_mod_attr_ls[[l]] <- rbind.data.frame(gene_filtered_mod_attr_df,transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                  
                )
                
              }
              
            }
            
            if(nrow(gene_filtered_mod_attr_df) < 1){
              
              if( nrow(transcript_filtered_mod_attr_df) == 1 ){
                
                mrna_count <<- mrna_count + 1
                
                transcript_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                transcript_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                
                transcript_filtered_mod_attr_df$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                
                #transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                
                # gsub("ID=","transID=",x)
                
                #)
                
                #transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                
                # gsub("Parent=","Origin=",x)
                
                #)
                
                transcript_filtered_mod_attr_df$attributes <- 
                  
                  sapply(transcript_filtered_mod_attr_df$attributes, function(x)
                    
                    gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x)    
                    
                  )
                
                transcript_filtered_mod_attr_df$attributes <- sapply(
                  
                  transcript_filtered_mod_attr_df$attributes,function(x)
                    paste0(x,";geneid=",paste0(filtered_mod_attr_ls[[l]]$Final_Mod_Attr,
                                               collapse = ";"))  
                  
                )
                
                parent_id <- gsub("ID=","",str_split(str_extract(
                  transcript_filtered_mod_attr_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])   
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                           gsub("Parent=(.)*\\|","",x),x)  
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "\\|Parent="),
                           gsub("\\|Parent=",";Parent=",x),x)
                    
                  )
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                    
                    
                    ifelse(str_detect(x,pattern = "ID=.+;"),
                           gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                          pattern = ";")[[1]][1],"",x),paste0(";",x)    
                           
                    )
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$attributes <- sapply(
                    
                    Non_gene_id$attributes,function(x)
                      paste0(x,";geneid=",paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),
                                                 collapse = ";"))  
                    
                  )
                  
                  Non_gene_id$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                final_mod_attr_ls[[l]] <- rbind.data.frame(transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                  
                )
                
              }
              
              else if(nrow(transcript_filtered_mod_attr_df) < 1){
                
                #parent_id <- gsub("ID=","",str_split(str_extract(gene_filtered_mod_attr_df$attributes,pattern = "ID=.+;"))[[1]][1])   
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  gsub(paste0(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],";"),"",x))
                
                if(nrow(Non_gene_id) >= 1){
                  
                  #Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  #gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                    
                    
                    ifelse(str_detect(x,pattern = "ID=.+;"),
                           gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                          pattern = ";")[[1]][1],"",x),paste0(";",x)    
                           
                    )
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                  
                  #Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  Non_gene_id$attributes <- sapply(
                    
                    Non_gene_id$attributes,function(x)
                      paste0(x,";geneid=",paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),
                                                 collapse = ";"))  
                    
                  )
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                final_mod_attr_ls[[l]] <- rbind.data.frame(Non_gene_id) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                  
                )
                
              }
              
            }  
            
          }
          
          final_mod_attr_df <- do.call("rbind.data.frame",final_mod_attr_ls)
          
          filtered_strand_ls[[s]] <- final_mod_attr_df %>% select(chrloc,source,type,start,end,score,
                                                                  strand,phase,Final_Mod_Attr,attributes) %>% 
            distinct(.keep_all = TRUE)
          
        }
        
        if((length(unique(part_gene[part_gene$type == "ncRNA_gene" | 
                                    part_gene$type == "miRNA",]$type)) >= 1 & nrow(part_gene[part_gene$type == 'tRNA',]) == 0) | 
           (nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene",]) >= 1 & 
            nrow(part_gene[part_gene$type == "ncRNA_gene",]) > 1)){
          
          #part_df <- part_df %>% arrange(rank)
          
          final_gb_ls <- list()
          
          if (nrow(part_gene[part_gene$source == "Genbank",]) >= 1){
            
            filtered_genbank_ls <- part_gene %>% filter(source == "Genbank")
            
            filtered_gb_ls <- list()
            
            for (h in 1 : length(unique(filtered_genbank_ls$modified_attributes))) {
              
              
              filtered_gb_ls[[h]] <- filtered_genbank_ls %>% filter(Final_Mod_Attr == unique(filtered_genbank_ls$modified_attributes)[h]) %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                
              )
              
              
              if(nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]) >= 1){
                
                gene_count <<- gene_count + 1
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes,function(x) 
                  gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XBXT10g",
                                                                                                 paste0(rep(0,6-nchar(gene_count)),
                                                                                                        collapse = ""),gene_count),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes,function(x) 
                  paste0(x,";geneid=",paste0(filtered_gb_ls[[h]]$Final_Mod_Attr),collapse = ";"))
                
              }
              
              if(nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]) >= 1){
                
                mrna_count <<- mrna_count + 1
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-mRNA",
                                                                                                 paste0(rep(0,6-nchar(mrna_count)),
                                                                                                        collapse = ""),mrna_count),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;")),
                       paste0("Parent=XBXT10g",
                              paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  paste0(x,";geneid=",paste0(filtered_gb_ls[[h]]$Final_Mod_Attr),collapse = ";"))
                
                parent_id <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))  
                
              }
              
              else if (nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]) < 1){
                
                
                parent_id <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))  
                
              }
              
              Non_gene_id <- filtered_gb_ls[[h]] %>% filter(!(type == "gene" | type == "mRNA"))
              
              Non_gene_f_c <- Non_gene_id %>% 
                group_by(type) %>% summarise(count = n())
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                exon_i_count <- exon_count + 1  
                exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                cds_i_count <- cds_count + 1    
                cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                match_i_count <- match_count + 1
                match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                V_gene_segment_i_count <- V_gene_segment_count + 1
                V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                C_gene_segment_i_count <- C_gene_segment_count + 1
                C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                or_i_count <- or_count + 1
                or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                D_loop_i_count <- D_loop_count + 1
                D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                
              }
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(filtered_gb_ls[[l]]$Final_Mod_Attr,
                                               collapse = ";"))  
                  
                )
                
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x)    
                         
                  )
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              filtered_gene_df <- filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",] %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
              filtered_trans_df <- filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",] %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
              
              final_gb_ls[[h]] <- rbind.data.frame(filtered_gene_df,filtered_trans_df,
              Non_gene_id) %>% select(
                                                     
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes                                       
                                                     
                )
              
            }
            
          }
          
          final_gb_df <- do.call("rbind.data.frame",final_gb_ls)
          
          part_gene <- part_gene %>% filter(!(source %in% "Genbank"))
          
          part_gene$transcript <- sapply(part_gene$attributes,function(x) 
            
            ifelse(str_detect(x,pattern = "Parent=.+"),
                   
                   gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1]),
                   
                   gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
            
          )
          
          part_gene <- part_gene %>% arrange(transcript)
          
          filtered_gene_ls <- list()
          filtered_M_ls <- list()
          
          for(k in 1: length(unique(part_gene$transcript))){
            
            gene_count <<- gene_count + 1
            
            filtered_gene_ls[[k]] <- part_gene %>% filter(transcript %in% unique(part_gene$transcript)[k]) %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            )
            
            filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | 
                                    filtered_gene_ls[[k]]$type == "pseudogene" | 
                                    filtered_gene_ls[[k]]$type == "gene",]$attributes <- 
              sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | 
                                             filtered_gene_ls[[k]]$type == "pseudogene" | 
                                             filtered_gene_ls[[k]]$type == "gene",]$attributes,function(x) 
                                               gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                                                    paste0("ID=XBXT10g",
                                                           paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x))
            
            
            if(nrow(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                          filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                          filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                          filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                          filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                          filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]) >= 1){                                                                                                                                              
              
              mrna_count <<- mrna_count + 1
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | filtered_gene_ls[[k]]$type == "miRNA" | 
                                      filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                      filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                      filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" |
                                      filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                      filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes,function(x) 
                                                                                                                                                                                 gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-rna",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x))
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                      filtered_gene_ls[[k]]$type == "miRNA" | 
                                      filtered_gene_ls[[k]]$type == "rRNA" | 
                                      filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                      filtered_gene_ls[[k]]$type == "snoRNA" | 
                                      filtered_gene_ls[[k]]$type == "scRNA" | 
                                      filtered_gene_ls[[k]]$type == "Y_RNA" | 
                                      filtered_gene_ls[[k]]$type == "ncRNA" | 
                                      filtered_gene_ls[[k]]$type == "guide_RNA" |
                                      filtered_gene_ls[[k]]$type == "transcript" | 
                                      filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "miRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "rRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "snoRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "Y_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "guide_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "transcript" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes,function(x) 
                                                                                                                                                                                 gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=XBXT10g",
                                                                                                                                                                                                                                                                    paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x))
              
              
              parent_id <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                          filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                          filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                          filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                          filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                          filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
              
            }
            
            else if (nrow(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]) < 1){
              
              parent_id <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
              
              
            }
            
            Non_gene_id <- filtered_gene_ls[[k]] %>% filter(!(type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | type == "snRNA" | 
                                                                type == "scRNA" | type == "miRNA" | type == "transcript" | type == "primary_transcript" |
                                                                type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "ncRNA_gene" | type == "mRNA" | type == "gene" | type == "pseudogene"
                                                              
            ))
            
            Non_gene_f_c <- Non_gene_id %>% 
              group_by(type) %>% summarise(count = n())
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              exon_i_count <- exon_count + 1  
              exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              cds_i_count <- cds_count + 1    
              cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              match_i_count <- match_count + 1
              match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              V_gene_segment_i_count <- V_gene_segment_count + 1
              V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              C_gene_segment_i_count <- C_gene_segment_count + 1
              C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              or_i_count <- or_count + 1
              or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              D_loop_i_count <- D_loop_count + 1
              D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
              
            }
            
            if(nrow(Non_gene_id) >= 1){
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x)
                
              )
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                       gsub("Parent=(.)*\\|","",x),x)  
                
              )
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                ifelse(str_detect(x,pattern = "\\|Parent="),
                       gsub("\\|Parent=",";Parent=",x),x)
                
              )
              
              Non_gene_id$type_count <- NA
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                
              }
              
              Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                
                ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                
              )
              
              Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
              
              Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                
                paste0("ID=XB-",x)
                
              )
              
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                
                
                ifelse(str_detect(x,pattern = "ID=.+;"),
                       gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                      pattern = ";")[[1]][1],"",x),paste0(";",x)   
                       
                )
                
              )
              
              Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
              
            }
            
            Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
            
            filtered_gene_df <- filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",] %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
            
            filtered_trans_df <- filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | filtered_gene_ls[[k]]$type == "miRNA" | 
                                                         filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                         filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                         filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" |
                                                         filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                         filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
            
            filtered_M_ls[[k]] <- rbind.data.frame(filtered_gene_df,filtered_trans_df,Non_gene_id) %>% select(
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
            )
            
          }
          
          Final_mod_d <- do.call("rbind.data.frame",filtered_M_ls)
          
          if(nrow(final_gb_df) >= 1){
            
            Final_mod_d <- rbind.data.frame(final_gb_df,Final_mod_d)
            
          }
          
          filtered_strand_ls[[s]] <- Final_mod_d %>% select(chrloc,source,type,start,end,score,
                                                            strand,phase,Final_Mod_Attr,attributes) %>% 
            distinct(.keep_all = TRUE)
          
        }
        
      } 
      
      
      else if(((nrow(part_gene[part_gene$type == "gene",]) >= 2) &
               (length(unique(part_gene[part_gene$type == "gene",]$type) %in% 
                       c("gene")) == 1 ) &
               nrow(part_gene[part_gene$type == "tRNA",]) == 0) |
              ((nrow(part_gene[part_gene$type == "pseudogene",]) >= 2) & 
               (length(unique(part_gene[part_gene$type == "pseudogene",]$type) %in%  
                       c("pseudogene")) == 1))){ 
        
        #else if(nrow(part_gene_id) > 1 & nrow(part_mRNA_id) >= 0){
        part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
        
        part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
                                               type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                               type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                               type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "primary_transcript" |
                                               type == "transcript" | type == "pseudogenic_transcript" ) 
        
        Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
          select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
          distinct(.keep_all = TRUE)
        
        Non_gene_f_c <- Non_gene_id %>% 
          group_by(type) %>% summarise(count = n())
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
          
          exon_i_count <- exon_count + 1  
          exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
          
          cds_i_count <- cds_count + 1   
          cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
          
          match_i_count <- match_count + 1
          match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
          
          V_gene_segment_i_count <- V_gene_segment_count + 1
          V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
          
          C_gene_segment_i_count <- C_gene_segment_count + 1
          C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
          
          or_i_count <- or_count + 1
          or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
          
          D_loop_i_count <- D_loop_count + 1
          D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
        }
        
        if(nrow(part_gene_id) >= 1){
          
          gene_count <<- gene_count + 1
          
          part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
            select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
            distinct(.keep_all = TRUE)
          
          part_gene_id_df$start <- min(part_gene$start)
          
          part_gene_id_df$end <- max(part_gene$end)
          
          part_gene_id_df$attributes <- sapply(part_gene_id_df$attributes,function(x) 
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x)
            
          )
          
          part_gene_id_df$attributes <- sapply(part_gene_id_df$attributes,function(x) 
            
            paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")) 
            
          )
          
          #part_gene_id_df$modified_attributes <- paste0(unique(part_gene$Final_mod_attr),collapse = ";")
          
        }
        
        if(nrow(part_mRNA_id) >= 1){
          
          mrna_count <<- mrna_count + 1
          
          part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
            select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
            distinct(.keep_all = TRUE)
          
        }
        
        #### finding the start and end positions for customary gene/transcript model
        
        if(nrow(part_mRNA_id_df) == 1){
          
          part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x)
            
          )
          
          part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
            
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x)
            
          )
          
          part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
            
            paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")) 
            
          )
          
          part_mRNA_id_df$start <- min(part_gene_id_df$start)
          
          part_mRNA_id_df$end <- max(part_gene_id_df$end)
          
          #part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$Final_Mod_attr),collapse = ";")    
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
          if(nrow(Non_gene_id) >= 1){
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
              gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
              gsub("Parent=(.)*\\|","",x),x)  
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "\\|Parent="),
                     gsub("\\|Parent=",";Parent=",x),x)
              
            )
            
            Non_gene_id$type_count <- NA
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
              
            }
            
            Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
              
              ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
              
            )
            
            Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
            
            Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
              
              paste0("ID=XB-",x)
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(!str_detect(x,pattern = "ID=(.)+\\|(.)+;") & 
                       str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),ifelse(
                                      str_detect(x,pattern = "ID=(.)+\\|(.)+;"),
                                      gsub("\\|",";",x),paste0(";",x)))
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
              
              
              ifelse(str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),x    
                     
              )
              
            )
            
            Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
            
            Non_gene_id$attributes <- sapply(
              
              Non_gene_id$attributes,function(x)
                paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),
                                           collapse = ";"))  
              
            )
            
          }
          
          Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
            select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% 
            distinct(.keep_all = TRUE) 
          
        }
        
        #}  
        
        else if(nrow(part_mRNA_id_df) < 1){
          
          #mrna_count <- mrna_count + 1
          
          #part_mRNA_id_df <- part_gene_id_df 
          
          #part_mRNA_id_df$type <- "transcript"
          
          #part_mRNA_id_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
          # mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
          # mrna_count,";")
          
          #part_mRNA_id_df$start <- part_gene_id_df$start
          
          #part_mRNA_id_df$end <- part_gene_id_df$end 
          
          #part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$Final_mod_attr),collapse = ";")
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
          if(nrow(Non_gene_id) >= 1){
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
              gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                     gsub("Parent=(.)*\\|","",x),x)  
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "\\|Parent="),
                     gsub("\\|Parent=",";Parent=",x),x)
              
            )
            
            Non_gene_id$type_count <- NA
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
              
            }
            
            Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
              
              ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
              
            )
            
            Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
            
            Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
              
              paste0("ID=XB-",x)
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
              
              
              ifelse(str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),paste0(";",x)    
                     
              )
              
            )
            
            Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x) 
              
              paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")) 
              
            )
            
            #Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
            
          }
          
          Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id)
          
        }
        
      }
      
      else if((nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | 
                              part_gene$type == "ncRNA_gene",]) == 0) & 
              (nrow(part_gene[part_gene$type == "mRNA" | 
                              part_gene$type == "tRNA" | 
                              part_gene$type == "rRNA" | part_gene$type == "lnc_RNA" | 
                              part_gene$type == "snoRNA" | part_gene$type == "snRNA" | 
                              part_gene$type == "scRNA" | part_gene$type == "miRNA" |  
                              part_gene$type == "Y_RNA" | part_gene$type == "ncRNA" | 
                              part_gene$type == "guide_RNA" | 
                              part_gene$type == "primary_transcript" |
                              part_gene$type == "transcript" | 
                              part_gene$type == "pseudogenic_transcript",]) > 1)){
        
        #if(nrow(part_mRNA_id) >= 1){
        
        part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
        
        part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
                                               type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                               type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                               type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "primary_transcript" |
                                               type == "transcript" | type == "pseudogenic_transcript" ) 
        
        Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
          select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
          distinct(.keep_all = TRUE)
        
        Non_gene_f_c <- Non_gene_id %>% 
          group_by(type) %>% summarise(count = n())
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
          
          exon_i_count <- exon_count + 1  
          exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
          
          cds_i_count <- cds_count + 1   
          cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
          
          match_i_count <- match_count + 1
          match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
          
          V_gene_segment_i_count <- V_gene_segment_count + 1
          V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
          
          C_gene_segment_i_count <- C_gene_segment_count + 1
          C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
          
          or_i_count <- or_count + 1
          or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
          
          D_loop_i_count <- D_loop_count + 1
          D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
        }
        
        mrna_count <<- mrna_count + 1
        
        part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
          select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
          distinct(.keep_all = TRUE)
        
        part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
          
          gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x)
          
        )
        
        part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
          
          paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")) 
          
        )
        
        part_mRNA_id_df$start <- min(part_gene$start)
        
        part_mRNA_id_df$end <- max(part_gene$end)
        
        parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
        
        if(nrow(Non_gene_id) >= 1){
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                   gsub("Parent=(.)*\\|","",x),x)  
            
          )
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            ifelse(str_detect(x,pattern = "\\|Parent="),
                   gsub("\\|Parent=",";Parent=",x),x)
            
          )
          
          Non_gene_id$type_count <- NA
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
            
          }
          
          Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
            
            ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
            
          )
          
          Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
          
          Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
            
            paste0("ID=XB-",x)
            
          )
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
            
            
            ifelse(str_detect(x,pattern = "ID=.+;"),
                   gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                  pattern = ";")[[1]][1],"",x),paste0(";",x)    
                   
            )
            
          )
          
          Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
          
          Non_gene_id$attributes <- sapply(
            
            Non_gene_id$attributes,function(x)
              paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),
                                         collapse = ";"))  
            
          )
          
          #Non_gene_id$modified_attributes <- paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")    
          
        }
        
        Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) 
        
      }
      
    }  
    
  }

  if(unique(part_gene$strand) == "-"){
    
    if(nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | 
                      part_gene$type == "ncRNA_gene",]) < 2 &
       nrow(part_gene[part_gene$type == "mRNA" | part_gene$type == "tRNA" | 
                      part_gene$type == "rRNA" | part_gene$type == "lnc_RNA" | part_gene$type == "snoRNA" | 
                      part_gene$type == "snRNA" | part_gene$type == "scRNA" | part_gene$type == "miRNA" |  
                      part_gene$type == "Y_RNA" | part_gene$type == "ncRNA" | part_gene$type == "guide_RNA" | 
                      part_gene$type == "primary_transcript" |
                      part_gene$type == "transcript" | part_gene$type == "pseudogenic_transcript",]) < 2){
      
      #if(nrow(part_gene_id) <= 1 & nrow(part_mRNA_id) <= 1){
      
      part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
      
      part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
                                             type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                             type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                             type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "primary_transcript" |
                                             type == "transcript" | type == "pseudogenic_transcript" ) 
      
      Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
        select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
        distinct(.keep_all = TRUE)
      
      Non_gene_f_c <- Non_gene_id %>% 
        group_by(type) %>% summarise(count = n())
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
        
        exon_i_count <- exon_count + 1  
        exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
        
        cds_i_count <- cds_count + 1   
        cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
        
        match_i_count <- match_count + 1
        match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
        
        V_gene_segment_i_count <- V_gene_segment_count + 1
        V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
        
        C_gene_segment_i_count <- C_gene_segment_count + 1
        C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
        
        or_i_count <- or_count + 1
        or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
      }
      
      if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
        
        D_loop_i_count <- D_loop_count + 1
        D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
      }
      
      if(nrow(part_gene_id) == 1 & nrow(part_mRNA_id) <= 1){  
        
        if(nrow(part_gene_id) == 1){
          
          gene_count <<- gene_count + 1
          
          part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x)
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                 paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x  
                 
            ))
          
        }
        
        if(nrow(part_mRNA_id) == 1){
          
          mrna_count <<- mrna_count + 1
          
          part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x)
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                 paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x 
                 
            ))
          
          part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x)
            
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                 paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x  
                 
            ))
          
          #part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
          
          # paste0(x,";geneid=",paste0(unique(part_gene$Final_mod_attr),collapse = ";")) 
          
          #)
          
        }
        
        #### finding the start and end positions for customary gene/transcript model
        
        if(nrow(part_mRNA_id) == 1){
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
        }
        
        else if(nrow(part_mRNA_id) < 1 & nrow(part_gene_id) == 1){
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
        }
        
        if(nrow(Non_gene_id) >= 1){
          
          #if(length(parent_id) == 1){
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
          
          #}
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                   gsub("Parent=(.)*\\|","",x),x)  
            
          )
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            ifelse(str_detect(x,pattern = "\\|Parent="),
                   gsub("\\|Parent=",";Parent=",x),x)
            
          )
          
          Non_gene_id$type_count <- NA
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
            
          }
          
          Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
            
            ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
            
          )
          
          Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
          
          Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
            
            paste0("ID=XB-",x)
            
          )
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            
            ifelse(str_detect(x,pattern = "ID=.+;"),
                   gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                  pattern = ";")[[1]][1],"",x),paste0(";",x))
            
          )
          
          Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
          
        }
        
        Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_gene_id <- part_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_mRNA_id <- part_mRNA_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id,part_mRNA_id,Non_gene_id)
        
      }
      
      if(nrow(part_gene_id) == 0 & nrow(part_mRNA_id) <= 1){
        
        if(nrow(part_mRNA_id) == 1){
          
          mrna_count <<- mrna_count + 1
          
          part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x)
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                 paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x 
                 
            ))
          
          #part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x)
          
          # gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
          #  paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x  
          
          # ))
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
          if(nrow(Non_gene_id) >= 1){
            
            #if(length(parent_id) == 1){
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
              gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
            
            #}
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                     gsub("Parent=(.)*\\|","",x),x)  
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "\\|Parent="),
                     gsub("\\|Parent=",";Parent=",x),x)
              
            )
            
            Non_gene_id$type_count <- NA
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
              
            }
            
            Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
              
              ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
              
            )
            
            Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
            
            Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
              
              paste0("ID=XB-",x)
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              
              ifelse(str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),paste0(";",x))
              
            )
            
            Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
            
          }
          
        }
        
        if(nrow(part_mRNA_id) == 0){
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            gsub(paste0(str_split(str_extract(x,pattern = "Parent=.+;"),
                                  pattern = ";")[[1]][1],";"),"",x)
            
          )
          
          if(nrow(Non_gene_id) >= 1){
            
            Non_gene_id$type_count <- NA
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
              
            }
            
            Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
              
              ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
              
            )
            
            Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
            
            Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
              
              paste0("ID=XB-",x)
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              
              ifelse(str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),paste0(";",x))
              
            )
            
            Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
            
          }
          
        }
        
        part_gene_id <- part_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_mRNA_id <- part_mRNA_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id,part_mRNA_id,Non_gene_id)
        
      }
      
    }  
    
    else if ((nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | part_gene$type == "ncRNA_gene",]) > 1) | (
      nrow(part_gene[part_gene$type == "mRNA" | part_gene$type == "tRNA" | 
                     part_gene$type == "rRNA" | part_gene$type == "lnc_RNA" | part_gene$type == "snoRNA" | 
                     part_gene$type == "snRNA" | part_gene$type == "scRNA" | part_gene$type == "miRNA" |  
                     part_gene$type == "Y_RNA" | part_gene$type == "ncRNA" | part_gene$type == "guide_RNA" | 
                     part_gene$type == "primary_transcript" |
                     part_gene$type == "transcript" | part_gene$type == "pseudogenic_transcript",]) > 1)){
      
      if(nrow(part_gene[part_gene$type == "gene",]) >= 1 & 
         nrow(part_gene[(part_gene$type == "pseudogene"),]) >= 1 & 
         nrow(part_gene[part_gene$type == "ncRNA_gene" | 
                        part_gene$type == "tRNA",]) == 0){
        
        Unified_ls <- list()
        
        part_gene_No <- part_gene %>% filter(type == "gene" | type == "pseudogene") %>% 
          group_by(type) %>% summarise( Count = n(), source_t = paste0(source,collapse = ";"),
                                        mod_attr = paste0(unique(modified_attributes),collapse = ";") )
        
        for(i in 1 : length(part_gene_No$Count)){
          
          if(part_gene_No$Count[i] == 1){
            
            gene_count <<- gene_count + 1
            
            source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
            
            part_gene_Mod_attr <- unique(str_split(part_gene_No$mod_attr[i],pattern = ";")[[1]])
            
            part_gene_type <- part_gene_No[i,]$type    
            
            part_gene_M <- part_gene %>% filter((source %in% source_tot) & (modified_attributes %in% part_gene_Mod_attr))
            
            part_trans_M <- part_gene_M %>% filter(type == "mRNA" | type == "tRNA" | 
                                                     type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | 
                                                     type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                     type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                     type == "transcript" | type == "pseudogenic_transcript" | 
                                                     type == "primary_transcript")
            
            part_trans_M_df <- part_trans_M[which.max(part_trans_M$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df <- part_gene_M %>% filter(type == part_gene_type) %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes,Overlap_Value
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df <- part_gene_M_df[which.max(part_gene_M_df$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            ) %>% distinct(.keep_all = TRUE)
            
            part_gene_M_df$type <- part_gene_type
            
            part_gene_M_df$start <- min(part_gene_M$start)
            part_gene_M_df$end <- max(part_gene_M$end)
            part_gene_M_df$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
            
            part_gene_M_df$attributes <- sapply(part_gene_M_df$attributes,function(x) 
              gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                   paste0("ID=XBXT10g",
                          paste0(rep(0,6-nchar(gene_count)),collapse = ""),
                          gene_count),x))
            
            part_gene_M_df$attributes <- sapply(
              
              part_gene_M_df$attributes,function(x)
                paste0(x,";geneid=",paste0(unique(part_gene_M$Final_Mod_Attr),
                                           collapse = ";"))  
              
            )
            
            Non_gene_id <- part_gene_M  %>% filter(!(type == part_gene_type | 
                                                       type %in% part_trans_M$type)) %>% select(
                                                         chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                                                       ) %>% distinct(.keep_all = TRUE)
            
            Non_gene_f_c <- Non_gene_id %>% 
              group_by(type) %>% summarise(count = n())
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              exon_i_count <- exon_count + 1  
              exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              cds_i_count <- cds_count + 1    
              cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              match_i_count <- match_count + 1
              match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              V_gene_segment_i_count <- V_gene_segment_count + 1
              V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              C_gene_segment_i_count <- C_gene_segment_count + 1
              C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              or_i_count <- or_count + 1
              or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              D_loop_i_count <- D_loop_count + 1
              D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
              
            }
            
            if(nrow(part_trans_M_df) == 1){
              
              mrna_count <<- mrna_count + 1
              
              part_trans_M_df$start <- min(part_gene_M$start)
              part_trans_M_df$end <- max(part_gene_M$end)
              part_trans_M_df$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
              
              part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                     paste0("Parent=XBXT10g",
                            paste0(rep(0,6-nchar(gene_count)),collapse = ""),
                            gene_count),x))
              
              part_trans_M_df$attributes <- sapply(part_trans_M_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                     paste0("ID=XB-mRNA",
                            paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
                            mrna_count),x))
              
              part_trans_M_df$attributes <- sapply(
                
                part_trans_M_df$attributes,function(x)
                  paste0(x,";geneid=",paste0(unique(part_gene_M$Final_Mod_Attr),
                                             collapse = ";"))  
                
              )
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_trans_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
              
              if(nrow(Non_gene_id) >= 1){  
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                       paste0("Parent=",parent_id),x))
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x))
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(unique(part_gene_M$Final_Mod_Attr),
                                               collapse = ";"))  
                  
                )
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M_df <- part_gene_M_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_trans_M_df <- part_trans_M_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                       strand,phase,Final_Mod_Attr,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            else if(nrow(part_trans_M_df) != 1){
              
              #mrna_count <- mrna_count + 1
              
              #part_trans_M_df <- part_gene_M_df
              
              #part_trans_M_df$type <- "transcript"
              
              #part_trans_M_df$start <- min(part_gene_M$start)
              #part_trans_M_df$end <- max(part_gene_M$end)
              #part_trans_M_df$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
              
              #part_trans_M_df$attributes <- 
              # paste0("ID=XB-mRNA",
              #paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
              #mrna_count,";",
              #"Parent=XBXT10g",
              #paste0(rep(0,6-nchar(gene_count)),collapse = ""),
              #gene_count,";")
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_gene_M_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x))
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id$Final_Mod_Attr <- paste0(unique(part_gene_M$Final_Mod_Attr),collapse = ";")
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(unique(part_gene_M$Final_Mod_Attr),
                                               collapse = ";"))  
                  
                )
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M_df <- part_gene_M_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_trans_M_df <- part_trans_M_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M <- rbind.data.frame(part_gene_M_df,part_trans_M_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                       strand,phase,Final_Mod_Attr,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            Unified_ls[[i]] <- part_gene_M %>% select(chrloc,source,type,start,end,score,
                                                      strand,phase,Final_Mod_Attr,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }   
          
          else if(part_gene_No$Count[i] > 1){
            
            gene_count <<- gene_count + 1
            
            source_tot <- unique(str_split(part_gene_No$source_t[i],pattern = ";")[[1]])
            
            part_gene_type <- part_gene_No[i,]$type    
            
            part_gene_Mod_attr <- unique(part_gene[
              part_gene$type == part_gene_type,]$modified_attributes)
            
            part_gene_Mod <- part_gene %>% filter((source %in% source_tot ) & (modified_attributes %in% part_gene_Mod_attr))
            
            part_gene_id <- part_gene_Mod %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
            
            part_mRNA_id <- part_gene_Mod %>% filter(type == "mRNA" | type == "tRNA"| 
                                                       type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                                       type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                                       type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | 
                                                       type == "transcript" | type == "pseudogenic_transcript" | 
                                                       type == "primary_transcript")
            
            Non_gene_id <-  part_gene_Mod %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
              select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) 
            
            Non_gene_f_c <- Non_gene_id %>% 
              group_by(type) %>% summarise(count = n())
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              exon_i_count <- exon_count + 1
              exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              cds_i_count <- cds_count + 1    
              cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              match_i_count <- match_count + 1
              match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              V_gene_segment_i_count <- V_gene_segment_count + 1
              V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              C_gene_segment_i_count <- C_gene_segment_count + 1
              C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              or_i_count <- or_count + 1
              or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              D_loop_i_count <- D_loop_count + 1
              D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
              
            }
            
            
            #### finding the start and end positions for customary gene/transcript model
            
            part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
              select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) 
            
            part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
              select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) 
            
            part_gene_id_df$type <- part_gene_type
            
            #part_gene_id$attributes <- sapply(part_gene_id$attributes,function(x) 
            
            # gsub("ID=","gID=",x)
            
            #)
            
            part_gene_id_df$start <- min(part_gene_Mod$start)
            
            part_gene_id_df$end <- max(part_gene_Mod$end)
            
            part_gene_id_df$Final_Mod_Attr <- paste0(unique(part_gene_Mod$Final_Mod_Attr),collapse = ";")
            
            part_gene_id_df$attributes <- 
              sapply(part_gene_id_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                     paste0("ID=XBXT10g",
                            paste0(rep(0,6-nchar(gene_count)),collapse = ""),
                            gene_count),x))
            
            part_gene_id_df$attributes <- sapply(
              
              part_gene_id_df$attributes,function(x)
                paste0(x,";geneid=",paste0(unique(part_gene_Mod$Final_Mod_Attr),
                                           collapse = ";"))  
              
            )
            
            if(nrow(part_mRNA_id_df) >= 1){
              
              #part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
              
              # gsub("ID=","transID=",x)
              
              #)
              
              #part_mRNA_id$attributes <- sapply(part_mRNA_id$attributes,function(x) 
              
              # gsub("Parent=","Origin=",x)
              
              #)
              
              mrna_count <<- mrna_count + 1
              
              part_mRNA_id_type <- part_mRNA_id_df$type
              
              part_mRNA_id_df$start <- min(part_gene_Mod$start)
              
              part_mRNA_id_df$end <- max(part_gene_Mod$end)
              
              #part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
              
              part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                     paste0("ID=XB-mRNA",
                            paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
                            mrna_count),x))  
              
              part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
                gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                     paste0("Parent=XBXT10g",
                            paste0(rep(0,6-nchar(gene_count)),collapse = ""),
                            gene_count),x))
              
              part_gene_id_df$attributes <- sapply(
                
                part_gene_id_df$attributes,function(x)
                  paste0(x,";geneid=",paste0(unique(part_gene_Mod$Final_Mod_Attr),
                                             collapse = ";"))  
                
              )
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                       paste0("Parent=",parent_id),x))
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x))
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id$Final_Mod_Attr <- paste0(unique(part_gene_Mod$Final_Mod_Attr),collapse = ";")
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(unique(part_gene_Mod$Final_Mod_Attr),
                                               collapse = ";"))  
                  
                )
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                       strand,phase,Final_Mod_Attr,attributes) %>% 
                distinct(.keep_all = TRUE)
              
              
            }
            
            else if(nrow(part_mRNA_id_df) != 1){
              
              #mrna_count <- mrna_count + 1
              
              #part_mRNA_id_df <- part_gene_id_df
              
              #part_mRNA_id_df$type <- "transcript"
              
              #part_mRNA_id_df$start <- min(part_gene_Mod$start)
              #part_mRNA_id_df$end <- max(part_gene_Mod$end)
              #part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene_Mod$modified_attributes),collapse = ";")
              
              #part_mRNA_id_df$attributes <- 
              #paste0("ID=XB-mRNA",
              #paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count,";","Parent=XBXT10g",
              #paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count,";")
              
              parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1]) 
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                       paste0("Parent=",parent_id),x))
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x))
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
                Non_gene_id$Final_Mod_Attr <- paste0(unique(part_gene_Mod$Final_Mod_Attr),collapse = ";")
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(unique(part_gene_Mod$Final_Mod_Attr),
                                               collapse = ";"))  
                  
                )
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              part_gene_M <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
                select(chrloc,source,type,start,end,score,                                                                           
                       strand,phase,Final_Mod_Attr,attributes) %>% 
                distinct(.keep_all = TRUE)
              
            }
            
            Unified_ls[[i]] <- part_gene_M %>%
              select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% 
              distinct(.keep_all = TRUE)
            
          }
          
        }
        
        filtered_strand_ls[[s]] <- do.call("rbind.data.frame",Unified_ls) 
        
      }
      
      else if(((nrow(part_gene[part_gene$type == "ncRNA_gene",]) > 1) &
               (length(unique(part_gene[part_gene$type == "ncRNA_gene",]$type) %in% 
                       c("ncRNA_gene")) == 1 )) | (nrow(part_gene[part_gene$type == "gene" | 
                                                                  part_gene$type == "pseudogene",]) >= 1 & 
                                                   nrow(part_gene[part_gene$type == "ncRNA_gene",]) >= 1) |
              (nrow(part_gene[part_gene$type == "tRNA",]) >= 1)) {
        
        #other_ls[[length(other_ls)+ 1]] <- part_gene
        
        #part_df <- part_gene %>% left_join(rank_data,by = c("type"))
        
        if(length(unique(part_gene[part_gene$type == "tRNA",]$type)) >= 1 | 
           ((length(unique(part_gene[part_gene$type == "tRNA",]$type)) >= 1) & 
            (length(unique(part_gene[part_gene$type == "ncRNA_gene" | 
                                     part_gene$type == "miRNA",]$type)) >= 1))){
          
          part_gene$ID <- sapply(part_gene$attributes,function(x) gsub("ID=","",
                                                                       str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
          
          part_gene$Parent <- sapply(part_gene$attributes,function(x) gsub("Parent=","",
                                                                           str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))
          
          filtered_rna_gene_ids <- part_gene %>% filter(type == "mRNA" |
                                                          type == "tRNA"| type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" |
                                                          type == "snRNA" | type == "scRNA" | type == "miRNA" | type == "Y_RNA" |
                                                          type == "ncRNA" | type == "guide_RNA" | type == "transcript" | 
                                                          type == "primary_transcript" |type == "pseudogenic_transcript") %>% select(ID,Parent)
          
          part_d <- part_gene %>% left_join(filtered_rna_gene_ids,by = c("Parent" = "ID"))
          
          part_d <- part_d %>% mutate(
            
            transcript = ifelse(is.na(Parent.y),Parent,Parent.y),
            final_transcript = ifelse(is.na(transcript),ID,transcript)
            
            
          ) %>% arrange(final_transcript) %>% select(chrloc,source,type,start,end,
                                                     score,strand,phase,modified_attributes,attributes,Final_Mod_Attr,final_transcript,Overlap_Value)
          
          filtered_mod_attr_ls <- list()
          final_mod_attr_ls <- list()
          
          for (l in 1 : length(unique(part_d$final_transcript))) {
            
            filtered_mod_attr_ls[[l]] <- part_d %>% filter(final_transcript %in% unique(part_d$final_transcript)[l]) 
            
            gene_filtered_mod_attr <- filtered_mod_attr_ls[[l]] %>% filter(type == "gene" | 
                                                                             type == "pseudogene" | type == "ncRNA_gene")
            
            transcript_filtered_mod_attr <- filtered_mod_attr_ls[[l]] %>% filter(type == "mRNA" |
                                                                                   type == "tRNA"| type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" |
                                                                                   type == "snRNA" | type == "scRNA" | type == "miRNA" | type == "Y_RNA" |
                                                                                   type == "ncRNA" | type == "guide_RNA" | type == "transcript" | 
                                                                                   type == "primary_transcript" |type == "pseudogenic_transcript")
            
            gene_filtered_mod_attr_df <- gene_filtered_mod_attr[which.max(gene_filtered_mod_attr$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            )
            
            transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr[which.max(transcript_filtered_mod_attr$Overlap_Value),] %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            )
            
            Non_gene_id <- filtered_mod_attr_ls[[l]] %>% filter(!(type %in% gene_filtered_mod_attr$type | type %in% transcript_filtered_mod_attr$type)) %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            )
            
            Non_gene_f_c <- Non_gene_id %>% 
              group_by(type) %>% summarise(count = n())
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              exon_i_count <- exon_count + 1 
              exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              cds_i_count <- cds_count + 1   
              cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              match_i_count <- match_count + 1
              match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              V_gene_segment_i_count <- V_gene_segment_count + 1
              V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              C_gene_segment_i_count <- C_gene_segment_count + 1
              C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              or_i_count <- or_count + 1
              or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              D_loop_i_count <- D_loop_count + 1
              D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
              
            }
            
            if(nrow(gene_filtered_mod_attr_df) == 1){
              
              gene_count <<- gene_count + 1
              
              gene_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
              gene_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
              
              gene_filtered_mod_attr_df$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
              
              #gene_filtered_mod_attr$attributes <- sapply(gene_filtered_mod_attr$attributes, function(x)
              
              # gsub("ID=","gID=",x)
              
              #)
              
              gene_filtered_mod_attr_df$attributes <- sapply(
                
                gene_filtered_mod_attr_df$attributes,function(x)
                  gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                       paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),
                       x))
              
              gene_filtered_mod_attr_df$attributes <- sapply(
                
                gene_filtered_mod_attr_df$attributes,function(x)
                  paste0(x,";geneid=",paste0(filtered_mod_attr_ls[[l]]$Final_Mod_Attr,
                                             collapse = ";"))  
                
              )
              
              if(nrow(transcript_filtered_mod_attr_df) == 1 ){
                
                mrna_count <<- mrna_count + 1
                
                transcript_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                transcript_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                
                transcript_filtered_mod_attr_df$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                
                transcript_filtered_mod_attr_df$attributes <-  sapply(
                  
                  transcript_filtered_mod_attr_df$attributes,function(x)
                    gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                         paste0("ID=XB-RNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),
                         x))
                
                transcript_filtered_mod_attr_df$attributes <-  sapply(
                  
                  transcript_filtered_mod_attr_df$attributes,function(x)
                    gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],
                         paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),
                         x))
                
                transcript_filtered_mod_attr_df$attributes <- sapply(
                  
                  transcript_filtered_mod_attr_df$attributes,function(x)
                    paste0(x,";geneid=",paste0(filtered_mod_attr_ls[[l]]$Final_Mod_Attr,
                                               collapse = ";"))  
                  
                )
                
                parent_id <- gsub("ID=","",str_split(str_extract(
                  transcript_filtered_mod_attr_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])   
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                           gsub("Parent=(.)*\\|","",x),x)  
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "\\|Parent="),
                           gsub("\\|Parent=",";Parent=",x),x)
                    
                  )
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                    
                    
                    ifelse(str_detect(x,pattern = "ID=.+;"),
                           gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                          pattern = ";")[[1]][1],"",x),paste0(";",x)    
                           
                    )
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                  
                  Non_gene_id$attributes <- sapply(
                    
                    Non_gene_id$attributes,function(x)
                      paste0(x,";geneid=",paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),
                                                 collapse = ";"))  
                    
                  )
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                gene_filtered_mod_attr_df <- gene_filtered_mod_attr_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                final_mod_attr_ls[[l]] <- rbind.data.frame(gene_filtered_mod_attr_df,transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                  
                )
                
              }
              
              else if(nrow(transcript_filtered_mod_attr_df) < 1){
                
                parent_id <- gsub("ID=","",str_split(str_extract(gene_filtered_mod_attr_df$attributes,pattern = "ID=.+;"))[[1]][1])   
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                           gsub("Parent=(.)*\\|","",x),x)  
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "\\|Parent="),
                           gsub("\\|Parent=",";Parent=",x),x)
                    
                  )
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                    
                    
                    ifelse(str_detect(x,pattern = "ID=.+;"),
                           gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                          pattern = ";")[[1]][1],"",x),paste0(";",x)   
                           
                    )
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                  
                  Non_gene_id$attributes <- sapply(
                    
                    Non_gene_id$attributes,function(x)
                      paste0(x,";geneid=",paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),
                                                 collapse = ";"))  
                    
                  )
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                gene_filtered_mod_attr_df <- gene_filtered_mod_attr_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                final_mod_attr_ls[[l]] <- rbind.data.frame(gene_filtered_mod_attr_df,transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                  
                )
                
              }
              
            }
            
            if(nrow(gene_filtered_mod_attr_df) < 1){
              
              if( nrow(transcript_filtered_mod_attr_df) == 1 ){
                
                mrna_count <<- mrna_count + 1
                
                transcript_filtered_mod_attr_df$start <- min(filtered_mod_attr_ls[[l]]$start)
                transcript_filtered_mod_attr_df$end <- max(filtered_mod_attr_ls[[l]]$end)
                
                transcript_filtered_mod_attr_df$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                
                #transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                
                # gsub("ID=","transID=",x)
                
                #)
                
                #transcript_filtered_mod_attr$attributes <- sapply(transcript_filtered_mod_attr$attributes, function(x)
                
                # gsub("Parent=","Origin=",x)
                
                #)
                
                transcript_filtered_mod_attr_df$attributes <- 
                  
                  sapply(transcript_filtered_mod_attr_df$attributes, function(x)
                    
                    gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x)    
                    
                  )
                
                transcript_filtered_mod_attr_df$attributes <- sapply(
                  
                  transcript_filtered_mod_attr_df$attributes,function(x)
                    paste0(x,";geneid=",paste0(filtered_mod_attr_ls[[l]]$Final_Mod_Attr,
                                               collapse = ";"))  
                  
                )
                
                parent_id <- gsub("ID=","",str_split(str_extract(
                  transcript_filtered_mod_attr_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])   
                
                if(nrow(Non_gene_id) >= 1){
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                           gsub("Parent=(.)*\\|","",x),x)  
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                    
                    ifelse(str_detect(x,pattern = "\\|Parent="),
                           gsub("\\|Parent=",";Parent=",x),x)
                    
                  )
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                    
                    
                    ifelse(str_detect(x,pattern = "ID=.+;"),
                           gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                          pattern = ";")[[1]][1],"",x),paste0(";",x)    
                           
                    )
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$attributes <- sapply(
                    
                    Non_gene_id$attributes,function(x)
                      paste0(x,";geneid=",paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),
                                                 collapse = ";"))  
                    
                  )
                  
                  Non_gene_id$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                transcript_filtered_mod_attr_df <- transcript_filtered_mod_attr_df  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                final_mod_attr_ls[[l]] <- rbind.data.frame(transcript_filtered_mod_attr_df,Non_gene_id) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                  
                )
                
              }
              
              else if(nrow(transcript_filtered_mod_attr_df) < 1){
                
                #parent_id <- gsub("ID=","",str_split(str_extract(gene_filtered_mod_attr_df$attributes,pattern = "ID=.+;"))[[1]][1])   
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  gsub(paste0(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],";"),"",x))
                
                if(nrow(Non_gene_id) >= 1){
                  
                  #Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  #gsub("Parent=",paste0("Parent=",parent_id,";Origin="),x))
                  
                  Non_gene_id$type_count <- NA
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                    
                  }
                  
                  if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                    
                    Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                    
                  }
                  
                  Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                    
                    ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                    
                  )
                  
                  Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                  
                  Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                    
                    paste0("ID=XB-",x)
                    
                  )
                  
                  Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                    
                    
                    ifelse(str_detect(x,pattern = "ID=.+;"),
                           gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                          pattern = ";")[[1]][1],"",x),paste0(";",x)    
                           
                    )
                    
                  )
                  
                  Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                  
                  Non_gene_id$Final_Mod_Attr <- paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),collapse = ";")
                  
                  #Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)
                  
                  Non_gene_id$attributes <- sapply(
                    
                    Non_gene_id$attributes,function(x)
                      paste0(x,";geneid=",paste0(unique(filtered_mod_attr_ls[[l]]$Final_Mod_Attr),
                                                 collapse = ";"))  
                    
                  )
                  
                }
                
                Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
                
                final_mod_attr_ls[[l]] <- rbind.data.frame(Non_gene_id) %>% select(
                  
                  chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                  
                )
                
              }
              
            }  
            
          }
          
          final_mod_attr_df <- do.call("rbind.data.frame",final_mod_attr_ls)
          
          filtered_strand_ls[[s]] <- final_mod_attr_df %>% select(chrloc,source,type,start,end,score,
                                                                  strand,phase,Final_Mod_Attr,attributes) %>% 
            distinct(.keep_all = TRUE)
          
        }
        
        if((length(unique(part_gene[part_gene$type == "ncRNA_gene" | 
                                    part_gene$type == "miRNA",]$type)) >= 1 & nrow(part_gene[part_gene$type == 'tRNA',]) == 0) | 
           (nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene",]) >= 1 & 
            nrow(part_gene[part_gene$type == "ncRNA_gene",]) > 1)){
          
          #part_df <- part_df %>% arrange(rank)
          
          final_gb_ls <- list()
          
          if (nrow(part_gene[part_gene$source == "Genbank",]) >= 1){
            
            filtered_genbank_ls <- part_gene %>% filter(source == "Genbank")
            
            filtered_gb_ls <- list()
            
            for (h in 1 : length(unique(filtered_genbank_ls$modified_attributes))) {
              
              
              filtered_gb_ls[[h]] <- filtered_genbank_ls %>% filter(Final_Mod_Attr == unique(filtered_genbank_ls$modified_attributes)[h]) %>% select(
                
                chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
                
              )
              
              
              if(nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]) >= 1){
                
                gene_count <<- gene_count + 1
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes,function(x) 
                  gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XBXT10g",
                                                                                                 paste0(rep(0,6-nchar(gene_count)),
                                                                                                        collapse = ""),gene_count),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes,function(x) 
                  paste0(x,";geneid=",paste0(filtered_gb_ls[[h]]$Final_Mod_Attr),collapse = ";"))
                
              }
              
              if(nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]) >= 1){
                
                mrna_count <<- mrna_count + 1
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-mRNA",
                                                                                                 paste0(rep(0,6-nchar(mrna_count)),
                                                                                                        collapse = ""),mrna_count),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;")),
                       paste0("Parent=XBXT10g",
                              paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x))
                
                filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes,function(x) 
                  paste0(x,";geneid=",paste0(filtered_gb_ls[[h]]$Final_Mod_Attr),collapse = ";"))
                
                parent_id <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))  
                
              }
              
              else if (nrow(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",]) < 1){
                
                
                parent_id <- sapply(filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))  
                
              }
              
              Non_gene_id <- filtered_gb_ls[[h]] %>% filter(!(type == "gene" | type == "mRNA"))
              
              Non_gene_f_c <- Non_gene_id %>% 
                group_by(type) %>% summarise(count = n())
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                exon_i_count <- exon_count + 1  
                exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                cds_i_count <- cds_count + 1    
                cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                match_i_count <- match_count + 1
                match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                V_gene_segment_i_count <- V_gene_segment_count + 1
                V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                C_gene_segment_i_count <- C_gene_segment_count + 1
                C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                or_i_count <- or_count + 1
                or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                D_loop_i_count <- D_loop_count + 1
                D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
                
              }
              
              if(nrow(Non_gene_id) >= 1){
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x)
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                         gsub("Parent=(.)*\\|","",x),x)  
                  
                )
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                  
                  ifelse(str_detect(x,pattern = "\\|Parent="),
                         gsub("\\|Parent=",";Parent=",x),x)
                  
                )
                
                Non_gene_id$type_count <- NA
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                  
                }
                
                if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                  
                  Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                  
                }
                
                Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                  
                  ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                  
                )
                
                Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
                
                Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                  
                  paste0("ID=XB-",x)
                  
                )
                
                Non_gene_id$attributes <- sapply(
                  
                  Non_gene_id$attributes,function(x)
                    paste0(x,";geneid=",paste0(filtered_gb_ls[[l]]$Final_Mod_Attr,
                                               collapse = ";"))  
                  
                )
                
                
                Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                  
                  
                  ifelse(str_detect(x,pattern = "ID=.+;"),
                         gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                        pattern = ";")[[1]][1],"",x),paste0(";",x)    
                         
                  )
                  
                )
                
                Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
                
              }
              
              Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
              
              filtered_gene_df <- filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "gene",] %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
              filtered_trans_df <- filtered_gb_ls[[h]][filtered_gb_ls[[h]]$type == "mRNA",] %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes)
              
              final_gb_ls[[h]] <- rbind.data.frame(filtered_gene_df,filtered_trans_df,
                                                   Non_gene_id) %>% select(
                                                     
                                                     chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes                                       
                                                     
                                                   )
              
            }
            
          }
          
          final_gb_df <- do.call("rbind.data.frame",final_gb_ls)
          
          part_gene <- part_gene %>% filter(!(source %in% "Genbank"))
          
          part_gene$transcript <- sapply(part_gene$attributes,function(x) 
            
            ifelse(str_detect(x,pattern = "Parent=.+"),
                   
                   gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1]),
                   
                   gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
            
          )
          
          part_gene <- part_gene %>% arrange(transcript)
          
          filtered_gene_ls <- list()
          filtered_M_ls <- list()
          
          for(k in 1: length(unique(part_gene$transcript))){
            
            gene_count <<- gene_count + 1
            
            filtered_gene_ls[[k]] <- part_gene %>% filter(transcript %in% unique(part_gene$transcript)[k]) %>% select(
              
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
              
            )
            
            filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | 
                                    filtered_gene_ls[[k]]$type == "pseudogene" | 
                                    filtered_gene_ls[[k]]$type == "gene",]$attributes <- 
              sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | 
                                             filtered_gene_ls[[k]]$type == "pseudogene" | 
                                             filtered_gene_ls[[k]]$type == "gene",]$attributes,function(x) 
                                               gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],
                                                    paste0("ID=XBXT10g",
                                                           paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x))
            
            
            if(nrow(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                          filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                          filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                          filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                          filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                          filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]) >= 1){                                                                                                                                              
              
              mrna_count <<- mrna_count + 1
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | filtered_gene_ls[[k]]$type == "miRNA" | 
                                      filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                      filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                      filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" |
                                      filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                      filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes,function(x) 
                                                                                                                                                                                 gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-rna",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x))
              
              filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                      filtered_gene_ls[[k]]$type == "miRNA" | 
                                      filtered_gene_ls[[k]]$type == "rRNA" | 
                                      filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                      filtered_gene_ls[[k]]$type == "snoRNA" | 
                                      filtered_gene_ls[[k]]$type == "scRNA" | 
                                      filtered_gene_ls[[k]]$type == "Y_RNA" | 
                                      filtered_gene_ls[[k]]$type == "ncRNA" | 
                                      filtered_gene_ls[[k]]$type == "guide_RNA" |
                                      filtered_gene_ls[[k]]$type == "transcript" | 
                                      filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "miRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "rRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "snoRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "scRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "Y_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "guide_RNA" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "transcript" | 
                                                                                                                                                                               filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes,function(x) 
                                                                                                                                                                                 gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=XBXT10g",
                                                                                                                                                                                                                                                                    paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x))
              
              
              parent_id <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                          filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                          filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                          filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                          filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                          filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
              
            }
            
            else if (nrow(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | 
                                                filtered_gene_ls[[k]]$type == "miRNA" | filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" | 
                                                filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]) < 1){
              
              parent_id <- sapply(filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",]$attributes, function(x) gsub("ID=","",str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
              
              
            }
            
            Non_gene_id <- filtered_gene_ls[[k]] %>% filter(!(type == "rRNA" | type == "lnc_RNA" | type == "snoRNA" | type == "snRNA" | 
                                                                type == "scRNA" | type == "miRNA" | type == "transcript" | type == "primary_transcript" |
                                                                type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "ncRNA_gene" | type == "mRNA" | type == "gene" | type == "pseudogene"
                                                              
            ))
            
            Non_gene_f_c <- Non_gene_id %>% 
              group_by(type) %>% summarise(count = n())
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              exon_i_count <- exon_count + 1  
              exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              cds_i_count <- cds_count + 1    
              cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              match_i_count <- match_count + 1
              match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              V_gene_segment_i_count <- V_gene_segment_count + 1
              V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              C_gene_segment_i_count <- C_gene_segment_count + 1
              C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              or_i_count <- or_count + 1
              or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              D_loop_i_count <- D_loop_count + 1
              D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
              
            }
            
            if(nrow(Non_gene_id) >= 1){
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x)
                
              )
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                       gsub("Parent=(.)*\\|","",x),x)  
                
              )
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
                
                ifelse(str_detect(x,pattern = "\\|Parent="),
                       gsub("\\|Parent=",";Parent=",x),x)
                
              )
              
              Non_gene_id$type_count <- NA
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
                
              }
              
              if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
                
                Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
                
              }
              
              Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
                
                ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
                
              )
              
              Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
              
              Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
                
                paste0("ID=XB-",x)
                
              )
              
              
              Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
                
                
                ifelse(str_detect(x,pattern = "ID=.+;"),
                       gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                      pattern = ";")[[1]][1],"",x),paste0(";",x)   
                       
                )
                
              )
              
              Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
              
            }
            
            Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
            
            filtered_gene_df <- filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "ncRNA_gene" | filtered_gene_ls[[k]]$type == "pseudogene" | filtered_gene_ls[[k]]$type == "gene",] %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
            
            filtered_trans_df <- filtered_gene_ls[[k]][filtered_gene_ls[[k]]$type == "snRNA" | filtered_gene_ls[[k]]$type == "miRNA" | 
                                                         filtered_gene_ls[[k]]$type == "rRNA" | filtered_gene_ls[[k]]$type == "lnc_RNA" | 
                                                         filtered_gene_ls[[k]]$type == "snoRNA" | filtered_gene_ls[[k]]$type == "scRNA" | 
                                                         filtered_gene_ls[[k]]$type == "Y_RNA" | filtered_gene_ls[[k]]$type == "ncRNA" |
                                                         filtered_gene_ls[[k]]$type == "guide_RNA" | filtered_gene_ls[[k]]$type == "transcript" | 
                                                         filtered_gene_ls[[k]]$type == "primary_transcript" | filtered_gene_ls[[k]]$type == "mRNA",]  %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
            
            filtered_M_ls[[k]] <- rbind.data.frame(filtered_gene_df,filtered_trans_df,Non_gene_id) %>% select(
              chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes
            )
            
          }
          
          Final_mod_d <- do.call("rbind.data.frame",filtered_M_ls)
          
          if(nrow(final_gb_df) >= 1){
            
            Final_mod_d <- rbind.data.frame(final_gb_df,Final_mod_d)
            
          }
          
          filtered_strand_ls[[s]] <- Final_mod_d %>% select(chrloc,source,type,start,end,score,
                                                            strand,phase,Final_Mod_Attr,attributes) %>% 
            distinct(.keep_all = TRUE)
          
        }
        
      } 
      
      
      else if(((nrow(part_gene[part_gene$type == "gene",]) >= 2) &
               (length(unique(part_gene[part_gene$type == "gene",]$type) %in% 
                       c("gene")) == 1 ) &
               nrow(part_gene[part_gene$type == "tRNA",]) == 0) |
              ((nrow(part_gene[part_gene$type == "pseudogene",]) >= 2) & 
               (length(unique(part_gene[part_gene$type == "pseudogene",]$type) %in%  
                       c("pseudogene")) == 1))){ 
        
        #else if(nrow(part_gene_id) > 1 & nrow(part_mRNA_id) >= 0){
        part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
        
        part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
                                               type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                               type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                               type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "primary_transcript" |
                                               type == "transcript" | type == "pseudogenic_transcript" ) 
        
        Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
          select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
          distinct(.keep_all = TRUE)
        
        Non_gene_f_c <- Non_gene_id %>% 
          group_by(type) %>% summarise(count = n())
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
          
          exon_i_count <- exon_count + 1  
          exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
          
          cds_i_count <- cds_count + 1   
          cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
          
          match_i_count <- match_count + 1
          match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
          
          V_gene_segment_i_count <- V_gene_segment_count + 1
          V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
          
          C_gene_segment_i_count <- C_gene_segment_count + 1
          C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
          
          or_i_count <- or_count + 1
          or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
          
          D_loop_i_count <- D_loop_count + 1
          D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
        }
        
        if(nrow(part_gene_id) >= 1){
          
          gene_count <<- gene_count + 1
          
          part_gene_id_df <- part_gene_id[which.max(part_gene_id$Overlap_Value),] %>% 
            select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
            distinct(.keep_all = TRUE)
          
          part_gene_id_df$start <- min(part_gene$start)
          
          part_gene_id_df$end <- max(part_gene$end)
          
          part_gene_id_df$attributes <- sapply(part_gene_id_df$attributes,function(x) 
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x)
            
          )
          
          part_gene_id_df$attributes <- sapply(part_gene_id_df$attributes,function(x) 
            
            paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")) 
            
          )
          
          #part_gene_id_df$modified_attributes <- paste0(unique(part_gene$Final_mod_attr),collapse = ";")
          
        }
        
        if(nrow(part_mRNA_id) >= 1){
          
          mrna_count <<- mrna_count + 1
          
          part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
            select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
            distinct(.keep_all = TRUE)
          
        }
        
        #### finding the start and end positions for customary gene/transcript model
        
        if(nrow(part_mRNA_id_df) == 1){
          
          part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
            
            gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x)
            
          )
          
          part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
            
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=XBXT10g",paste0(rep(0,6-nchar(gene_count)),collapse = ""),gene_count),x)
            
          )
          
          part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
            
            paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")) 
            
          )
          
          part_mRNA_id_df$start <- min(part_gene_id_df$start)
          
          part_mRNA_id_df$end <- max(part_gene_id_df$end)
          
          #part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$Final_Mod_attr),collapse = ";")    
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
          if(nrow(Non_gene_id) >= 1){
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
              gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                     gsub("Parent=(.)*\\|","",x),x)  
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "\\|Parent="),
                     gsub("\\|Parent=",";Parent=",x),x)
              
            )
            
            Non_gene_id$type_count <- NA
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
              
            }
            
            Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
              
              ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
              
            )
            
            Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
            
            Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
              
              paste0("ID=XB-",x)
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(!str_detect(x,pattern = "ID=(.)+\\|(.)+;") & 
                       str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),ifelse(
                                      str_detect(x,pattern = "ID=(.)+\\|(.)+;"),
                                      gsub("\\|",";",x),paste0(";",x)))
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
              
              
              ifelse(str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),x    
                     
              )
              
            )
            
            Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
            
            Non_gene_id$attributes <- sapply(
              
              Non_gene_id$attributes,function(x)
                paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),
                                           collapse = ";"))  
              
            )
            
          }
          
          Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) %>% 
            select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% 
            distinct(.keep_all = TRUE) 
          
        }
        
        #}  
        
        else if(nrow(part_mRNA_id_df) < 1){
          
          #mrna_count <- mrna_count + 1
          
          #part_mRNA_id_df <- part_gene_id_df 
          
          #part_mRNA_id_df$type <- "transcript"
          
          #part_mRNA_id_df$attributes <- paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
          # mrna_count,";","Parent=XBXT10g",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),
          # mrna_count,";")
          
          #part_mRNA_id_df$start <- part_gene_id_df$start
          
          #part_mRNA_id_df$end <- part_gene_id_df$end 
          
          #part_mRNA_id_df$modified_attributes <- paste0(unique(part_gene$Final_mod_attr),collapse = ";")
          
          parent_id <- gsub("ID=","",str_split(str_extract(part_gene_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
          
          if(nrow(Non_gene_id) >= 1){
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
              gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                     gsub("Parent=(.)*\\|","",x),x)  
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
              
              ifelse(str_detect(x,pattern = "\\|Parent="),
                     gsub("\\|Parent=",";Parent=",x),x)
              
            )
            
            Non_gene_id$type_count <- NA
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
              
            }
            
            if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
              
              Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
              
            }
            
            Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
              
              ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
              
            )
            
            Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
            
            Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
              
              paste0("ID=XB-",x)
              
            )
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
              
              
              ifelse(str_detect(x,pattern = "ID=.+;"),
                     gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                    pattern = ";")[[1]][1],"",x),paste0(";",x)    
                     
              )
              
            )
            
            Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
            
            Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x) 
              
              paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")) 
              
            )
            
            #Non_gene_id$modified_attributes <- paste0(unique(part_gene$modified_attributes),collapse = ";")    
            
          }
          
          Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
          
          filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id)
          
        }
        
      }
      
      else if((nrow(part_gene[part_gene$type == "gene" | part_gene$type == "pseudogene" | 
                              part_gene$type == "ncRNA_gene",]) == 0) & 
              (nrow(part_gene[part_gene$type == "mRNA" | 
                              part_gene$type == "tRNA" | 
                              part_gene$type == "rRNA" | part_gene$type == "lnc_RNA" | 
                              part_gene$type == "snoRNA" | part_gene$type == "snRNA" | 
                              part_gene$type == "scRNA" | part_gene$type == "miRNA" |  
                              part_gene$type == "Y_RNA" | part_gene$type == "ncRNA" | 
                              part_gene$type == "guide_RNA" | 
                              part_gene$type == "primary_transcript" |
                              part_gene$type == "transcript" | 
                              part_gene$type == "pseudogenic_transcript",]) > 1)){
        
        #if(nrow(part_mRNA_id) >= 1){
        
        part_gene_id <- part_gene %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene")
        
        part_mRNA_id <- part_gene %>% filter(type == "mRNA" | type == "tRNA" | 
                                               type == "rRNA" |type == "lnc_RNA" | type == "snoRNA" | 
                                               type == "snRNA" | type == "scRNA" | type == "miRNA" |  
                                               type == "Y_RNA" | type == "ncRNA" | type == "guide_RNA" | type == "primary_transcript" |
                                               type == "transcript" | type == "pseudogenic_transcript" ) 
        
        Non_gene_id <-  part_gene %>% filter(!(type %in% part_gene_id$type | type %in% part_mRNA_id$type)) %>% 
          select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
          distinct(.keep_all = TRUE)
        
        Non_gene_f_c <- Non_gene_id %>% 
          group_by(type) %>% summarise(count = n())
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
          
          exon_i_count <- exon_count + 1  
          exon_count <<- exon_count + Non_gene_f_c[Non_gene_f_c$type == "exon",]$count
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
          
          cds_i_count <- cds_count + 1   
          cds_count <<- cds_count + Non_gene_f_c[Non_gene_f_c$type == "CDS",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
          
          match_i_count <- match_count + 1
          match_count <<- match_count + Non_gene_f_c[Non_gene_f_c$type == "match",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
          
          V_gene_segment_i_count <- V_gene_segment_count + 1
          V_gene_segment_count <<- V_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
          
          C_gene_segment_i_count <- C_gene_segment_count + 1
          C_gene_segment_count <<- C_gene_segment_count + Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
          
          or_i_count <- or_count + 1
          or_count <<- or_count + Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]$count   
        }
        
        if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
          
          D_loop_i_count <- D_loop_count + 1
          D_loop_count <<- D_loop_count + Non_gene_f_c[Non_gene_f_c$type == "D_loop",]$count   
        }
        
        mrna_count <<- mrna_count + 1
        
        part_mRNA_id_df <- part_mRNA_id[which.max(part_mRNA_id$Overlap_Value),] %>%
          select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>%
          distinct(.keep_all = TRUE)
        
        part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
          
          gsub(str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1],paste0("ID=XB-mRNA",paste0(rep(0,6-nchar(mrna_count)),collapse = ""),mrna_count),x)
          
        )
        
        part_mRNA_id_df$attributes <- sapply(part_mRNA_id_df$attributes,function(x) 
          
          paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")) 
          
        )
        
        part_mRNA_id_df$start <- min(part_gene$start)
        
        part_mRNA_id_df$end <- max(part_gene$end)
        
        parent_id <- gsub("ID=","",str_split(str_extract(part_mRNA_id_df$attributes,pattern = "ID=.+;"),pattern = ";")[[1]][1])
        
        if(nrow(Non_gene_id) >= 1){
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x) 
            gsub(str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1],paste0("Parent=",parent_id),x))
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            ifelse(str_detect(x,pattern = "Parent=(.)+\\|Parent=(.)+;"),
                   gsub("Parent=(.)*\\|","",x),x)  
            
          )
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes, function(x)
            
            ifelse(str_detect(x,pattern = "\\|Parent="),
                   gsub("\\|Parent=",";Parent=",x),x)
            
          )
          
          Non_gene_id$type_count <- NA
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "exon",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "exon",]$type_count <- exon_i_count : exon_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "CDS",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "CDS",]$type_count <- cds_i_count : cds_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "match",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "match",]$type_count <- match_i_count : match_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "V_gene_segment",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "V_gene_segment",]$type_count <- V_gene_segment_i_count : V_gene_segment_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "C_gene_segment",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "C_gene_segment",]$type_count <- C_gene_segment_i_count : C_gene_segment_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "origin_of_replication",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "origin_of_replication",]$type_count <- or_i_count : or_count
            
          }
          
          if(nrow(Non_gene_f_c[Non_gene_f_c$type == "D_loop",]) >= 1){
            
            Non_gene_id[Non_gene_id$type == "D_loop",]$type_count <- D_loop_i_count : D_loop_count
            
          }
          
          Non_gene_id$type_mod <- sapply(Non_gene_id$type_count, function(x)
            
            ifelse(nchar(x) <= 6,paste0(paste0(rep(0,6 - nchar(x)),collapse = ""),x),x)
            
          )
          
          Non_gene_id$type_f <- paste0(casefold(Non_gene_id$type),Non_gene_id$type_mod)
          
          Non_gene_id$type_final <- sapply(Non_gene_id$type_f, function(x)
            
            paste0("ID=XB-",x)
            
          )
          
          Non_gene_id$attributes <- sapply(Non_gene_id$attributes,function(x)
            
            
            ifelse(str_detect(x,pattern = "ID=.+;"),
                   gsub(str_split(str_extract(x,pattern = "ID=.+;"),
                                  pattern = ";")[[1]][1],"",x),paste0(";",x)    
                   
            )
            
          )
          
          Non_gene_id$attributes <- paste0(Non_gene_id$type_final,Non_gene_id$attributes)
          
          Non_gene_id$attributes <- sapply(
            
            Non_gene_id$attributes,function(x)
              paste0(x,";geneid=",paste0(unique(part_gene$Final_Mod_Attr),
                                         collapse = ";"))  
            
          )
          
          #Non_gene_id$modified_attributes <- paste0(unique(part_gene$Final_Mod_Attr),collapse = ";")    
          
        }
        
        Non_gene_id <- Non_gene_id %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_gene_id_df <- part_gene_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        part_mRNA_id_df <- part_mRNA_id_df %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Mod_Attr,attributes) %>% distinct(.keep_all = TRUE)
        
        filtered_strand_ls[[s]] <- rbind.data.frame(part_gene_id_df,part_mRNA_id_df,Non_gene_id) 
        
      }
      
    }  
    
  }
  
 }

Final_gene_ls[[u]] <- do.call("rbind.data.frame",filtered_strand_ls)

}

Final_gene_df <- do.call("rbind.data.frame",Final_gene_ls)
return(Final_gene_df)

}

filtered_Final_ls <- list()
Unified_Final_ls <- list()

for (i in 1: length(unique(Com_Filtered_d$chrloc))){
  
  filtered_Final_ls[[i]] <- Com_Filtered_d %>% filter(chrloc == unique(Com_Filtered_d$chrloc)[i])
  Unified_Final_ls[[i]] <- Sort_Unified_Data(filtered_Final_ls[[i]])
  
}

Unified_final <- do.call("rbind.data.frame",Unified_Final_ls)

######################################################

error_file <- read.table("GFF3_error.txt",header = TRUE,sep = "\t")

names(error_file) <- c("description")  

#error_file$description <- sapply(error_file$description,function(x) 
#str_trim(x,c("both")))

error_file$description <- sapply(error_file$description,function(x) 

  str_trim(str_split(str_extract(x,pattern = "(.)*\\|[[:space:]]Parent="),pattern = "\\|[[:space:]]Parent=")[[1]][1],c("both"))

)

Missing_P <- Unified_final %>% filter(ID %in% error_file$description)  

###############################################################

#Com_Filtered_data <- Com_Filtered_d %>% filter(!(is.na(chrloc) & is.na(strand) & is.na(start) & is.na(end) & is.na(type) & is.na(score) & is.na(strand) & is.na(phase) & is.na(modified_attributes) & is.na(attributes) & is.na(Final_Mod_Attr)))

######### FInal SOrting of Data ##################

#UD_s <- str_sort(unique(Unified_final$chrloc),numeric = TRUE)

#rank_d <- data.frame(chr = UD_s, rank = 1:77)

#Sorted_Unified_Data <- Unified_final %>% left_join(rank_d,by = c("chrloc" = "chr"))

Sort_Unified_F <- Unified_final %>% 
select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
distinct(.keep_all = TRUE)

Sort_Unified_F$source <- "Xenbase"

write.table(Sort_Unified_F,file = "Final_Modified_NCBI_Ens_XGC_Unified_Model.gff3",
sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing into final file

############################################################

Unified_final <- read.delim("Final_Modified_NCBI_Ens_XGC_Unified_Model.gff3",sep = "\t",header = FALSE)

names(Unified_final) <- c("chrloc","source","type","start","end","score","strand","phase","attributes")

Unified_final$ID <- sapply(Unified_final$attributes,function(x) gsub("ID=","",
                                                                 
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

Unified_final$Parent <- sapply(Unified_final$attributes,function(x) gsub("Parent=","",
                                                                     
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))  

#filtered_d <- list()
#filtered_rep_ids <- list()

#for (j in 1: length(filtered_id_list)) {
  
  #filtered_d[[j]] <- 
    
   # lapply(Filtered_Gene_Data$mod_attr_split,function(x)
    
    #    Filtered_Gene_Data[filtered_id_list[j] %in% x,]
    
    #)
  
  #filtered_rep_ids[[length(filtered_rep_ids)+1]] <- lapply(filtered_d,function(x)
    
  # ifelse(nrow(x) > 1,x[!which.max(x$mod_attr_c),],NA)  
    
  #)
  
#}

#removal_rep_data <- unique(unlist(filtered_rep_ids))

#removal_rep_ids <- unique(removal_rep_data$modified_attributes)
  
#Filtered_f_data <- Filtered_df %>% filter(!(modified_attributes %in% removal_rep_ids))

######## sorting data acc to chr location #############

#Test_Chr1_Unified_File <- Filtered_df %>% filter(chrloc == "Chr1") %>% 
 # select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
  #distinct(.keep_all = TRUE)

#Test_Chr1_Unified_File$source <- "Xenbase"

#UD_s <- str_sort(unique(Filtered_f_data$chrloc),numeric = TRUE)

#rank_d <- data.frame(chr = UD_s, rank = 1:77)

#Sorted_Unified_Data <- Filtered_df %>% left_join(rank_d,by = c("chrloc" = "chr"))

#Sort_Unified_F <- Sorted_Unified_Data %>% arrange(rank) %>% 
 # select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
#  distinct(.keep_all = TRUE)

#Sort_Unified_F$source <- "Xenbase"

#write.table(Sort_Unified_F,file = "Final_NCBI_Ens_XGC_Unified_Model.gff3",
# sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing into final file

#write.table(Test_Chr1_Unified_File,file = "Test_Chr1_NCBI_Ens_XGC_Unified_File.gff3",
 #           sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file


     ################ inspecting duplicates ###############

#filtered_rep_Unified_d <- Final_Unified_d %>% filter(type == "gene" | type == "pseudogene" | type == "ncRNA_gene") %>%
 # group_by(modified_attributes) %>%
#  summarise(
    
 #   mod_attr_c = n()
    
#  )

#filtered_df <- filtered_rep_Unified_d[filtered_rep_Unified_d$mod_attr_c > 1,]

#filtered_d_2 <- filtered_df[filtered_df$mod_attr_c == 2,]

#filtered_d_2 <- filtered_df[filtered_df$mod_attr_c > 2,]

####################################################################

g_elf1 <- Unified_final %>% filter(Final_Mod_Attr == "elf1")

g_LOC <- Unified_final %>% filter(Final_Mod_Attr == "LOC100145764;ENSXETG00000033287;smyd5;LOC116406512")

g_LOC_1 <- Unified_final %>% filter(Final_Mod_Attr == "smyd5")

################ checking for models ################################

g_trnaq_cug <- Unified_final %>% filter(Final_Mod_Attr == "trnaq-cug")

g_trnaq_cug$ID <- sapply(g_trnaq_cug$attributes,function(x) gsub("ID=","",
                                                                 
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_trnaq_cug$Parent <- sapply(g_trnaq_cug$attributes,function(x) gsub("Parent=","",
                                                                     
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))  

##################################################################  

g_trnat_ugu <- Unified_final %>% filter(Final_Mod_Attr == "trnat-ugu")

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

g_mod_abhd17b <- Unified_final %>% filter(Final_Mod_Attr == "abhd17b")

g_mod_abhd17b_XENTR_v10001772 <- Unified_final %>% filter(Final_Mod_Attr == "abhd17b;XENTR_v10001772")

g_mod_abhd17b_XENTR_v10001772$ID <- sapply(g_mod_abhd17b_XENTR_v10001772$attributes,function(x) gsub("ID=","",
                                                                                                     
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_mod_abhd17b_XENTR_v10001772$Parent <- sapply(g_mod_abhd17b_XENTR_v10001772$attributes,function(x) gsub("Parent=","",
                                                                                                         
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

#################################################################

g_trnac_gca <- Unified_final %>% filter(Final_Mod_Attr == "trnac-gca")

g_trnac_gca$ID <- sapply(g_trnac_gca$attributes,function(x) gsub("ID=","",
                                                                 
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_trnac_gca$Parent <- sapply(g_trnac_gca$attributes,function(x) gsub("Parent=","",

str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

##############################################################

g_U4 <- Unified_final %>% filter(Final_Mod_Attr == "U4")

g_U4$ID <- sapply(g_U4$attributes,function(x) gsub("ID=","",
                                                   
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_U4$Parent <- sapply(g_U4$attributes,function(x) gsub("Parent=","",
                                                       
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

##################################################################

g_jade1 <- Unified_final %>% filter(Final_Mod_Attr == "jade1")

#Unified_Final_D$source <- "Xenbase"

#write.table(Unified_Final_D,file = "NCBI_Ens_XGC_Unified_Final.gff3",
#sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

#write.table(Test_Chr1_Unified_File,file = "Test_Chr1_NCBI_Ens_XGC_Unified_File.gff3",
 #           sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

######### sorting data ########

#Unified_D <- read.delim("NCBI_Ens_XGC_Unified_Final.gff3",sep = "\t",header = FALSE)

#names(Unified_D) <- c("chrloc","source","type","start","end","score","strand","phase","attributes")

#Sort_Unified_F <- Sorted_Unified_Data %>% arrange(rank) %>% select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% distinct(.keep_all = TRUE)

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

Missing_P <- Unified_final %>% filter(ID %in% error_file$description)  

#####################################################################

#Test_Chr1_Unified_File <- Sorted_Unified_Data %>% filter(chrloc == "Chr1") %>% 
 # select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
#  distinct(.keep_all = TRUE)

#Test_Chr1_Unified_File$source <- "Xenbase"

#Test_Chr1_Unified_File$ID <- sapply(Test_Chr1_Unified_File$attributes,function(x) gsub("ID=","",
#str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))
#Test_Chr1_Unified_File$Parent <- sapply(Test_Chr1_Unified_File$attributes,function(x) gsub("Parent=","",
#str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

#error_file <- read.table("GFF3_error_uniq.txt",header = TRUE,sep = "\t")

#names(error_file) <- c("description")  

#error_file$description  <- sapply(error_file$description, function(x)
  
 # gsub("Parent=","",str_split(str_extract(x,pattern = "Parent=.+"),pattern = ";")[[1]][1])
  
#)

#Missing_P <- Test_Chr1_Unified_File %>% filter(Parent %in% error_file$description)  

#Test_Unified_Gene_Models <- Unified_M_Df %>% filter(modified_attributes == "aacs") %>%
# select(chrloc,source,type,start,end,score,strand,phase,attributes) %>% 
#distinct(.keep_all = TRUE)

#Test_Unified_Gene_Models$source <- "Xenbase"

#write.table(Test_Unified_Gene_Models,file = "Test_NCBI_Ens_XGC_Unified_Models.gff3",
#sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

#### checking the custom models of some complicated genes overall and chromosome in particular

#abv1 <- Unified_final %>% filter(Final_Mod_Attr == "b3gnt4" | Final_Mod_Attr == "pycr3")

#g_pycr3 <- Unified_final %>% filter(Final_Mod_Attr == "pycr3")

#abv <- Unified_final %>% filter(Final_Mod_Attr == "b3gnt4;pycr3" | Final_Mod_Attr == "b3gnt4;b3gnt4" | Final_Mod_Attr == "pycr3;pycr3") %>% 
 # distinct(.keep_all = TRUE)

#abv$ID <- sapply(abv$attributes,function(x) gsub("ID=","",
                                                 
#str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

#abv$Parent <- sapply(abv$attributes,function(x) gsub("Parent=","",
                                                     
#str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

#abg1 <- Unified_final %>% filter(Final_Mod_Attr == "5S_rRNA") %>% distinct(.keep_all = TRUE)

######################################################################

g_aacs <- Unified_final %>% filter(Final_Mod_Attr == "aacs")

#abg3_o <- Final_Unified_df %>% filter(modified_attributes == "trnat-ugu")

g_aacs$ID <- sapply(g_aacs$attributes,function(x) gsub("ID=","",
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

g_aacs$Parent <- sapply(g_aacs$attributes,function(x) gsub("Parent=","",
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

####### upon observation, the complicated models are properly assigned according to 
### their respective strands,chromosomes and start and end.

#################### long ncRNA genes ###################

#ncRNA_D <- Final_Unified_d %>% filter(modified_attributes == "U4" | modified_attributes == "xtr-mir-101a-1" | modified_attributes == "xtr-mir-146b" | modified_attributes == "smyd5" | modified_attributes == "RNAaseP_nuc" | modified_attributes == "LOC116406437" | modified_attributes == "LOC101734477")

#ncRNA_D_o <- Final_grped_df %>% filter(modified_attributes == "U4" | modified_attributes == "xtr-mir-101a-1" | modified_attributes == "xtr-mir-146b" | modified_attributes == "smyd5" | modified_attributes == "RNAaseP_nuc" | modified_attributes == "LOC116406437" | modified_attributes == "LOC101734477")

g_hra2 <- Unified_final %>% filter(Final_Mod_Attr == "htra2;ENSXETG00000049421")

g_bbc3_hra3 <- Unified_final %>% filter(Final_Mod_Attr == "bbc3;htra2;ENSXETG00000049421")

g_entire <- Unified_final %>% filter(Final_Mod_Attr == "LOC101732307;dok1;ENSXETG00000021845;ENSXETG00000048397;ENSXETG00000046822")

g_entire_d <- Unified_final %>% filter(Final_Mod_Attr == "bbc3;htra2;ENSXETG00000049421;htra2;ENSXETG00000049421")

#g_bbc3_hra3$ID <- sapply(g_bbc3_hra3$attributes,function(x) gsub("ID=","",
#str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

#g_bbc3_hra3$Parent <- sapply(g_bbc3_hra3$attributes,function(x) gsub("Parent=","",
#str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

#g_hra2$ID <- sapply(g_hra2$attributes,function(x) gsub("ID=","",
#str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

#g_hra2$Parent <- sapply(g_hra2$attributes,function(x) gsub("Parent=","",
#str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

####################################################################

abg <- Unified_final[1:500,]
g_chr1 <- Unified_final[Unified_final$chrloc == "Chr1",][1:500,]
g_chr2 <- Unified_final[Unified_final$chrloc == "Chr2",][1:500,]

part_gene$ID <- sapply(part_gene$attributes,function(x) gsub("ID=","",
str_split(str_extract(x,pattern = "ID=.+;"),pattern = ";")[[1]][1]))

part_gene$Parent <- sapply(part_gene$attributes,function(x) gsub("Parent=","",
str_split(str_extract(x,pattern = "Parent=.+;"),pattern = ";")[[1]][1]))

##############################################################

###### Nodal5 genes ########

g_nodal5_ENSXETG00000017442_entire <- Unified_final %>% filter(Final_Mod_Attr == "nodal5;nodal5.2;LOC100490072;nodal5.2;LOC101734231;LOC116409495")

unique(g_nodal5_ENSXETG00000017442_entire$modified_attributes)

g_nodal_3 <- Com_Filtered_d %>% filter(start == 130796558 & end == 130822828)
g_nodal_4 <- Com_Filtered_d %>% filter(start == 130796558 & end == 130822795)
g_nodal_5 <- Com_Filtered_d %>% filter(start == 130796244 & end == 130822795)

g_nodal5 <- Unified_final %>% filter(Final_Mod_Attr == "nodal5")

####################################################

Unified_final <- Unified_final %>% mutate(
  
  mod_attr_split = str_split(Final_Mod_Attr,pattern = ";")
  
)

Unified_final$mod_attr_c <- lapply(Unified_final$mod_attr_split, function(x)
  
  length(unique(x))
  
)

filtered_gene_mod_attr_count <- Unified_final[Unified_final$type == "gene" & Unified_final$mod_attr_c > 1,] %>% select(
  
  chrloc,start,end,Final_Mod_Attr,mod_attr_c
  
) %>% distinct(.keep_all = TRUE) %>% summarise(count = n())

filtered_gene_mod_attr_single_count <- Unified_final[Unified_final$type == "gene" & Unified_final$mod_attr_c == 1,] %>% select(
  
  chrloc,start,end,Final_Mod_Attr,mod_attr_c
  
) %>% distinct(.keep_all = TRUE)

filtered_gene_mod_attr_count$mod_attr_c <- unlist(filtered_gene_mod_attr_count$mod_attr_c)

filtered_ncrnagene_mod_attr_count <- Unified_final[Unified_final$type == "ncRNA_gene",] %>% select(
  
  chrloc,start,end,Final_Mod_Attr,mod_attr_c
  
) %>% distinct(.keep_all = TRUE) %>% summarise(count = n())

filtered_pseudogene_mod_attr_count <- Unified_final[Unified_final$type == "pseudogene",] %>% select(
  
  chrloc,start,end,Final_Mod_Attr,mod_attr_c
  
) %>% distinct(.keep_all = TRUE) %>% summarise(count = n())

#filtered_gene_LOC_mod_attr_count <- Unified_final[Unified_final$type == "gene" & str_detect(Unified_final$Final_Mod_Attr,pattern = "LOC") & Unified_final$mod_attr_c > 1,] %>% select(
  
 # chrloc,start,end,Final_Mod_Attr,mod_attr_c
  
#) %>% distinct(.keep_all = TRUE) %>% summarise(count = n())

write.csv(filtered_gene_mod_attr_count,file = "Complictaed_gene_count.csv",quote = FALSE,row.names = FALSE)

model_count <- data.frame(type = c("single","complicated","ncRNA_gene","pseudogene"),count = c(filtered_gene_mod_attr_single_count$count, filtered_gene_mod_attr_count$count,filtered_ncrnagene_mod_attr_count$count,filtered_pseudogene_mod_attr_count$count))

model_count$type <- as.factor(model_count$type)
model_count$count <- as.factor(model_count$count)

plotly::plot_ly(data = model_count,x = ~type,y = ~count,type = "bar",color = ~type,colors = c("seagreen","maroon")) %>% layout(
  
  yaxis = list(showgrid = FALSE,tickformat = "digits")

) 