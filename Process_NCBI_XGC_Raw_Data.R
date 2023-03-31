
#########################################################################################

########## Processing NCBI and XGC Datasets to be infed into fjoin algorithm ############

#########################################################################################

##loading packages

library(dplyr)
library(stringr)
library(tidyr)

########################################################################################

## loading in raw datasets ##

########################################################################################

##Reading in NCBI Data from datasets

NCBI_Data <- read.csv('genomic_NCBI.gff',header = FALSE,sep = "\t",comment.char = "#")

names(NCBI_Data) <- c('chrloc','source','type','start','end','score','strand','phase','attributes')


##Reading in XGC Data from datasets

XGC_Data <- read.csv('genomic_XGC.gff',header = FALSE,sep = "\t",comment.char = "#")

names(XGC_Data) <- c('chrloc','source','type','start','end','score','strand','phase','attributes')

#####################################################################################

########## Analyzing the datasets to extract gene ids for NCBI and XGC Data #########

#####################################################################################

## removing region type from NCBI and XGC datasets ##

NCBI_Data <- NCBI_Data %>% filter(type != "region" & type != "cDNA_match")
XGC_Data <- XGC_Data %>% filter(type != "region")


## changing Accession IDS to chromosome locations of column 1 through assembly report data

assem_file <- read.csv('https://ftp.ncbi.nlm.nih.gov/genomes/genbank/vertebrate_other/Xenopus_tropicalis/latest_assembly_versions/GCA_000004195.4_UCB_Xtro_10.0/GCA_000004195.4_UCB_Xtro_10.0_assembly_report.txt')

names(assem_file) <- c('Description')

assem_file <- data.frame(Desc = assem_file$Description[35:201])

assem_data <- assem_file %>% separate(Desc,into=c('Sequence_Name','Sequence-Role','Assigned-Molecule','Assigned-Molecule-Location/Type','GenBank_Accn','Relationship','RefSeq_Accn','Assembly-Unit'	,'Sequence-Length'	,'UCSC-style-name'),sep = '\t')

NCBI_assem_df <- assem_data[,c("Sequence_Name","RefSeq_Accn")]
UCB_assem_df <-  assem_data[,c("Sequence_Name","GenBank_Accn")]

NCBI_Data <- NCBI_Data %>% left_join(NCBI_assem_df,by = c('chrloc' = 'RefSeq_Accn')) %>% select('Sequence_Name','source','type','start','end','score','strand','phase','attributes') 

XGC_Data <- XGC_Data %>% left_join(UCB_assem_df,by = c('chrloc' = 'GenBank_Accn')) %>% select('Sequence_Name','source','type','start','end','score','strand','phase','attributes') %>% distinct(.keep_all = TRUE) 

### Analyzing NCBI_Data merged with chromosome locations ##

unique(NCBI_Data$Sequence_Name)
unique(XGC_Data$Sequence_Name)

## renaming source name for NCBI Data ##

NCBI_Data$source <- "NCBI" 

## adding modified attributes to NCBI Data to be prepared for fjoin phase ##

NCBI_Data$modified_attributes <- sapply(NCBI_Data$attributes, function(x) gsub("gene=","",str_split(str_extract(x,pattern = "gene=.+;"),pattern = ";")[[1]][1]))

### checking for some abberant gene ids hidden in attributes and adding it to gene id info

NCBI_M_Data <- NCBI_Data %>% mutate(
  
  Mod_Attr = ifelse(is.na(modified_attributes),
  ifelse(str_detect(attributes,pattern = "gene=.+"),
  gsub("gene=","",str_extract(attributes,pattern = "gene=.+")),modified_attributes),
  modified_attributes)
  
  
)

### extracting ids for those features with non-existent gene id information ###

Processed_NCBI_Data <- NCBI_M_Data %>% mutate(
  

  Modified_Attr = ifelse(is.na(Mod_Attr),
   strsplit(str_extract(attributes,pattern = "ID=.+;"),split = ";"),Mod_Attr                      
  )
                        
) %>% select(Sequence_Name,source,type,start,end,score,strand,phase,Modified_Attr,attributes) %>% distinct(.keep_all = TRUE)

### preparing NCBI Datas ids info for non-existent gene IDs models

Processed_NCBI_Data$Modified_Attr <- sapply(Processed_NCBI_Data$Modified_Attr, function(x) gsub("ID=","",x[1]))
  
##checking if all gene ids have been extracted and any missing gene ids

Processed_NCBI_Data %>% filter(Modified_Attr == "rpl34p")

Processed_NCBI_Data %>% filter(is.na(Modified_Attr) | Modified_Attr == "NA")

## renaming Processed NCBI Data and XGC Data to be further processed and to extract gene id information

names(Processed_NCBI_Data) <- c('chrloc','source','type','start','end','score','strand','phase','modified_attributes','attributes')

names(XGC_Data) <- c('chrloc','source','type','start','end','score','strand','phase','attributes')

##########################################################################################

### Outcome: gene details were extracted to the attributes and added as column for NCBI Data.

#########################################################################################

## writting final processed NCBI Data as gff files ###

write.table(Processed_NCBI_Data,file = "NCBI.gff3",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)

#########################################################################################

## merging gene id info to XGC Data from NCBI gene id data and adding it to XGC Data ###

## Merge Operation:

# Initially, NCBI and XGC Data are overlapped along with gene id information from NCBI 
# Data through foverlap.
# After the above overlap, XGC Data is further merged with the overlap data containing
# gene id information.

#########################################################################################

## Preparation of NCBI and XGC Data

# preparing NCBI and XGC data for foverlap.Gene id info is contained in the gene feature
# of each models in NCBI Data.So only gene features and other features which doesnt 
# comprise of gene id info are filtered from NCBI and merged into NCBI through foverlap.

#########################################################################################

NCBI_Mod_df <- Processed_NCBI_Data %>% 
              filter(type == "gene" | type == "pseudogene" | type == "D_loop"| 
              type == "match" | type == "orgin_of_replication" | 
              chrloc == "MT") %>% select(chrloc,strand,start,end,modified_attributes)

XGC_Mod_df <- XGC_Data[,c('chrloc','strand','start','end')]

########## foverlap algorithm ##############

require(data.table)
setDT(NCBI_Mod_df)
setDT(XGC_Mod_df)
setkey(XGC_Mod_df)

overlaps <- foverlaps(NCBI_Mod_df,XGC_Mod_df,type = "any",nomatch = NA)

overlaps_D <- overlaps %>% filter(!is.na(start) & !is.na(end)) %>% select(chrloc,strand,start,end,modified_attributes) %>% distinct(.keep_all = TRUE)

###### Note: overlap data has already been merged with modified attributes which are the 
### desired gene ids

## merging overlap data to XGC Data to merge them with desired gene ids.

Overlap_NCBI_XGC_Data <- XGC_Data %>% left_join(overlaps_D,by = c("chrloc","strand","start","end")) %>% select(chrloc,source,type,start,end,score,strand,phase,modified_attributes,attributes) %>% distinct(.keep_all = TRUE)

###Arranging Data and manipulating multiple overlapped gene models leading to as much modified attributes consistent
### with XGC Data.Processing such info into modified attributes by grouping them

Final_Overlap_NCBI_XGC_D <- Overlap_NCBI_XGC_Data %>% group_by(chrloc,source,type,start,end,score,strand,phase,attributes) %>% mutate(
  
  Modified_Attr = paste0(unique(modified_attributes),collapse = ";"),
  Modified_Attr_Count = length(unique(modified_attributes))
  
) %>% ungroup() %>% select(chrloc,source,type,start,end,score,strand,phase,Modified_Attr,Modified_Attr_Count,attributes) %>% distinct(.keep_all = TRUE)

#Final_D <- Final_Overlap_NCBI_XGC_D %>% mutate(
  
 # Attr = ifelse( Modified_Attr_Count > 2,
   #              paste0(str_split(Modified_Attr,pattern = ";")[[1]][1:2],collapse = ";"),
  #               
    #            ifelse(
     #             Modified_Attr == 2,
      #            str_split(Modified_Attr,pattern = ";")[[1]][1],
                  
       #           ifelse(
        #            Modified_Attr_Count == 1,
         #           Modified_Attr,NA
          #          
           #       )
                   
            #    )
            
             # )
  
         #  )

#Final_D <- Final_Overlap_NCBI_XGC_D %>% mutate(
  
  
  
 # Attr = ifelse(Modified_Attr_Count > 2,
   #             str_split(Modified_Attr,pattern = ";")[[1]],
  #              ifelse(Modified_Attr_Count == 2,
    #                   str_split(Modified_Attr,pattern = ";")[[1]],
     #                  ifelse(Modified_Attr_Count == 1,
      #                        Modified_Attr,
       #                       NA
        #                
         #              )
          #        
 #               )       
                       
                       
          #)
  
    
#  )
  
  
#Final_Df <- Final_D %>% mutate(
  
  
 # Final_Attr =    ifelse(is.na(Attr) & Modified_Attr_Count  == 2,      
  #                str_split(Modified_Attr,pattern = ";")[[1]][1],
 #                 ifelse(is.na(Attr) & Modified_Attr_Count == 1,      
   #               Modified_Attr,Attr))
#)


###checking for unique XGC Models

Unique_XGC_Models <- Overlap_NCBI_XGC_Data %>% filter(is.na(modified_attributes)) %>% distinct(.keep_all = TRUE) #39219 unique XGC Models

###checking for multiple gene ids and eliminating unecessary identifiers

Mulitple_Mod_Attr_Overlap_NCBI_XGC_D <- Final_Overlap_NCBI_XGC_D %>% filter(Modified_Attr_Count == 2) %>% select(Modified_Attr,Modified_Attr_Count) %>% distinct(.keep_all = TRUE)
##4180 is the number of features with multiple gene ids.

Mulitple_Mod_Attr_Overlap_NCBI_XGC_D_2m <- Final_Overlap_NCBI_XGC_D %>% filter(Modified_Attr_Count > 2) %>% select(Modified_Attr,Modified_Attr_Count) %>% distinct(.keep_all = TRUE)
#508

#### extracting identifiers for unique XGC Models ####

Final_Processed_XGC_Data <- Final_Overlap_NCBI_XGC_D %>% mutate(
  
  Mod_Attr = ifelse(Modified_Attr == "NA",gsub("locus_tag=","",str_extract(attributes,pattern = "locus_tag=XENTR_v[0-9]+")),Modified_Attr)
  
) %>% select(chrloc,source,type,start,end,score,strand,phase,Mod_Attr,Modified_Attr_Count,attributes) %>% distinct(.keep_all = TRUE)

Final_Processed_XGC_Data$Mod_Attr <- sapply(Final_Processed_XGC_Data$Mod_Attr,function(x) 
  
  ifelse(str_detect(x,pattern = ".+;.+"),paste0(str_split(x,pattern = ";")[[1]],collapse = ";"),x))

Final_Processed_XGC_Df <- Final_Processed_XGC_Data %>% mutate(
  
  Final_Attr = ifelse(
    
    Modified_Attr_Count == 2,
    str_split(Mod_Attr,pattern = ";"),
    Mod_Attr
  
  )

)

Final_Processed_XGC_Df$Final_Attr <- sapply(Final_Processed_XGC_Df$Final_Attr,function(x)
  
  ifelse(length(x) > 1,x[1],x))


Processed_XGC_final_df <- Final_Processed_XGC_Df %>% mutate(
  
  
  modified_attributes = ifelse( Modified_Attr_Count > 2,
                                ifelse(grepl("^LOC",Final_Attr),
                                  NA,
                                  str_split(Final_Attr,pattern = ";")
  
    
      ),Final_Attr    
  
    ) 
  
)
  

Processed_XGC_final_df$modified_attributes <- sapply(Processed_XGC_final_df$modified_attributes,
                                                     
function(x) ifelse(length(x) > 1,unique(as.vector(x[1])),x))
  

Processed_XGC_final_df <- Processed_XGC_final_df %>% mutate(
  
  #$modified_attributes <- sapply(Processed_XGC_final_df$modified_attributes, function(x)
  
  Mod_Attribute =  ifelse(is.na(modified_attributes),str_split(Final_Attr,";"),
                        modified_attributes)
  )
  
Processed_XGC_final_df$Mod_Attribute <- sapply(Processed_XGC_final_df$Mod_Attribute,function(x) 
           ifelse(length(x) > 1,x[which(!str_detect(x,"^LOC"))],x))

Processed_XGC_final_df$Mod_Attribute <- sapply(Processed_XGC_final_df$Mod_Attribute,function(x)
  
  ifelse( length(x) > 1,unique(as.vector(x)),x 
          
  )
  
)

Processed_XGC_final <- Processed_XGC_final_df %>% mutate(
  
  
  Final_Att = ifelse(is.na(Mod_Attribute),str_split(Final_Attr,pattern = ";"),Mod_Attribute)
  
  
  )


Processed_XGC_final$Final_Att <- sapply(Processed_XGC_final$Final_Att,function(x)
  
  
  ifelse(length(x) > 1,unique(as.vector(x[1])),x)
  
  
) 


Processed_XGC_df <- Processed_XGC_final %>% select(chrloc,source,type,start,end,score,strand,phase,Final_Att,attributes) %>% distinct(.keep_all = TRUE)

##correcting format for start and end

#Final_Processed_XGC_Data$start <- sapply(Final_Processed_XGC_Data$start,function(x) format(x,scientific = F))  
#Final_Processed_XGC_Data$end <- sapply(Final_Processed_XGC_Data$end,function(x) format(x,scientific = F))  

##checking for existence of  any missing gene id information

abg <- Processed_XGC_df %>% filter(is.na(Final_Att) | Final_Att == "NA")

######################################################################################

#### writing the Final Processed XGC Data into a gff file to be infed into fjoin #####

######################################################################################

write.table(Processed_XGC_df,file = "Processed_XGC.gff3",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE) #writing processed XGC Data into gff file

#########################################################################################