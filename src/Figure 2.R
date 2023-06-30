#######################################################################################################
#Library load
#######################################################################################################

.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")
#remove.packages("ggtree")

pacman::p_load(readr,ggpubr,readxl,devtools,taxize,rotl,ape,treeio,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)

#######################################################################################################
# Load Roary data 
#######################################################################################################

Roary_raf <- read_csv("E:/Excel/Metagenomic data/BBC/Roary/Lactococcus raffinolactis/gene_presence_absence.csv")

Roary <- c()
Roary[["raf"]] <- Roary_raf 



Roary_summary_raf <- read_delim("E:/Excel/Metagenomic data/BBC/Roary/Lactococcus raffinolactis/summary_statistics.txt", 
                                delim = "\t", escape_double = FALSE, 
                                col_names = FALSE, trim_ws = TRUE)

Roary_summary <- c()
Roary_summary[["raf"]] <- Roary_summary_raf 




for (name in names(Roary_summary)){
  strain_specific <- 
    length(which(Roary[[name]]$`No. isolates`==1))
  Roary_summary[[name]] <- rbind(Roary_summary[[name]], c("Strain specific genes", "(1% of strains)", strain_specific))
  
  
  Roary_summary[[name]]$X3[which(Roary_summary[[name]]$X1=="Shell genes")] <- as.numeric(Roary_summary[[name]]$X3[which(Roary_summary[[name]]$X1=="Shell genes")]) -strain_specific
  
}



#myPalette <- brewer.pal(5, "Set2")

#png(filename="E:/STORE N GO/R/Plots/BBC/pangenome_lc_raffinolactis.png",width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)

#pie(c(as.numeric(Roary_summary$X3[1]),as.numeric(Roary_summary$X3[3]),as.numeric(Roary_summary$X3[[6]])), labels =levels(as.factor(Roary_summary$X1[c(1,3,6)])), border="white", col=myPalette,main = "")

#graphics.off()

#######################################################################################################
# Load MAG data and incorprate with reference data 
#######################################################################################################

high_med_quaility_MAGs_summary <- read_csv("E:/Excel/Metagenomic data/metawrap/high_med_quaility_MAGs_summary.csv")

high_med_quaility_MAGs_summary[grep("Lactococcus raffinolactis",high_med_quaility_MAGs_summary$Species),"Assembly id"]

metabolism_list <- c()

setwd("E:/Excel/Metagenomic data/dram/mag annotation_metawrap/")
temp = list.files(pattern="Lactococcus_raffinolactis*")



#getElement(high_med_quaility_MAGs_summary[grep("Lactococcus raffinolactis",high_med_quaility_MAGs_summary$species),],"id")
temp <- c(temp, paste(getElement(high_med_quaility_MAGs_summary[grep("Lactococcus raffinolactis",high_med_quaility_MAGs_summary$species),],"id"), "_mag_distill",sep = ""))

temp <- paste(temp, "/metabolism_summary.xlsx",sep="")



temp <- temp[-c(which(temp== "V8_S69_merge_host_removed_bin.2_mag_distill/metabolism_summary.xlsx"))]
temp <- temp[-c(which(temp== "V7_S68_merge_host_removed_bin.2_mag_distill/metabolism_summary.xlsx"))]






# create function to read multiple sheets per excel file
#######################################################################################################
# Import metabolism data 
#######################################################################################################
read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  sapply(sheets, function(f) as.data.frame(readxl::read_excel(filename, sheet = f)), 
         simplify = FALSE)
}



#######################################################################################################
# Load health terms 
#######################################################################################################
metabolism_list <- c()
metabolism_list [["raf"]] <- lapply(temp,read_excel_allsheets)


names(metabolism_list[["raf"]]) <- gsub("_mag_distill/metabolism_summary.xlsx","",temp) 

dram_total <- c()
dram_total_wide <- c()

for (species in names(metabolism_list)){
  
  
  
  for (i in names(metabolism_list[[species]])){
    for (j in names(
      metabolism_list[[species]][[i]])){
      data <- metabolism_list[[species]][[i]][[j]] 
      data$category <- j
      data$sample_id <- i
      colnames(data)[6] <- "gene_count"
      metabolism_list[[species]][[i]][[j]]<- data[which(data[,"gene_count"]>0),]
      
      
      dram_total[[species]] <- rbind(dram_total[[species]], metabolism_list[[species]][[i]][[j]])
      
      dram_total[[species]] <- dram_total[[species]][!duplicated(dram_total[[species]]),]
      data.frame(genes=levels(as.factor(dram_total[[species]]$gene_id))) %>% 
        filter(grepl("K0", genes)) %>% 
        mutate(gene_number=paste("genes",c(1:nrow(.)),sep=""), .before=genes) #%>% 
       # write.table(paste("E:/Excel/Metagenomic data/BBC/KEGG Mapper/ko_gene_list_dram_",species,".txt",sep=""),col.names = FALSE,row.names = FALSE, quote = FALSE)
    }
    
  }
  
  
  
  #dram_total[[species]]$wider_value <-  paste(dram_total[[species]]$gene_id, dram_total[[species]]$gene_description,sep="_")
  dram_total_wide[[species]] <-   dram_total[[species]] %>%  pivot_wider(names_from = "sample_id", values_from = "gene_count")
  
}


##############################
amino_acid <- c()
biosynthesis <- c()
category <- "Organic Nitrogen" 
header <- "Amino Acid"

category <- "carbon utilization"
data <- c()
category_data <- data.frame(categories=as.character(levels(as.factor(dram_total[[species]]$category))), raf=as.numeric(NA),cas=as.numeric(NA) )





for (species in names(dram_total)){
  
  for (category in levels(as.factor(dram_total[[species]]$category))){
    
    data[[species]][[category]] <- dram_total_wide[[species]][which(dram_total_wide[[species]]$category==category),]
    
    
    category_data[which(category_data$categories==category), which(colnames(category_data)==species)] <-   length(names(table( data[[species]][[category]]$gene_id)))
    
    
  }
}




##############################

category_data$categories[which(category_data$categories=="carbon utilization")] <- "carbon utilisation"

pc <- category_data[-c(which(category_data$categories== "carbon utilization (Woodcroft)")),] %>% pivot_longer(!categories, names_to = "species", values_to = "count") %>% 
  mutate( species=gsub("raf", "Lactococcus\n raffinolactis",species)) %>%
  mutate( species=gsub("cas", "Lacticaseibacillus\n casei",species)) %>% 
  ggplot(aes(x = species, y = count, fill = categories)) +
  geom_col(width = .5, position = "fill") +
  #geom_col(aes(x=0,y=0))+
  #coord_polar("y") +
  theme_bw()+
  labs(title="Pan- genome annotations - functional categories", y="Annotations",fill="Gene type")+
  scale_fill_viridis_d() +
  theme(plot.title = element_text(face = "bold", size = 25,hjust = 0.5),
        aspect.ratio=2,axis.title.x = element_blank(),
        axis.title = element_text(face = "bold", size = 17.5),
        axis.text = element_text(size = 12.5),
        legend.title = element_text(face = "bold", size = 17.5),
        legend.text = element_text(size = 12.5),
        axis.text.y = element_text(size=15),
        #axis.text.x = element_text(size=15, face="italic", vjust = -0.15, hjust=.5),
        # axis.text.x = element_text(size=15, face="italic", vjust = .5,hjust=-.04,angle = -45),
        axis.text.x = element_text(size=15, face="italic"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(1.5, 'cm'), #change legend key size
        legend.key.height = unit(1.5, 'cm'), #change legend key height
        legend.key.width = unit(1.5, 'cm'))



#check this plz
print("check this")
sum(category_data$raf)


#######################################################################################################
#update this section
#######################################################################################################
cazy <- c()
cazy_wide <- c()
peptidases <- c()
SCFA <- c()
for (species in names(dram_total)){
  
  cazy[[species]] <-  dram_total[[species]][c(which(dram_total[[species]]$header=="CAZY")),]
  # 
  peptidases[[species]] <- dram_total[[species]][c(which(dram_total[[species]]$header=="Peptidase")),]
  # 
  # 
  #t <- peptidases[grep("Catlytic type:",peptidases$gene_description),]
  # 
  SCFA[[species]] <- dram_total[[species]][c(which(dram_total[[species]]$header=="SCFA  and alcohol conversions")),]
  # #dplyr::select(cazy,gene_id,gene_count,sample_id) %>%  pivot_wider(names_from = gene_id,values_from = gene_count)
  # 
  cazy_wide[[species]] <- 
    cazy[[species]] %>%  pivot_wider(names_from = sample_id,values_from = gene_count)
  # 
  # 
  cazy_wide[[species]]$prevalence <- length(colnames(cazy_wide[[species]]))-
    rowSums(is.na( cazy_wide[[species]][7:length(colnames( cazy_wide[[species]]))]))
  # 
  # 
  
  for (item in levels(as.factor(cazy_wide[[species]]$module))){
    
    print(species)
    print(paste(item," - ",nrow(cazy[[species]][which(cazy_wide[[species]]$module==item),]),sep=""))
  }
}
#
cazy_of_interest <- c("GH36", "GH42","GH109", "GH29",  "GH20" ) #note actual one]

#View(cazy_wide[which(cazy_wide$gene_id %in% cazy_of_interest),])
#o_nitrogen <-  dram_total[c(which(dram_total$category=="Organic Nitrogen")),]

species <- "raf"
for (species in names(dram_total)){
  
  t <- peptidases[[species]][grep("Catlytic type:",peptidases[[species]]$gene_description),]
  
  t$gene_description <- gsub(";.*","",t$gene_description)
  
  print(species)
  print(levels(as.factor(t$gene_id)))
  print(length(levels(as.factor(t$gene_id))))
  print(table(t$gene_id))
  print(length(t$gene_description))
  #print(table(t$gene_description))
  
  
}



#######################################################################################################
# Now that i have the data in a usable format 
#######################################################################################################
KEGG <- c()
KEGG[["raf"]]<- read_csv("E:/Excel/Metagenomic data/BBC/KEGG Mapper/Kegg_complete_raf_modified.csv")

#######################################################################################################
#Library load
#######################################################################################################

.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")
#remove.packages("ggtree")

pacman::p_load(readr,ggpubr,readxl,devtools,taxize,rotl,ape,treeio,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)

#######################################################################################################
# Load prokka data 
#######################################################################################################
setwd("E:/Excel/Metagenomic data/BBC/Prokka")

temp <- c()
temp[["raf"]] = list.files(pattern="*.tsv")

#myfiles = lapply(temp,read_table2, comment="--")


myfiles <- c()
myfiles[["raf"]] = lapply(temp[["raf"]],read_delim,  delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

#######################################################################################################
# Load drep data 
#######################################################################################################


dereplicated_genomes <- read_csv("E:/Excel/Metagenomic data/drep/dereplicated_genomes.txt", 
                                 col_names = FALSE)
dereplicated_genomes$X1 <- gsub(".fa","",dereplicated_genomes$X1)

drep_clusters <- read_csv("E:/Excel/Metagenomic data/drep/data_tables/Cdb.csv")
# remove the .fa ending gsub interfers with the acetobacter_fabarum naming system 
drep_clusters$genome <- substr(drep_clusters$genome,1,nchar(drep_clusters$genome)-3)


drep_clusters$assembly_id <- drep_clusters$genome


pos_genomic <- grep("_genomic",drep_clusters$assembly_id)


drep_clusters$assembly_id[pos_genomic] <- sub('^([^_]+_[^_]+_[^_]+_[^_]+).*', '\\1',drep_clusters$assembly_id[pos_genomic])


#drep_clusters$genome <-  gsub("_ASM.*","",drep_clusters$genome)

d <- unlist(strsplit(drep_clusters$assembly_id[pos_genomic],"_"))

d1 <-
  paste(d[seq(3,length(d),by=4)],"_" ,d[seq(4,length(d),by=4)],sep="")

drep_clusters$assembly_id[pos_genomic] <- d1



drep_clusters <- drep_clusters[which(drep_clusters$assembly_id %in%  gsub(".tsv|_ASM.*","",temp[[species]])),]


drep_clusters$assembly_id <- 
  gsub("_merge_host_removed","",drep_clusters$assembly_id)




#######################################################################################################
# Get list of genes 
#######################################################################################################

prokka_gene_list <- c()

for (species in names(myfiles)){
  
  names(myfiles[[species]]) <- gsub("_merge_host_removed|.tsv","",temp[[species]])
  
  
  
  t <- c()
  gene_list <- data.frame(gene=character(),
                          product=character(),
                          ftype=character(),
                          sample=character())
  
  
  for (i in names(myfiles[[species]])){
    
    t <- unique(dplyr::select(myfiles[[species]][[i]],gene,ftype, product))
    
    t <- t[which(t$ftype=="CDS"),]
    # <- data.frame(genes=levels(as.factor((myfiles[[i]]$gene)))
    
    
    t$sample <- i
    
    gene_list <- rbind(gene_list,t)
    
    
    
  }
  gene_list$value <- 1
  
  
  
  gene_list_wide <- gene_list %>%  pivot_wider(names_from = sample,values_from = value,  values_fill = 0 )
  prokka_gene_list[[species]] <-  gene_list_wide
}


View(  dram_total[[species]])

View(  prokka_gene_list[[species]])
#######################################################################################################
# Load health terms 
#######################################################################################################


Health_search_terms <- read_excel("E:/Excel/Metagenomic data/Humann3/Health_search-terms_john.xlsx",sheet = 7)

health_genes <- c()
for (species in names(myfiles)){
  
  health_genes[[species]] <-  prokka_gene_list[[species]][which(prokka_gene_list[[species]]$gene %in% Health_search_terms$Gene),]
  
  
  health_genes[[species]] <- merge(Health_search_terms,   health_genes[[species]], by.x="Gene",by.y="gene")
  
  
  # calculate prevalance 
  
  
  health_genes[[species]]$prevalence <- 
    rowSums(  health_genes[[species]][6:length(colnames(  health_genes[[species]]))], na.rm=TRUE)
  
  
  # move prevalence column 
  health_genes[[species]]  <-  health_genes[[species]] %>% relocate(prevalence, .before=colnames(health_genes[[species]])[6])
  
  health_genes[[species]]$pangenome <- "Accessory"
  
  health_genes[[species]]$pangenome[which(health_genes[[species]]$prevalence==length(names(myfiles[[species]])))] <- "Core"
  health_genes[[species]]$pangenome[which(health_genes[[species]]$prevalence==1)] <- "Strain specific"
  
  
  #write.csv(  health_genes[[species]],paste("E:/Excel/Metagenomic data/BBC/Phagscs/",species,"phagscs_gene_list.csv",sep = ""),row.names = FALSE )
  
  
  
  
  
}

# Nmber of phagsc across refereces 

colSums(health_genes[["raf"]][,c(7:length(health_genes[["raf"]])-1)])




table(health_genes[["raf"]]$group)
data <- 
 health_genes[["raf"]]

colnames(data) <- gsub("_ASM.*","",colnames(data)  )


colnames(data)<- gsub("_bin.*","",colnames(data)  )

data[is.na(
  data)] <- 0
data <- data[-c(which(data$product=="Negative regulator of genetic competence ClpC/MecB")),]

col_details <- data.frame(Species=as.character("Lactococcus raffinolactis"), sample_id= as.character(names(myfiles[["raf"]]))) %>% remove_rownames %>% column_to_rownames(var="sample_id")

rownames(col_details) <- gsub("_ASM.*","",rownames(col_details)  )


rownames(col_details)<- gsub("_bin.*","",rownames(col_details)  )





library(RColorBrewer)

col_details$Source <- "Kefir"
col_details$Source[grep("GCF",rownames(col_details) )] <- "Reference genome"


drep_clusters$assembly_id <- gsub("_bin.*","",drep_clusters$assembly_id) 








library(RColorBrewer)


row_details <- dplyr::select(data, Gene,pangenome,pangenome,group) %>%remove_rownames %>% column_to_rownames(var="Gene")




colnames(row_details) <- c("Lactococcus raffinolactis","Gene type")
col_details$Source <- "Kefir"
col_details$Source[grep("GCF",rownames(col_details) )] <- "Reference genome"



rownames(row_details)[which(row_details$`Gene type`=="modulation")]

data <- 
  data[-c(which(data$Gene %in% rownames(row_details)[which(row_details$`Gene type`=="modulation")] )),]


row_details <- row_details[-c(which(row_details$`Gene type`=="modulation")),]


row_details$`Gene type`<- str_to_title(row_details$`Gene type`)


graph_data <- dplyr::select(data, Gene,rownames(col_details)) %>%remove_rownames %>% column_to_rownames(var="Gene")










col_details_v2 <- 
  merge(
    drep_clusters, col_details,by.x="assembly_id",by.y=0)



col_details_v2 <- dplyr::select(col_details_v2, assembly_id,secondary_cluster,Source) %>% column_to_rownames("assembly_id")

col_details_v2$secondary_cluster <- gsub("7_1","1",col_details_v2$secondary_cluster)
col_details_v2$secondary_cluster <- gsub("8_1","2",col_details_v2$secondary_cluster)

colnames(col_details_v2)[which(colnames(col_details_v2)=="secondary_cluster")] <- "Sub cluster"

colnames(row_details)[which(colnames(row_details)=="Lactococcus raffinolactis")] <- "Pan- genome"
pd <- 
  pheatmap::pheatmap(graph_data, annotation_col =col_details_v2, annotation_row= row_details, scale = "none",clustering_method="ward.D",
                     border_color= "black", #main = "Prevalence of health promoting genes",
                     #cellwidth = 99, 
                     #cellheight = 30,
                     #cluster_rows=FALSE,cluster_cols=FALSE,
                     color =brewer.pal(3,"Set1")[2:1],  fontsize = 15,angle_col = 45,legend_breaks = c(0,1),legend_labels = c("Absence","Presence"))





########################################################################################################################
#Research objectives in this script
# See if there is a differnece between sterile and non sterile conditions
# Identify pathogenic microbes
# see how household environment affects there prevalence 
################################################################################################################################################################################################################################################



# See if there is a differnece between sterile and non sterile conditions
#############################################################################################################################################################################################################################################

#Libraries used
########################################################################################################################

.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")
pacman::p_load(readxl,readr,reshape2,dplyr, gplots,Heatplus,vegan,RColorBrewer,tidyr,gtools,stringr,tidyverse,ComplexHeatmap,magick,viridis)
pacman::p_load(readxl,devtools,taxize,rotl,ape,treeio,ggtree,DECIPHER,ggdendro,ggplot2,tidyr,optmatch,rentrez,plyr,dplyr,RColorBrewer,stringr,scales)
library(vegan)
library(ggplot2)
library(grid)


########################################################################################################################
#Import metacache data
########################################################################################################################
#metacache <- read_csv("Q:/H2020 Master/Citizen Science Project/Results/04_short_read_profiling/04_Metacache/04_metacache_total_species_profile_v2.csv")
metacache <- read_csv("E:/Excel/Metagenomic data/Overall taxonomic profiling/combined/kraken_custom_species.csv")

metacache <- 
  metacache %>% dplyr::select(-...1) %>% column_to_rownames("name") %>% t() %>% as.data.frame() %>% rownames_to_column("sample_id")


global_mk_metadata <- read_csv("Q:/H2020 Master/Citizen Science Project/Citizen science metadata/sample_metadata/global_milk_kefir_metadata_v1.csv")



metacache <- 
  
  metacache[which(metacache$sample_id %in% global_mk_metadata$Sample[-c(which(global_mk_metadata$`kefir type`== "Medium control"))]),]



#metacache [is.na(metacache )] <- 0

#metacache <- 


# maxab_species <- apply(metacache %>%remove_rownames %>%  column_to_rownames("sample_id"),2, max, na.rm=TRUE)
# 
# n1_species <-names(which(maxab_species < .1))
# 
# metacache <- metacache[,-c(which(colnames(metacache)%in% n1_species)),]


colour_pal= c("#332288", "#6699CC", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#CC6677", "#AA4466", "#882255", "#AA4499")

metacache_long <- 
  metacache %>% pivot_longer(!sample_id, names_to = "species",values_to="RA") 


metacache_long$species_updated <- "Other"
metacache_long$species_updated[which(metacache_long$species=="Lactococcus raffinolactis")] <- "Lactococcus raffinolactis"

#metacache_long$RA <- 
#metacache_long$RA*100
pa <- 
ggplot(metacache_long, aes(x= sample_id , y=RA, fill=species_updated))+
  #geom_col()+
  geom_bar(position="fill",stat="identity")+
  theme_bw()+
  # scale_fill_manual(values = colour_pal) +
  theme(#legend.position = "none",
    axis.title = element_text(face = "bold", size = 17.5),
   # axis.text = element_text(size = 1.5),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 20.5),
    legend.text = element_text(size = 15.5),
    strip.text = element_text(face="bold",size = 20),
    legend.key.size = unit(1, 'cm'), #change legend key size
    legend.key.height = unit(1, 'cm'), #change legend key height
    legend.key.width = unit(1, 'cm'),
   axis.text.x = element_text(size = 20),
    
   # axis.text.x = element_text(size = 10, colour = "black", vjust = 0.5, hjust = 1, face= "bold"), 
    #axis.title.y = element_text(size = 1, face = "bold"),legend.title = element_text(size = 16, face = "bold"), 
    #legend.text = element_text(size = 8, face = "bold", colour = "black"),
    axis.text.y = element_blank(),
   #axis.text.y = element_text(size=12),
   panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
   panel.border = element_blank(),
   panel.background = element_blank())+ 
  scale_y_continuous(expand = c(0,0)) + 
  labs(x = "Sample", y = "Relative Abundance (%)", fill="Species")+
  coord_flip()+
  scale_y_continuous(labels=c("0.00"="0%","0.25" = "25%", "0.50" = "50%", "0.75" = "75%",
                              "1.00" = "100%"))+
  scale_fill_brewer(palette="Set1")















library(ggplotify)


png(filename="Q:/H2020 Master/Liam WALSH/BBC/Figure 2_v2.png", width = 7864, height=5200,res =300,pointsize = 15)

ggarrange(pa,
as.ggplot(
pd), labels = c("A.","B."),nrow=2,ncol=1, heights=c(.4,.8), widths = c(.5,.8))

graphics.off()



