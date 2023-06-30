.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")


library(tidyverse)
library(readxl)
library(dplyr)
library(Hmisc)
## here are the codes I was using
aaron <- read_tsv("E:/Excel/Metagenomic data/BBC/bbc/aaron_codes.txt", show_col_types = FALSE) %>%
  dplyr::select(2, 5, 3) %>%
  setNames(c("subject_id", "treatment", "stage")) %>%
  mutate_all(function(x) {gsub("-.*", "", str_to_lower(x))})

## There are 58 metabolomics samples.
## The metadata variables are named differently:
## - CLASS is the same as treatment. Isabel coded treatments as:
##    - `1` = "yakult"
##    - `2` = "kefir"
##    - `3` = "inulin"
## - pre/post is the same as stage. Isabel coded groups are coded as:
##    - `0` = "pre"
##    - `1` = "post"
## - ID is the same as subject_id. We coded subject_id the same.
## The command below will recode Isabel's metadata variables to match ours.

isabel <- read_tsv("E:/Excel/Metagenomic data/BBC/bbc/isabel_codes.txt", show_col_types = FALSE) %>%
  mutate(treatment = recode(
    CLASS,
    `1` = "yakult",
    `2` = "kefir",
    `3` = "inulin"
  )) %>%
  mutate(stage = recode(
    `pre/post`,
    `0` = "pre",
    `1` = "post"
  )) %>%
  mutate(subject_id = ID)

## import the NMR results
nmr <- read_excel("E:/Excel/Metagenomic data/BBC/bbc/nmr_results.xlsx", skip = 4) %>%
  dplyr::select(-c(1:7)) %>%
  inner_join(isabel, ., by = c("CLASS", "pre/post", "ID")) %>%
  mutate(sample_id = paste0(treatment, "_", stage, "_", subject_id)) %>%
  dplyr::select(ncol(.), 7:(ncol(.) - 1)) %>%
  arrange(sample_id)

# import scfa results

scfa_data <- read_xlsx(
  "E:/STORE N GO/Writings/BBC/bbc_study/isabel_metabolomics/Results SCFA urine BBC.xlsx",sheet = "Sheet5")


colnames(scfa_data)[1] <- "sample_id"


scfa_data$sample_id <- 
     gsub("yoghurt drink_|Yoghurt drink_", "yakult_", 
                      gsub("Inulin","inulin", 
                      gsub("Kefir","kefir",
                      gsub("Pre", "pre", 
                      gsub("Post","post",scfa_data$sample_id)))))


nmr <- merge(scfa_data, nmr,by="sample_id")# nothing significant in scfa's 
# is there anything obvious in the quantity of metabolites between colonization resistant and susceptible 
nmr_kefir <- nmr[grep("kefir",nmr$sample_id),]
metadata <-  data.frame(sample_id=nmr$sample_id) %>%
  separate(sample_id, c("treatment", "stage", "V4"), "_", remove = F) %>%
  #dplyr::select("treatment", "stage", "V4") %>%
  mutate(treatment = recode(treatment, "fmp" = "FMP", "inulin" = "Inulin", "kefir" = "Kefir")) %>%
  mutate(stage = recode(stage, "pre" = "Pre", "post" = "Post")) %>%
  mutate(stage = factor(stage, levels = c("Pre", "Post")))




########################################################################################################################
#PCOA
########################################################################################################################
data <- scfa_data
rows <- data$sample_id
  data <- dplyr::select(data,-c("sample_id"))
  rownames(data) <- rows
species_dist <- vegdist(as.matrix(data), method="bray")
species_pcoa <- pcoa(species_dist)

points <- species_pcoa$vectors[, 1:2] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  merge(., metadata, by = "sample_id")

pcoa_values <- data.frame(relative_eig = species_pcoa$values$Relative_eig) %>%
  dplyr::mutate(axis = paste0("Axis.", row_number())) %>%
  dplyr::mutate(relative_eig = 100 * relative_eig) %>%
  dplyr::select(axis, relative_eig)

axis_1 <- paste0("Axis 1 [", round(pcoa_values[1, 2], 2), "%]")
axis_2 <- paste0("Axis 2 [", round(pcoa_values[2, 2], 2), "%]")

xmin <- 0.95 * min(points$Axis.1)
ymax <- 0.95 * max(points$Axis.2)
permanova_tab1 <- c()
permanova_tab <- c()
rval <-c()
pval <-c()
for (category in levels(as.factor(metadata$treatment))){
  subset_meta1 <- metadata %>%
    filter(treatment == category ) %>%
    column_to_rownames("sample_id")
  
  sample_ids1 <- rownames(subset_meta1)
  
  subset_data1 <- data[which(rownames(data) %in% sample_ids1),] 
  
  
  # test to see if there is dissimilairty in the microbiome due to the adoption of lc.raffinolactis
  # test <- subset_meta1[which(subset_meta1$treatment==i),]
  # test$raffino <- "Absent"
  # test$raffino[which(rownames( test) %in% c("kefir_post_12_S24", "kefir_post_19_S38","kefir_post_23_S46", "kefir_post_28_S55"))] = "Present"
  # 
  # 
  # 
  # 
  # aov_tab <- adonis(subset_data1 ~ raffino, data = test)$aov.tab
  
  aov_tab <- adonis2(subset_data1 ~ stage, data =   subset_meta1) %>%
  as.data.frame() %>%
    drop_na() %>%
    dplyr::select(`R2`, `Pr(>F)`)
  
  rval <- paste0("R^2 = ", format(100 * round(aov_tab[1, 1], 3), nsmall = 1), "%")
  pval <- paste0("p = ", format(round(aov_tab[1, 2], 3), nsmall = 3))
  
  permanova_tab <- data.frame(Axis.1 = xmin, Axis.2 = ymax, label = paste0(rval, "<br>", pval), name=category)
  
  permanova_tab1 <- rbind(permanova_tab,permanova_tab1)
}

chrom_order <- c("FMP","Inulin", "Kefir")

permanova_tab1 <- 
  permanova_tab1[order(match(permanova_tab1$name, chrom_order)),]

points$stage <- gsub("Pre", "Before",points$stage )

#################

#plot PCOA-SCFA no longer in use 
#################
# p <- ggplot(points, aes(x = Axis.1, y = Axis.2 )) +
#   geom_point(aes(fill=stage),color="black",pch = 21, size = 8) +
#   
#   facet_wrap(~ treatment) +
#   stat_ellipse(aes(color=stage))+
#   # stat_ellipse(geom = "polygon",
#   #              alpha = 0,
#   #              type = "norm"
#   #              )+
#   theme_bw()+
#   labs( x = axis_1, y = axis_2,title = "Bray-Curtis PCoA plot of SCFA")+
#   theme(#plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
#     #box.color = "black",size=25, lineheight = 2),
#     plot.title = element_text(face = "bold", size = 25,hjust = 0.5),
#     axis.title = element_text(face = "bold", size = 12.5),
#     axis.text = element_text(size = 20),
#     legend.title = element_text(face = "bold", size = 12.5),
#     legend.position = "right",
#     strip.text = element_text(face="bold",size = 20),
#     legend.key.size = unit(1, 'cm'), #change legend key size
#     legend.key.height = unit(1, 'cm'), #change legend key height
#     legend.key.width = unit(1, 'cm') #change legend key width
#   ) +
#   scale_fill_manual(name="Stage",values = c("Before" = "#899DA4", "Post" = "#C93312"))+
#   scale_colour_manual(name="Stage",values = c("Before" = "#899DA4", "Post" = "#C93312"))

# Add text labels r  

#f_labels <- data.frame(treatment =  permanova_tab1$name, label = permanova_tab1$label)
# p <- 
#   p+
#   ggtext::geom_richtext(x =   permanova_tab1$Axis.1[1], y = permanova_tab1$Axis.2[1]*1.45 ,size=7, aes(label = label),
#                         
#                         data =  f_labels,
#                         #aes(x = Axis.1, y = Axis.2, label = label),
#                         fill = NA,
#                         hjust = 0,
#                         label.colour = NA,
#                         label.padding = grid::unit(rep(0, 4), "pt")
#   )


#png(filename=paste("F:/STORE N GO/R/Plots/BBC/Functional_output/",i,".png",sep=""), width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
#plot(p)
#graphics.off()


########################################################################################################################
#Explore if increase in metabolites between interventions is significant
########################################################################################################################

#nmr_raf <- nmr
#nmr_raf$raffionlactis <- "Absent"

#nmr_raf$raffionlactis[which(nmr_raf$sample_id %in% c("kefir_post_12","kefir_post_19","kefir_post_23","kefir_post_28"))] <- "Present"


isabel$sample_id <- paste(isabel$treatment,isabel$stage,isabel$subject_id,sep="_")
#nmr_raf$stage <- "Post"

#nmr_raf$stage[grep ("pre", nmr_raf$sample_id)] <- "Pre"

#nmr_raf <-

#merge(nmr_raf,dplyr::select(isabel,treatment,stage,subject_id,sample_id),by="sample_id",all.x=TRUE, sort=FALSE)


Metabolite_data <- nmr %>% pivot_longer(!sample_id,names_to = "nmr", values_to = "scfa")


Metabolite_data <- merge(Metabolite_data,  dplyr::select(isabel,treatment,stage,subject_id,sample_id),by="sample_id",all.x=TRUE, sort=FALSE)


id <- 12
data <- c()
metabolite <- c()
metabolite <-"Acetic (uM)"

growth <- c()
growth_metabolite <- data.frame(sample_id=as.character(),
                                nmr=as.character(),
                                growth=as.numeric(),
                                treatment=as.character(), 
                                subject_id=as.character())

t <- c()


for (id in unique(Metabolite_data$subject_id)){
  
  data <- subset(x=Metabolite_data,
         subset =Metabolite_data$subject_id==id )
  
  
  
  for (metabolite in levels(as.factor(data$nmr))){
    
    # get the growth of each metabolote per sample
  growth <-   subset(x=data,
           subset = nmr ==metabolite & stage=="post") %>%  dplyr::select(scfa) -
    subset(x=data,
           subset = nmr ==metabolite & stage=="pre")%>%  dplyr::select(scfa)
  
  
  # save the details of each participant including grwoth
  t <- data.frame(sample_id=as.character(unique(paste(data$treatment,data$subject_id,sep="_"))),
                  nmr=as.character(metabolite),
                  growth=as.numeric(growth),
                  treatment=as.character(unique(paste(data$treatment))),
                  subject_id=as.character(id))
  
# save to a bigger data frame
  growth_metabolite <- rbind(growth_metabolite,t)
  }
  
}


growth_metabolite$type <- "H NMR Metabolic profiling analysis"


growth_metabolite$type[which(growth_metabolite$nmr %in% colnames(scfa_data))] <- "SCFACs analysis"






# plot metabolite growth 

p_scfa <- 
  ggplot(subset(x=growth_metabolite,
                  type=="SCFACs analysis"), aes(x=nmr, y=growth, fill=treatment)) +
  geom_boxplot(color="black")+
  geom_point(
    color="black",
    size=2,
    pch = 21, 
    alpha = 0.5, 
    position = position_jitterdodge(jitter.width = .25))+
  #geom_jitter()+
  facet_wrap(~ nmr , scales = "free")+
  theme_bw()+
  labs(y="Change", x="SCFA",fill="Intervention")+ 
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold.italic"),
        strip.text = element_text(face="bold",size = 20),
        axis.title = element_text(face = "bold", size = 25.5),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #
        legend.text = element_text(size = 15)
  )

p_nmr <- 
  ggplot(subset(x=growth_metabolite,
                type=="H NMR Metabolic profiling analysis"), aes(x=nmr, y=growth, fill=treatment)) +
  geom_boxplot(color="black")+
  geom_point(
    color="black",
    size=2,
    pch = 21, 
    alpha = 0.5, 
    position = position_jitterdodge(jitter.width = .25))+
  facet_wrap(~ nmr , scales = "free")+
  theme_bw()+
  labs(y="Change", x="nmr",fill="Intervention")+ 
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold.italic"),
        strip.text = element_text(face="bold",size = 20),
        axis.title = element_text(face = "bold", size = 25.5),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #
        legend.text = element_text(size = 15)
  )


#png(filename="E:/STORE N GO/R/Plots/BBC/Supplementary figures/Figure S1A.png", width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
#plot(p_scfa)
#graphics.off()


#png(filename="E:/STORE N GO/R/Plots/BBC/Supplementary figures/Figure S1B.png", width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
#plot(p_nmr)
#graphics.off()




metabolite <- "Acetic (uM)"
metabolite <- "D-Glucose"

# wilcoxon to see if growth differs between interventions per each metabolite

growth_metabolite$raf<- "Absent"


growth_metabolite$raf[which(growth_metabolite$sample_id %in% c("kefir_12","kefir_19","kefir_23","kefir_28"))] <- "Present"



stats <-  as.data.frame(matrix(nrow=length(levels(as.factor(growth_metabolite$nmr))), ncol=1 , dimnames=list(c(levels(as.factor(growth_metabolite$nmr))),c("p_value"))))

stat <- c()
  for (metabolite in levels(as.factor(growth_metabolite$nmr))){
    
  
  x <- subset(x=growth_metabolite,  #[which(growth_metabolite$treatment=="kefir"),]
              subset=raf=="Present"&
                nmr==metabolite)

y <- subset(x=growth_metabolite,  #[which(growth_metabolite$treatment=="kefir"),]
            subset=raf=="Absent" &
              nmr==metabolite) 


stat <- wilcox.test(x$growth,y$growth, alternative = "two.sided")

stats[which(rownames(stats)==metabolite),"p_value"] <- stat$p.value

   }
    
    
    
stats$fdr <-   
p.adjust(stats$p_value, method = "fdr")


########################################################################################################################
#Explore if increase in metabolites between interventions is significant
########################################################################################################################

growth_metabolite <-subset(x=growth_metabolite[which(growth_metabolite$nmr %in% c("Acetic (uM)", "N,N-Dimethylglycine","Betaine","Succinic acid" )),],
       subset = treatment=="kefir" )


growth_metabolite <- merge(growth_metabolite ,stats,by.x="nmr",by.y=0,all.x=TRUE)



p_metabolite <- 
ggplot(growth_metabolite, aes(x=nmr, y=growth, fill=as.factor(raf))) +
  geom_boxplot(color="black", outlier.shape = NA)+
  geom_point(
    color="black",
    size=2,
    pch = 21,
    alpha = 0.5,
    position = position_jitterdodge(jitter.width = .25))+
  facet_wrap(~ nmr +type , scales = "free")+
  theme_bw()+
  labs(y="Change", x="",fill="Lactococcus raffinolactis")+ 
  theme(axis.text.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(face = "bold.italic",size=20),
        strip.text = element_text(face="bold",size = 20),
        axis.title = element_text(face = "bold", size = 17.5),
        axis.text = element_text(size = 12.5),
        legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #
        legend.text = element_text(size = 17.5)
        )+
  scale_fill_manual(values = c("Absent" = "#899DA4", "Present" = "#C93312"))

# Gt the top value for each metabolite of interest 
text <- growth_metabolite %>% group_by(nmr) %>%  filter(growth==max(growth))%>% as.data.frame


p_metabolite <- 
  p_metabolite +
  geom_text(data=text ,aes( x=.55 ,y=growth,label=paste( "p value = ", round(p_value,digits = 2),sep="")),size=7) #,"\n", "FDR = ", round(text$fdr)



########################################################################################################################
#Load and modify bracken species 
########################################################################################################################
bracken_species <- read_delim("E:/Excel/Metagenomic data/BBC/kraken2/Bracken/merged_bracken_report.report", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
positions_to_remove <- grep (".report_num",colnames(bracken_species))
bracken_species <- bracken_species[,-c(positions_to_remove)]

species_names <- bracken_species$name
bracken_species <- bracken_species[,-c(1,2,3)]
rownames(bracken_species) <- species_names


bracken_species <- as.data.frame(t(bracken_species))
id <- rownames(bracken_species) 
colnames(bracken_species) <- gsub(" ","_",colnames(bracken_species) )#modify here

colnames(bracken_species) <- paste("s__",colnames(bracken_species),sep="")

#bracken_species[,which(duplicated(colnames(bracken_species)))] # check if there is any duplicates 

bracken_species <- sapply(bracken_species, as.numeric)*100 # multiply by 100 to get the percentage 

#maxab_species <- apply(bracken_species,1, max, na.rm=TRUE)

#bracken_species <- bracken_species[-c(which(maxab_species < 1)),]
#bracken_species_rows <- bracken_species_rows[-c(which(maxab_species < 1))]
bracken_species <- as.data.frame(bracken_species)


id <- gsub("_S.*","",id)

#rownames(bracken_species) <- id
bracken_species$sample_id <-  id 
bracken_species$sample_id <- gsub("fmp_","yakult_",bracken_species$sample_id )
bracken_species <- bracken_species[-c(which(duplicated(bracken_species$sample_id))),] 
bracken_species <- bracken_species %>%
  relocate(sample_id, .before=colnames(bracken_species)[1]) 



bracken_species <- 
bracken_species[order(bracken_species$sample_id, nmr$sample_id),]


identical(bracken_species$sample_id, nmr$sample_id)
bracken_species[-c(which(bracken_species$sample_id %in% c("kefir_post_12","kefir_post_19","kefir_post_23","kefir_post_28"))),which(colnames(bracken_species)=="s__Lactococcus_raffinolactis")] <- 0.000



Metabolite_data_wide <- Metabolite_data %>% pivot_wider(names_from = nmr,values_from = scfa )





#bracken_species <- bracken_species %>%
  #relocate(sample_id, .before="s__Lactococcus_virus_ASCC465") 
bracken_species_test <- bracken_species %>% dplyr::select(-c(sample_id)) %>% as.matrix()

## run the correlations on raw metabolite data 
output <- Hmisc::rcorr(
  as.matrix(bracken_species[,-1]),
  as.matrix(nmr[,-1]),
  type = "spearman"
)






## run the correlations on growth metabolite data 
output_growth <- Hmisc::rcorr(
  as.matrix(bracken_species[,-1]),
  as.matrix(Metabolite_data_wide [,-c(1,2,3,4)]),
  type = "spearman"
)



## function to tidy the outputs
tidy_corr <- function(x, value) {
  
  x %>%
    as.data.frame() %>%
    rownames_to_column("species") %>%
    filter(species %in% colnames(bracken_species[,-1])) %>%
    dplyr::select(species, all_of(colnames(nmr[,-1]))) %>%
    pivot_longer(!species, names_to = "metabolite", values_to = value)
  
}



## get the r-values
r_values <- output$r %>%
  tidy_corr(., "r")


r_values_growth <- output_growth$r %>%
  tidy_corr(., "r")

## get the p-values and do p-value correction
p_values <- output$P %>%
  tidy_corr(., "p") %>%
  mutate(fdr = p.adjust(p, method = "fdr"))


p_values_growth <- output_growth$P %>%
  tidy_corr(., "p") %>%
  mutate(fdr = p.adjust(p, method = "fdr"))



## combine the results
combo <- inner_join(r_values, p_values, by = c("species", "metabolite")) %>%
  arrange(desc(r * r))




combo_growth <- inner_join(r_values_growth, p_values_growth, by = c("species", "metabolite")) %>%
  arrange(desc(r * r))

combo <- subset(x=combo,
                subset=species=="s__Lactococcus_raffinolactis") #&
                  #r >=0.2)

combo$type <- "H NMR Metabolic profiling analysis"


combo$type[which(combo$metabolite %in% colnames(scfa_data))] <- "SCFACs analysis"


# plot correlations for lc. raffinolactis
library(ggplot2)


p_cor <- 
ggplot(combo ,aes(r,reorder(metabolite,r))) +geom_bar(stat="identity",fill="coral1",color="black",width = 1,size=1.5)+
  #geom_text(aes(label=paste("p =", round(p,2))),nudge_x = ifelse(combo$r > 0, 0.01, -0.01))+
  #facet_wrap(~type,ncol=1, scales='free_y',margins=TRUE)+
  facet_grid(type ~ . , scales = "free_y", space = "free_y")+
  theme_bw()+
  labs(y="Metabolite", x="R")+ 
  theme(
        strip.text = element_text(face="bold",size = 20),
        axis.title = element_text(face = "bold", size = 17.5),
        axis.text = element_text(size = 12.5)
        
  )



range(
combo[which(combo$species=="s__Lactococcus_raffinolactis"),"r"])
range(
  combo_growth[which(combo_growth$species=="s__Lactococcus_raffinolactis"),"r"])



library(ggpubr)
p <- ggarrange(p_cor,
          p_metabolite,ncol=1, labels = c("A.", "B."))



png(filename="E:/STORE N GO/R/Plots/BBC/Figure 5.png", width = 7864, height=5200,res =300,pointsize = 15)


plot(p)
graphics.off()



#write.table(combo, "E:/Excel/Metagenomic data/BBC/correlations/metabolites/rcorr_output.tsv", sep = "\t", quote = F, row.names = F)

#write.table(combo_growth, "E:/Excel/Metagenomic data/BBC/correlations/metabolites/rcorr_output_growth.tsv", sep = "\t", quote = F, row.names = F)








