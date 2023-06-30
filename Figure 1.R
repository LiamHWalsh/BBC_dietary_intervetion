########################################################################################################################
#Libraries used 
########################################################################################################################

.libPaths("E:/STORE N GO/R/R-4.0.2/win-library/4.0")
library(tidyverse)
library(data.table)
library(vegan)
library(ape)
library(cowplot)
library(Maaslin2)
library(dplyr)
library(ggtext)
########################################################################################################################


bracken_species <- read_delim("E:/Excel/Metagenomic data/BBC/kraken2/Bracken/merged_bracken_report.report", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)

positions_to_remove <- grep (".report_num",colnames(bracken_species))
bracken_species <- bracken_species[,-c(positions_to_remove)]


species_names <- bracken_species$name
species_names <-gsub("\\.","_",species_names)
species_names <-gsub(" ","_",species_names)
species_names <-gsub("\\("," ",species_names)
species_names <-gsub("\\)","",species_names)
species_names <-gsub("-"," ",species_names)
species_names <-gsub("\\'","",species_names)
species_names <-gsub("\\[","",species_names)
species_names <-gsub("\\]","",species_names)
species_names <-gsub("/","",species_names)
species_names <-gsub(":","",species_names)

bracken_species <- bracken_species[,-c(1,2,3)]
rownames(bracken_species) <- species_names


bracken_species <- as.data.frame(t(bracken_species))
id <- rownames(bracken_species) 
colnames(bracken_species) <- gsub(" ","_",colnames(bracken_species) )#modify here

colnames(bracken_species) <- paste("s__",colnames(bracken_species),sep="")

#bracken_species[,which(duplicated(colnames(bracken_species)))] # check if there is any duplicates 

bracken_species <- sapply(bracken_species, as.numeric)*100 # multiply by 100 to get the percentage 

#maxab_species <- apply(bracken_species,2, max, na.rm=TRUE)

#maxab_species <- as.data.frame(maxab_species)

#species_to_remove <- rownames(maxab_species)[which(maxab_species$maxab_species<1)]

#bracken_species <- bracken_species[,-c(which(maxab_species< 1))]
#bracken_species_rows <- bracken_species_rows[-c(which(maxab_species < 1))]
bracken_species <- as.data.frame(bracken_species)

id <- gsub("_bracken.species.report_frac","",id)
bracken_species$sample_id <-  id 



########################################################################################################################
# tidy input
########################################################################################################################


metadata <- bracken_species %>%
  relocate(sample_id, .before="s__Lactococcus_virus_ASCC465") %>% 
  separate(sample_id, c("treatment", "stage", "V4", "V5"), "_", remove = F) %>%
  dplyr::select(1:3) %>%
  mutate(treatment = recode(treatment, "fmp" = "FMP", "inulin" = "Inulin", "kefir" = "Kefir")) %>%
  mutate(stage = recode(stage, "pre" = "Pre", "post" = "Post")) %>%
  mutate(stage = factor(stage, levels = c("Pre", "Post")))



maaslin2_data <-  bracken_species %>%
  column_to_rownames("sample_id")

maaslin2_metadata <- metadata %>%
  column_to_rownames("sample_id")
########
# PCoA #
########

treatments <- levels(as.factor(metadata$treatment))

plot_list <- c()

i<- "FMP"
permanova_tab1 <-c()
#for (i in treatments) {
  
  subset_meta <- metadata %>%
    #filter(treatment == i) %>%
    column_to_rownames("sample_id")
  
  sample_ids <- rownames(subset_meta)
  
  subset_data <- bracken_species %>%
    filter(sample_id %in% sample_ids) %>%
    column_to_rownames("sample_id")
  
  species_dist <- vegdist(as.matrix(subset_data), method="bray")
  
  species_pcoa <- pcoa(species_dist)
  
  points <- species_pcoa$vectors[, 1:2] %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    merge(., metadata, by = "sample_id")
  
  pcoa_values <- data.frame(relative_eig = species_pcoa$values$Relative_eig) %>%
    mutate(axis = paste0("Axis.", row_number(.))) %>%
    mutate(relative_eig = 100 * relative_eig) %>%
    dplyr::select(axis, relative_eig)
  
  axis_1 <- paste0("Axis 1 [", round(pcoa_values[1, 2], 2), "%]")
  axis_2 <- paste0("Axis 2 [", round(pcoa_values[2, 2], 2), "%]")
  

  
    xmin <- 0.95 * min(points$Axis.1)
    ymax <- 0.95 * max(points$Axis.2)
    
    for (i in treatments) {
      subset_meta1 <- metadata %>%
        filter(treatment == i) %>%
        column_to_rownames("sample_id")
      
      sample_ids1 <- rownames(subset_meta1)
      
      subset_data1 <- bracken_species %>%
        filter(sample_id %in% sample_ids1) %>%
        column_to_rownames("sample_id")
      
      
      # test to see if there is dissimilairty in the microbiome due to the adoption of lc.raffinolactis
      # test <- subset_meta1[which(subset_meta1$treatment==i),]
      # test$raffino <- "Absent"
      # test$raffino[which(rownames( test) %in% c("kefir_post_12_S24", "kefir_post_19_S38","kefir_post_23_S46", "kefir_post_28_S55"))] = "Present"
      # 
      # 
      # 
      # 
      # aov_tab <- adonis(subset_data1 ~ raffino, data = test)$aov.tab
      
    aov_tab <- adonis(subset_data1 ~ stage, data = subset_meta1[which(subset_meta1$treatment==i),])$aov.tab %>%
      as.data.frame() %>%
      drop_na() %>%
      dplyr::select(`R2`, `Pr(>F)`)
    
    rval <- paste0("R^2 = ", format(100 * round(aov_tab[1, 1], 3), nsmall = 1), "%")
    pval <- paste0("p = ", format(round(aov_tab[1, 2], 3), nsmall = 3))
    
    permanova_tab <- data.frame(Axis.1 = xmin, Axis.2 = ymax, label = paste0(rval, "<br>", pval), name=i)
    
    permanova_tab1 <- rbind(permanova_tab,permanova_tab1)
    }
    
 chrom_order <- c("FMP","Inulin", "Kefir")
    
 permanova_tab1 <- 
    permanova_tab1[order(match(permanova_tab1$name, chrom_order)),]
   
 points$stage <- gsub("Pre", "Before",points$stage )
  p <- ggplot(points, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(fill=stage),color="black",pch = 21, size = 8) +
    
    facet_wrap(~ treatment) +
    stat_ellipse(aes(color=stage))+
    # stat_ellipse(geom = "polygon",
    #              alpha = 0,
    #              type = "norm"
    #              )+
    theme_bw()+
    labs( x = axis_1, y = axis_2,title = "Bray-Curtis PCoA plot based on Bracken output")+
    
    theme(#plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
      #box.color = "black",size=25, lineheight = 2),
      plot.title = element_text(face = "bold", size = 25,hjust = 0.5),
      axis.title = element_text(face = "bold", size = 17.5),
      axis.text = element_text(size = 25.5),
      legend.title = element_text(face = "bold", size = 25.5),
      legend.text = element_text(size = 20.5),
      legend.position = "right",
      strip.text = element_text(face="bold",size = 20),
      legend.key.size = unit(1.5, 'cm'), #change legend key size
      #legend.key.height = unit(1.5, 'cm'), #change legend key height
      #legend.key.width = unit(1.5, 'cm') #change legend key width
    ) +
     scale_fill_manual(name="Stage",values = c("Before" = "#899DA4", "Post" = "#C93312"))+
    scale_colour_manual(name="Stage",values = c("Before" = "#899DA4", "Post" = "#C93312"))+
  guides(colour="none",fill= guide_legend(override.aes = list(size=21)))


  # Add text labels r  
    
    f_labels <- data.frame(treatment =  permanova_tab1$name, label = permanova_tab1$label)
  p <- 
p+
  ggtext::geom_richtext(x = -0.4006096, y = 0.3705523 ,size=7, aes(label = label),
    
    data =  f_labels,
    #aes(x = Axis.1, y = Axis.2, label = label),
    fill = NA,
    hjust = 0,
    label.colour = NA,
    label.padding = grid::unit(rep(0, 4), "pt")
  )
  
  p

  #colnames(maaslin2_data)


###################
# alpha diversity #
###################

  indexs <- "shannon"
for (indexs in c("shannon","simpson")){
alpha_div <- vegan::diversity(maaslin2_data[], indexs) %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  rename(shannon=2) %>%
  merge(., metadata, by="sample_id")

alpha_div$stage <- gsub("Pre", "Before",alpha_div$stage)


treatments <- alpha_div %>%
  dplyr::select(treatment, stage) %>%
  distinct() %>%
  group_by(treatment, stage) %>%
  tally() %>%
  ungroup() %>%
  group_by(treatment) %>%
  tally() %>%
  ungroup() %>%
  filter(n > 1) %>%
  .$treatment

alpha_ann_total <- c()

alpha_wilcox <- alpha_div %>%
  #filter(treatment %in% products) %>%
  group_by(treatment) %>%
  rstatix::wilcox_test(shannon ~ stage) %>%
  ungroup() %>%
  dplyr::select(treatment, p)

#Get the impact of lc.raffinolactis presence in the kefir group
alpha_div_raffino <- 
alpha_div %>%
  filter(treatment %in% "Kefir") %>%
  group_by(treatment) %>%
  mutate(colonisation_raffino = "Absent")

alpha_div_raffino$colonisation_raffino[which(alpha_div_raffino$sample_id %in% c("kefir_post_12_S24", "kefir_post_19_S38","kefir_post_23_S46", "kefir_post_28_S55"))] = "Present"

alpha_div_raffino %>% 
  rstatix::wilcox_test(shannon ~ colonisation_raffino) %>%
  ungroup() %>%
  dplyr::select(treatment, p)




for (products in treatments){
alpha_ann <- alpha_div %>%
  filter(treatment %in% products) %>% 
  filter(shannon == max(shannon)) %>%
 mutate(shannon =  shannon) %>%
  dplyr::select(treatment, shannon) %>%
  mutate(p = paste0("p = ", alpha_wilcox[which(alpha_wilcox$treatment==treatment), "p"]))

alpha_ann_total <- rbind(alpha_ann_total,alpha_ann)
}

# if(indexs=="shannon"){
#   print("index = shannon")
# text_y <- max(alpha_ann_total$shannon) 
#   
#   }else{
#   print("index = simpson")
#   
#     }


alpha_ann_total$y<- 
  max(alpha_ann_total$shannon)+.05  




alpha_div$title <- str_to_title(indexs)

plot_list[[indexs]] <- 
  
  
  
  ggplot(alpha_div, aes(x = treatment, y = shannon)) +
  geom_boxplot(aes(fill=stage), outlier.shape = NA) +
  geom_point(
    aes(x = treatment, colour = stage, fill = stage),
    color="black",
    size=2,
    pch = 21, 
    alpha = 0.5, 
    position = position_jitterdodge(jitter.width = .25)
  ) +
  geom_text(data = alpha_ann_total, aes(x = treatment, y =  y[1], label = p),size=7) +
  facet_wrap(~title)+
  theme_bw() +
  theme( plot.title = element_text(face = "bold", size = 25,hjust = 0.5),
    #plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
                                               # box.color = "black",fill="lightgray",size=25, lineheight = 2),
         axis.title = element_text(face = "bold", size = 17.5),
         axis.text = element_text(size = 25.5),
         legend.position = "right",
    legend.title = element_text(face = "bold", size = 25.5),
    legend.text = element_text(size = 20.5),
    strip.text = element_text(face="bold",size = 20),
    legend.key.size = unit(3, 'cm'), #change legend key size
    legend.key.height = unit(3, 'cm'), #change legend key height
    legend.key.width = unit(3, 'cm') #change legend key width
  ) +
  labs(x = "Treatment", y = "Alpha diversity measure", colour = "Stage", fill = "Stage")+#,title="Alpha diversity of groups before/after treatment")+
  scale_colour_manual(values = c("Before" = "#899DA4", "Post" = "#C93312")) +
  scale_fill_manual(values = c("Before" = "#899DA4", "Post" = "#C93312"))
#graphics.off()


}

################
# Significantly altered taxon 
################

#bracken_species$s__Lactococcus_raffinolactis[-c(which(bracken_species$sample_id %in% c("kefir_post_12_S24", "kefir_post_19_S38","kefir_post_23_S46", "kefir_post_28_S55" )))] <-0 

bracken_species_raffino <- dplyr::select(bracken_species,s__Lactococcus_raffinolactis,sample_id)
bracken_species_raffino <- merge(bracken_species_raffino,metadata,by="sample_id")


  bracken_species_raffino$species <- "Lactococcus raffinolactis"



bracken_species_raffino$stage <- gsub("Pre", "Before",bracken_species_raffino$stage  )

dodge <- position_dodge(width = 0.5)
bracken_species_raffino$title <- "Kefir"
pc <- ggplot(bracken_species_raffino[which(bracken_species_raffino$treatment=="Kefir"),], aes(species, as.numeric(s__Lactococcus_raffinolactis), fill=stage))+
         geom_violin(position = dodge)+
  geom_boxplot(colour="white",width=.035, position = dodge)+
  facet_wrap(~title)+
  theme_bw()+
  theme( plot.title = element_text(face = "bold", size = 25,hjust = 0.5),
    #plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
                                                #box.color = "black",fill="lightgray",size=25, lineheight = 2),
         axis.title = element_text(face = "bold", size = 17.5),
         axis.text.x = element_text(face = "italic",size = 25.5),
         axis.text.y = element_text(size = 25.5),
         legend.position = "right",
         legend.title = element_text(face = "bold", size = 25.5),
         legend.text = element_text(size = 20.5),
    strip.text = element_text(face="bold",size = 20),
    legend.key.size = unit(2, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm') #change legend key width
  )+
  labs(x="Taxon",y="Relative abundance (%)", colour = "Stage", fill = "Stage",title = "Significantly altered taxon")+ 
  scale_colour_manual(values = c("Before" = "#899DA4", "Post" = "#C93312")) +
  scale_fill_manual(values = c("Before" = "#899DA4", "Post" = "#C93312"))



################
# Plot all three plots together 
################
library(ggpubr)
pb <- ggarrange( plot_list[["shannon"]],  plot_list[["simpson"]]+ rremove("ylab"),common.legend = TRUE,legend="right")
  
pb <- annotate_figure(pb,
top = text_grob("Alpha diversity of groups before/after treatment", face = "bold", size = 25,hjust = 0.5))

Figure <- 
ggarrange(p,
          
          ggarrange(pb,pc,nrow=1, labels=c("B.", "C.")),ncol=1,nrow=2, labels="A.")


png(filename="E:/STORE N GO/R/Plots/BBC/Figure 1_v2.png", width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
Figure
graphics.off()

