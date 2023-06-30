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
my.files <- c()


temp=list.files("E:/Excel/Metagenomic data/BBC/Humann3/")
names(my.files ) <- temp

for (file in temp[c(5)]){
  data <- read_delim(paste("E:/Excel/Metagenomic data/BBC/Humann3/",file,sep=""), 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
  data <- data[,-c(which(colnames(data)=="inulin_post_4_S1__Abundance-CPM"))]
  
  
  colnames(data)[1] <- "functional_feature"
  
  positions_genus <- grep("g__", getElement(data,"functional_feature"))
  positions_unclassified <- grep("unclassified", getElement(data,"functional_feature"))#
  
  
  
  
  # total_humann3_pathabundance<- read_delim("E:/Excel/Metagenomic data/BBC/Humann3/total_humann3_pathabundance-cpm.tsv", 
  #                                               delim = "\t", escape_double = FALSE, 
  #                                               trim_ws = TRUE)
  # 
  # total_humann3_pathabundance <- total_humann3_pathabundance[,-c(which(colnames(total_humann3_pathabundance)=="inulin_post_4_S1__Abundance-CPM"))]
  
  
  # positions_genus <- grep("g__", total_humann3_pathabundance$`# Pathway`)
  # positions_unclassified <- grep("unclassified", total_humann3_pathabundance$`# Pathway`)
  
  positions_to_remove <- c(positions_genus,positions_unclassified)
  
  data <- data[-c(positions_to_remove),]
  if (file!="total_humann3_genefamilies_go-cpm.tsv"){
    data$functional_feature<- sub('.*:', '', data$functional_feature)}
  
  
  
  data$functional_feature<-trimws(data$functional_feature)
  rows <-data$functional_feature 
  data <- dplyr::select(data,-c(functional_feature))
  #data <- data[-c(1:2),]
  data<- t(data)
  
  colnames(data) <- rows
  
  #colnames(data) <- data[1,]
  #rownames_data <- row.names(data)
  #rownames_data <- rownames_data[-c(1)]
  #data <- data[-c(1),]

  data<- as.data.frame(data)
  
  #data <- data %>% mutate_if(is.character,as.numeric)
  #functions_name <- cbind(functions_name,column_names)
  #rownames(functions_name) <- rownames_data
  
  
  rownames(data) <-gsub("_S.*","",rownames(data))
  my.files[[file]] <- data



########################################################################################################################
# tidy input
########################################################################################################################


metadata <-  data.frame(sample_id=rownames(data)) %>%
  separate(sample_id, c("treatment", "stage", "V4"), "_", remove = F) %>%
  #dplyr::select("treatment", "stage", "V4") %>%
  mutate(treatment = recode(treatment, "fmp" = "FMP", "inulin" = "Inulin", "kefir" = "Kefir")) %>%
  mutate(stage = recode(stage, "pre" = "Pre", "post" = "Post")) %>%
  mutate(stage = factor(stage, levels = c("Pre", "Post")))

  maaslin2_data <- data
  
  maaslin2_metadata <- metadata %>%
    column_to_rownames("sample_id")
  

# 
# 
# ########################################################################################################################
# # investiagted gaba shunt 
# ########################################################################################################################
# 
# 
# gaba <- data %>% rownames_to_column("sample_id")
# 
# gaba <- gaba[,c(1,grep("GABA",colnames(gaba)))]
# 
# gaba <- merge(gaba, metadata, by = "sample_id")
# 
# 
# 
# ggplot(gaba, aes(x = treatment, y = `GABA shunt`)) +
#   geom_boxplot(aes(fill=stage), outlier.shape = NA) +
#   geom_point(
#     aes(x = treatment, colour = stage, fill = stage),
#     color="black",
#     size=2,
#     pch = 21, 
#     alpha = 0.5, 
#     position = position_jitterdodge(jitter.width = .25)
#   ) +
#   #geom_text(data = alpha_ann_total, aes(x = treatment, y =  y[1], label = p),size=7) +
#   #facet_wrap(~title)+
#   theme_bw() +
#   theme( plot.title = element_text(face = "bold", size = 25,hjust = 0.5),
#          #plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
#          # box.color = "black",fill="lightgray",size=25, lineheight = 2),
#          axis.title = element_text(face = "bold", size = 17.5),
#          axis.text = element_text(size = 12.5),
#          legend.position = "right",
#          legend.title = element_text(face = "bold", size = 17.5),
#          legend.text = element_text(size = 12.5),
#          strip.text = element_text(face="bold",size = 20),
#          legend.key.size = unit(1.5, 'cm'), #change legend key size
#          legend.key.height = unit(1.5, 'cm'), #change legend key height
#          legend.key.width = unit(1.5, 'cm') #change legend key width
#   ) +
#   labs(x = "Treatment", y = "Copies Per Million", colour = "Stage", fill = "Stage")+#,title="Alpha diversity of groups before/after treatment")+
#   scale_colour_manual(values = c("Before" = "#899DA4", "Post" = "#C93312")) +
#   scale_fill_manual(values = c("Before" = "#899DA4", "Post" = "#C93312"))
# 
# 
# gaba %>%
#   #filter(treatment %in% products) %>%
#   group_by(treatment) %>%
#   rstatix::wilcox_test(`GABA shunt` ~ stage) %>%
#   ungroup() %>%
#   dplyr::select(treatment, p)
# 
# 
# 
# 
# #graphics.off()
# 
# 
# alpha_div_raffino$colonisation_raffino[which(alpha_div_raffino$sample_id %in% c("kefir_post_12_S24", "kefir_post_19_S38","kefir_post_23_S46", "kefir_post_28_S55"))] = "Present"
# 
# 
# gaba$raff <- "Absent"
# 
#  gaba_reduced <- gaba[which(gaba$sample_id  %in% c("kefir_post_12", "kefir_pre_12", "kefir_post_19","kefir_pre_19","kefir_post_23","kefir_pre_23", "kefir_post_28", "kefir_pre_28")),]
#  
#  ggplot(gaba_reduced, aes(x = treatment, y = `GABA shunt`)) +
#    geom_boxplot(aes(fill=raff), outlier.shape = NA) +
#    geom_point(
#      aes(x = treatment, colour = raff, fill = raff),
#      color="black",
#      size=2,
#      pch = 21, 
#      alpha = 0.5, 
#      position = position_jitterdodge(jitter.width = .25)
#    ) +
#    #geom_text(data = alpha_ann_total, aes(x = treatment, y =  y[1], label = p),size=7) +
#    #facet_wrap(~title)+
#    theme_bw() +
#    theme( plot.title = element_text(face = "bold", size = 25,hjust = 0.5),
#           #plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
#           # box.color = "black",fill="lightgray",size=25, lineheight = 2),
#           axis.title = element_text(face = "bold", size = 17.5),
#           axis.text = element_text(size = 12.5),
#           legend.position = "right",
#           legend.title = element_text(face = "bold", size = 17.5),
#           legend.text = element_text(size = 12.5),
#           strip.text = element_text(face="bold",size = 20),
#           legend.key.size = unit(1.5, 'cm'), #change legend key size
#           legend.key.height = unit(1.5, 'cm'), #change legend key height
#           legend.key.width = unit(1.5, 'cm') #change legend key width
#    ) +
#    labs(x = "Treatment", y = "Copies Per Million", colour = "Stage", fill = "Stage")+#,title="Alpha diversity of groups before/after treatment")+
#    scale_colour_manual(values = c("Before" = "#899DA4", "Post" = "#C93312")) +
#    scale_fill_manual(values = c("Before" = "#899DA4", "Post" = "#C93312"))
# ########################################################################################################################
# # lc.raffinolactis vs non lc. raffinolactis 
# ########################################################################################################################
# 
#  my.files[[1]] <-  my.files[[1]] %>%  rownames_to_column("sample_id")
#  
#  gaba_EC <- my.files[[1]][,c(1,grep("2.6.1.19",colnames(my.files[[1]])))]
#  
#  gaba_EC <- merge(gaba_EC, metadata, by = "sample_id")
#  
# 
#  
#  gaba_EC_reduced <- gaba_EC[which(gaba_EC$sample_id  %in% c("kefir_post_12", "kefir_pre_12", "kefir_post_19","kefir_pre_19","kefir_post_23","kefir_pre_23", "kefir_post_28", "kefir_pre_28")),]
#  
# 
#  ggplot(gaba_EC_reduced, aes(x = treatment, y = `2.6.1.19`)) +
#    geom_boxplot(aes(fill=stage), outlier.shape = NA) +
#    geom_point(
#      aes(x = treatment, colour = stage, fill = stage),
#      color="black",
#      size=2,
#      pch = 21, 
#      alpha = 0.5, 
#      position = position_jitterdodge(jitter.width = .25)
#    ) 
#  
#  
#  gaba_EC_reduced %>%
#    #filter(treatment %in% products) %>%
#    group_by(treatment) %>%
#    rstatix::wilcox_test(`2.6.1.19` ~ stage) %>%
#    ungroup() %>%
#    dplyr::select(treatment, p)
 
 
 
# PCoA #
########

treatments <- levels(as.factor(metadata$treatment))

plot_list <- c()

permanova_tab1 <-c()

  
  permanova_tab1 <-c()
  
subset_meta <- metadata %>%
  #filter(treatment == i) %>%
  column_to_rownames("sample_id")

sample_ids <- rownames(subset_meta)

subset_data <- data 

species_dist <- vegdist(as.matrix(subset_data), method="bray")

species_pcoa <- pcoa(species_dist)

points <- species_pcoa$vectors[, 1:2] %>%
  as.data.frame() %>%
  rownames_to_column("sample_id") %>%
  merge(., metadata, by = "sample_id")

pcoa_values <- data.frame(relative_eig = species_pcoa$values$Relative_eig) %>%
  mutate(axis = paste0("Axis.", row_number())) %>%
  mutate(relative_eig = 100 * relative_eig) %>%
  select(axis, relative_eig)

axis_1 <- paste0("Axis 1 [", round(pcoa_values[1, 2], 2), "%]")
axis_2 <- paste0("Axis 2 [", round(pcoa_values[2, 2], 2), "%]")

xmin <- 0.95 * min(points$Axis.1)
ymax <- 0.95 * max(points$Axis.2)

for (category in treatments){
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
  
  aov_tab <- adonis(subset_data1 ~ stage, data =   subset_meta1)$aov.tab %>%
    as.data.frame() %>%
    drop_na() %>%
    select(`R2`, `Pr(>F)`)
  
  rval <- paste0("R^2 = ", format(100 * round(aov_tab[1, 1], 3), nsmall = 1), "%")
  pval <- paste0("p = ", format(round(aov_tab[1, 2], 3), nsmall = 3))
  
  permanova_tab <- data.frame(Axis.1 = xmin, Axis.2 = ymax, label = paste0(rval, "<br>", pval), name=category)
  
  permanova_tab1 <- rbind(permanova_tab,permanova_tab1)
}

chrom_order <- c("FMP","Inulin", "Kefir")

permanova_tab1 <- 
  permanova_tab1[order(match(permanova_tab1$name, chrom_order)),]

points$stage <- gsub("Pre", "Before",points$stage )


p <- ggplot(points, aes(x = Axis.1, y = Axis.2 )) +
  geom_point(aes(fill=stage),color="black",pch = 21, size = 8) +
  
  
  facet_wrap(~ treatment) +
  stat_ellipse(aes(color=stage))+
  # stat_ellipse(geom = "polygon",
  #              alpha = 0,
  #              type = "norm"
  #              )+
  theme_bw()+
  labs( x = axis_1, y = axis_2,title = "Bray-Curtis PCoA plot based on HUMAnN3 output")+
  theme(#plot.title = ggtext::element_textbox_simple(halign  = 0.5,linetype = 1, # turn on border
    #box.color = "black",size=25, lineheight = 2),
    plot.title = element_text(face = "bold", size = 25,hjust = 0.5),
    axis.title = element_text(face = "bold", size = 12.5),
    axis.text = element_text(size = 17.5),
    legend.title = element_text(face = "bold", size = 12.5),
    legend.position = "right",
    strip.text = element_text(face="bold",size = 20),
    legend.key.size = unit(2, 'cm'), #change legend key size
    legend.key.height = unit(2, 'cm'), #change legend key height
    legend.key.width = unit(2, 'cm') #change legend key width
  ) +
  scale_fill_manual(name="Stage",values = c("Before" = "#899DA4", "Post" = "#C93312"))+
  scale_colour_manual(name="Stage",values = c("Before" = "#899DA4", "Post" = "#C93312"))

# Add text labels r  

f_labels <- data.frame(treatment =  permanova_tab1$name, label = permanova_tab1$label)
p <- 
  p+
  ggtext::geom_richtext(x =   permanova_tab1$Axis.1[1], y = permanova_tab1$Axis.2[1]*1.45 , aes(label = label),
                        
                        data =  f_labels,
                        #aes(x = Axis.1, y = Axis.2, label = label),
                        fill = NA,
                        hjust = 0,
                        label.colour = NA,
                        size=10,
                        label.padding = grid::unit(rep(0, 4), "pt")
  )


png(filename=paste("E:/STORE N GO/R/Plots/BBC/Functional_output/",file,".png",sep=""), width = 7864, height=5200,res =300,pointsize = 15) #, width=2000, height=1950)
plot(p)
graphics.off()
}



########################################################################################################################
# run MaAsLin2 #
########################################################################################################################
treatments <- levels(as.factor(metadata$treatment))


'%!in%' <- function(x,y)!('%in%'(x,y))
for (i in treatments){
  
  
  sub_meta <- metadata %>%
    filter(treatment == i) %>%
    column_to_rownames("sample_id")
  
  sub_data <- maaslin2_data %>%
    filter(rownames(.) %in% rownames(sub_meta))
  #colnames(sub_data) <- gsub("\\(",".",colnames(sub_data))
  #colnames(sub_data) <- gsub("\\)",".",colnames(sub_data))
  unlink(paste0("E:/Excel/Metagenomic data/BBC/Maaslin2/Humann3/",file,"/"),recursive = TRUE,force = TRUE)
  dir.create(paste0("E:/Excel/Metagenomic data/BBC/Maaslin2/Humann3/",file,"/", i),recursive = TRUE)
  out_dir <-paste0("E:/Excel/Metagenomic data/BBC/Maaslin2/Humann3/",file,"/", i)
  
  Maaslin2(
    sub_data,
    sub_meta,
    out_dir,
    fixed_effects = c('stage'),
    random_effects = c('group')
  )
  # 
  # all_results <- fread(paste0(out_dir, "/all_results.tsv"))
  # 
  # 
  # 
  # res <- all_results %>%
  #   filter(qval < 0.25)
  # 
  # if (nrow(res) == 0) {
  #   
  #   res <- all_results %>%
  #     filter(pval < 0.25)
  #   
  # }
  # 
  # 
  # spp <- res %>%
  #   .$feature
  # 
  # #spp <- gsub("[[:punct:]]","", spp)
  # 
  # 
  # 
  # # spp <- gsub("sp..","sp.", spp)
  # # spp <- gsub("sp...","sp.", spp)
  # # all_results$feature<- gsub("sp..","sp.",all_results$feature)
  # # all_results$feature<- gsub("sp...","sp.",all_results$feature)
  # 
  # 
  # for (sp in spp) {
  #   
  #   
  #   stats <- all_results %>%
  #     filter(feature == sp) %>%
  #     mutate(feature=gsub("\\.","",feature)) %>% 
  #     dplyr::select(feature, coef, qval)
  #   # if(sp=="s__Calothrix_sp__336.3"){
  #   #  sp <- "s__Calothrix_sp__336/3"
  #   # }
  #   
  #   
  #   mpa <-  sub_data %>%
  #     rownames_to_column("sample_id") %>%
  #     merge(., metadata, by = "sample_id") %>%
  #     rename_with(~ gsub("-",".",.x)) %>%
  #     dplyr::select(sample_id, stage, all_of(sp)) %>%
  #     dplyr::rename("abundance"=3 )
  #   
  #   
  #   coef <- round(stats[1, "coef"], 3)
  #   qval <- round(stats[1, "qval"], 3)
  #   subtitle <- paste0("qval: ", qval, "; ", "coeE: ", coef)
  #   
  #   p <- ggplot(mpa, aes(x = stage, y = abundance)) +
  #     geom_boxplot(aes(colour = stage), outlier.shape = NA) +
  #     geom_point(
  #       aes(x = stage, colour = stage, fill = stage),
  #       pch = 21, 
  #       alpha = 0.5, 
  #       position = position_jitterdodge(jitter.width = 0.25)
  #     ) +
  #     theme_classic() +
  #     theme(
  #       plot.title = element_text(face = "bold", size = 12.5, hjust = 0.5),
  #       plot.subtitle = element_text(size = 12.5, hjust = 0.5),
  #       axis.title = element_text(face = "bold", size = 12.5),
  #       axis.text = element_text(size = 12.5),
  #       legend.position = "none"
  #     ) +
  #     labs(x = "Stage", y = "Shannon index", title = sp, subtitle = subtitle) +
  #     scale_colour_manual(values = c("Pre" = "#899DA4", "Post" = "#C93312")) +
  #     scale_fill_manual(values = c("Pre" = "#899DA4", "Post" = "#C93312"))
  #   
  #   ggsave(plot = p, file = paste0(out_dir, "/", sp, ".png"), height = 4, width = 5, dpi = 300)
  #   
}

#

}
