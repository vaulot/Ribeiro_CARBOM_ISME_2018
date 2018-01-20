library("mapproj")
library("maps")
# library("mapdata")
library("ggmap")

library("treemap") # To create tree maps...

library("ggplot2")
library("gplots") # to convert colors to hex and also for heat plots
# library("heatmaply")   # for interactive heatmaps 

library("scales") # For formatting the scales in graphic
library("gridExtra") # necessary to arrange graphics into grid
library("plotrix") # needed to draw circles without problem with character size in pdf files
library("RColorBrewer") # for color options

library("akima") # To interpolate the grid for contour plots

library("mapplots") # for pie charts on a map
# library("reshape2")  # Needed to reshape the column with the classes
# library("plyr")   # To normalize data by groups...
library("dplyr")  # To filter data
library("tidyr")  # To melt the data using a more refined way
library("stringr") #Ë To process strings

library("vegan") # For ecological computations....

library("igraph") # For networks graphs

setwd("C:/Users/vaulot/Google Drive/Papers/2016 Ribeiro CARBOM ISME/R")

# =======================================================
# ======================  Samples  ======================
# =======================================================

# Read the sample data files
  samples<- read.table("carbom_samples.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Read the environmental data file
  envdata<- read.table("carbom_envdata.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Merge the environmental data into the sample table
  envdata$station=as.character(envdata$station)  # neeeded to do the merge
  samples <- left_join(samples, select(envdata, -longitude, -latitude, -bottom_depth), by= c("station","depth"))

# Read the flow cytometry data file
  fcm <- read.table("carbom_fcm.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)


# Create label for sample
  samples$sample_label<-paste("TR",samples$transect,"_St",samples$station,"_",samples$depth,"m_",samples$sample_illumina, sep="")

  fractions<-c("Pico","Nano")

# =======================================================
#   Create map of stations - Fig. 1 
#-----------------------------------------------------------
  
# --- FIG. 1A - map of region with normal map
pdf("./pdf_figures_2.0/Fig_1_map_global_bw 2.0.pdf",width = 6, height = 6, useDingbats=FALSE)
map(database = "world", xlim = c(-80,-30), ylim = c(-40,10),fill=TRUE, boundary=TRUE, col="grey80")
points(samples$longitude, samples$latitude, pch=19, col="black", cex=0.8)
# text(samples$longitude+0.25, samples$latitude+0.25, labels=samples$station, cex=0.7)
map.axes(las=1) # las=1 means all labels are horizontal
dev.off()

# --- map of stations with Google map
map_1 <- get_map(location = c(-42, -25), zoom = 6,  scale = 2, color = "color", maptype = "satellite")
map_stations <- ggmap(map_1) + 
                geom_point(data=samples, aes(x=longitude, y=latitude), color="white", size=3.5) +
                geom_text(data=samples, aes(x=longitude+0.25,y=latitude+0.25, label=station), 
                          color="white", size=5, hjust=0)
ggsave( plot=map_stations, 
        filename="./pdf_figures_2.0/Fig_1_map_local 2.0.pdf",
        width = 20, height = 20, scale=1.2, units="cm", useDingbats=FALSE)

# --- FIG. 1A - map of region with Google map
map_2 <- get_map(location = c(-45, -20), zoom = 4,  scale = 2, color = "color",  maptype = "roadmap")
map_wide <- ggmap(map_2) + 
            geom_point(data=samples, aes(x=longitude, y=latitude), color="black", size=0.5)
map_wide
ggsave( plot=map_wide, 
        filename="./pdf_figures_2.0/Fig_1_map_wide terrain 2.0.pdf",
        width = 20, height = 20, scale=1.2, units="cm", useDingbats=FALSE)


# =======================================================
#   Do contour plots of fcm 
#-----------------------------------------------------------
  fcm_tr <- filter(fcm, (transect == 2)) 
  ggplot(data = fcm_tr, aes(x = longitude, y = -depth, z = picoeukaryotes)) +
  geom_tile(aes(fill = picoeukaryotes)) +
  stat_contour(bins=10) +
  ggtitle("Picoeuks") +
  xlab("Longitude") +
  ylab("Depth") +
  scale_fill_continuous(name = "picoeukaryotes",low = "white", high = "blue") +
  theme(plot.title = element_text(size = 25, face = "bold"),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title.x = element_text(size = 20, vjust = -0.5),
        axis.title.y = element_text(size = 20, vjust = 0.2),
        legend.text = element_text(size = 10))

  fld <- with(fcm_tr, interp(x = longitude, y = -depth, z = picoeukaryotes, duplicate = "mean", nx= 1000, ny=1000))
  filled.contour(x = fld$x, y = fld$y,z= fld$z, color.palette =  colorRampPalette(c("blue", "white", "red")),
               xlab = "Longitude",
               ylab = "Depth",
               main = "Picoeuks",
               key.title = title(main = "picoeuks", cex.main = 1))
# =======================================================
#   Reformat the table for the paper - Table 1  
#-----------------------------------------------------------

sample_table <- samples %>% filter (Select_18S_nifH=="Yes") %>% 
                            group_by (station_order, transect, station, depth, latitude, longitude) %>% 
                            summarise(fraction_number = n())
sample_table_pico <- samples %>%  filter ((Select_18S_nifH=="Yes") & (fraction == "Pico")) %>% 
                                  transmute (station_order, 
                                             total_18S_pico = total_18S, 
                                             total_nifH_pico = total_nifH, 
                                             sample_pico = sample_illumina )
sample_table_nano <- samples %>%  filter ((Select_18S_nifH=="Yes") & (fraction == "Nano")) %>% 
                                  transmute (station_order, 
                                             total_18S_nano = total_18S, 
                                             total_nifH_nano = total_nifH, 
                                             sample_nano = sample_illumina ) 
sample_table <- left_join(sample_table, sample_table_pico)
sample_table <- left_join(sample_table, sample_table_nano)




# =======================================================
# ==========================  Colors  ======================
# =======================================================

# colorRampPalette is in the RColorBrewer package.  This creates a colour palette that shades from light yellow to red in RGB space with 100 unique colours
  scaleyellowred <- colorRampPalette(c("lightyellow", "red"), interpolate="spline", space = "rgb")(30)

# Read the 18S classes files
  classes<- read.table("classes_18S.txt", header=TRUE, sep="\t", 
                       na.strings="NA", dec=".", strip.white=TRUE, comment.char = "")
  #  classes$color_hex<-col2hex(classes$color)
  classes_ordered <- rev(as.character(classes$class))
# This named vector is necessary for the ggplot2 below
  class_color<-structure(as.character(classes$color_hex),.Names=as.character(classes$class))

# Read the nifH color files
  classes_nifH<- read.table("classes_nifH.txt", header=TRUE, sep="\t", 
                            na.strings="NA", dec=".", strip.white=TRUE, comment.char = "")
  # classes_nifH$color_hex<-col2hex(classes_nifH$color)
# This named vector is necessary for the ggplot2 below
  class_nifH_color<-structure(as.character(classes_nifH$color),.Names=as.character(classes_nifH$class)) 
  taxo_nifH_ordered <- rev(as.character(classes_nifH$class))
  # classes_nifH$color_otu_hex<-col2hex(classes_nifH$color_otu)
# This named vector is necessary for the ggplot2 below
  class_nifH_color_otu<-structure(as.character(classes_nifH$color_otu_hex),.Names=as.character(classes_nifH$class)) 

  
# =====================================================================================================================================================================
#                                  ==========================  18S  ======================
# =====================================================================================================================================================================

# Read the otu table
  otu_table<- read.table("carbom_18S_otu_0.02.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

# Merge with classes for colors
  otu_table <-merge(otu_table, subset(classes, select = class:color_hex), by.x="Final_Taxo4", by.y="class") 

# Create label for otu
  otu_table$otu_label <- paste(otu_table$otu,"_",otu_table$Final_Taxo7, sep="")

# Remove unwanted OTU from Opitokontha and Sptrepto
  otu_total_total <- sum(otu_table$otu_total)
  otu_table <- filter(otu_table, Final_Taxo2 != "Opisthokonta")				   # Remove  otus from Opisthokonta
  otu_table <- filter(otu_table, Final_Taxo3 != "Streptophyta")  				 # Remove  otus from Streptophyta

# =======================================================
# 18S Get some basic statistics
#-----------------------------------------------------------

  otu_total_total_clean <- sum(otu_table$otu_total)
  pct_otu_contam <- (otu_total_total - otu_total_total_clean)/otu_total_total  # Compute fraction in Opistonkons and Sptreptos

  otu_stats_1 <- otu_table %>% group_by(trophic_mode) %>%summarise(n(), sum(otu_total)/otu_total_total_clean)
  otu_stats_2 <- otu_table %>% filter(trophic_mode!= "auto") %>% group_by(Final_Taxo4) %>%summarise(n(), sum(otu_total)/otu_total_total_clean)

# =======================================================
# 18S Melt the otu table and filter to keep the autotrophs
#-----------------------------------------------------------

# Melt the otu table (use the new tidyr library)
  otu_melt <- otu_table %>% gather(sample_id, sequence_number, X10n:tri03p, factor_key = TRUE)  

# Merge otu table and sample
  otu_melt <-merge(otu_melt, samples, by.x="sample_id", by.y="sample_id") 
  

# Only keep selected stations
  otu_melt <- filter(otu_melt, Select_18S_nifH == "Yes")					 # Keep selected stations

# Write melted table
  write.table(otu_melt, file="carbom_18S_otu_0.02.otu_melt.txt", sep="\t", row.names = FALSE) 

# Only keep the autotrophic groups
  otu_melt_auto <- filter(otu_melt, trophic_mode=="auto")	
  total_sequence_corrected <- otu_melt_auto %>% 
                              group_by (sample_illumina) %>% 
                              summarize (total_18S_auto=sum(sequence_number)) %>%
                              select(sample_illumina, total_18S_auto)
  otu_melt_auto <-merge(otu_melt_auto, total_sequence_corrected, by.x="sample_illumina", by.y="sample_illumina") 

# =======================================================
# 18S Compute matrixes of abundances for OTUs for autotrophs
#-----------------------------------------------------------

# Reshape the abundance matrix in condensed form 
  # Keep only necessary columns
  otu_matrix <- otu_melt_auto %>% transmute(sample_label, otu_label, sequence_number, Final_Taxo2, Final_Taxo3, Final_Taxo4, Final_Taxo7, Final_Taxo8)
  # Reshape in wide form, sort by taxonomy and then remove the taxonomy information
  otu_matrix <- otu_matrix %>% 
                spread(sample_label, sequence_number) %>% 
                arrange(Final_Taxo2, Final_Taxo3, Final_Taxo4, otu_label ) %>% 
                select (-Final_Taxo2, -Final_Taxo3, -Final_Taxo4, -Final_Taxo7, -Final_Taxo8)
  
# Provide row.names using the labels for the samples  
  row.names(otu_matrix)=otu_matrix$otu_label
  
# Remove the column that contains the labels
  otu_matrix=select(otu_matrix, -otu_label)  

# Transpose the matrix
  otu_matrix=t(otu_matrix)  

    
# Normalize the data per sample
  otu_matrix_relative <- otu_matrix/rowSums(otu_matrix)
  
# determine the maximum relative abundance for each column
  maxab <- apply(otu_matrix_relative, 2, max)  # 1 indicates row and 2 columns
  
# remove the otus with less than 20 % as their maximum relative abundance -> for heatmaps
  otu_low <- names(which(maxab < 0.20))
  otu_abund <- names(which(maxab >= 0.20))
  otu_matrix_relative_abund <- otu_matrix_relative[, -which(colnames(otu_matrix_relative) %in% otu_low)]
  otu_matrix_abund <- otu_matrix[, -which(colnames(otu_matrix_relative) %in% otu_low)]

# Do a small table with the information from major OTUs
  otu_table_paper <- otu_table %>% filter(otu_label %in% otu_abund) %>% transmute(otu, otu_total, Division=Final_Taxo3, Class=Final_Taxo4, Genus=Final_Taxo7, Species=Final_Taxo8, BLAST_best_Accession, BLAST_best_identity, BLAST_best_Description) %>% arrange(otu)
  write.table(otu_table_paper, file="otu_18S_table_paper.txt", sep = "\t", row.names = FALSE, quote=FALSE)

# Number of OTU per class
  otu_table_class <- otu_table %>% group_by(Final_Taxo3, Final_Taxo4) %>% summarize(number_of_otu = n(), number_of_reads = sum(otu_total)) 
  
# =======================================================
# 18S -  Fig. 2A - Treemap - based on Class - Using the Mean (and not the sum of all sequences)
#-----------------------------------------------------------

# First aggregate the different otus by Class for each sample
  otu_melt_agg_class<-otu_melt %>%
                      group_by(sample_illumina,fraction, level, total_18S,Final_Taxo4, color_hex) %>%
                      summarize (sequence_number=sum(sequence_number) )
# Compute mean for each class, each fraction
  otu_melt_agg_mean<-ddply(otu_melt_agg_class, 
                           c("fraction","Final_Taxo4", "color_hex"), 
                           summarize, sequence_pct=mean(sequence_number/total_18S) )
# Compute mean for each class
  otu_melt_agg_mean_2_fractions<-ddply(otu_melt_agg_class, 
                                       c("Final_Taxo4", "color_hex"), 
                                       summarize, sequence_pct=mean(sequence_number/total_18S) )

# Do the tree maps separate for pico and nano
treemap(otu_melt_agg_mean[(otu_melt_agg_mean$fraction=="Nano"),], 
        index=c("Final_Taxo4"),vSize="sequence_pct",vColor="color_hex",type="color",
        title="18S Nano mean of sequences",asp=1, lowerbound.cex.labels= 1, fontsize.labels = 6)
treemap(otu_melt_agg_mean[(otu_melt_agg_mean$fraction=="Pico"),], 
        index=c("Final_Taxo4"),vSize="sequence_pct",vColor="color_hex",type="color",
        title="18S Pico mean of sequences",asp=1)

# FIG 2A - Tree map for pico + nano
pdf("./pdf_figures_2.0/Fig_2A Treemap 18S classes 2.0.pdf",width = 6, height = 6, useDingbats=FALSE)
treemap(otu_melt_agg_mean_2_fractions, 
        index=c("Final_Taxo4"),vSize="sequence_pct",vColor="color_hex",type="color",
        title="18S mean of sequences",asp=1, lowerbound.cex.labels= 1, fontsize.labels = 8)
dev.off()

pdf("./pdf_figures_2.0/Fig_2A Treemap 18S classes details 2.0.pdf",width = 6, height = 6, useDingbats=FALSE)
treemap(otu_melt_agg_mean_2_fractions, 
        index=c("Final_Taxo4"),vSize="sequence_pct",vColor="color_hex",type="color",
        title="18S mean of sequences",asp=1, lowerbound.cex.labels= 0, fontsize.labels = 10)
dev.off()


# The next treemap put side by the Pico and Nano....
treemap(otu_melt_agg_mean, index=c("fraction", "Final_Taxo4"),
        vSize="sequence_pct",vColor="color_hex",type="color",
        title="18S Pico mean of sequences",asp=1)

# =======================================================
# 18S - Treemap - based on otu - Using the Mean (and not the sum of all sequences)
#-----------------------------------------------------------

# First compute mean per fraction
  otu_melt_agg_mean <-ddply(otu_melt, c("fraction", "otu_label", "color_hex"), summarize, sequence_pct=mean(sequence_number/total_18S) )

# Compute summary for all fractions
  otu_melt_agg_mean_no_fraction <-ddply(otu_melt, c("otu_label", "color_hex", "Final_Taxo4"), summarize, sequence_pct=mean(sequence_number/total_18S) )

# Do the tree maps per fraction
  treemap(otu_melt_agg_mean[(otu_melt_agg_mean$fraction=="Nano"),], 
          index=c("otu_label"),vSize="sequence_pct",
          vColor="color_hex",type="color",
          title="18S Nano mean of sequences",
          asp=1, lowerbound.cex.labels= 1, fontsize.labels = 6)
  treemap(otu_melt_agg_mean[(otu_melt_agg_mean$fraction=="Pico"),], 
          index=c("otu_label"),vSize="sequence_pct",vColor="color_hex",type="color",
          title="18S Pico mean of sequences",asp=1)
  
# Do the tree maps for both fractions
  treemap(otu_melt_agg_mean_no_fraction, 
          index=c("otu_label"),
          vSize="sequence_pct",
          vColor="color_hex",
          type="color",
          title="18S mean of sequences",
          asp=1)
  
# =======================================================
# 18S - FIG 2B - BarPlot of major otus - Using the Mean (and not the sum of all sequences)
#-----------------------------------------------------------
    size_factor=10
    otu_to_plot <- otu_melt_agg_mean_no_fraction %>% 
                    filter(sequence_pct>0.005) %>% 
                    arrange(-sequence_pct) %>% 
                    mutate (otu_label=str_replace_all(otu_label, "_", " ")) %>%
                    mutate (otu_label=str_replace_all(otu_label, " X", ""))
    plot1<-ggplot(otu_to_plot,  
                  aes(y = sequence_pct, 
                      x = reorder(otu_label,sequence_pct) , 
                      fill=Final_Taxo4) ) +
          geom_bar(stat = "identity") +  coord_flip()+  
          theme_bw(size_factor) + 
          ggtitle("OTUs") + 
          theme(axis.text.y=element_text(angle = 0, hjust = 0)) + 
                xlab("")+ylab("Mean percent of metabarcodes")+  
          scale_y_continuous(labels = percent_format(), breaks=seq(0,0.14,0.02)) +
          scale_fill_manual(values = class_color)
    plot1
    ggsave( plot=plot1, 
          filename="./pdf_figures_2.0/Fig_2B_otu_abnundance 2.0.pdf",
          width = 20, height = 15, scale=1, units="cm", useDingbats=FALSE)  

  
# =======================================================
# 18S - Bar Plots of all 
#-----------------------------------------------------------

  size_factor=10
  otu_melt<-otu_melt[order(otu_melt$transect, otu_melt$station, otu_melt$depth),]
  
  
  for (one_fraction in fractions) {
    otu_melt_fraction <- filter(otu_melt, fraction==one_fraction)	
#   write.table(otu_melt, paste(one_fraction,"otu_melt.tsv"), sep="\t")
    plot1<-ggplot(otu_melt_fraction,  
                  aes(x = reorder(sample_label,-station_order), 
                      y = sequence_number/total_18S , 
                      fill=factor(Final_Taxo4, levels = classes_ordered)) ) +
          geom_bar(stat = "identity") +  coord_flip()+
          theme_bw(size_factor) + 
          ggtitle(paste ("CARBOM 18S - All -", one_fraction)) +
          theme(axis.text.x=element_text(angle = 0, hjust = 1)) + 
                xlab("")+ylab("Percent of metabarcodes")+  
          scale_y_continuous(labels = percent_format()) +
          scale_fill_manual(values = class_color, name="Taxonomy", 
                            guide=guide_legend(reverse = TRUE, nrow = 6, ncol=7, 
                                               keywidth=0.7, keyheight=0.7, 
                                               label.hjust = 0, label.vjust = 0, 
                                               label.position = "right"))+
          theme(legend.position="top", 
                legend.text=element_text(size=size_factor),
                legend.title=element_text(size=size_factor))
    ggsave(plot=plot1, filename=paste("Carbom_18S_all",one_fraction," 1.1.pdf",sep=""),
           width = 20, height = 20, scale=1.2, units="cm")
}

# =======================================================
# 18S - For autotrophs determine the classes which do not contribute at least at 20% somewhere
#-----------------------------------------------------------  
  
# First aggregate the different otus by Class for each sample
  class_melt_auto<-ddply(otu_melt_auto, c("sample_label","Final_Taxo4"), summarize, sequence_number=sum(sequence_number) )

# Reshape the abundance matrix in condensed form
  class_matrix <- dcast(class_melt_auto, sample_label ~ Final_Taxo4 , value.var="sequence_number")  

# Provide row.names using the labels for the samples  
  row.names(class_matrix)=class_matrix$sample_label
  
# Remove the column that contains the labels
  class_matrix=select(class_matrix, -sample_label)  
  
# Normalize the data per sample
  class_matrix_relative <- class_matrix/rowSums(class_matrix)
  
# determine the maximum relative abundance for each column
  maxab_class <- apply(class_matrix_relative, 2, max)  # 1 indicates row and 2 columns
  
# determine the classes with less than 20 % as their maximum relative abundance 
  class_low <- names(which(maxab_class < 0.20))
  
# Low abundance classes are changed to "Others"
  otu_melt_auto$Final_Taxo3 <- as.character(otu_melt_auto$Final_Taxo3)  # Necessary to add new factors !
  otu_melt_auto$Final_Taxo4 <- as.character(otu_melt_auto$Final_Taxo4)
  otu_melt_auto$Final_Taxo3[otu_melt_auto$Final_Taxo4 %in% class_low]<- "Others_auto"
  otu_melt_auto$Final_Taxo4[otu_melt_auto$Final_Taxo4 %in% class_low]<- "Others_auto"
  otu_melt_auto$Final_Taxo3 <- as.factor(otu_melt_auto$Final_Taxo3)     # Necessary to add new factors !
  otu_melt_auto$Final_Taxo4 <- as.factor(otu_melt_auto$Final_Taxo4)

  
# =======================================================
# 18S - Fig. S1 - Bar Plots of autotrophs
#-----------------------------------------------------------

  otu_melt_auto_agg_class <- ddply(otu_melt_auto, 
                                   c("sample_illumina","sample_label", "station_order", 
                                     "fraction", "level", "total_18S_auto","Final_Taxo3","Final_Taxo4"), 
                                   summarize, sequence_number=sum(sequence_number) )
  
  for (one_fraction in fractions) {
    otu_melt_fraction <- filter(otu_melt_auto_agg_class, fraction==one_fraction)	
	otu_melt_fraction <- arrange (otu_melt_fraction, desc(station_order), Final_Taxo3, Final_Taxo4)
#   write.table(otu_melt, paste(one_fraction,"otu_melt.tsv"), sep="\t")
    plot1<-ggplot(otu_melt_fraction,  
                  aes(x = reorder(sample_label,-station_order), y = sequence_number/total_18S_auto , 
                      fill=factor(Final_Taxo4, levels  = classes_ordered)) )+
      geom_bar(stat = "identity") +  coord_flip()+
      theme_bw(size_factor) + 
      ggtitle(paste ("CARBOM 18S - Autotrophs -", one_fraction)) +
      theme(axis.text.x=element_text(angle = 0, hjust = 0.5, vjust=0.5)) +  # to align left the y legend (vjust = 0)
      theme(axis.text.y=element_text(angle = 0, hjust = 0, vjust=0.5)) +  # to align left the y legend (vjust = 0)
      xlab("")+
      ylab("Percent of metabarcodes")+  
      scale_y_continuous(labels = percent_format()) +
      scale_fill_manual(values = class_color, name="Taxonomy", 
                        guide=guide_legend(reverse=TRUE, nrow = 4, ncol=4, 
                                           keywidth=0.7, keyheight=0.7, 
                                           label.hjust = 0, label.vjust = 0, label.position = "right"))+
      theme(legend.position="top", legend.text=element_text(size=size_factor),legend.title=element_text(size=size_factor))
    ggsave( plot=plot1, filename=paste("./pdf_figures_2.0/Fig. S1 - Carbom_18S_autotrophs_",one_fraction," 2.0.pdf",sep=""), 
            width = 20, height = 20, scale=0.8, units="cm")
}


# =======================================================
# 18S - Fig. 3 - Pie charts for transect 3 for Pico and Nano
#-----------------------------------------------------------  
  xlim=c(50,450)
  ylim=c(-120,20)
  pico_nano_factor <- 2.5
  nano_divider <- 5
  pico_divider <- nano_divider * pico_nano_factor
  
for (one_fraction in fractions) {
  otu_melt_fraction <- filter(otu_melt_auto, (fraction==one_fraction) & (transect==2))
  
# make xyz structure that is used to draw piecharts
  xyz <- make.xyz(otu_melt_fraction$transect_distance,-otu_melt_fraction$depth,otu_melt_fraction$sequence_number,otu_melt_fraction$Final_Taxo4)
  
# set the color in the same order than in the xyz array
  pie_color <- left_join(data.frame(class=colnames(xyz$z)), classes)
  
# Next 2 lines if for pie size proportional to cell concentration measure by flow cytometry 
  pie_size <- left_join(data.frame(transect_distance=xyz$x, depth=-xyz$y), filter(samples, (fraction==one_fraction) & (transect==2)) )
  if (one_fraction== "Pico") pie_size$pie_size<- sqrt(pie_size$picoeuks/pico_divider) else pie_size$pie_size<- sqrt(pie_size$nanoeuks/nano_divider)
  
# Fixed pie size 
  pie_size$pie_size <- 12
  
  pdf(paste("./pdf_figures_2.0/Fig_3 Pies 18S ", one_fraction," fixed size 2.0.pdf", spe=""),width = 10, height = 5, useDingbats=FALSE)
  plot(xyz$x, xyz$y, xlim=xlim, ylim=ylim , xlab = "Distance from coast (km)", ylab="Depth (m)", main=one_fraction, frame.plot = TRUE, asp=1.5)
  draw.pie(xyz$x, xyz$y, xyz$z, radius = pie_size$pie_size, col=as.character(pie_color$color_hex))
  # legend.pie(50,-120,labels=colnames(xyz$z), radius=20, bty="n", col=as.character(pie_color$color_hex), cex=0.6, label.dist=1.1)
  dev.off()
}

# =======================================================
#   Compute absolute abundance of two major OTUs per mL 
#-----------------------------------------------------------


otu_mL <-as.data.frame(otu_matrix_relative_abund)
otu_mL$sample_label=row.names(otu_matrix_relative_abund)
otu_mL <- transmute (otu_mL,sample_label,Otu001_Ostreococcus, Otu002_UCYN_A1_host  )
otu_mL <- left_join(samples, otu_mL)

otu_mL_pico <- otu_mL %>% filter ((Select_18S_nifH=="Yes") & (fraction == "Pico")) %>% transmute (station_order, otu001_pico=Otu001_Ostreococcus*picoeuks, otu002_pico = Otu002_UCYN_A1_host*picoeuks )
otu_mL_nano <- otu_mL %>% filter ((Select_18S_nifH=="Yes") & (fraction == "Nano")) %>% transmute (station_order, otu001_nano=Otu001_Ostreococcus*nanoeuks, otu002_nano = Otu002_UCYN_A1_host*nanoeuks )

otu_mL <- otu_mL %>% filter (Select_18S_nifH=="Yes") %>% group_by (station_order, transect, station, depth, latitude, longitude, level, transect_distance) %>% summarise(fraction_number = n())
 
otu_mL<- left_join(otu_mL, otu_mL_pico)
otu_mL<- left_join(otu_mL, otu_mL_nano)
otu_mL[is.na(otu_mL)] <- 0
otu_mL <- mutate(otu_mL, otu001_mL = otu001_pico + otu001_nano, otu002_mL = otu002_pico + otu002_nano)  

  
  
# =======================================================
# 18S  - Similarity Matrix for samples  : Compute it for different metrix 
#-----------------------------------------------------------
  
# Compute dissimil between samples using Bray-Curtis on non normalized matrix (can lead to problems - see paper in Nature Microbio reviews)
  sample_dissimil_bray=vegdist(otu_matrix_abund,method = "bray", na.rm = F)
  
# Compute dissimil between otus using Bray-Curtis
  otu_dissimil_bray=vegdist(t(otu_matrix_abund),method = "bray", na.rm = F)
  
# Change from dissimil to simil and transform to matrix
  otu_simil_bray <- as.matrix(1 - otu_dissimil_bray)

# =======================================================
# 18S - Compute Shannon diversity based on otu 
#-----------------------------------------------------------

  otu_diversity<-as.data.frame(vegan::diversity(otu_matrix, "shannon"))

# =======================================================
#  18S - Fig. 4 - Heat maps based on tutorial  : #  http://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/  
#-----------------------------------------------------------

  
# Do average linkage hierarchical clustering of the rows. Other options are 'aver', 'complete' or 'single'
  sample_cluster=hclust(sample_dissimil_bray, method = "complete")
  plot(sample_cluster, cex=1)

# Do average linkage hierarchical clustering of the otus. Other options are 'complete' or 'single'
  otu_cluster <- hclust(otu_dissimil_bray, "complete")
  plot(otu_cluster, cex=1)

# Do the heatmap with the samples and otu clustered
  heatmap(as.matrix(otu_matrix_relative_abund), Rowv = as.dendrogram(sample_cluster), Colv = as.dendrogram(otu_cluster), col = scaleyellowred, margins = c(12, 4))
  
  
# Do the heatmap with the samples only clustered
  heatmap(as.matrix(otu_matrix_relative_abund), Rowv = as.dendrogram(sample_cluster), Colv = NA, col = scaleyellowred, margins = c(20, 4))


# Fig. 4 - Do the heatmap with the samples only clustered Better to provide the scale
  pdf("./pdf_figures_2.0/Fig_4 Heatmap 18S  2.0.pdf",width = 8, height = 10, useDingbats=FALSE)
  heatmap.2(as.matrix(otu_matrix_relative_abund), 
            Rowv = as.dendrogram(sample_cluster), Colv = FALSE, 
            dendrogram="row", col = scaleyellowred, trace="none", 
            key = TRUE, density.info="none", cexRow=0.8, cexCol=0.8, margins = c(15, 10), 
            lmat = rbind(c(0,3),c(2,1),c(0,4)), lwid = c(1,4), lhei = c(0.2,4,0.6), srtCol=45, 
            key.title= NA , key.xlab = "Relative abundance")
  dev.off()
  
# =======================================================
#  18S - MSDS based using vegan on Bray distance of relative abundance, 
#--------------------------------------------------------
# Compute MSDS
  otu_matrix_relative_abund.mds <- metaMDS(otu_matrix_relative_abund, trace = FALSE) 
  otu_matrix_relative_abund.mds
  plot(otu_matrix_relative_abund.mds, type = "t", display="species")

# Create environmental table to put on msds
  envdata_table <- select(filter(samples,Select_18S_nifH == "Yes"), sample_label, depth, transect_distance, phosphates, silicates, nitrates, temperature, fluorescence, salinity, fraction)
  row.names(envdata_table) <- envdata_table$sample_label
  envdata_table <- select(envdata_table, -sample_label)

# Fit environmental variables including Pico and Nano
  otu_matrix_relative_abund.ef <- envfit(otu_matrix_relative_abund.mds, envdata_table, permu = 999, na.rm = TRUE)
  otu_matrix_relative_abund.ef
  plot(otu_matrix_relative_abund.ef, p.max = 1)

# =====================================================================================================================================================================
#                                  ==========================  nifH  ======================
# =====================================================================================================================================================================
  
# Threshold for total sequence number.  Below this number, nifH is assumed to be absent  
  nifH_sequence_threshold <- 2000
  
# Read the nifH otu file
  otu_nifH<- read.table("carbom_nifH_otu_0.02.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
# Create label for otu
  otu_nifH$otu_label <- paste(otu_nifH$otu,"_nifH_",otu_nifH$Taxo2, sep="")
# Merge with classes for colors
  otu_nifH <-merge(otu_nifH, select(classes_nifH,class:color_otu_hex), by.x="Taxo", by.y="class") 
# Only keep cyanobacteria
  otu_nifH <- otu_nifH %>% filter (Taxo1=="Cyanobacteria")
  
# Melt the otu table (use the new tidyr library)
  otu_nifH_melt <- otu_nifH %>% gather(sample_id, sequence_number, X10n:tri03p, factor_key = TRUE)  

# Merge otu table and sample
  otu_nifH_melt <-merge(otu_nifH_melt, samples, by.x="sample_id", by.y="sample_id") 
  
# Recompute number of sequences using only cyanobacteria  and merge
  total_nifH_cyano <- otu_nifH_melt %>% group_by(sample_id) %>% summarize(total_nifH_cyano=sum(sequence_number))
  otu_nifH_melt <- left_join(otu_nifH_melt, total_nifH_cyano)
  
# If less than 2000 sequences (10% of the mean) set to 0
  otu_nifH_melt$sequence_number[ otu_nifH_melt$total_nifH_cyano < nifH_sequence_threshold ] <- 0 
  
# Keep only selected samples
  otu_nifH_melt <- filter(otu_nifH_melt, Select_18S_nifH == "Yes")					 # Remove  unselected stations

# Write melted table
  write.table(otu_nifH_melt, file="carbom_nifH_otu_0.02.otu_melt.txt", sep="\t", row.names = FALSE)

# =======================================================
# nifH - Compute matrixes of abundances for OTUs
#-----------------------------------------------------------

# Keep only necessary columns
  otu_nifH_matrix <- otu_nifH_melt %>% transmute(sample_label, otu_label, sequence_number, Taxo)
# Reshape in wide form, sort by taxonomy and then remove the taxonomy information
  otu_nifH_matrix <- otu_nifH_matrix %>% 
                spread(sample_label, sequence_number) %>% 
                arrange(Taxo, otu_label ) %>% 
                select (-Taxo)


# Provide row.names using the labels for the samples  
  row.names(otu_nifH_matrix)=otu_nifH_matrix$otu_label
  
# Remove the column that contains the labels
  otu_nifH_matrix=select(otu_nifH_matrix, -otu_label)  
  
# Transpose the matrix
  otu_nifH_matrix=t(otu_nifH_matrix)   

# Normalize the data per sample
  otu_nifH_matrix_relative <- as.matrix(otu_nifH_matrix/rowSums(otu_nifH_matrix))
  
# Since some samples have no nifH, it induces division by 0 and therefore naN which needs to be removed  
  otu_nifH_matrix_relative[is.nan(otu_nifH_matrix_relative)] <- 0
  
# determine the maximum relative abundance for each column
  maxab_nifH <- apply(otu_nifH_matrix_relative, 2, max)  # 1 indicates row and 2 columns
  sample_max_nifH <- apply(otu_nifH_matrix_relative, 1, max)  # 1 indicates row and 2 columns
  length(sample_max_nifH[sample_max_nifH >= 0.90])/length(sample_max_nifH)  # Just to compute
  
# remove the otus with less than 20 % as their maximum relative abundance
  otu_nifH_low <- names(which(maxab_nifH < 0.20))
  otu_nifH_abund <- names(which(maxab_nifH >= 0.20))
  otu_nifH_matrix_relative_abund <- otu_nifH_matrix_relative[, -which(colnames(otu_nifH_matrix_relative) %in% otu_nifH_low)]
  otu_nifH_matrix_abund <- otu_nifH_matrix[, -which(colnames(otu_nifH_matrix_relative) %in% otu_nifH_low)]
  
# Do a small table with the information from major OTUs
  otu_nifH_table_paper <- otu_nifH %>% filter(otu_label %in% otu_nifH_abund) %>% transmute(otu, otu_total, Class=Taxo1, Genus=Taxo2, BLAST_best_Accession, BLAST_best_identity, BLAST_best_Description) %>% arrange(otu)
  write.table(otu_nifH_table_paper, file="otu_nifH_table_paper.txt", sep = "\t", row.names = FALSE, quote=FALSE)

  
  
# =======================================================
# nifH - Treemap - based on otu - Separate Pico and Nano
#-----------------------------------------------------------
# First compute mean per fraction and removing the samples with no nifH
  otu_nifH_melt_agg_mean <-otu_nifH_melt %>% filter(total_nifH_cyano >= nifH_sequence_threshold )%>% 
                group_by (fraction, otu_label, color_otu_hex)%>% 
                summarise( sequence_pct=mean(sequence_number/total_nifH_cyano) )

# Do the tree maps
  treemap(otu_nifH_melt_agg_mean[(otu_nifH_melt_agg_mean$fraction=="Nano"),], 
          index=c("otu_label"),vSize="sequence_pct",vColor="color_otu_hex",type="color",
          title="nifH Nano mean of sequences",asp=1)
  treemap(otu_nifH_melt_agg_mean[(otu_nifH_melt_agg_mean$fraction=="Pico"),], 
          index=c("otu_label"),vSize="sequence_pct",vColor="color_otu_hex",type="color",
          title="nifH Pico mean of sequences",asp=1)

# =======================================================
# nifH - Fig. 5A - Treemap - based on otu - Pico and Nano together 
#-----------------------------------------------------------
# First compute mean per fraction and removing the samples with no nifH
  otu_nifH_melt_agg_mean <-otu_nifH_melt %>% filter(total_nifH_cyano >= nifH_sequence_threshold )%>% 
                group_by (otu_label, color_otu_hex)%>% 
                summarise( sequence_pct=mean(sequence_number/total_nifH_cyano) )

# Do the tree maps
  pdf("./pdf_figures_2.0/Fig_5 Treemap nifH 2.0.pdf",width = 6, height = 6, useDingbats=FALSE)
  treemap(otu_nifH_melt_agg_mean, 
          index=c("otu_label"),vSize="sequence_pct",vColor="color_otu_hex",type="color",
          title="nifH mean of sequences",asp=1)
  dev.off()
# =======================================================
# nifH - Fig. S1 - Bar Plots of Taxo 
#-----------------------------------------------------------
# Low abundance otus are changed to "Others"
  otu_nifH_melt$Taxo <- as.character(otu_nifH_melt$Taxo)  # Necessary to add new factors !
# otu_nifH_melt$Taxo[otu_nifH_melt$otu_label %in% otu_nifH_low] <- "Others"
  otu_nifH_melt$Taxo <- as.factor(otu_nifH_melt$Taxo)

  size_factor=10
  otu_nifH_melt<-arrange(otu_nifH_melt, transect, station, depth)

# Aggregate by taxonomy  
  otu_nifH_melt_agg_taxo <- ddply(otu_nifH_melt, c("sample_illumina","sample_label", "station_order", "fraction", "level", "total_nifH_cyano","Taxo"), summarize, sequence_number=sum(sequence_number) )
  
  for (one_fraction in fractions) {
    otu_nifH_melt_fraction <- filter(otu_nifH_melt_agg_taxo, fraction==one_fraction)	
    plot1<-ggplot(otu_nifH_melt_fraction,  
                  aes(x = reorder(sample_label,-station_order), 
                      y = sequence_number/total_nifH_cyano , 
                      fill=factor(Taxo, levels=taxo_nifH_ordered)) )+
      geom_bar(stat = "identity") +  coord_flip()+
      theme_bw(size_factor) + ggtitle(paste ("CARBOM nifH - ", one_fraction)) +
      theme(axis.text.x=element_text(angle = 0, hjust = 0.5, vjust=0.5)) +  
      theme(axis.text.y=element_text(angle = 0, hjust = 0, vjust=0.5)) +  # to align left the y legend (hjust = 0)
      xlab("") + ylab("Percent of metabarcodes")+  
      scale_y_continuous(labels = percent_format()) +
      scale_fill_manual(values = class_nifH_color_otu, name="Taxo", 
                        guide=guide_legend(reverse = TRUE, nrow = 4, ncol = 3, 
                                           keywidth=0.7, keyheight=0.7, 
                                           label.hjust = 0, label.vjust = 0, 
                                           label.position = "right")) +
      theme(legend.position="top", 
            legend.text=element_text(size=size_factor),
            legend.title=element_text(size=size_factor))
      ggsave(plot=plot1, filename=paste("Fig S1 nifH_",one_fraction," 2.0.pdf",sep=""),
            width = 20, height = 20, scale=0.8, units="cm")
}



# =======================================================
# nifH - Similarity Matrix for samples  : Compute it for different metrix 
#-----------------------------------------------------------
# Remove the rows with no nifH
  otu_nifH_matrix_relative_abund <- otu_nifH_matrix_relative_abund[rowSums(otu_nifH_matrix_relative_abund)!=0, ]
  
# Compute dissimil between samples using Bray-Curtis
  sample_nifH_dissimil=vegdist(otu_nifH_matrix_relative_abund,method = "bray", na.rm = F)
# Replace NaN by 0 for samples without nifH
  sample_nifH_dissimil[is.nan(sample_nifH_dissimil)] <- 0
# Compute dissimil between otus using Bray-Curtis - need to transpose the matrix
  otu_nifH_dissimil=vegdist(t(otu_nifH_matrix_relative_abund),method = "bray", na.rm = F)
  
  
# =======================================================
#  nifH -Heat maps based on tutorial  : #  http://www.molecularecologist.com/2013/08/making-heatmaps-with-r-for-microbiome-analysis/  
#-----------------------------------------------------------

# Do average linkage hierarchical clustering of the rows. Other options are 'aver', 'complete' or 'single'
  sample_nifH_cluster=hclust(sample_nifH_dissimil, method = "complete")
# Do average linkage hierarchical clustering of the otus. Other options are 'complete' or 'single'
  otu_nifH_cluster <- hclust(otu_nifH_dissimil, "complete")
	
# Do the heatmap with the otu clustered
  heatmap(as.matrix(otu_nifH_matrix_relative_abund), Rowv = as.dendrogram(sample_nifH_cluster), Colv = as.dendrogram(otu_nifH_cluster), col = scaleyellowred, margins = c(12, 4))
# Do the heatmap with the samples only clustered
  heatmap(as.matrix(otu_nifH_matrix_relative_abund), Rowv = as.dendrogram(sample_nifH_cluster), Colv = NA, col = scaleyellowred, margins = c(15, 6), cexRow=0.5, cexCol=0.7, keep.dendro=FALSE)

#  pdf("./pdf_figures_2.0/Fig_5B Heatmap nifH  2.0.pdf",width = 8, height = 10, useDingbats=FALSE)
  heatmap.2(as.matrix(otu_nifH_matrix_relative_abund), 
            Rowv = as.dendrogram(sample_nifH_cluster), Colv = FALSE, dendrogram="row", 
            col = scaleyellowred, trace="none", key = TRUE, density.info="none", 
            cexRow=1, cexCol=1, margins = c(20, 10), 
            lmat = rbind(c(0,3),c(2,1),c(0,4)), lwid = c(1,4), lhei = c(0.2,4,0.6), srtCol=45, 
            key.title= NA , key.xlab = "Relative abundance")
#  dev.off()
  
# The next one creates interactive heatmap, nice légend, but does not work because cannot save as pdf
#  heatmaply(as.matrix(otu_nifH_matrix_relative_abund),colors = scaleyellowred, dendrogram = "row", row_dend_left=FALSE, margins = c(200, 200, NA, 0))


# =======================================================
#  nifH - MSDS based using vegan on Bray distance of relative abundance, 
#--------------------------------------------------------
# Compute MSDS
  otu_nifH_matrix_relative_abund.mds <- metaMDS(otu_nifH_matrix_relative_abund, trace = FALSE) 
  otu_nifH_matrix_relative_abund.mds
  plot(otu_nifH_matrix_relative_abund.mds, type = "t", display="species")

# Create environmental table to put on msds
  envdata_table <- select(filter(samples,Select_18S_nifH == "Yes"), sample_label, depth, transect_distance, phosphates, silicates, nitrates, temperature, fluorescence, salinity, fraction)
  row.names(envdata_table) <- envdata_table$sample_label
  envdata_table <- select(envdata_table, -sample_label)
  envdata_table_nifH <- envdata_table[rownames(otu_nifH_matrix_relative_abund),]  # Only select rows for which there are some nifH data

# Fit environmental variables including Pico and Nano
  otu_nifH_matrix_relative_abund.ef <- envfit(otu_nifH_matrix_relative_abund, envdata_table_nifH, permu = 999, na.rm = TRUE)
   otu_nifH_matrix_relative_abund.ef
  plot( otu_nifH_matrix_relative_abund.ef, p.max = 1)

# =====================================================================================================================================================================
#                            ===  18S and nifH  - Network analysis OTU level====
# =====================================================================================================================================================================

  
# --- Combine 18S and nifH into a single matrix (use the the abund - < 0.20 ) ---
  
  otu_combined_matrix_abund <- merge (otu_matrix_abund, otu_nifH_matrix_abund, by="row.names")
  row.names(otu_combined_matrix_abund)=otu_combined_matrix_abund$Row.names
  otu_combined_matrix_abund=select(otu_combined_matrix_abund, -Row.names)   
# remove samples for which no nifH (8 samples have no nifH) - NOT TO BE DONE because done before
#  sample_no_nifH <- names(which(sample_max_nifH == 0)) 
#  otu_combined_matrix_abund <- otu_combined_matrix_abund[-which(rownames(otu_combined_matrix_abund) %in% sample_no_nifH), ]
# --- Export combined matrix for use with sparCC
  write.table(t(otu_combined_matrix_abund), file="carbom_18S_nifH_combined_abund.txt", sep="\t", row.names = TRUE, quote=FALSE)

# # --- NOT USED Combine 18S and nifH into a single matrix (use the the abund - < 0.20 ) - For relative abundance and not absolute --- 
#   otu_combined_matrix_relative_abund <- merge (otu_matrix_relative_abund, otu_nifH_matrix_relative_abund, by="row.names")
#   row.names(otu_combined_matrix_relative_abund)=otu_combined_matrix_relative_abund$Row.names
#   otu_combined_matrix_relative_abund=select(otu_combined_matrix_relative_abund, -Row.names)   
# # --- Export combined matrix for use with sparCC
#   write.table(t(otu_combined_matrix_relative_abund), file="carbom_18S_nifH_combined_abund_relative.txt", sep="\t", row.names = TRUE, quote=FALSE)

  
#  --- Compute dissimil between samples using Bray-Curtis
  otu_combined_dissimil_bray=vegdist(t(otu_combined_matrix_abund),method = "bray", na.rm = F)
  otu_combined_simil_bray <- as.matrix(1 - otu_combined_dissimil_bray)
#  --- Compute Spearman rank correlation - does not vary with relative or absolute abundance
  otu_combined_spearman <- cor(otu_combined_matrix_abund, method = "spearman")

# =======================================================
#  Network graphs - Using iGraph - 
#  Note : The matrix contains all smaples even those without Cyanos nifH
#-----------------------------------------------------------

# --- Bray --- 
  
# Create graph based on combined matrix using the bray distance
  otu_graph_bray <- graph.adjacency(otu_combined_simil_bray, weighted=TRUE, mode="upper")
# Remove the links that correspond to simil <0.25
  otu_graph_bray <- delete_edges( otu_graph_bray, E(otu_graph_bray)[weight<0.25])
  
  otu_graph <- otu_graph_bray

# --- Spearman ---   
  
# Create graph based on spearman correlation for otu plus nifH
  otu_graph_spearman <- graph.adjacency(otu_combined_spearman , weighted=TRUE, mode="upper")
# Remove the links that correspond to correl <0.35 (excludes all negative correlations)
  otu_graph_spearman <- delete_edges( otu_graph_spearman, E(otu_graph_spearman)[weight<0.35])
# Remove the links that correspond to self nodes
  otu_graph_spearman <- simplify(otu_graph_spearman)

  otu_graph <- otu_graph_spearman
  
# --- Fig. 6 - SparCC ---
  
# Use Python 2.. and not 3...  DOS Command to copy : 
# "C:\Program Files\Anaconda\python.exe" SparCC.py carbom_18S_nifH_combined_abund.txt -i 10 --cor_file=carbom_18S_nifH_combined_abund_sparcc.txt > sparcc.log
# "C:\Program Files\Anaconda\python.exe" MakeBootstraps.py carbom_18S_nifH_combined_abund.txt -n 100 -p ./bootstrap/ -t carbom_permutated_#.txt
# FOR /L %C IN (0,1,99) DO ("C:\Program Files\Anaconda\python.exe" SparCC.py ./bootstrap/carbom_permutated_%C.txt -i 10 --cor_file=./bootstrap/boot_cor_%C.txt >> sparcc.log)
# "C:\Program Files\Anaconda\python.exe" PseudoPvals.py carbom_18S_nifH_combined_abund_sparcc.txt ./bootstrap/boot_cor_#.txt 100 -o carbom_18S_nifH_combined_abund_sparcc.pvals_two_sided.txt -t two_sided >> sparcc.log

# Read SparCC correlation 
  sparcc <- read.table("carbom_18S_nifH_combined_abund_sparcc.txt",header=TRUE,sep="\t")
  row.names(sparcc) <- sparcc[,1]
  sparcc=as.matrix(sparcc[,2:ncol(sparcc)])
  
# Read SparCC pseudo_pvalues 
  sparcc_pvalue <- read.table("carbom_18S_nifH_combined_abund_sparcc.pvals_two_sided.txt",header=TRUE,sep="\t")
  row.names(sparcc_pvalue) <- sparcc_pvalue[,1]
  sparcc_pvalue=as.matrix(sparcc_pvalue[,2:ncol(sparcc_pvalue)])
  
# Remove the links that correspond to correl <0.20 (excludes all negative correlations) and for which pseudo p-value is > 0.05
  sparcc [sparcc<0.20] <- 0
  sparcc [sparcc_pvalue>0.05] <- 0  # This does not remove any edge... so all edges are supported for p <=0.05

# create the graph
  otu_graph_sparcc=graph.adjacency(sparcc, weighted=TRUE, mode="upper")
  
# Remove the links that correspond to self nodes
  otu_graph_sparcc <- simplify(otu_graph_sparcc)
  otu_graph <- otu_graph_sparcc

  

# --- Draw graph ---  
    
# Graph layout
  otu_graph$layout <- layout_with_fr #The Fruchterman-Reingold layout algorithm
# Make the width of the edges proportionnel to the simil
  E(otu_graph)$width<-50^E(otu_graph)$weight
  Junk <- E(otu_graph)

# Make the vertex size and color related to origin of otu
  otu_graph_vertex <- data.frame(V(otu_graph)$name)
  otu_graph_vertex <- rename(otu_graph_vertex, otu_label=V.otu_graph..name)

# igraph replace "-" and " " by "."  - Not necessary anymore because the names have been fixed
  otu_graph_vertex$otu_label <- as.character(gsub( ".", "_", otu_graph_vertex$otu_label,fixed=TRUE))

# Colors and shape for 18S
  foo <- select(otu_table, otu_label, otu_fraction_of_total, color)
  otu_graph_vertex_18S <-inner_join(otu_graph_vertex, foo) 

  otu_graph_vertex_18S$size <- sqrt(otu_graph_vertex_18S$otu_fraction_of_total)*100
  otu_graph_vertex_18S$shape <- "circle"
  otu_graph_vertex_18S$label.font <- 1 # plain

# Colors and shape for nifH
  foo <- select(otu_nifH, otu_label, otu_fraction_of_total, color)  
  otu_graph_vertex_nifH <-inner_join(otu_graph_vertex,foo ) 
  
  otu_graph_vertex_nifH$size <- sqrt(otu_graph_vertex_nifH$otu_fraction_of_total)*30
  otu_graph_vertex_nifH$shape <- "circle"  # change to rectangle for final plot....
  otu_graph_vertex_nifH$label.font <- 2 # bold
  
  otu_graph_vertex <-merge(otu_graph_vertex_18S, otu_graph_vertex_nifH, all=T, sort=F) 

  V(otu_graph)$color <- otu_graph_vertex$color
  V(otu_graph)$shape <- otu_graph_vertex$shape
  V(otu_graph)$size  <- otu_graph_vertex$size
  V(otu_graph)$label.font <- otu_graph_vertex$label.font

  tkplot(otu_graph) 

  # This file is used by Gephi
  # export Gephi output as pdf and rework in Corel Draw
  write_graph(otu_graph, "otu_graph_sparcc.gml", format = "gml") 
  

# =======================================================
# 18S - Compare composition of samples with and without cyanos
#-----------------------------------------------------------
  
  total_nifH_cyano <- total_nifH_cyano %>% mutate(sample_illumina= str_replace(sample_id, "X",""))
  otu_melt_agg_class <- left_join(otu_melt_agg_class, total_nifH_cyano)
  nifH_sequence_threshold
  
  # Compute mean for each class, each fraction
  otu_melt_agg_mean_with_cyanos <- otu_melt_agg_class %>%
                                   filter(total_nifH_cyano>nifH_sequence_threshold) %>%
                                   group_by(fraction,Final_Taxo4,color_hex ) %>% 
                                   summarize(sequence_pct=mean(sequence_number/total_18S) )
  otu_melt_agg_mean_without_cyanos <- otu_melt_agg_class %>%
                                   filter(total_nifH_cyano<nifH_sequence_threshold) %>%
                                   group_by(fraction,Final_Taxo4,color_hex ) %>% 
                                   summarize(sequence_pct=mean(sequence_number/total_18S) )

# Do the tree maps
treemap(otu_melt_agg_mean_with_cyanos, 
        index=c("Final_Taxo4"),vSize="sequence_pct",vColor="color_hex",type="color",
        title="18S with cyanos mean of sequences",asp=1, lowerbound.cex.labels= 1, fontsize.labels = 6)
treemap(otu_melt_agg_mean_without_cyanos, 
        index=c("Final_Taxo4"),vSize="sequence_pct",vColor="color_hex",type="color",
        title="18S without cyanos mean of sequences",asp=1, lowerbound.cex.labels= 1, fontsize.labels = 6)


  
  
  
  
  
  
  
  
  
    
# =======================================================
# ======  NOT USED  ======================
# =======================================================
	
	
# =======================================================
# FIG : Ven diagram for share OTUs
otu_venn<-ddply(otu_melt, c("otu","fraction","level"), summarize, occurence=sum(sequence_number>0) )
otu_venn_level_fraction <-dcast(otu_venn, otu ~ fraction + level, sum)
otu_venn_level <-dcast(otu_venn, otu ~ level, sum)
otu_venn_fraction <-dcast(otu_venn, otu ~ fraction, sum)

otu_venn_fraction[(otu_venn_fraction$Nano>0)&(otu_venn_fraction$Pico>0),"summary"]<-"Nano+Pico"
otu_venn_fraction[(otu_venn_fraction$Nano>0)&(otu_venn_fraction$Pico==0),"summary"]<-"Nano"
otu_venn_fraction[(otu_venn_fraction$Nano==0)&(otu_venn_fraction$Pico>0),"summary"]<-"Pico"
otu_venn_fraction_plot<-ddply(otu_venn_fraction, c("summary"), summarize, otu_number=length(summary) )

otu_venn_level[(otu_venn_level$Surf>0)&(otu_venn_level$Deep>0)&(otu_venn_level$Tricho>0),"summary"]<-"Surf+Deep+Tricho"
otu_venn_level[(otu_venn_level$Surf>0)&(otu_venn_level$Deep>0)&(otu_venn_level$Tricho==0),"summary"]<-"Surf+Deep"
otu_venn_level[(otu_venn_level$Surf>0)&(otu_venn_level$Deep==0)&(otu_venn_level$Tricho==0),"summary"]<-"Surf"
otu_venn_level[(otu_venn_level$Surf==0)&(otu_venn_level$Deep>0)&(otu_venn_level$Tricho==0),"summary"]<-"Deep"
otu_venn_level[(otu_venn_level$Surf>0)&(otu_venn_level$Deep==0)&(otu_venn_level$Tricho>0),"summary"]<-"Surf+Tricho"
otu_venn_level[(otu_venn_level$Surf==0)&(otu_venn_level$Deep>0)&(otu_venn_level$Tricho>0),"summary"]<-"Deep+Tricho"
otu_venn_level[(otu_venn_level$Surf==0)&(otu_venn_level$Deep==0)&(otu_venn_level$Tricho>0),"summary"]<-"Tricho"
otu_venn_level_plot<-ddply(otu_venn_level, c("summary"), summarize, otu_number=length(summary) )

# =======================================================
# FIG : BLAST results
otu_BLAST<- read.table("carbom_16S_otu_0.02 with BLAST.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
otu_BLAST <- otu_BLAST[otu_BLAST$taxo1!="Bacteria",]																    # Remove  otus from Bacteria
otu_BLAST <- otu_BLAST[otu_BLAST$Final_Taxo4!="Embryophyceae",]															# Remove  otus from higher plants
otu_BLAST_melt <- melt(otu_BLAST, id.vars=c( "otu","taxo1","taxo2","taxo3","Final_Taxo4","taxo5","taxo6","taxo7","taxo8","taxo9","taxo10","otu_total"),measure.vars=c("BLAST_id_pct","BLAST_cult_id_pct"),variable.name = "BLAST_target", value.name = "BLAST_id_pct")

ggplot(otu_BLAST_melt, aes(Final_Taxo4, BLAST_id_pct)) + geom_boxplot(aes(fill = BLAST_target))+ theme(axis.text.x=element_text(angle = 45, hjust = 1)) 
ggplot(otu_BLAST_melt, aes(Final_Taxo4, BLAST_id_pct)) + geom_point(aes(color = BLAST_target), size=6)+ theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggplot(otu_BLAST_melt, aes(otu_total, BLAST_id_pct)) + geom_point(aes(color = Final_Taxo4, shape=BLAST_target), size=6)+ scale_x_log10(limits=c(100,200000))

# =======================================================
# FIG : Phylogenetic Trees

library("ggtree")

# Read the tree
tree<-read.tree(file = "C:/Daniel/Papers/Ribeiro CARBOM PLoS 2015/data 16S NGS/trees/Haptos/Hapto_RefSeqs_BenteTree_1.7 FastTree.newick")

# Export the tips ID and the names
tree.tips<-data.frame(tip=c(1:length(tree.haptos$tip.label)),label=tree.haptos$tip.label)

# Export to Excel and then add the size of the OTUs
tree.tips<- read.table("C:/Daniel/Papers/Ribeiro CARBOM PLoS 2015/data 16S NGS/carbom_16S_otu_0.02_tree_hapto_labels.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

tree.graph<-ggtree(tree.haptos)+ geom_tiplab(size=1, align=FALSE)
tree.inset <- lapply(tree.tips$size_relative, function(circle_size) { 
	if (circle_size>0) {ggplot()+ scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1))+ geom_point(aes(x=0.5, y=0.5), size=30*circle_size, shape=21, fill="red") + theme_inset() }
	else {ggplot()+ scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,1))+ geom_point(aes(x=10, y=10), size=0, shape=21, fill="red") + theme_inset() }
	}) 
names(tree.inset) <- tree.tips$node 
inset(tree.graph, tree.inset, width=.03, height=.1, vjust=0, hjust=-0.02, x = "branch")



bars <- nodebar(tree.tips, cols=4)
inset(tree.graph, bars)
tree.graph

# Tests
ggtree(tree.haptos)+ geom_tiplab(size=3, align=FALSE)
ggtree(tree.haptos, layout="circular")+ geom_tiplab(size=3, aes(angle=angle), align=TRUE)
ggtree(tree.haptos) + geom_text2(aes(subset==isTip, label=node), hjust=-.3)


