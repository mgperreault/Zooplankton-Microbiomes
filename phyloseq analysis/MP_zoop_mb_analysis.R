###### investigate the final 16S library - learning from Chris's "plots in phyloseq.R" script

getwd()

#loading packages

library("tidyverse")
library("ggplot2")
library("phyloseq")
library("vegan")
library("cowplot")
library("gridExtra")
library("ggpubr")
library("dplyr")


remotes::install_github("microbiome/microbiome")
library("microbiome")

devtools::install_github("microbiome/mia")
library("mia")


# ps.prune is the final phyloseq object from pipeline


############### Reading in data - ASV count table, metadata, taxa table
# clear workspace 
rm(list = ls())



# importing ASV count table
ASV_table <- read.csv("otu_table.csv")
str(ASV_table)

#importing taxa table
taxa_table <- as.matrix(read.csv("tax_table.csv"))
str(taxa_table)

#importing metadata
metadata <- read.csv("sam_data.csv")

# set row names
row.names(metadata)<-metadata$X 

# remove junk column
metadata <- metadata[-1]


# In the "Organism" column, change "NA" to "Water"
metadata$Organism <- metadata$Organism %>% replace_na("Water")

# Make a data frame with a column for the basin type
metadata$Basin <- ifelse(metadata$Lake == "Cooney" & "Blue" & "Virginia", "North",
                         ifelse(metadata$Lake == "Eastern Brook" & "Serene" & "South", "Middle"))

metadata$Tree_Line <- ifelse(metadata$Elevation.ft < 10000, "Below", "Above")


#factors
# format metadata (have to redo this if loading back in)
make.fac<-c("Original_ID", "Sequencing_ID", "Sample_Type", "Time.Point", "Organism", "Lake", "Plate", "Plate_name", "Well", "sample_control")
metadata[make.fac]<-lapply(metadata[make.fac], factor) # make all these factors



#################### assemble the phyloseq object
PS.fin<-
  read_phyloseq(
    otu.file = "otu_table.csv",
    taxonomy.file = "tax_table.csv",
    metadata.file = "sam_data.csv",
    type = c("simple"),
    sep = ","
  )


#replace sample data with reformatted metadata
metadata_no_rarefy <- sample_data(PS.fin)
str(metadata_no_rarefy)

# check that levels exist
levels(get_variable(PS.fin, "Lake"))
levels(get_variable(PS.fin, "Organism"))
###



##################### Now on to plotting

###### some richness plots
#richness by Lake
richness.plot <- plot_richness(PS.fin, x="Lake", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "figures/richness.plot.location.pdf", height=4, width=10)
dev.off() 

#richness by organism (or sample type i.e, water)
richness.plot<-plot_richness(PS.fin, x="Organism", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "figures/richness.plot.organism.pdf", height=4, width=12)
dev.off() 

#richness by Time.point
richness.plot<-plot_richness(PS.fin, x="Time.Point", measures=c("Observed", "Shannon")) + theme_bw()
richness.plot
dev.copy(pdf, "figures/richness.plot.Time.Point.pdf", height=4, width=10)
dev.off() 
###### 




############ Read counts and rarefaction curves

# Make a data frame with a column for the read counts of each sample
sample_sum_df <- as.data.frame(sample_sums(PS.fin))
colnames(sample_sum_df) <- "read.sum"
sample_sum_df$Sequencing_ID <- rownames(sample_sum_df)

# merge in the reads
run.metaD <- merge(metadata, sample_sum_df, by="Sequencing_ID", all.y=TRUE)
write.csv(run.metaD, "output/run.metaD.csv")


# Histogram of sample read counts
hist.depth<-ggplot(sample_sum_df, aes(x = read.sum)) + 
  geom_histogram(color = "black", fill = "gold2", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) + theme_classic() + geom_vline(xintercept=1000, lty=2)
hist.depth
dev.copy(pdf, "figures/hist.depth.pdf", height=4, width=5)
dev.off() 


# zoom in a bit on the mode of the histogram
zoom.hist.depth<-ggplot(sample_sum_df, aes(x = read.sum)) + 
  geom_histogram(color = "black", fill = "gold2", binwidth = 1000) + coord_cartesian(xlim=c(0, 5000)) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank()) + theme_classic() + geom_vline(xintercept=1000, lty=2)
zoom.hist.depth
dev.copy(pdf, "figures/zoom.hist.depth.pdf", height=4, width=5)
dev.off() 

# a different view
read_depth_violin <- ggplot(sample_sum_df, aes(x=1, y=read.sum)) +
  geom_violin() + ggtitle("Distribution of sample sequencing depth") + scale_y_log10()
read_depth_violin
dev.copy(pdf, "figures/read_depth_violin.pdf", height=4, width=5)
dev.off()



# Calculating Good's coverage
summarize(read.sum = sum(read.sum),
            n_sings = sum(value == 1),
            goods = 1 - (n_sings / read.sum))



#### more plots
pdf(file="figures/read.by.species.pdf", height=4, width=10)
boxplot(run.metaD$read.sum~run.metaD$Organism)
dev.off() 

pdf(file="figures/read.by.sample.pdf", height=4, width=5)
boxplot(run.metaD$read.sum~run.metaD$sample_control)
dev.off() 

pdf(file="figures/log.reads.sample.pdf", height=4, width=7)
ggplot(run.metaD, aes(x=sample_control, y=log(read.sum), color=Organism)) + geom_boxplot()
dev.off() 

pdf(file="figures/reads.sample.pdf", height=4, width=7)
ggplot(run.metaD, aes(x=sample_control, y=read.sum, color=Organism)) + geom_boxplot()
dev.off() 
##### 


############# Rarefaction


# first plotting a rarefying curve
count_table_filt <- as.data.frame(otu_table(PS.fin))
rarecurve((count_table_filt), step=50, cex=0.5, xlim=c(0,100000), ylab = "ASVs", label=FALSE)
abline(v = 5000, lty = "dotted", col="red", lwd=2)

dev.copy(pdf, "figures/rare.raw.pdf", height=4, width=5)
dev.off() 

# there's some very low sequencing depths here, we need to get rid of them 

# looking at the number of reads in a new dataframe (df2) that mirrors the first df we made with our phyloseq object 

df2 <- as.data.frame(sample_data(PS.fin)) # Put sample_data into a ggplot-friendly data.frame
df2$LibrarySize <- sample_sums(PS.fin) # this is the number of reads corresponding to each sample 
df2 <- df2[order(df2$LibrarySize),] # ordering the df by the number of reads from lowest to highest - quite a lot of variation here 
df2$Index <- seq(nrow(df2))

# looking at LibrarySize before and after getting rid of NAs 
df2$LibrarySize
sort(rowSums(otu_table(t(PS.fin))))


# 500 read sequencing depth
physeq.rare500 <- rarefy_even_depth(PS.fin, rngseed=1e4, sample.size = 500)
sort(rowSums(otu_table(t(physeq.rare500))))
# 44 samples removed because they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 03_16S107_16S110_16S111_16S112_16S
# 3902OTUs were removed because they are no longer present in any sample after random sub-sampling



 # 750 read sequencing depth
physeq.rare750 <- rarefy_even_depth(PS.fin, rngseed=1e4, sample.size = 750)
sort(rowSums(otu_table(t(physeq.rare750))))
# 64 samples removed because they contained fewer reads than `sample.size`.
# Up to first five removed samples are: 102_16S103_16S107_16S110_16S111_16S
# 3531OTUs were removed because they are no longer present in any sample after random sub-sampling



# 1000 sequencing depth
physeq.rare1k <- rarefy_even_depth(PS.fin, rngseed=1e4, sample.size = 1e3)
sort(rowSums(otu_table(t(physeq.rare1k))))
# 78 samples removed because they contained fewer reads than `sample.size`.
# Up to first five removed samples are:08_16S102_16S103_16S107_16S110_16S
# 3313OTUs were removed because they are no longer present in any sample after random sub-sampling 

table(sample_data(physeq.rare500)[, "Lake"], exclude = NULL)


############ PCoA
sample_variables(PS.fin)

# make colors for sites
library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


###
T1<- subset_samples(PS.fin, Time.Point=="1")
ORD.T1 <- ordinate(T1, method='MDS', distance='bray')

NMDS.ord.T1<-plot_ordination(
  physeq = T1,                                                   
  ordination = ORD.T1) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time1") +
  theme_classic()   

NMDS.ord.T1
dev.copy(pdf, "figures/NMDS.ord.T1.pdf", height=6, width=7)
dev.off() 

###


###
T2<- subset_samples(PS.fin, Time.Point=="2")
ORD.T2 <- ordinate(T2, method='MDS', distance='bray')

NMDS.ord.T2<-plot_ordination(
  physeq = T2,                                                   
  ordination = ORD.T2) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time2") +
  theme_classic()     

NMDS.ord.T2
dev.copy(pdf, "figures/NMDS.ord.T2.pdf", height=6, width=7)
dev.off() 

####
T3<- subset_samples(PS.fin, Time.Point=="3")
ORD.T3 <- ordinate(T3, method='MDS', distance='bray')

NMDS.ord.T3<-plot_ordination(
  physeq = T3,                                                   
  ordination = ORD.T3) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_shape_manual(values=c(19,17,15,3,7,8,2,0)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time3") +
  theme_classic()     

NMDS.ord.T3
dev.copy(pdf, "figures/NMDS.ord.T3.pdf", height=6, width=7)
dev.off() 


###
T4<- subset_samples(PS.fin, Time.Point=="4")
ORD.T4 <- ordinate(T4, method='MDS', distance='bray')

NMDS.ord.T4<-plot_ordination(
  physeq = T4,                                                   
  ordination = ORD.T4) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time4") +
  theme_classic()     

NMDS.ord.T4
dev.copy(pdf, "figures/NMDS.ord.T4.pdf", height=6, width=7)
dev.off() 


###
T5<- subset_samples(PS.fin, Time.Point=="5")
ORD.T5 <- ordinate(T5, method='MDS', distance='bray')

NMDS.ord.T5<-plot_ordination(
  physeq = T5,                                                   
  ordination = ORD.T5) +                                                
  geom_point(aes(color = Lake, shape=Organism), size = 3) +    
  stat_ellipse(level=0.9, linetype = 2, aes(color=Lake)) +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") + ggtitle("Time5") +
  theme_classic()     

NMDS.ord.T5
dev.copy(pdf, "figures/NMDS.ord.T5.pdf", height=6, width=7)
dev.off() 


#### testing other NMDS plots



str(PS.fin)

asv_rel_abund <- inner_join(run.metaD, otu_table, by="Sequencing_ID") %>%
  inner_join(., taxa_table, by = "ASV") %>%
  group_by(Sequencing_ID) %>%
  mutate(rel_abund = count / sum(count)) %>%
  ungroup() %>%
  select(-count) %>%
  pivot_longer(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
               names_to="level",
               values_to="taxon")
  
  
  asv_rel_abund %>%
    filter(level=="Phylum") %>%
    group_by(Sequencing_ID, taxon) %>%
    summarize(rel_abund = sum(rel_abund))
  
  # Throwing out rare species in ASV table
  
  hist(log10(colSums(ASV_table[,-1])))
  tots<-colSums(ASV_table)
  subset(tots, tots>10)  
  
  
  
  ##### New metadata in physeq object #####
  
 rm(list = ls())
  
  
metadata <- read.csv('metadata_mp.csv')
ASV_table <- read.csv('otu_table.csv')
tax_table <- read.csv('tax_table.csv')  


# modify metadata
str(metadata)
table(metadata$organism)

# set row names
row.names(metadata)<-metadata$X 

# remove junk column
metadata <- metadata[-1]

# format metadata (have to redo this if loading back in)
make.fac <- as.factor(c("sample_ID", "sequencing_ID", "sample_Type", "organism_water", "time_point", "lake"))



# make new column for cladocerans, copepoda, water


metadata$phy_group <- ifelse(metadata$organism == "Calanoid" | metadata$organism == "Cyclopoid" | metadata$organism =="Large Calanoid", "copepoda", ifelse(metadata$organism == "Water", "water", "cladocera"))

levels(as.factor(metadata$phy_group))

metadata$phy_group <- as.factor(metadata$phy_group)


# making new physeq object
PS.fin<-
  read_phyloseq(
    otu.file = "otu_table.csv",
    taxonomy.file = "tax_table.csv",
    metadata.file = "metadata_mp.csv",
    type = c("simple"),
    sep = ","
  ) 
str(PS.fin)



# change sample names to match 
sample_names(PS.fin) <- metadata$sequencing_ID
rownames(metadata) <- metadata$sequencing_ID
  
sample_data(PS.fin) <- metadata

# now we have final physeq object with correct metadata 


# Make a data frame with a column for the read counts of each sample

sample_sum_df <- data.frame(read.sum = sample_sums(PS.fin))

# make row names the sample names

colnames(sample_sum_df) <- "read.sum"
sample_sum_df$sequencing_ID <- rownames(sample_sum_df)

# merge in the reads
run.metaD <- merge(metadata, sample_sum_df, by="sequencing_ID", all.y=TRUE)
write.csv(run.metaD, "output/run.metaD.final.csv")



  
  
  
  
  
  
