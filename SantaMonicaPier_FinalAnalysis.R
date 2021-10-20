save.image("FinalAnalysis.Rdata")
load("FinalAnalysis.Rdata")
setwd("~/SMP_R_working_directory")
#### Load required libraries ####
options(scipen=999)
library(reshape2)
library(ggplot2)
library(vegan)
library(edgeR)
library(UpSetR)
library(decontam)
library(tidyverse)

#### Curated Table: Normalize, Blanks, Taxonomy, Filter, Env ####
# read in the ASV table and the taxonomy table
table = read.table("SMP2018-2019_averaged.txt", header = TRUE)
taxonomy = read.table("taxonomy-2018_2019.txt", header = TRUE)
#Check to see if the lengths match
nrow(table) == nrow(taxonomy)
#slap the taxonomy on to the end of the ASV table
ASV_table = plyr::join(table, taxonomy, by = "ASV.ID", type = "left")
nrow(ASV_table) == nrow(table)
dim(ASV_table)
#ASV_table<-read.csv("Qiime1-format-ASVtable_B_pr2111.csv",header=TRUE)
emptyvector = c()
for (i in 1:nrow(ASV_table)) {
  emptyvector = c(emptyvector, paste("ASV_",i, sep = ""))
}
rownames(ASV_table) = ASV_table$ASV.ID = emptyvector#Here I want to give the rownames the ASV IDs
dim(ASV_table) #Check the dimensions. There are 29643rows and 42 columns
columname= c("ASV.ID","blank1","blank2","blank3","Day1-2018", "Day2-2018", "Day3-2018", "Day4-2018", "Day5-2018", "Day6-2018", "Day7-2018", "Day8-2018", "Day9-2018", "Day10-2018", "Day11-2018", "Day12-2018", "Day13-2018", "Day14-2018", "Day15-2018", "Day1-2019", "Day2-2019", "Day3-2019", "Day4-2019", "Day5-2019", "Day6-2019", "Day7-2019", "Day8-2019", "Day9-2019", "Day10-2019", "Day11-2019", "Day12-2019", "Day13-2019", "Day14-2019", "Day15-2019", "Day16-2019", "Day17-2019", "Day18-2019", "Day19-2019", "Day20-2019", "Day21-2019", "Day22-2019", "taxonomy")
colnames(ASV_table) = columname
# Filter out rare ASVs 
ASV_total<-apply(ASV_table[2:41],1,sum) # Tally the abundances of each ASV (rows) across all samples (columns)
ASV_table.no1 = ASV_table[ ASV_total >100,] # Keep only ASVs with 100 or more reads across all samples
ASVs_lost = nrow(ASV_table)-nrow(ASV_table.no1);ASVs_lost
# 28540 ASVs lost

### Blank Treatment ###
#use the decontam R package to handle blanks
#Requires two pieces of information:
# 1)Numeric count table with samples as rows (no strings!) 2) vector identifying the blank samples
# 1) is satisfied with t(ASV_table[2:41])
# 2) Create the vector. The first four samples are the blanks
vector_for_decontamination = c(rep(TRUE,3),rep(FALSE, 37))
# Use decontam to produce a data.frame of the probability of an ASV being an contaminant sequence
contam_df = decontam::isContaminant(t(ASV_table.no1[2:41]), neg = vector_for_decontamination)
contam_df %>% filter(contaminant == "TRUE") #note use of the filter() function. Its handy
contam_asvs = row.names(contam_df[contam_df$contaminant == "TRUE", ])
ASV_table.no1[rownames(ASV_table.no1) %in% contam_asvs,]
ASV_table_noContam = ASV_table.no1[!row.names(ASV_table.no1) %in% contam_asvs, ]
contam_loss = nrow(ASV_table.no1)-nrow(ASV_table_noContam);contam_loss
# lost 1 contaminant ASV

ASV_table_noContam = ASV_table_noContam %>% select(-blank1,-blank2,-blank3)
#write.csv(ASV_table_noContam, row.names = F, file = "ASV_preTMM.csv")
#head(ASV_table_noContam)
#nrow(ASV_table_noContam)

#### Normalize using edgeR ###
ListDGE = DGEList(ASV_table_noContam[,2:38])
#ListDGE
ListDGE = calcNormFactors(ListDGE, method = "TMM") 
#ListDGE
TMMNorm_ASV_table = cpm(ListDGE)
#head(TMMNorm_ASV_table)
TMMNorm_ASV_table = as.data.frame(TMMNorm_ASV_table)
TMMNorm_ASV_table$ASV.ID = row.names(TMMNorm_ASV_table) #slap the ASV IDs back on the table for joining with the other sample info
Joined<-plyr::join(TMMNorm_ASV_table, ASV_table_noContam[c(1,39)], by="ASV.ID", type="left", match="first")
# the join function joins two data.frames together. Here we are joining the ASV counts (TMMNorm...) and the taxonomy information by the ASV.IDs. type, refers to using the TMMNorm, which is in the left of the function input, as the guide. Match refers to matching the first case of an asv Id.
dim(Joined)
head(Joined)
## Import the full physical and chemical parameters table
Env_table_full = read_csv("Env_Table1.csv", col_types = "ccddddddd") %>% 
  mutate(Date = lubridate::mdy(Date), 
         Year = lubridate::year(Date)) %>% # extract the year from the date using lubridate
  select(Year,Sample,Nitrogen,Phosphorus,Temperature,Salinity,Chlorophyll,-DO,-DA,-Date) # reorder table and remove the fields not used for plotting
### Make discreet Env tables
# 2018
Env_df2018 = Env_table_full %>% 
  filter(Year == 2018) %>% 
  select(-Year,-Sample) %>% 
  t() %>% 
  as.data.frame() 
colnames(Env_df2018) = Env_table_full %>% filter(Year == 2018) %>% pull(Sample)
Env_df2018 = Env_df2018 %>% 
  mutate(ASV.ID = rownames(.)) %>% 
  select(ASV.ID,everything())

# 2019
Env_df2019 = Env_table_full %>% 
  filter(Year == 2019) %>% 
  select(-Year,-Sample) %>% 
  t() %>% 
  as.data.frame() 
Env_df2019 = Env_df2019 %>% 
  mutate(V24 = rownames(.)) %>% 
  select(everything(),V24,-V2)
Env_df_combined = cbind(Env_df2018,Env_df2019)
colnames(Env_df_combined) = columname[c(1,5:42)]
ASVs_W_Envs = rbind(Joined,Env_df_combined)

#Rename function designed to deal with PR2 databases

pr2_rename_taxa<-function(df){
  library(reshape2)
  split<-colsplit(df$taxonomy, ";", c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8"))
  split[ split == "" ] = NA
  split$Taxa<-"Other/unknown"
  split$Taxa[split$Level1 == "No blast hit"]="No blast hit"
  split$Taxa[split$Level1 == "Unassigned"]="Unassigned"
  split$Taxa[split$Level1 == "None"]="None"
  split$Taxa[split$Level1 == "Nitrogen"]="Nitrogen"
  split$Taxa[split$Level1 == "Phosphorus"]="Phosphorus"
  split$Taxa[split$Level1 == "Temperature"]="Temperature"
  split$Taxa[split$Level1 == "Salinity"]="Salinity"
  split$Taxa[split$Level1 == "Chlorophyll"]="Chlorophyll"
  split$Taxa[split$Level2=="Amoebozoa"]="Amoebozoa"
  split$Taxa[split$Level2=="Apusozoa"]="Other/unknown"
  split$Taxa[split$Level2=="Eukaryota"]="Other/unknown"
  split$Taxa[split$Level2=="Stramenopiles"]="Stramenopiles-Other"
  split$Taxa[split$Level2=="Alveolata"]="Alveolates-Other"
  split$Taxa[split$Level2=="Opisthokonta"]="Opisthokonts-Other"
  split$Taxa[split$Level2=="Archaeplastida"]="Archaeplastids"
  split$Taxa[split$Level2=="Excavata"]="Excavates"
  split$Taxa[split$Level2=="Rhizaria"]="Rhizaria-Other"
  split$Taxa[split$Level2=="Hacrobia"]="Hacrobia-other"
  split$Taxa[split$Level3=="Haptophyta"]="Haptophytes"
  split$Taxa[split$Level3=="Fungi"]="Fungi"
  split$Taxa[split$Level3=="Metazoa"]="Metazoa"
  split$Taxa[split$Level3=="Foraminifera"]="Foraminifera"
  split$Taxa[split$Level3=="Dinophyta"]="Dinoflagellates-other"
  split$Taxa[split$Level4=="Syndiniales"]="Syndiniales"
  split$Taxa[split$Level3=="Cryptophyta"]="Cryptophytes"
  split$Taxa[split$Level3=="Ciliophora"]="Ciliates"
  split$Taxa[split$Level3=="Cercozoa"]="Cercozoa"
  split$Taxa[split$Level3=="Radiolaria"]="Radiolaria"
  split$Taxa[split$Level3=="Ochrophyta"]="Unclassified-Ochrophyte"
  split$Taxa[split$Level3=="Choanoflagellida"]="Choanoflagellida"
  split$Taxa[split$Level4=="Pelagophyceae"]="Pelagophytes"
  split$Taxa[split$Level4=="Bacillariophyta"]="Diatoms"
  split$Taxa[split$Level3=="Telonemia"]="Telonemia"
  split$Taxa[split$Level4=="MAST-10"]="MAST"
  split$Taxa[split$Level4=="MAST-2"]="MAST"
  split$Taxa[split$Level4=="MAST-3"]="MAST"
  split$Taxa[split$Level4=="MAST-6"]="MAST"
  split$Taxa[split$Level4=="MAST-8"]="MAST"
  split$Taxa[split$Level4=="MAST-1"]="MAST"
  split$Taxa[split$Level4=="MAST-12"]="MAST"
  split$Taxa[split$Level4=="MAST-25"]="MAST"
  split$Taxa[split$Level4=="MAST-4"]="MAST"
  split$Taxa[split$Level4=="MAST-7"]="MAST"
  split$Taxa[split$Level4=="Dinophyceae"]="Dinophyceae"
  return(split)
} 

NT = pr2_rename_taxa(ASVs_W_Envs)
DBin = data.frame(ASVs_W_Envs, NT)
DBin_noMetazoa = DBin[DBin$Taxa != "Metazoa",] # Remove metazoa from further analysis
dim(DBin_noMetazoa)
colnames(DBin_noMetazoa)
# Breakdown #
nrow(ASV_table) # start 29643
nrow(ASV_table) - nrow(ASV_table.no1) # lost 28540
nrow(ASV_table.no1) - nrow(ASV_table_noContam) # lost 1
nrow(Joined) - nrow(Joined_env) # gain of 5 rows here
nrow(Joined_env) - nrow(DBin) # 0 lost between steps
nrow(DBin) - nrow(DBin_noMetazoa) # lost 30 metazoan ASVs

# Fill in the taxonomy in the blank columns
# Want to pull the names off, rotate them, fill them, rotate them back, and put them back on the count table.
#Pull the names off
DBin_names = DBin_noMetazoa[,40:48]

#Rotate and fill them in
DBin_names_rotated = as.data.frame(t(DBin_names))
DBin_names_rotated = DBin_names_rotated %>% fill(1:length(DBin_names_rotated))

#rerotate
DBin_names = as.data.frame(t(DBin_names_rotated))
columname2= c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8","Tax")
colnames(DBin_names)=columname2
DBin_names$ASV.ID = rownames(DBin_names)

# Reconnect
DBin_2018_2019_joined = plyr::join(DBin_noMetazoa[,1:38],DBin_names[c(1:10)], by = "ASV.ID", match = "first")

#write.csv(DBin_2018_2019_joined, file = "DBin_2018_2019_joined.csv")
#save(DBin_2018_2019_joined, file = "DBin_2018_2019_joined.rda")
## DBin_2018_2019_joined is the datafile that will be used across all analysis - 
## Both years joined, Environmental parameters added, Treatment for blanks, TMM normalized, taxonomy added and filled to level 8, metazoa filtered.
#2018 ASV table
Table_2018 = DBin_2018_2019_joined[c(1:16,39:47)]
ASV_total = apply(Table_2018[2:16],1,sum)
Table_2018.asvfilt = Table_2018[ASV_total > 5,]
nrow(Table_2018) - nrow(Table_2018.asvfilt)# there were 1347 ASVs lost in making this 2018 table because they belonged to the 2019 dataset
#save(Table_2018, file = "Table_2018.rda")

#2019 ASV table
Table_2019 = DBin_2018_2019_joined[c(1,17:47)]
ASV_total = apply(Table_2019[2:23],1,sum)
Table_2019.asvfilt = Table_2019[ASV_total > 5,]
nrow(Table_2019)-nrow(Table_2019.asvfilt)# lost 14,341 ASVs at this step because they were from the 2018 dataset
#save(Table_2019, file = "Table_2019.rda")


##### Figure 3:taxa barplot #####
DBin_2018_2019_joined_dummy = DBin_2018_2019_joined[1:1072,]
DBin_2018_2019_joined_dummy = DBin_2018_2019_joined_dummy %>% 
  mutate(Day16.2018 = 0, Day17.2018 = 0, Day18.2018 = 0, Day19.2018 = 0, Day20.2018 = 0, Day21.2018 = 0, Day22.2018 = 0) %>% 
  select(c(colnames(DBin_2018_2019_joined_dummy[1:16]),Day16.2018:Day22.2018),everything())
data.m = reshape2::melt(DBin_2018_2019_joined_dummy[1:17307,])
data.agg<-aggregate(data.m$value, by=list(Taxa=data.m$Tax,Samples=data.m$variable),sum) #sum sequences by taxonomic group

# Make the order and colors
tax_order=c("Ciliates","Alveolates-Other","Dinophyceae","Syndiniales","Amoebozoa","Excavates","Rhizaria-other","Radiolaria","Cercozoa","Foramnifera","Archaeplastids","Cryptophytes",    "Unclassified-Ochrophyte","Diatoms","Pelagophytes","Stramenopiles-Other","MAST",    "Hacrobia-other","Haptophytes","Telonemia","Choanoflagellida","Fungi","Opisthokonts-Other","Other/unknown")
tax_color = c("#7f0000","#b30000","#e31a1c","#ec7014",   "#8c510a","#dfc27d",    "#08306b","#08519c","#6baed6","#7fcdbb",     "#e7298a","#ffff99",      "#e5f5e0","#a1d99b","#41ab5d","#006d2c","#00441b",      "#8c510a","#f6e8c3","coral1",   "#d8daeb","#b2abd2","#542788","#2d004b")
#The lengths must be the same

names(tax_color)<-tax_order
data.agg = data.agg %>% separate(Samples, into = c("Day","Year"), sep = "\\.")
data.agg$Taxa<-factor(data.agg$Taxa, levels=(tax_order)) #factoring to keep the order

### 2019 plots - which is more informative?
#Bar plot of community composition
Bars_2019 = ggplot(data.agg[order(data.agg$Taxa),] %>% filter(Year == "2019"), aes(y=x,fill=Taxa,x=factor(Day, levels = unique(data.agg$Day)))) +
  geom_bar(position = "fill", stat = "identity", color="black",aes(fill=Taxa))+
  scale_fill_manual(values=tax_color)+
  labs(title="", x="",y="Relative abundance of reads")+
  ggthemes::theme_fivethirtyeight() +
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1.2))

### 2018 plots - which is more informative?
#Bar plot of community composition
Bars_2018 = ggplot(data.agg[order(data.agg$Taxa),] %>% filter(Year == "2018"), aes(y=x,fill=tax,x=factor(Day, levels = unique(data.agg$Day)))) +
    geom_bar(position = "fill", stat = "identity", color="black",aes(fill=Taxa), show.legend = F)+
    scale_fill_manual(values=tax_color)+
    labs(title="", x="",y="Relative abundance of reads")+
    ggthemes::theme_fivethirtyeight() +
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1.2))
 
library(patchwork)
#Joined Relativeplots
Bars_2018 + Bars_2019 + plot_layout(ncol = 1, widths = .5) 



#### Generate Figure 4: Dominant diatom and dino taxa ####

#load("SMP-2018_2019-RE-REanalysis.Rdata")
Diatom_dino_table = DBin_2018_2019_joined_dummy[DBin_2018_2019_joined_dummy$Level4 == "Bacillariophyta",]
melted = Diatom_dino_table %>% reshape2::melt()

#Get the taxa that contain at least 10% of the total reads in 2018 
Dom_tax_2018 = aggregate(melted$value, by=list(Taxa=melted$Level7,Samples=melted$variable),sum) %>% #sum sequences by taxonomic group
  separate(Samples, into = c("Day", "Year")) %>% #split the Samples column into two in order to group by year
  filter(Year == "2018") %>% # only view 2018
  group_by(Taxa) %>% # group by taxa; prepare to sum reads across all samples.
  summarise(total = sum(x)) %>% # number of reads per for each taxa.
  mutate(ratio = total/sum(total)) %>% # make a new column containing the percent of total diatom reads
  filter(ratio > .1) %>% # filter reads below 10% of the total amount of reads
  pull(Taxa) # pull of a vector of the taxa that meet this criterion

# Repeat for 2019. Get the taxa that contained at least 10% of the total reads in 2019
# note: this is done separately because of the difference in read counts between 2019 and 2018; the scales are off.

Dom_tax_2019 = aggregate(melted$value, by=list(Taxa=melted$Level7,Samples=melted$variable),sum) %>% #sum sequences by taxonomic group
  separate(Samples, into = c("Day", "Year")) %>% 
  filter(Year == "2019") %>% 
  group_by(Taxa) %>% 
  summarise(total = sum(x)) %>% 
  mutate(ratio = total/sum(total)) %>% 
  filter(ratio > .1) %>% 
  pull(Taxa)

Dom_tax_table_2018 = aggregate(melted$value, by=list(Taxa=melted$Level7,Samples=melted$variable),sum) %>% #sum sequences by taxonomic group
  separate(Samples, into = c("Day", "Year")) %>% 
  filter(Year == "2018") %>% 
  mutate(ratio = x/sum(x)) %>% 
  filter(Taxa %in% c(Dom_tax_2018,Dom_tax_2019))

Dom_tax_table_2019 = aggregate(melted$value, by=list(Taxa=melted$Level7,Samples=melted$variable),sum) %>% #sum sequences by taxonomic group
  separate(Samples, into = c("Day", "Year")) %>% 
  filter(Year == "2019") %>% 
  mutate(ratio = x/sum(x)) %>% 
  filter(Taxa %in% c(Dom_tax_2018,Dom_tax_2019))

### Box Plots Diatoms ###

Diatom_bars_2018 = Dom_tax_table_2018 %>%
  ggplot(aes(y=(ratio),fill=Taxa,x=factor(Day, levels = c("Day1","Day2","Day3","Day4","Day5","Day6","Day7","Day8","Day9","Day10","Day11","Day12","Day13","Day14","Day15","Day16","Day17","Day18","Day19","Day20","Day21","Day22","Day23"))))+
  geom_bar(position = "stack", stat = "identity", color="black",aes(fill=Taxa),show.legend = F)+
  scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10),
                     labels = c(0,2.5,5.0,7.5,10)) +
  scale_fill_manual(name = "",
                    values = c("#227C9D","#17C3B2","#FFCB77","#faf0dc","#FE6D73","#07004D")) +
  labs(title="", x="",y="Proportion of reads")+
  ggthemes::theme_fivethirtyeight() +
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1,color="black"))

Diatom_bars_2019 = Dom_tax_table_2019 %>%
  ggplot(aes(y=(ratio),fill=Taxa,x=factor(Day, levels = c("Day1","Day2","Day3","Day4","Day5","Day6","Day7","Day8","Day9","Day10","Day11","Day12","Day13","Day14","Day15","Day16","Day17","Day18","Day19","Day20","Day21","Day22","Day23"))))+
  geom_bar(position = "stack", stat = "identity", color="black",aes(fill=Taxa))+
  scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10),
                     labels = c(0,2.5,5.0,7.5,10)) +
  scale_fill_manual(name = "",
                    values = c("#227C9D","#17C3B2","#FFCB77","#faf0dc","#FE6D73","#07004D")) +
  labs(title="", x="",y="Proportion of reads")+
  ggthemes::theme_fivethirtyeight() +
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1,color="black"))

### Repeat for Dinoflagellates ###

Dino_table = DBin_2018_2019_joined_dummy[DBin_2018_2019_joined_dummy$Level4 == "Dinophyceae",]
melted_dinos = Dino_table %>% reshape2::melt()

# Get the taxa that contain at least 10% of the total reads for their respective assemblage
# in 2018 
Dom_dino_tax_2018 = aggregate(melted_dinos$value, by=list(Taxa=melted_dinos$Level7,Samples=melted_dinos$variable),sum) %>% #sum sequences by taxonomic group
  separate(Samples, into = c("Day", "Year")) %>% 
  filter(Year == "2018") %>% # only view 2018
  group_by(Taxa) %>% # group by taxa; prepare to sum reads across all samples.
  summarise(total = sum(x)) %>% # number of reads per for each taxa.
  mutate(ratio = total/sum(total)) %>% # make a new column containing the percent of total diatom reads
  filter(ratio > .1) %>% # filter reads below 10% of the total amount of reads
  pull(Taxa) # pull of a vector of the taxa that meet this criterion

# Repeat for 2019. Get the taxa that contained at least 10% of the total reads in 2019
# note: this is done separately because of the difference in read counts between 2019 and 2018; the scales are off.

Dom_dino_tax_2019 = aggregate(melted_dinos$value, by=list(Taxa=melted_dinos$Level7,Samples=melted_dinos$variable),sum) %>% #sum sequences by taxonomic group
  separate(Samples, into = c("Day", "Year")) %>% 
  filter(Year == "2019") %>% 
  group_by(Taxa) %>% 
  summarise(total = sum(x)) %>% 
  mutate(ratio = total/sum(total)) %>% 
  filter(ratio > .1) %>% 
  pull(Taxa)

Dom_dino_tax_table_2018 = aggregate(melted_dinos$value, by=list(Taxa=melted_dinos$Level7,Samples=melted_dinos$variable),sum) %>% #sum sequences by taxonomic group
  separate(Samples, into = c("Day", "Year")) %>% 
  filter(Year == "2018") %>% 
  mutate(ratio = x/sum(x)) %>% 
  filter(Taxa %in% c(Dom_dino_tax_2018,Dom_dino_tax_2019))

Dom_dino_tax_table_2019 = aggregate(melted_dinos$value, by=list(Taxa=melted_dinos$Level7,Samples=melted_dinos$variable),sum) %>% #sum sequences by taxonomic group
  separate(Samples, into = c("Day", "Year")) %>% 
  filter(Year == "2019") %>% 
  mutate(ratio = x/sum(x)) %>% 
  filter(Taxa %in% c(Dom_dino_tax_2018,Dom_dino_tax_2019))

## Stacked box plots ###

Dino_bars_2018 = Dom_dino_tax_table_2018 %>%
  ggplot(aes(y=(ratio),fill=Taxa,x=factor(Day, levels = c("Day1","Day2","Day3","Day4","Day5","Day6","Day7","Day8","Day9","Day10","Day11","Day12","Day13","Day14","Day15","Day16","Day17","Day18","Day19","Day20","Day21","Day22","Day23"))))+
  geom_bar(position = "stack", stat = "identity", color="black",aes(fill=Taxa),show.legend = F)+
  scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10),
                     labels = c(0,2.5,5.0,7.5,10)) +
  scale_fill_manual(name = "",
                    values = c("#227C9D","#17C3B2","#FFCB77","#FE6D73","#faf0dc")) +
  labs(title="", x="",y="Proportion of reads")+
  ggthemes::theme_fivethirtyeight() +
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1,color="black"))

Dino_bars_2019 = Dom_dino_tax_table_2019 %>%
  ggplot(aes(y=(ratio),fill=Taxa,x=factor(Day, levels = c("Day1","Day2","Day3","Day4","Day5","Day6","Day7","Day8","Day9","Day10","Day11","Day12","Day13","Day14","Day15","Day16","Day17","Day18","Day19","Day20","Day21","Day22","Day23"))))+
  geom_bar(position = "stack", stat = "identity", color="black",aes(fill=Taxa))+
  scale_y_continuous(breaks = c(0,0.025,0.05,0.075,0.10),
                     labels = c(0,2.5,5.0,7.5,10)) +
  scale_fill_manual(name = "",
                    values = c("#227C9D","#17C3B2","#FFCB77","#FE6D73","#faf0dc")) +
  labs(title="", x="",y="Proportion of reads")+
  ggthemes::theme_fivethirtyeight() +
  theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1,color="black"))

#Combine plots
#Combine Diatom boxplots
Diatom_bars_2018 + Diatom_bars_2019 + plot_layout(ncol = 1, widths = .5) 
#combine dinoflagellate boxplots
Dino_bars_2018 + Dino_bars_2019 + plot_layout(ncol = 1, widths = .5) 

#### Figure 5 Dominant diatoms and dinoflagellates #####
##2. Diatoms(generally) (done)
Dtom_Year1_totals = apply(DBin_noMetazoa[DBin_noMetazoa$Level4 == "Bacillariophyta",2:16],1,sum)
Dtom_Year2_totals = apply(DBin_noMetazoa[DBin_noMetazoa$Level4 == "Bacillariophyta",17:38],1,sum)
Dtom_totalsDF = data.frame(Dtom_Year1_totals,Dtom_Year2_totals)
#Filtering
Dtom_Year1_proportion = Dtom_Year1_totals/sum(Dtom_Year1_totals)
Dtom_Year2_proportion = Dtom_Year2_totals/sum(Dtom_Year2_totals)
Dtom_totalsDF_proportion = data.frame(Dtom_Year1_proportion,Dtom_Year2_proportion)
Dtom_filt_proportion = apply(Dtom_totalsDF_proportion,1,sum) # calculate the total 'proportion' across the two years (this is basically to get some decimal because I will later set a cutoff at .01 using this vector to do it)
Dtom_totalsDF_filt100 = Dtom_totalsDF[Dtom_filt_proportion > .01, ] # filter out/ identify ASVs that account for abundances that are more than 1% of all ASVs
#Plot formatting
Dtom_totalsDF_filt100$ASV.ID = rownames(Dtom_totalsDF_filt100)
Dtom_totalsDF_filt100 = Dtom_totalsDF_filt100 %>% select(ASV.ID,"Year1"=Dtom_Year1_totals,"Year2"=Dtom_Year2_totals)
Dtom_Tidy = Dtom_totalsDF_filt100 %>% mutate(CombTotal = apply(Dtom_totalsDF_filt100[2:3],1, sum)) %>% gather("Year", "Count",2:3) %>% arrange(ASV.ID) %>% mutate(Perc=case_when(Year=="Year1"~round(Count/CombTotal*100),TRUE~-round(Count/CombTotal*100)), signal=case_when(Year=="Year1"~1,TRUE~-1))
Dtom_Tidy$ASV.ID = factor(Dtom_Tidy$ASV.ID, levels = unique(Dtom_Tidy$ASV.ID), ordered = F)
#Pyramid plot
Diatom_filt = ggplot(Dtom_Tidy)+
  geom_bar(aes(x=reorder(ASV.ID, -Perc),y=Perc, fill= Year), stat = 'identity')+ 
  coord_flip()+
  ggthemes::theme_fivethirtyeight() +
  scale_y_reverse()

#make upset plots for presence absence
#DF formatting for upset
Dtom_totalsDF$bin.Year1 = ifelse(Dtom_totalsDF$Dtom_Year1_totals > 5,1,0)
Dtom_totalsDF$bin.Year2 = ifelse(Dtom_totalsDF$Dtom_Year2_totals > 5,1,0)
#Plot upset
Dtom_Rawupset = UpSetR::upset(Dtom_totalsDF, sets = c("bin.Year2","bin.Year1"), keep.order = T, order.by = "freq", mainbar.y.label = "Number of ASVs Exclusive to Set/sample", text.scale = 1.5, sets.x.label = "Total ASVs")

### Repeat for dinoflagellates ###
Dino_Year1_totals = apply(DBin_noMetazoa[DBin_noMetazoa$Level4 == "Dinophyceae",2:16],1,sum)
Dino_Year2_totals = apply(DBin_noMetazoa[DBin_noMetazoa$Level4 == "Dinophyceae",17:38],1,sum)
Dino_totalsDF = data.frame(Dino_Year1_totals,Dino_Year2_totals)
#Filtering
Dino_Year1_proportion = Dino_Year1_totals/sum(Dino_Year1_totals)
Dino_Year2_proportion = Dino_Year2_totals/sum(Dino_Year2_totals)
Dino_totalsDF_proportion = data.frame(Dino_Year1_proportion,Dino_Year2_proportion)
Dino_filt_proportion = apply(Dino_totalsDF_proportion,1,sum) # calculate the total 'proportion' across the two years (this is basically to get some decimal because I will later set a cutoff at .01 using this vector to do it)
Dino_totalsDF_filt100 = Dino_totalsDF[Dino_filt_proportion > .01, ] # filter out/ identify ASVs that account for abundances that are more than 1% of all ASVs
#Plot formatting
Dino_totalsDF_filt100$ASV.ID = rownames(Dino_totalsDF_filt100)
Dino_totalsDF_filt100 = Dino_totalsDF_filt100 %>% select(ASV.ID,"Year1"=Dino_Year1_totals,"Year2"=Dino_Year2_totals)
Dino_Tidy = Dino_totalsDF_filt100 %>% mutate(CombTotal = apply(Dino_totalsDF_filt100[2:3],1, sum)) %>% gather("Year", "Count",2:3) %>% arrange(ASV.ID) %>% mutate(Perc=case_when(Year=="Year1"~round(Count/CombTotal*100),TRUE~-round(Count/CombTotal*100)), signal=case_when(Year=="Year1"~1,TRUE~-1))
Dino_Tidy$ASV.ID = factor(Dino_Tidy$ASV.ID, levels = unique(Dino_Tidy$ASV.ID), ordered = TRUE)
#levels(Dino_Tidy$ASV.ID)
Dino_pyramid = ggplot(Dino_Tidy)+
  geom_bar(aes(x=reorder(ASV.ID, -Perc),y=Perc, fill= Year), stat = 'identity') +
  #  scale_fill_manual(values = c("green","blue")) +
  coord_flip() +
  ggthemes::theme_fivethirtyeight() +
  scale_y_reverse()

#make upset plot for presence/absence
#Raw Community (done)
Dino_totalsDF$bin.Year1 = ifelse(Dino_totalsDF$Dino_Year1_totals > 5,1,0)
Dino_totalsDF$bin.Year2 = ifelse(Dino_totalsDF$Dino_Year2_totals > 5,1,0)
Dino_Rawupset = UpSetR::upset(Dino_totalsDF, sets = c("bin.Year2","bin.Year1"), keep.order = F, mainbar.y.label = "Number of ASVs (Dinoflagellates)", text.scale = 1.5, sets.x.label = "Total ASVs")



### Ordination and Alpha diversity ####

DBin_nums_only = DBin_2018_2019_joined[1:1072,2:38]
DBin_relabund = vegan::decostand(DBin_nums_only, MARGIN = 2, method = "total")
DBin_relabund.t = t(DBin_relabund)

cluster<-hclust(dist(DBin_relabund.t), method="average")
plot(cluster, main = "SMP 2018 - 2019 Cluster Dendrogram", ylab = "similarity")

######## Alpha diversity 
#Calculate inverse simpson index
invsimp<-vegan::diversity(DBin_nums_only,index="invsimpson",2)
#Calculate the number of ASVs
ASV_count <-colSums(DBin_nums_only>0) #to evaluate species richness
#Plot
alpha<-data.frame(invsimp,ASV_count) #combine measurements into one dataframe
alpha$samples<-row.names(alpha)
alpha$samples= factor(alpha$samples, levels = alpha$samples)##****** Trick to keep the order in the plot
alpha.m<-reshape2::melt(alpha)
Alpha_plot = ggplot(alpha.m, aes(x=samples, y=value, group=1)) +
  geom_line() +
  geom_point(size=1) +
  facet_grid(variable~.,scales="free") +
  theme_bw() +
  ggthemes::theme_fivethirtyeight()+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5,size=8),axis.text.y=element_text(size=12),legend.position = "top")

