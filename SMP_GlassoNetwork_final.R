setwd("~/SMP_R_working_directory")
#save.image(file = "SMP2019_netword.RData")
#load("~/SMP_R_working_directory/2019/SMP2019_anal-21Jun.RData")

#Load required libraries
library(reshape2)
library(tidyverse)
library(ggplot2)
library(vegan)
#library(plyr)
library(SpiecEasi)
#library(intergraph)
library(igraph)
library(limma)
library(viridis)
library(decontam)

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
columname= c("ASV.ID","blank1","blank2","blank3","Day1-2018", "Day2-2018", "Day3-2018", "Day4-2018", "Day5-2018", "Day6-2018", "Day7-2018", "Day8-2018", "Day9-2018", "Day10-2018", "Day11-2018", "Day12-2018", "Day13-2018", "Day14-2018", "Day15-2018", "Day1-2019", "Day2-2019", "Day3-2019", "Day4-2019", "Day5-2019", "Day6-2019", "Day7-2019", "Day8-2019", "Day9-2019", "Day10-2019", "Day11-2019", "Day12-2019", "Day13-2019", "Day14-2019", "Day15-2019", "Day16-2019", "Day17-2019", "Day18-2019", "Day19-2019", "Day20-2019", "Day21-2019", "Day22-2019", "taxonomy")
colnames(ASV_table) = columname
# remove global singletons
ASV_total<-apply(ASV_table[2:41],1,sum)
ASV_table.no1 = ASV_table[ ASV_total >1,]
singles_lost = nrow(ASV_table)-nrow(ASV_table.no1);singles_lost
# 11693 singletons lost

### Blank Treatment ###
#use the decontam R package to handle blanks
#Requires two pieces of information:
# 1)Numeric count table with samples as rows (no strings!) 2) vector identifying the blank samples
# 1) is satisfied with t(ASV_table[2:41])
# 2) Create the vector. The first four samples are the blanks
vector_for_decontamination = c(rep(TRUE,3),rep(FALSE, 37))
# Use decontam to produce a data.frame of the probability of an ASV being an contaminant sequence
contam_df = isContaminant(t(ASV_table.no1[2:41]), neg = vector_for_decontamination)
contam_df %>% filter(contaminant == "TRUE") #note use of the filter() function. Its handy
contam_asvs = row.names(contam_df[contam_df$contaminant == "TRUE", ])
ASV_table.no1[rownames(ASV_table.no1) %in% contam_asvs,]
ASV_table_noContam = ASV_table.no1[!row.names(ASV_table.no1) %in% contam_asvs, ]
contam_loss = nrow(ASV_table.no1)-nrow(ASV_table_noContam) 
# lost 3 contaminant ASVs

ASV_table_noContam = ASV_table_noContam %>% select(-blank1,-blank2,-blank3)

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

NT = pr2_rename_taxa(ASV_table_noContam)
DBin = data.frame(ASV_table_noContam, NT)
DBin_noMetazoa = DBin[DBin$Taxa != "Metazoa",]
#head(DBin_noMetazoa)
# I want to fill in all of the na values with the classification from one level up.
# For this i will use the fill function from tidyverse
#colnames(Table_2018)
Table_2018 = DBin_noMetazoa[,c(1:16,40:48)]
ASV_total = apply(Table_2018[2:16],1,sum)
DBin_noMetazoa = Table_2018[ASV_total > 1,]
nrow(Table_2018)-nrow(DBin_noMetazoa) #lost 1426 ASVs from 2019
colnames(DBin_noMetazoa)

# want to pull the names off, rotate them, fill them, rotate them back, and put them back on the count table.
#Pull the names off
DBin_names = DBin_noMetazoa[,17:25]

#Rotate and fill them in
DBin_names_rotated = as.data.frame(t(DBin_names))
DBin_names_rotated = DBin_names_rotated %>% fill(1:length(DBin_names_rotated))

#rerotate
DBin_names = as.data.frame(t(DBin_names_rotated))
columname2= c("Level1","Level2","Level3","Level4","Level5","Level6", "Level7","Level8","Tax")
colnames(DBin_names)=columname2
DBin_names$ASV.ID = rownames(DBin_names)
head(DBin_names)
colnames(DBin_noMetazoa)

DBin_joined = plyr::join(DBin_noMetazoa[,1:16],DBin_names[c(1:10)], by = "ASV.ID", match = "first")
head(DBin_joined)
#Add a new column in which Level 8 and the Tax columns are joined. I will spit these into two later.
DBin_tax_joined = DBin_joined %>% unite(newTax, Tax, Level8, sep = "__")

Species_Totals = plyr::ddply(DBin_tax_joined, "newTax", summarise, Day1 = sum(Day1.2018), Day2 = sum(Day2.2018), Day3 = sum(Day3.2018), Day4 = sum(Day4.2018), Day5 = sum(Day5.2018), Day6 = sum(Day6.2018), Day7 = sum(Day7.2018), Day8 = sum(Day8.2018), Day9 = sum(Day9.2018), Day10 = sum(Day10.2018), Day11 = sum(Day11.2018), Day12 = sum(Day12.2018), Day13 = sum(Day13.2018), Day14 = sum(Day14.2018), Day15 = sum(Day15.2018))
rownames(Species_Totals) = Species_Totals$newTax
Species_Totals = Species_Totals %>% select(-newTax)
Species_Totals.t = t(Species_Totals)

ASV_total = apply(Species_Totals.t,2, sum)
Networ_taxa_DF_WideMat = data.matrix(Species_Totals.t[,ASV_total>10])


###### Create and quality trim network object ######
#Recommend doing the spiec easi step on a server!! My macbook pro couldn't handle it.
GL <- spiec.easi(Networ_taxa_DF_WideMat, method = 'glasso', lambda.min.ratio=1e-2, nlambda=20, pulsar.params = list(rep.num=50))
### lines for glasso - make graph ###
GL
gl.cor  <- cov2cor(as.matrix(getOptCov(GL))) # If we dont do this step, there is no 'weight' information to add to the adjacency matrix
colnames(gl.cor) <- rownames(gl.cor) <- colnames(Networ_taxa_DF_WideMat) #here i may be able to give the verticies the taxa names if i feed it the vector of names from the levels.

weighted.adj.mat <- gl.cor*getRefit(GL)
colnames(weighted.adj.mat) <- rownames(weighted.adj.mat) <- colnames(Networ_taxa_DF_WideMat)
#Make igraph object from the adjacency matrix
adj_matrix = as.matrix(weighted.adj.mat)
grph <- adj2igraph(as.matrix(weighted.adj.mat))
#nodes:561  links:2582 
### Labels for each node are the rownames of the adjacency matrix
Net_label = rownames(weighted.adj.mat)
# The rownames have two pieces of information, species level and coarse level resolution
# Use split to add both pieces of information to the network

N = as.data.frame(str_split_fixed(Net_label, "__", 2)) # Create a dataframe where the original netlabel containing group and species joined by __ are split up the group and species. These will be assigned to each node
colnames(N) = c("group","species") #Name the columns of this new dataframe appropriately
Net_label1 = N$species
Net_label2 = N$group
### @@ Note to self: there is no need to save the two above into new variables. Just assign them straight up. V(grph)$Taxgroup = N$group

# Add the labels to the nodes(vertexes) of the grph object
V(grph)$Taxgroup = Net_label2
V(grph)$Taxspecies = Net_label1
V(grph)$ID_Numb = V(grph)$name
V(grph)$name = V(grph)$Taxgroup
# Below is my attempt to add information about abundance information. 
# This method will be particularly useful for those taxa that have episodic 
# blooms because their abundances should be concentrated around specific periods of a few days.
Network_longDF = as.data.frame.matrix(t(Networ_taxa_DF_WideMat))
Network_longDF$RelAb = apply(Network_longDF[,1:14],1,sum)/sum(apply(Network_longDF[,1:14],1,sum))
# Below is my attempt to create some sort of variable that is proportional to the taxa's (node's) relative abundance.
y = sqrt(Network_longDF$RelAb)
#y = (y+(min(y)+1)^2/10)
V(grph)$RelAbund = y

#check the attributes of the nodes/vertex
vertex_attr(grph, index = V(grph)) #check the attributes of the nodes/vertex

# Color positive and negative associations (links)
E(grph)[weight > 0.0]$color<-"#0570b0" #now color the edges based on their values positive is steelblue
E(grph)[weight < 0.0]$color<-"#fe9929"  #now color the edges based on their values that are less than the .5 (strongly positive) grey

# Remove edges with very low weight (lost 433 edges) #
weight_threshold <- 0.01
grph_sparse <- delete.edges(grph,which(abs(E(grph)$weight)<0.01))
# Remove nodes with no links
grph_sparse = delete.vertices(grph_sparse,which(igraph::degree(grph_sparse)<1))

plot(grph_sparse,
     vertex.label=NA,
     vertex.size=degree(grph_sparse),
     layout=layout_with_fr(grph_sparse))
title(main="Full Sparse 2018", sub="317-1099")

#note# Its difficult to see any patterns. There are two groups of hubs that are positioned on opposite sides of the ring.
#note# Separating the network into positive and negative networks may help me see the patterns therin

###### separate out positive and negative associations #####
# Create negative network
grph.neg = delete.edges(grph_sparse,which(E(grph_sparse)$weight>0))
grph.neg = delete.vertices(grph.neg,which(igraph::degree(grph.neg)<1))
#Nodes:197 Links:768
dd.grph.neg = degree.distribution(grph.neg) # calculate the degree distribution of the negative network
plot(0:(length(dd.grph.neg)-1), dd.grph.neg, type='b', lty = 1,col = "darkblue",
     ylab="Frequency", xlab="Degree", main="Negative Degree Distributions")

which(igraph::degree(grph.neg)>20) # find out which taxonomic groups have high degrees
plot(grph.neg,
     #vertex.label=V(grph)$Taxgroup ,
     #vertex.label=V(grph.pos)$Taxgroup,
     vertex.label=NA,
     #vertex.color=,
     #edge.color="black",
     #vertex.size = V(grph.pos)$RelAbund/8,
     vertex.size = degree(grph.neg),
     layout=layout_with_fr(grph.neg))

# Negative association network is a big hairball

#Create Positive network 
grph.pos <- delete.edges(grph_sparse,which(E(grph_sparse)$weight<0))
#V(grph.pos)$name = V(grph.pos)$Taxspecies
# r Remove unconnected nodes from the positive network
grph.pos <- delete.vertices(grph.pos,which(igraph::degree(grph.pos)<1))
l = layout_with_fr(grph.pos, niter = 5000)
plot(grph.pos,
     #vertex.label=V(grph)$Taxgroup ,
     #vertex.label=V(grph.pos)$Taxgroup,
     vertex.label=NA,
     #vertex.color=,
     #edge.color="black",
     #vertex.size = V(grph.pos)$RelAbund/8,
     vertex.size = igraph::degree(grph.pos),
     layout=l)
title(main = "pos 2018", sub = "277-727")

# Patterns show up in the positive network




##### Circos plots #####

library(circlize)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
lasso.list.g = adj_matrix %>% get_upper_tri() %>% reshape2::melt() %>% na.omit()

lasso.list.g_taxgroups = lasso.list.g %>% separate(Var1, into = c("Var1", "Var1a"), sep = "__") %>% separate(Var2, into = c("Var2", "Var2a"), sep = "__") %>% select(Var1,Var2,value)

#nrow(lasso.list.g_taxgroups) == nrow(lasso.list.g) # got some warning messages so this is a sanity check that I am not losing data between these two steps.
#******* step 3 remove all non-zeros
lasso.list.g_subset = subset(lasso.list.g_taxgroups, abs(value) > weight_threshold)
#nrow(lasso.list.g_subset) #Saniity check. The rows correspond to an edge
#******** step 4 add positive and neg
lasso.list.g_subset$sign = ifelse(lasso.list.g_subset$value>0,"pos","neg")

######

summary.g = lasso.list.g_subset %>% group_by(Var1,Var2,sign) %>% tally() %>% arrange(desc(n))
tail(summary.g)

sum(summary.g$n); grph_sparse #sanity check. The number of associations matches that of the total graph
# Deal with the fact that you have repeats of phylum1 - phylum2//phylum2 - phylum1
summary.g_neg = summary.g %>% filter(sign == "neg")
summary.g_pos = summary.g %>% filter(sign == "pos")
sum(summary.g_neg$n);grph.neg # 1110 edges in neg network. Sanity check
sum(summary.g_pos$n);grph.pos# 1622 edges in pos network. Sanity check
sum(sum(summary.g_pos$n),sum(summary.g_neg$n)) # 2732 edges in combined network. Sanity check

#optional script to view the proportion of network and some sanity numbers
summary.g_neg %>% 
  # filter(Var1 == "Stramenopiles_Diatoms" | Var2 == "Stramenopiles_Diatoms") %>%
  group_by(Var1,Var2) %>% 
  
  mutate(rat = n/sum(n), sanity = sum(n), sanity2=sum(rat)) %>%
  print(n=100)

# Color codes for taxonomy
taxcolors <- c("Ciliates"="#DD9898","Dinophyceae"="#7490E1","Alveolates-Other"="#D55198", "Syndiniales"="#D846E0" ,"Archaeplastids"="#83E754","Cryptophytes"="#DE9D4E", "Haptophytes"="#C7E080","Hacrobia-other"="#C7E080","Telonemia"="salmon","Choanoflagellida"="#7D8095","Fungi"="#79E1D6","Other/unknown"="#79E1D6","Cercozoa"="#BDE7B6","Radiolaria"="#DA94DC","Diatoms"="#D3DDDA","Unclassified-Ochrophyte"="#78976A","Stramenopiles-Other"="#E7CB9F","Pelagophytes"="#E7CB9F","MAST"="orange")
c(taxcolors)
##################### TIME TO PLOT #######################
# circlize library for circos plots - chordDiagram command
# I think there are other library options for circos plots in R but I haven't explored them yet
# FOR MORE INFO: https://jokergoo.github.io/circlize_book/book/the-chorddiagram-function.html

library(circlize)
### Full network
chordDiagram(summary.g, annotationTrack =c("grid"),
             annotationTrackHeight = c(0.05, 0.3),
             grid.col = taxcolors,
             link.visible = c(summary.g$sign=="pos"))
             #order = c("Alveolates_Ciliates","Alveolates_Dinophyceae","Alveolates_Other","Alveolates_Syndiniales","Archaeplastids","Cryptophytes","Hacrobia_Haptophytes","Hacrobia_other","Opisthokont_Choanoflagellida","Opisthokonts_Fungi","Opisthokonts_Other","Rhizaria_Cercozoa","Rhizaria_Radiolaria","Stramenopiles_Diatoms","Stramenopiles_Ochrophyta","Stramenopiles_Other","Stramenopiles_Pelagophytes"))
# Labeling (doesn't look very good for my plot but perhaps with some tweaking)
# I got this forloop from online FYI
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, facing = "downward", niceFacing = TRUE, col = "black",cex=0.6)
}

#Negative network
chordDiagram(summary.g_neg, annotationTrack =c("grid"),
             annotationTrackHeight = c(0.05, 0.3),
             grid.col = c(taxcolors),
             link.visible = c(summary.g$sign=="pos"))
# Labeling (doesn't look very good for my plot but perhaps with some tweaking)
# I got this forloop from online FYI
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, facing = "downward", niceFacing = TRUE, col = "black",cex=0.6)
}

## Positive network
chordDiagram(summary.g_pos, annotationTrack =c("grid"),annotationTrackHeight = c(0.05, 0.3),grid.col = c(taxcolors),link.visible = c(summary.g$sign=="pos"))
# Labeling (doesn't look very good for my plot but perhaps with some tweaking)
# I got this forloop from online FYI
for(si in get.all.sector.index()) {
  xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
  ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)
  circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 1, facing = "downward", niceFacing = TRUE, col = "black",cex=0.6)
}

# adding legends
legend("topleft",legend = "depth_order", col="red",lty=1,lwd=5)

