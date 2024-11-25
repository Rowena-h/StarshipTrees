#Written in R v4.3.1
library(tidyverse)    #2.0.0
library(ape)          #5.7-1
library(cowplot)      #1.1.3
library(ggnewscale)   #0.4.10
library(ggtree)       #3.9.1
library(ggtreeExtra)  #1.10.0
library(glottoTrees)  #0.1.10
library(phangorn)     #2.11.1
library(phytools)     #2.1-1
library(tgutil)       #0.1.15


#Directory paths
dir.starship <- "R:/GaeumannomycesStarships/"


###############################
## CURATED ELEMENT KMER TREE ##
###############################

#Read in metadata
metadata <- read.csv(paste0(dir.starship, "/metadata.csv")) %>%
  mutate(genus=word(species, 1))

#Read in tree
kmer.element.tree <- 
  read.tree(paste0(dir.starship, "mashtree/starship_tree_conservative.bootstrap.tre"))
kmer.element.tree$tip.label <- 
  sub("__.*", "",
      sub("gaeumannomyces.elements.id_", "",
          sub("Starships.id_", "", kmer.element.tree$tip.label)))

#Plot tree
gg.kmer.element.tree <- ggtree(kmer.element.tree, linewidth=NA, layout="fan") %<+% metadata +
  xlim(-0.3, 2)

# #Select nodes for clades
# gg.kmer.element.tree + 
#   geom_tree() +
#   geom_tiplab(aes(label=element),
#               size=2,
#               align=TRUE) +
#   geom_text(aes(label=node), size=3)

#Make dataframe of clade nodes
clades.df <- data.frame(
  clade=c("Gaeumannomyces", "Fusarium", "Fonsecaea", "Trichoderma",
          "Podospora", "Alternaria", "Aspergillus", "Macrophomina"),
  node=c(68, 60, 57, 65,
         87, 90, 95, 100)
)

#Add features
gg.kmer.element.tree2 <- gg.kmer.element.tree +
  geom_highlight(data=clades.df,
                 aes(node=node),
                 extend=1.3,
                 fill="#F2F2F2",
                 alpha=1) +
  geom_tree(linewidth=0.2,
            aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  scale_colour_manual(values="grey",
                      na.value="black") +
  new_scale_colour() +
  geom_tiplab(geom="label",
              aes(label=element, fill=genus, colour=genus),
              label.padding=unit(0, "pt"),
              label.size=0.5,
              size=1.5,
              align=TRUE,
              linesize=0.2,
              show.legend=FALSE) +
  geom_tiplab(geom="label",
              aes(label=element),
              fill="white",
              alpha=0.6,
              colour="black",
              label.padding=unit(0, "pt"),
              label.size=NA,
              size=1.5,
              align=TRUE,
              linetype=NULL,
              show.legend=FALSE) +
  scale_fill_manual(
    breaks=sort(unique(metadata$genus)),
    values=c('#1965B0', '#D1BBD7', '#81C4E7', '#E8ECFB', '#AA6F9E', '#E65518',
             '#E8601C', '#4EB265', '#5289C7', '#F1932D', '#F6C141', '#F7F056',
             '#994F88', '#EE8026', '#7BAFDE', '#AA6F9E', '#DC050C')
  ) +
  scale_colour_manual(
    breaks=sort(unique(metadata$genus)),
    values=c('#1965B0', '#D1BBD7', '#81C4E7', '#E8ECFB', '#AA6F9E', '#E65518',
             '#E8601C', '#4EB265', '#5289C7', '#F1932D', '#F6C141', '#F7F056',
             '#994F88', '#EE8026', '#7BAFDE', '#AA6F9E', '#DC050C')
  ) +
  geom_strip("Gt-CB1_s00036", "Gt-23d_s00099",
             barsize=0.3,
             offset=1.3) +
  geom_strip("Starship-1_Ffuj", "Starship-1_Foxy_Fo47",
             barsize=0.3,
             offset=1.3) +
  geom_strip("Starship-1_Fped", "Starship-1_Fnub",
             barsize=0.3,
             offset=1.3) +
  geom_strip("Starship-1_Tasp", "Starship-2_Tasp",
             barsize=0.3,
             offset=1.3) +
  geom_strip("Starship-2_Pans", "Enterprise_Pans",
             barsize=0.3,
             offset=1.3) +
  geom_strip("Argo_Mpha", "Phoenix_Mpha",
             barsize=0.3,
             offset=1.3) +
  geom_strip("Starship-1_Agai", "Starship-2_Aten",
             barsize=0.3,
             offset=1.3) +
  geom_strip("Starship-1_Anig", "Starship-1_Aory",
             barsize=0.3,
             offset=1.3) +
  geom_fruit(geom=geom_bar,
             aes(y=tip,
                 x=element.length),
             stat="identity",
             offset=1.8, 
             pwidth=0.5) +
  ggpreview(width=4, height=4)


##################################
## CURATED ELEMENT CAPTAIN TREE ##
##################################

#Read in tree
ml.cap.tree <- read.tree(paste0(dir.starship, "raxmlng/captains.raxml.support"))
ml.cap.tree$tip.label <- 
  sub("__.*", "",
      sub("gaeumannomyces.elements.id_", "",
          sub("Starships.id_", "", ml.cap.tree$tip.label)))

#Plot tree
gg.ml.cap.tree <- ggtree(ml.cap.tree, linewidth=NA, layout="fan") %<+%
  (metadata %>% select(tip=captain, element, genus)) +
  xlim(-0.1, 5) +
  geom_tree(linewidth=0.2,
            aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
            show.legend=FALSE) +
  scale_colour_manual(values="grey",
                      na.value="black") +
  new_scale_colour() +
  geom_tiplab(geom="label",
              aes(label=element, fill=genus, colour=genus),
              label.padding=unit(0, "pt"),
              label.size=0.5,
              size=1.5,
              align=TRUE,
              linesize=0.2,
              show.legend=FALSE) +
  geom_tiplab(geom="label",
              aes(label=element),
              fill="white",
              alpha=0.6,
              colour="black",
              label.padding=unit(0, "pt"),
              label.size=NA,
              size=1.5,
              align=TRUE,
              linetype=NULL,
              show.legend=FALSE) +
  scale_fill_manual(
    breaks=sort(unique(metadata$genus)),
    values=c('#1965B0', '#D1BBD7', '#81C4E7', '#E8ECFB', '#AA6F9E', '#E65518',
             '#E8601C', '#4EB265', '#5289C7', '#F1932D', '#F6C141', '#F7F056',
             '#994F88', '#EE8026', '#7BAFDE', '#AA6F9E', '#DC050C')
  ) +
  scale_colour_manual(
    breaks=sort(unique(metadata$genus)),
    values=c('#1965B0', '#D1BBD7', '#81C4E7', '#E8ECFB', '#AA6F9E', '#E65518',
             '#E8601C', '#4EB265', '#5289C7', '#F1932D', '#F6C141', '#F7F056',
             '#994F88', '#EE8026', '#7BAFDE', '#AA6F9E', '#DC050C')
  ) +
  guides(fill=guide_legend(override.aes=aes(label=""))) +
  theme(legend.position="none") +
  ggpreview(width=2.5, height=2.5)

#Write to file
pdf(file=paste0(dir.starship, "starship_captain_tree-", Sys.Date(), ".pdf"), 
    height=2.5, width=2.5)
gg.ml.cap.tree
dev.off()


########################
## CURATED TANGLEGRAM ##
########################

#Set arbitrary element to root on
arbitrary.root <- "Enterprise_Msp"
arbitrary.cap.root <- "melsp.1_278976_Enterprise-Captain"
kmer.element.tree.rooted <- root(kmer.element.tree, 
                                 arbitrary.root,
                                 resolve.root=TRUE,
                                 edgelabel=TRUE)
ml.cap.tree.rooted <- root(ml.cap.tree,
                           arbitrary.cap.root,
                           resolve.root=TRUE,
                           edgelabel=TRUE)

#Make tip labels match
kmer.element.tree.rooted$tip.label <- 
  metadata$element[match(kmer.element.tree.rooted$tip.label, metadata$tip)]
ml.cap.tree.rooted$tip.label <- 
  metadata$element[match(ml.cap.tree.rooted$tip.label, metadata$captain)]

#Make tanglegram
tanglegram <- cophylo(kmer.element.tree.rooted, ml.cap.tree.rooted)

#Extract tanglegram trees
kmer.element.tree.untangled <- tanglegram$trees[[1]]
ml.cap.tree.untangled <- tanglegram$trees[[2]]

#Plot trees
for (status in c("kmer.element", "ml.cap")) {
  
  tree.untangled <- get(paste0(status, ".tree.untangled"))
  other.untangled <- get(paste0(c("kmer.element", "ml.cap")[which(!c("kmer.element", "ml.cap") %in% status)], ".tree.untangled"))
  
  #Plot base tree
  gg.tree <- ggtree(tree.untangled, branch.length="none",
                    ladderize=FALSE, lwd=NA) %<+% 
    (metadata %>% select(element, everything())) +
    geom_tree(aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
              lwd=0.3,
              show.legend=FALSE) +
    scale_colour_manual(values="darkgrey", na.value="black")
  
  #Flip captain tree and add highlights and title
  if (status == "ml.cap") {
    
    gg.tree2 <- gg.tree +
      scale_x_reverse() +
      xlim(33, 0) +
      new_scale_colour() +
      geom_tiplab(geom="label",
                  aes(fill=genus, colour=genus),
                  label.padding=unit(0, "pt"),
                  label.size=0.5,
                  size=1.2,
                  offset=-13.9,
                  align=TRUE,
                  linesize=0.3) +
      geom_tiplab(geom="label",
                  fill="white",
                  alpha=0.6,
                  colour="black",
                  label.padding=unit(0, "pt"),
                  label.size=NA,
                  size=1.2,
                  offset=-13.9,
                  align=TRUE,
                  linetype=NULL) +
      ggtitle("Captain gene tree (ML)")
    
  } else {
    
    gg.tree2 <- gg.tree +
      xlim(0, 23) +
      new_scale_colour() +
      geom_tiplab(geom="label",
                  aes(fill=genus, colour=genus),
                  label.padding=unit(0, "pt"),
                  label.size=0.5,
                  size=1.2,
                  offset=9.9,
                  hjust=1,
                  align=TRUE,
                  linesize=0.3) +
      geom_tiplab(geom="label",
                  fill="white",
                  alpha=0.6,
                  colour="black",
                  label.padding=unit(0, "pt"),
                  label.size=NA,
                  size=1.2,
                  offset=9.9,
                  hjust=1,
                  align=TRUE,
                  linetype=NULL) +
      ggtitle("Element tree (kmer)")
    
  }
  
  #Add scales and theme
  gg.tree3 <- gg.tree2 +
    scale_fill_manual(
      breaks=sort(unique(metadata$genus)),
      values=c('#1965B0', '#D1BBD7', '#81C4E7', '#E8ECFB', '#AA6F9E', '#E65518',
               '#E8601C', '#4EB265', '#5289C7', '#F1932D', '#F6C141', '#F7F056',
               '#994F88', '#EE8026', '#7BAFDE', '#AA6F9E', '#DC050C')
    ) +
    scale_colour_manual(
      breaks=sort(unique(metadata$genus)),
      values=c('#1965B0', '#D1BBD7', '#81C4E7', '#E8ECFB', '#AA6F9E', '#E65518',
               '#E8601C', '#4EB265', '#5289C7', '#F1932D', '#F6C141', '#F7F056',
               '#994F88', '#EE8026', '#7BAFDE', '#AA6F9E', '#DC050C')
    ) +
    scale_y_continuous(expand=c(0, 2)) +
    theme_void() +
    theme(legend.position="none",
          plot.title=element_text(hjust=0.5, size=7, face="bold"))
  
  assign(paste0("gg.tree.", status), gg.tree3)
  
}

#Make dataframe for connecting lines
lines.df <- 
  data.frame(
    label1=gg.tree.kmer.element$data$label[gg.tree.kmer.element$data$isTip == TRUE],
    label2=gg.tree.ml.cap$data$label[match(gg.tree.kmer.element$data$label[gg.tree.kmer.element$data$isTip == TRUE],
                                           gg.tree.ml.cap$data$label[gg.tree.ml.cap$data$isTip == TRUE])],
    y1=gg.tree.kmer.element$data$y[gg.tree.kmer.element$data$isTip == TRUE],
    y2=gg.tree.ml.cap$data$y[match(gg.tree.kmer.element$data$label[gg.tree.kmer.element$data$isTip == TRUE],
                                   gg.tree.ml.cap$data$label[gg.tree.ml.cap$data$isTip == TRUE])]
  )

#Add genus for colouring and fix y to account for missing tips in captain tree
lines.df <- lines.df %>%
  mutate(genus=metadata$genus[match(label1, metadata$element)],
         y2=length(!which(is.na(gg.tree.kmer.element$data$tip)))/
           length(!which(is.na(gg.tree.ml.cap$data$tip))) * y2)

#Plot lines and labels
gg.lines <- ggplot(lines.df) +
  geom_segment(aes(x=0, y=y1, xend=1, yend=y2,
                   colour=genus),
               lwd=0.2,
               lty="dashed") +
  scale_y_continuous(expand=c(0, 2)) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_colour_manual(
    breaks=sort(unique(metadata$genus)),
    values=c('#1965B0', '#D1BBD7', '#81C4E7', '#E8ECFB', '#AA6F9E', '#E65518',
             '#E8601C', '#4EB265', '#5289C7', '#F1932D', '#F6C141', '#F7F056',
             '#994F88', '#EE8026', '#7BAFDE', '#AA6F9E', '#DC050C')
  ) +
  theme_void() +
  theme(legend.position="none")

#Combine into tanglegram
gg.tanglegram <- plot_grid(gg.tree.kmer.element, gg.lines, gg.tree.ml.cap,
                           ncol=3, rel_widths=c(1, 0.3, 1),
                           align="h", axis="bt")

ggpreview(gg.tanglegram, width=3, height=4)

#Write to file
pdf(file=paste0(dir.starship, "starship_fig1-", Sys.Date(), ".pdf"),
    height=4, width=7)
plot_grid(gg.kmer.element.tree2, gg.tanglegram,
          rel_widths=c(4, 3), labels="auto")
dev.off()

#Calculate topological distance of captain tree from element tree
kmer.element.tree.drop.rooted <- drop.tip(kmer.element.tree.rooted, c("Bdot_Voyager", "Mpha_Derelict"))
RF.dist(kmer.element.tree.drop.rooted, ml.cap.tree.rooted, normalize=TRUE)
#Number of bipartitions that differ
print(paste0(RF.dist(unroot(ml.cap.tree.rooted), unroot(kmer.element.tree.drop.rooted), normalize=TRUE) *
               length(bitsplits(unroot(kmer.element.tree.drop.rooted))[3]$freq), "/",
             length(bitsplits(unroot(kmer.element.tree.drop.rooted))[3]$freq)))

#Calculate percentage of genera that are monophyletic in captain and element trees
monophyly.df <- data.frame(genus=metadata %>% group_by(genus) %>% summarise(num=n()) %>% filter(num > 1) %>% pull(genus),
                           kmer.element=NA,
                           ml.cap=NA)

for (genus in unique(metadata$genus)) {
  
  monophyly.df$kmer.element[monophyly.df$genus == genus] <- 
    is.monophyletic(kmer.element.tree.rooted, metadata$element[metadata$genus == genus])
  monophyly.df$ml.cap[monophyly.df$genus == genus] <- 
    is.monophyletic(ml.cap.tree.rooted, metadata$element[metadata$genus == genus])
  
}

monophyly.df %>%
  summarise(kmer.element=sum(kmer.element, na.rm=TRUE)/n()*100,
            ml.cap=sum(ml.cap, na.rm=TRUE)/n()*100)


#######################################################
## BIG ELEMENT KMER TREE (FULL SUPPLEMENTARY FIGURE) ##
#######################################################

#Read in tree
kmer.element.big.tree <- 
  read.tree(paste0(dir.starship, "mashtree/starship_tree_big.bootstrap.tre"))
kmer.element.big.tree$tip.label <- 
  sub("__.*", "",
      sub("gaeumannomyces.elements.id_", "",
          sub("mycodb.final.starships.id_", "", kmer.element.big.tree$tip.label)))
kmer.element.big.tree$tip.label <- 
  gsub("Gt-3aA1", "Ga-3aA1",
       gsub("Gt-CB1", "Ga-CB1",
            gsub("Gt14LH10", "Gt-LH10", kmer.element.big.tree$tip.label)))

#Read in and format metadata
metadata.big <- read.csv(paste0(dir.starship, "metadata_Gluck-Thaler2024.tsv"), sep="\t")

metadata.tmp <- bind_rows(
  (metadata %>%
     filter(genus == "Gaeumannomyces") %>%
     mutate(genomeCode=sub("_.*", "", element),
            species=word(species, 2),
            lineage="{'clade': 'sordariomyceta',
            'kingdom': 'Fungi',
            'phylum': 'Ascomycota', 'subphylum': 'Pezizomycotina',
            'class': 'Sordariomycetes',
            'subclass': 'Sordariomycetidae',
            'order': 'Magnaporthales',
            'family': 'Magnaporthaceae',
            'subfamily': ''}") %>%
     select(genomeCode, genus, species, isolate="strain", lineage) %>%
     distinct()),
  metadata.big
) %>%
  mutate(family=sub("}", "", word(sub(",.*", "", sub(".* family: ", "", gsub("[']", "", lineage))), 1)),
         order=sub("}", "", word(sub(",.*", "", sub(".* order: ", "", gsub("[']", "", lineage))), 1)),
         class=sub("}", "", word(sub(",.*", "", sub(".* class: ", "", gsub("[']", "", lineage))), 1)))

#Captain classification from Gluck-Thaler et al. 2024
tyr.fams <- data.frame(
  family.name=c("Phoenix", "Hephaestus", "Tardis", "Serenity",
                "Prometheus", "Enterprise", "Galactica", "Moya",
                "Arwing", "Voyager", "Family 11"),
  familyID=c("fam01", "fam02", "fam03", "fam04",
             "fam05", "fam06", "fam07", "fam08",
             "fam09", "fam10", "fam11")
)

#Read in captain classifications for elements and add to metadata
tyr.class <- 
  readxl::read_xlsx(paste0(dir.starship, "Gluck-Thaler2024/SupplementaryTables1-15.xlsx"),
                    skip=1, sheet="S9") %>%
  mutate(family.name=tyr.fams$family.name[match(familyID, tyr.fams$familyID)])

#Add to classifications to metadata
metadata.all <- data.frame(tip=kmer.element.big.tree$tip.label) %>%
  mutate(genomeCode=sub("_.*", "", tip)) %>%
  left_join(metadata.tmp, by="genomeCode") %>%
  mutate(captain.family=tyr.class$family.name[match(tip, tyr.class$starshipID)])

#Plot base tree
gg.kmer.element.big.tree <- 
  ggtree(kmer.element.big.tree,
         aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
         linewidth=0.2, layout="fan") %<+% metadata.all +
  scale_colour_manual(values="grey",
                      na.value="black",
                      guide="none") +
  xlim(-0.4, NA) +
  new_scale_colour() +
  geom_tiplab(geom="label",
              aes(fill=class),
              label.size=NA,
              size=0.8,
              align=TRUE,
              linesize=0.2,
              show.legend=FALSE) +
  geom_tippoint(aes(colour=class, shape=captain.family),
                stroke=0.3,
                size=0.8) +
  scale_shape_manual(breaks=c("Phoenix", "Hephaestus", "Tardis", "Serenity",
                              "Prometheus", "Enterprise", "Galactica", "Moya",
                              "Arwing", "Voyager", "Family 11"),
                     values=c(15, 16, 17, 11, 7, 6, 14, 10, 3, 13, 8),
                     na.translate=FALSE) +
  scale_fill_manual(values=c('#7BAFDE','#D1BBD7', '#4EB265',
                             '#CAE0AB', '#F7F056', '#F4A736', '#DC050C')) +
  scale_colour_manual(values=c('#7BAFDE','#D1BBD7', '#4EB265',
                               '#CAE0AB', '#F7F056', '#F4A736', '#DC050C')) +
  guides(colour=guide_legend(override.aes=list(size=2), order=1),
         shape=guide_legend(override.aes=list(size=2), order=2, ncol=2)) +
  theme(legend.position=c(0.5, 0.5),
        legend.background=element_blank(),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(face="italic", margin=margin(l=0.1, b=0.1, t=0.1)))

#Plot skeleton tree to make labels for monophyletic clades
tmp.tree <- ggtree(kmer.element.big.tree) %<+% metadata.all

#Nodes for monophyletic clades
clade.nodes <- 
  c(860, 1079, 1078, 1075, 1068, 1062, 1089,
    1052, 1044, 1039, 846, 806, 837, 835, 832,
    841, 1092, 1108, 1188, 1215, 1218, 1208,
    1184, 1183, 1180, 1223, 1226, 1176, 1162,
    1172, 1153, 1158, 1125, 1143, 1118, 1149,
    1150, 779, 777, 771, 768, 655, 726, 712,
    709, 711, 650, 647, 759, 756, 642, 638,
    630, 621, 1229)

big.tree.data <- tmp.tree$data %>%
  mutate(genus.node=NA)

#Match nodes to tip labels
for (i in 1:length(clade.nodes)) {
  
  #Extract tips in clade
  tip <- getDescendants(kmer.element.big.tree, clade.nodes[i])
  tip.names <- na.omit(kmer.element.big.tree$tip.label[tip])
  
  #Check that clade is monophyletic as expected
  if (length(unique(tmp.tree$data$genus[tmp.tree$data$label %in% tip.names])) > 1) {
    
    print(paste0("Not monophyletic: ", clade.nodes[i]))
    
  }
  
  big.tree.data$genus.node[big.tree.data$label %in% tip.names] <-
    clade.nodes[i]
  
}

#Add consecutive numbers and number of elements to genus clade labels
big.tree.labels.data <- big.tree.data %>%
  filter(!is.na(genus.node)) %>%
  arrange(y) %>%
  mutate(genus.id=consecutive_id(genus.node)) %>%
  group_by(genus.id) %>%
  mutate(num.elements=n(),
         num.families=n_distinct(captain.family)) %>%
  filter(num.elements > 1) %>%
  slice(1) %>%
  group_by(genus, class) %>%
  mutate(num.clades=n(),
         genus.num=1:n(),
         genus.label=ifelse(
           num.clades > 1,
           paste0(genus, "-", genus.num, " (", num.elements, ")"),
           paste0(genus, " (", num.elements, ")")
         )) %>%
  ungroup() %>%
  select(genus.node, genus.id, genus.label, num.elements, num.families, class)

#Plot tree with clade labels
gg.kmer.element.big.genera.tree <- gg.kmer.element.big.tree +
  geom_cladelab(data=big.tree.labels.data,
                mapping=aes(node=genus.node, label=genus.label, colour=class),
                textcolour="black",
                barsize=0.4,
                offset=0.06,
                fontsize=2,
                offset.text=0,
                align=TRUE,
                angle="auto")

#Write to file
pdf(file=paste0(dir.starship, "starship_kmer_big_tree-", Sys.Date(), ".pdf"),
    height=10, width=10)
gg.kmer.element.big.genera.tree
dev.off()


#################################################
## BIG ELEMENT KMER TREE (MAIN SUMMARY FIGURE) ##
#################################################

#Plot tree with no tip labels
gg.kmer.element.big.summary.tree <- 
  ggtree(kmer.element.big.tree, linewidth=NA, layout="fan") %<+% metadata.all +
  xlim(-0.5, 1) +
  geom_hilight(data=big.tree.labels.data, 
               mapping=aes(node=genus.node, fill=class),
               alpha=0.2,
               align="right",
               show.legend=FALSE) +
  geom_cladelab(data=big.tree.labels.data, 
                mapping=aes(node=genus.node, label=genus.label, colour=class),
                textcolour="black",
                barsize=0.7,
                offset=0.01,
                fontsize=1.5,
                offset.text=0.03,
                align=TRUE,
                angle="auto", 
                show.legend=FALSE) +
  scale_fill_manual(values=c('#7BAFDE','#D1BBD7', '#F7F056', '#F4A736', '#DC050C')) +
  scale_colour_manual(values=c('#7BAFDE','#D1BBD7', '#F7F056', '#F4A736', '#DC050C')) +
  new_scale_colour() +
  geom_tree(aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
            linewidth=0.1, show.legend=FALSE) +
  scale_colour_manual(values="grey",
                      na.value="black") +
  new_scale_colour() +
  geom_tippoint(aes(colour=class),
                size=0.1) +
  scale_colour_manual(values=c('#7BAFDE','#D1BBD7', '#4EB265',
                               '#CAE0AB', '#F7F056', '#F4A736', '#DC050C')) +
  guides(colour=guide_legend(override.aes=list(size=0.4))) +
  theme(legend.position=c(0.5, 0.5),
        legend.background=element_blank(),
        legend.key.size=unit(1, "pt"),
        legend.title=element_blank(),
        legend.text=element_text(face="italic", size=6,
                                 margin=margin(l=2, b=0.1, t=0.1))) +
  ggpreview(width=4, height=4)

#Write to file
pdf(file=paste0(dir.starship, "starship_kmer_big_tree_summary-", Sys.Date(), ".pdf"),
    height=4, width=4)
gg.kmer.element.big.summary.tree
dev.off()


##################################################
## TRIMMED KMER TREE WITH FAMILY CLASSIFICATION ##
##################################################

#Get tips for each clade node
descendants.list <- list()

for (i in 1:length(big.tree.labels.data$genus.node)) {
  
  descendants.list[[big.tree.labels.data$genus.label[i]]] <- 
    na.omit(kmer.element.big.tree$tip.label[
      getDescendants(kmer.element.big.tree, big.tree.labels.data$genus.node[i])
    ])
  
}

#Replace support value node labels with consecutive numbering
kmer.element.big.tree.blank <- makeNodeLabel(kmer.element.big.tree, prefix="")
kmer.element.big.tree.blank$node.label <-
  as.character(as.numeric(kmer.element.big.tree.blank$node.label) +
                 length(kmer.element.big.tree.blank$tip.label))

#Prune tree to one branch per clade
kmer.element.clades.tree <- 
  keep_as_tip(kmer.element.big.tree.blank,
              as.character(big.tree.labels.data$genus.node))

#Truncate clade label
kmer.element.clades.tree.data <- big.tree.labels.data %>%
  mutate(tip.label=as.character(genus.node),
         genus.label=sub(" \\(.*", "", genus.label),
         tmp=ifelse(grepl("-", genus.label), sub(".*-", "", genus.label), NA),
         genus.label2=ifelse(
           grepl("-", genus.label),
           paste0(substr(genus.label, 1, 3), "-", tmp),
           substr(genus.label, 1, 3))) %>%
  select(tip.label, everything())

#Plot base tree
gg.kmer.element.big.family.tree <- 
  ggtree(kmer.element.clades.tree, linewidth=NA, branch.length="none",
         layout="fan") %<+% kmer.element.clades.tree.data

#Check nodes for flipping to match order of other figure
# gg.kmer.element.big.family.tree + 
#   geom_tree(linewidth=0.2) +
#   geom_tiplab(aes(label=genus.label2), size=3, offset=2) +
#   geom_text(aes(label=node), size=3)

gg.kmer.element.big.family.tree2 <- flip(gg.kmer.element.big.family.tree, 83, 108)

for (node in c(91, 92, 85, 86, 88, 129, 130, 133, 114, 71, 74, 66)) {
  
  gg.kmer.element.big.family.tree2 <- rotate(gg.kmer.element.big.family.tree2, node)
  
}

#Reorder tree
kmer.element.clades.ordered.tree <- rotateConstr(kmer.element.clades.tree,
                     (gg.kmer.element.big.family.tree2$data %>%
                        arrange(y) %>%
                        pull(label)))

#Extract nodes corresponding to tips for labelling
tip.nodes <- gg.kmer.element.big.family.tree2$data %>%
  filter(isTip) %>%
  pull(node)

gg.kmer.element.big.family.tree3 <-  
  ggtree(kmer.element.clades.ordered.tree, linewidth=NA, branch.length="none",
         layout="fan", open.angle=25, ladderize=FALSE) %<+% kmer.element.clades.tree.data +
  xlim(-10, 65) +
  geom_hilight(mapping=aes(subset=node %in% tip.nodes, fill=class),
               extend=17,
               alpha=0.2,
               show.legend=FALSE) +
  geom_tree(linewidth=0.2) +
  geom_tiplab(aes(label=genus.label2), size=2, offset=2) +
  geom_tippoint(aes(colour=class),
                size=0.5,
                show.legend=FALSE) +
  scale_fill_manual(values=c('#7BAFDE','#D1BBD7', '#F7F056', '#F4A736', '#DC050C')) +
  scale_colour_manual(values=c('#7BAFDE','#D1BBD7', '#F7F056', '#F4A736', '#DC050C'))

#Make dataframe of captain families per clade
grid.colour <- big.tree.data %>%
  filter(!is.na(genus.node)) %>%
  arrange(y) %>%
  mutate(genus.id=consecutive_id(genus.node),
         tip.label=as.character(genus.node)) %>%
  group_by(genus.id) %>%
  mutate(num.elements=n()) %>%
  filter(num.elements > 1) %>%
  group_by(tip.label, captain.family, num.elements) %>%
  summarise(count=n()) %>%
  mutate(prop=count/num.elements) %>%
  ungroup() %>%
  filter(!is.na(captain.family)) %>%
  complete(tip.label, captain.family, fill=list(prop=0)) %>%
  select(tip.label, captain.family, prop)

#Get angles for text
grid.text <- grid.colour %>%
  arrange(match(tip.label, (gg.kmer.element.big.family.tree3$data %>% arrange(y) %>% pull(label)))) %>%
  filter(prop > 0) %>%
  mutate(new.angle=gg.kmer.element.big.family.tree3$data$angle[match(tip.label, gg.kmer.element.big.family.tree3$data$label)]-90) %>%
  pull(new.angle)

gg.kmer.element.big.family.tree4 <- gg.kmer.element.big.family.tree3 +
  new_scale_fill() +
  geom_fruit(data=grid.colour,
             geom=geom_tile,
             aes(y=tip.label,
                 x=captain.family,
                 fill=prop),
             colour="dimgrey",
             offset=0.45, 
             pwidth=0.6,
             axis.params=list(axis="x",
                              title="Captain family",
                              text.angle=-90,
                              text.size=1.5,
                              line.size=0,
                              hjust=0,
                              vjust=0.5)) +
  geom_fruit(data=grid.colour %>% filter(prop > 0),
             geom=geom_text,
             aes(y=tip.label,
                 x=captain.family,
                 label=round(prop, 2)),
             size=0.8,
             angle=rev(grid.text),
             offset=-0.3,
             pwidth=0.6) +
  scale_fill_gradient(high="orange", low="white",
                       breaks=c(0, 0.5, 1)) +
  labs(fill="Proportion of captain\nfamilies in clade") +
  theme(legend.position=c(0.75, 0.85),
        legend.direction="horizontal",
        legend.title.position="top",
        legend.background=element_blank(),
        legend.key.size=unit(7, "pt"),
        legend.title=element_text(face="bold", size=6, hjust=1),
        legend.text=element_text(size=5, margin=margin(t=2))) +
        ggpreview(width=4, height=4)

#Write to file
pdf(file=paste0(dir.starship, "starship_kmer_big_tree_family_summary-", Sys.Date(), ".pdf"),
    height=4, width=4)
gg.kmer.element.big.family.tree4
dev.off()


##################################
## TAXONOMIC CLASS DISTRIBUTION ##
##################################

#Percentage of explored genomes with elements in Gluck-Thaler et al (2024)
class.elements.df <- metadata.tmp %>% 
  filter(class != "") %>%
  mutate(element=ifelse(genomeCode %in% metadata.all$genomeCode, "Y", "N")) %>%
  group_by(class) %>% summarise(elements=sum(element == "Y"), total=n()) %>%
  mutate(prop=elements/total) %>%
  arrange(desc(prop), desc(total), class) %>%
  filter(elements > 0)

gg.class.prop <- ggplot(class.elements.df, aes(x=prop, y=fct_reorder(class, prop), fill=class)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=paste0(round(prop*100), "%")),
            size=2, hjust=-0.3) +
  scale_fill_manual(values=c('#7BAFDE','#D1BBD7', '#4EB265',
                             '#CAE0AB', '#F7F056', '#F4A736', '#DC050C')) +
  scale_x_continuous(limits=c(0, 1),
                     expand=c(0, 0),
                     labels=scales::label_percent()) +
  labs(x="Percentage of explored genomes with element(s)\n(Gluck-Thaler et al. 2024)", y=NULL) +
  theme(legend.position="none",
        axis.title=element_text(size=6),
        axis.text.y=element_text(size=6, face="italic"),
        axis.text.x=element_text(size=6),
        panel.grid.major.y=element_blank(),
        plot.margin=margin(5.5, 10, 5.5, 5.5)) +
  ggpreview(width=3, height=2)

gg.class.total <- ggplot(class.elements.df, aes(x=total, y=fct_reorder(class, prop), fill=class)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=total),
            size=2, hjust=-0.3) +
  scale_fill_manual(values=c('#7BAFDE','#D1BBD7', '#4EB265',
                             '#CAE0AB', '#F7F056', '#F4A736', '#DC050C')) +
  scale_x_continuous(expand=expansion(mult=c(0, 0.15))) +
  labs(x="Total explored genomes\n(Gluck-Thaler et al. 2024)", y=NULL) +
  theme(legend.position="none",
        axis.title=element_text(size=6),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=6),
        panel.grid.major.y=element_blank(),
        plot.margin=margin(5.5, 10, 5.5, 5.5)) +
  ggpreview(width=3, height=2)

#Write to file
pdf(file=paste0(dir.starship, "gluckthaler2024_classes-", Sys.Date(), ".pdf"),
    height=2, width=5)
plot_grid(gg.class.prop, gg.class.total, rel_widths=c(1,0.75))
dev.off()


###############################################################
## CONSERVATIVE BIG ELEMENT KMER TREE (SUPPLEMENTARY FIGURE) ##
###############################################################

#Read in tree
kmer.element.big.cons.tree <- 
  read.tree(paste0(dir.starship, "mashtree/starship_tree_big_conservative.bootstrap.tre"))
kmer.element.big.cons.tree$tip.label <- 
  sub("__.*", "",
      sub("gaeumannomyces.elements.id_", "",
          sub("mycodb.final.starships.id_", "", kmer.element.big.cons.tree$tip.label)))
kmer.element.big.cons.tree$tip.label <- 
  gsub("Gt-3aA1", "Ga-3aA1",
       gsub("Gt-CB1", "Ga-CB1",
            gsub("Gt14LH10", "Gt-LH10", kmer.element.big.cons.tree$tip.label)))

#Read in and format metadata
metadata.big.cons <- read.csv(paste0(dir.starship, "metadata_Gluck-Thaler2024.tsv"), sep="\t")

metadata.tmp <- bind_rows(
  (metadata %>%
     filter(genus == "Gaeumannomyces") %>%
     mutate(genomeCode=sub("_.*", "", element),
            species=word(species, 2),
            lineage="{'clade': 'sordariomyceta',
            'kingdom': 'Fungi',
            'phylum': 'Ascomycota', 'subphylum': 'Pezizomycotina',
            'class': 'Sordariomycetes',
            'subclass': 'Sordariomycetidae',
            'order': 'Magnaporthales',
            'family': 'Magnaporthaceae',
            'subfamily': ''}") %>%
     select(genomeCode, genus, species, isolate="strain", lineage) %>%
     distinct()),
  metadata.big.cons
) %>%
  mutate(family=sub("}", "", word(sub(",.*", "", sub(".* family: ", "", gsub("[']", "", lineage))), 1)),
         order=sub("}", "", word(sub(",.*", "", sub(".* order: ", "", gsub("[']", "", lineage))), 1)),
         class=sub("}", "", word(sub(",.*", "", sub(".* class: ", "", gsub("[']", "", lineage))), 1)))

#Add to classifications to metadata
metadata.all.cons <- data.frame(tip=kmer.element.big.cons.tree$tip.label) %>%
  mutate(genomeCode=sub("_.*", "", tip)) %>%
  left_join(metadata.tmp, by="genomeCode") %>%
  mutate(captain.family=tyr.class$family.name[match(tip, tyr.class$starshipID)])

#Plot base tree
gg.kmer.element.big.cons.tree <- 
  ggtree(kmer.element.big.cons.tree,
         aes(colour=ifelse(as.numeric(label) < 70, "insig", NA)),
         linewidth=0.2, layout="fan") %<+% metadata.all.cons +
  scale_colour_manual(values="grey",
                      na.value="black",
                      guide="none") +
  xlim(-0.4, NA) +
  new_scale_colour() +
  geom_tiplab(geom="label",
              aes(fill=class),
              label.size=NA,
              size=0.9,
              align=TRUE,
              linesize=0.2,
              show.legend=FALSE) +
  geom_tippoint(aes(colour=class, shape=captain.family),
                stroke=0.3,
                size=0.8) +
  scale_shape_manual(breaks=c("Phoenix", "Hephaestus", "Tardis", "Serenity",
                              "Prometheus", "Enterprise", "Galactica", "Moya",
                              "Arwing", "Voyager", "Family 11"),
                     values=c(15, 16, 17, 11, 7, 6, 14, 10, 3, 13, 8),
                     na.translate=FALSE) +
  scale_fill_manual(values=c('#7BAFDE','#D1BBD7', '#CAE0AB',
                             '#F7F056', '#F4A736', '#DC050C')) +
  scale_colour_manual(values=c('#7BAFDE','#D1BBD7', '#CAE0AB',
                               '#F7F056', '#F4A736', '#DC050C')) +
  guides(colour=guide_legend(override.aes=list(size=2), order=1),
         shape=guide_legend(override.aes=list(size=2), order=2, ncol=2)) +
  theme(legend.position=c(0.5, 0.5),
        legend.background=element_blank(),
        legend.key.size=unit(10,"pt"),
        legend.title=element_blank(),
        legend.text=element_text(face="italic", margin=margin(l=0.1, b=0.1, t=0.1))) +
  ggpreview(width=8, height=8)

#Write to file
pdf(file=paste0(dir.starship, "starship_kmer_big_conservative_tree-", Sys.Date(), ".pdf"),
    height=8, width=8)
gg.kmer.element.big.cons.tree
dev.off()
