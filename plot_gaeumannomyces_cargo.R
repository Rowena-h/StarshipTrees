#Written in R v4.3.1
library(tidyverse)    #2.0.0
library(ape)          #5.7-1
library(ComplexUpset) #1.3.3
library(cowplot)      #1.1.3
library(ggforce)      #0.4.2
library(gggenomes)    #1.0.0
library(ggnewscale)   #0.4.10
library(ggpubr)       #0.6.0
library(ggrepel)      #0.9.5
library(ggtree)       #3.9.1
library(matrixStats)  #1.3.0
library(patchwork)    #1.2.0
library(phangorn)     #2.11.1
library(phytools)     #2.1-1
library(scales)       #1.3.0
library(tgutil)       #0.1.15
library(topGO)        #2.52.0


#Directory paths
dir.phylo <- "R:/GaeumannomycesGenomics/05_phylogenomics/"
dir.comp <- "R:/GaeumannomycesGenomics/07_comparative_genomics/"
dir.starship <- "R:/GaeumannomycesStarships/"


## Import some data from Hill et al. 2024 ##

#Get some variables from ideograms data
attach(paste0(dir.comp, "plot_ideograms_2023-12-19.RData"))
all.genes.df <- all.genes
all.tes.df <- all.genes.tes
detach()

#Load orthogroup data
load(paste0(dir.comp, "orthogroup-matrices-2024-05-16.RData"))

#Read in metadata
strains <- read.csv("R:/GaeumannomycesGenomics/strains", sep="\t", header=FALSE)
metadata <- read.csv(paste0(dir.phylo, "raxmlng/metadata.csv"))

#Fix orthogroup column names
colnames(orthogroups.count) <- strains$V4[match(colnames(orthogroups.count), strains$V1)]

#Get genes in orthogroups
genes.df <- orthogroups %>%
  rownames_to_column("orthogroup") %>%
  gather(strain, gene, -orthogroup) %>%
  separate_rows(sep=", ", gene) %>%
  mutate(new.strain=metadata$new.strain[match(strain, metadata$strain)])

#Read in Starship elements
elements.df <- 
  read.csv(paste0(dir.comp, "starfish/elementFinder/gaeumannomyces.elements.bed"),
           sep="\t", header=FALSE) %>%
  dplyr::rename(seqnames="V1", start="V2", end="V3", gene="V4", category="V5",
                strain="V6", elementID="V7", flank="V8") %>%
  mutate(strain=sub("_.*", "", seqnames),
         seqnames=sub(".*_", "", seqnames),
         gene=sub("-", "", sub("$", ".1", sub("_", "_EIv1_", gene))),
         elementID=str_replace_all(elementID, c("Gt14LH10"="Gt-LH10",
                                                "Gt-3aA1"="Ga-3aA1",
                                                "Gt-CB1"="Ga-CB1"))) %>%
  filter(!elementID %in% c("Gh-1B17_s00001", "Gh-2C17_s00015" , "Gt-4e_s00053",
                           "Gt-4e_s00062", "Gt-4e_s00063", "Gt-4e_s00064",
                           "Gt-LH10_s00079", "Gt-23d_s00102", "Gt-23d_s00109")) %>%
  left_join(genes.df, by=c("gene", "strain")) %>%
  mutate(new.strain=strains$V4[match(strain, strains$V1)])

#Filter for cargo genes
element.cargo.df <- elements.df %>%
  mutate(elementID=sub("Ga-3aA1_s00046,Ga-3aA1_s00047", "Ga-3aA1_s00047", elementID)) %>%
  filter(category != "insert",
         category != "flank",
         category != "cap")

#Number of cargo genes
range(element.cargo.df %>%
        group_by(elementID) %>%
        summarise(n()) %>%
        pull())

#Read in tree for context
kmer.element.big.tree <- read.tree("R:/GaeumannomycesStarships/mashtree/starship_tree_big.bootstrap.tre")
kmer.element.big.tree$tip.label <- sub("__.*", "", sub("gaeumannomyces.elements.id_", "", sub("mycodb.final.starships.id_", "", kmer.element.big.tree$tip.label)))
kmer.element.big.tree$tip.label <- gsub("Gt-3aA1", "Ga-3aA1", gsub("Gt-CB1", "Ga-CB1", gsub("Gt14LH10", "Gt-LH10", kmer.element.big.tree$tip.label)))

#Subset for Gaeumannomyces tips
tree.backbone <- extract.clade(kmer.element.big.tree, 1116)
tree.backbone <- drop.tip(
  tree.backbone,
  tree.backbone$tip.label[!tree.backbone$tip.label %in% elements.df$elementID]
)

#Plot base tree
gg.tree <- ggtree(tree.backbone, branch.length="none", lwd=0.2) +
  xlim(0, 10)

#Vector to order plots as in phylogeny
tree.levels <- gg.tree$data %>%
  filter(isTip) %>%
  arrange(desc(y)) %>%
  pull(label)


###################################
## INSIDE/OUTSIDE CARGO LOCATION ##
###################################

#Get list of non-cargo genes
other.genes <- genes.df$gene[which(!genes.df$gene %in% element.cargo.df$gene)]

#Check for missing genes
length(element.cargo.df$gene) + length(other.genes)
nrow(genes.df)
element.cargo.df$gene[which(!element.cargo.df$gene %in% genes.df$gene)]

#Assign genes as inside elements, outside or both
cargo.location <- list()
inside <- list()
outside <- list()
both <- list()

for (element in unique(element.cargo.df$elementID)) {
  
  strain <- sub("_.*", "", element)
  
  element.cargo.df.tmp <- element.cargo.df$gene[element.cargo.df$elementID == element]
  genes.tmp <- genes.df[intersect(which(genes.df$new.strain == strain),
                                  which(genes.df$gene != "")),]
  other.genes.tmp <- genes.tmp$gene[which(!genes.tmp$gene %in% element.cargo.df.tmp)]
  
  orthogroups.tmp <- unique(genes.df$orthogroup[intersect(which(genes.df$new.strain == strain),
                                                          which(genes.df$gene != ""))])
  
  #Set progress bar
  progress.bar <- txtProgressBar(1, length(orthogroups.tmp), initial=0, char="=", style=3)
  for (j in 1:length(orthogroups.tmp)) {
    
    #Update progress bar
    setTxtProgressBar(progress.bar, j)
    
    genes2.tmp <- genes.tmp$gene[which(genes.tmp$orthogroup %in% orthogroups.tmp[j])]
    
    if (all(!genes2.tmp %in% element.cargo.df.tmp)) {
      cargo.location[[element]][orthogroups.tmp[j]] <- "outside"
      outside[[element]][orthogroups.tmp[j]] <- orthogroups.tmp[j]
    } else if (all(genes2.tmp %in% element.cargo.df.tmp)) {
      cargo.location[[element]][orthogroups.tmp[j]] <- "inside"
      inside[[element]][orthogroups.tmp[j]] <- orthogroups.tmp[j]
    } else if (any(genes2.tmp %in% element.cargo.df.tmp) & 
               any(genes2.tmp %in% other.genes.tmp)) {
      cargo.location[[element]][orthogroups.tmp[j]] <- "both"
      both[[element]][orthogroups.tmp[j]] <- orthogroups.tmp[j]
    }
  }
  close(progress.bar)
  
}

#Combine results in dataframe
cargo.location.df <- do.call(rbind, Map(cbind, lapply(cargo.location, stack), element=names(cargo.location)))
#Convert to matrix
cargo.location.matrix <- outer(inside, both, function(x, y) as.numeric(Map(function(i, j) length(intersect(i, j)), x, y)))
rownames(cargo.location.matrix) <- rownames(cargo.location.matrix)
#Format for plotting
element.cargo.location.df <- cargo.location.matrix %>%
  as_tibble(rownames=NA) %>%
  rownames_to_column("inside") %>%
  pivot_longer(-inside, names_to="outside", values_to="value") %>%
  filter(outside != "Gt-19d1_s00091",
         outside != "Gt-8d_s00067") %>%
  mutate(value=ifelse(inside == outside, NA, value),
         insideID=as.numeric(sub(".*_s", "", inside)),
         outsideID=as.numeric(sub(".*_s", "", outside)),
         inside=fct_reorder(inside, insideID),
         outside=fct_reorder(outside, outsideID))

#Plot as grid
gg.location.grid <- ggplot(element.cargo.location.df, aes(x=inside, y=outside)) +
  geom_tile(aes(fill=value)) +
  geom_hline(data=data.frame(y=c(2.5, 6.5, 9.5, 13.5)),
             aes(yintercept=y),
             colour="grey",
             linewidth=0.3) +
  geom_vline(data=data.frame(x=c(2.5, 6.5, 9.5, 13.5)),
             aes(xintercept=x),
             colour="grey",
             linewidth=0.3) +
  geom_hline(data=data.frame(y=c(6.5)),
             aes(yintercept=y),
             linewidth=0.3) +
  geom_vline(data=data.frame(x=c(6.5)),
             aes(xintercept=x),
             linewidth=0.3) +
  geom_text(aes(label=value),
            size=1.8) +
  geom_text(data=data.frame(y=c(3.5, 12.5), x=0.5,
                            label=c("Ga", "GtB")),
            aes(x=x, y=y, label=label),
            vjust=-1, size=2, angle=90, fontface="italic") +
  geom_text(data=data.frame(x=c(3.5, 12.5), y=0.5,
                            label=c("Ga", "GtB")),
            aes(x=x, y=y, label=label),
            vjust=1.7, size=2, fontface="italic") +
  coord_cartesian(clip="off") +
  labs(x="Orthologous genes only inside element",
       y="Orthologous genes also outside element") +
  scale_x_discrete(position="top") +
  scale_y_discrete(position="right") +
  scale_fill_gradient(low="white", high="#a58337") +
  theme(legend.position="none",
        axis.title=element_text(size=7),
        axis.text.x.top=element_text(size=5, angle=90, vjust=0.5, hjust=0,
                                     margin=margin(b=-1)),
        axis.text.y.right=element_text(size=5,
                                       margin=margin(l=-1)),
        axis.ticks=element_blank(),
        plot.margin=margin(2, 2, 9, 10)) +
  ggpreview(width=3, height=3.2)

#Summarise total count of inside/outside genes per element
cargo.location.counts.df <- cargo.location.df %>%
  filter(values != "outside",
         element %in% element.cargo.location.df$outside) %>%
  group_by(element, values) %>%
  summarise(num.orthos=n()) %>%
  mutate(element=factor(element, levels=levels(element.cargo.location.df$outside)))

#Plot as bargraph
gg.location <- ggplot(cargo.location.counts.df, aes(x=num.orthos, y=element)) +
  geom_bar(stat="identity", aes(fill=values), width=0.8) +
  geom_segment(data=data.frame(y=c(2.5, 6.5, 9.5, 13.5)),
               aes(y=y, yend=y, x=-200, xend=150),
             colour="grey",
             linewidth=0.3) +
  geom_segment(data=data.frame(y=c(6.5)),
             aes(y=y, yend=y, x=-200, xend=150),
             linewidth=0.3) +
  scale_fill_manual(values=c("#4B854C", "#C3D6C3"),
                    labels=c("also outside", "only inside")) +
  coord_cartesian(clip="off", xlim=c(0, 150)) +
  labs(x="Total number\nof orthogroups", y=NULL, fill="Location of all\northogroup copies") +
  theme_minimal() +
  theme(axis.title.x=element_text(size=7), 
        axis.text.x=element_text(size=5),
        axis.text.y=element_text(size=5, margin=margin(r=-1)),
        axis.ticks.x=element_line(),
        legend.position=c(0.2, 1.2),
        legend.key.size=unit(7, "pt"),
        legend.title=element_text(size=7, face="bold"),
        legend.text=element_text(size=7),
        panel.grid.major.y=element_blank(),
        panel.background=element_blank()) +
  ggpreview(width=1.14, height=3)

#Write to file
pdf(file=paste0("R:/GaeumannomycesStarships/shared_orthos_grid-",
                Sys.Date(), ".pdf"),
    height=3.25, width=4.2)
(gg.location.grid | gg.location) + 
  plot_layout(widths=c(1, 0.3)) +
  plot_annotation(tag_levels="a") & 
  theme(plot.tag=element_text(face="bold"))
dev.off()


##############################
## SHARED CARGO UPSET PLOTS ##
##############################

## Main upset plot (intersections > 5) ##

#Make cargo presence-absence matrix
element.cargo.count <- element.cargo.df %>%
  filter(!is.na(orthogroup)) %>%
  dplyr::select(elementID, orthogroup) %>%
  distinct() %>%
  mutate(value=as.numeric(1)) %>%
  pivot_wider(names_from="elementID", values_fill=as.numeric(0)) %>%
  mutate(category=NA) %>%
  dplyr::select(orthogroup, category, rev(tree.levels))

#Categorise genes as accessory/specific to elements
progress.bar <- txtProgressBar(1, nrow(element.cargo.count), initial=0, char="=", style=3)
for (i in 1:nrow(element.cargo.count)) {
  
  #Update progress bar
  setTxtProgressBar(progress.bar, i)
  
  if (length(which(element.cargo.count[i,-c(1,2)] == 0)) == (length(colnames(element.cargo.count[i,-c(1,2)])) - 1)) {
    element.cargo.count$category[element.cargo.count$orthogroup == element.cargo.count$orthogroup[i]] <-
      "specific"
  }
  if (length(which(element.cargo.count[i,-c(1,2)] == 0)) == 0) {
    element.cargo.count$category[element.cargo.count$orthogroup == element.cargo.count$orthogroup[i]] <-
      "core"
  }
  if (length(which(element.cargo.count[i,-c(1,2)] == 0)) < (length(colnames(element.cargo.count[i,-c(1,2)])) - 1) &&
      length(which(element.cargo.count[i,-c(1,2)] == 0)) > 1) {
    element.cargo.count$category[element.cargo.count$orthogroup == element.cargo.count$orthogroup[i]] <-
      "accessory"
  }
}
close(progress.bar)

upset.text.position <- get_size_mode('exclusive_intersection')

#Plot upset plot for intersections with size > 5
gg.upset.summary <-
  upset(
    as.data.frame(element.cargo.count),
    colnames(element.cargo.count)[-c(1,2)],
    width_ratio=0.2,
    height_ratio=2.2,
    min_size=5,
    sort_sets=FALSE,
    sort_intersections_by=c("degree", "cardinality"),
    guides="over",
    matrix=
      (intersection_matrix(geom=geom_point(size=1)) +
         scale_colour_manual(
           values=c("TRUE"="black", "FALSE"="white",
                    "Gh"="#F1F1F1",
                    "Ga"="#d8d8d8",
                    "GtA"="darkgrey",
                    "GtB"="#888888"),
           breaks=c("Gh", "Ga", "GtA", "GtB"),
           labels=c(expression(italic("Gh")),
                    expression(italic("Ga")),
                    expression(paste(italic("Gt"), "A")),
                    expression(paste(italic("Gt"), "B"))))),
    stripes=
      upset_stripes(
        geom=geom_segment(linewidth=3.2),
        data=data.frame(set=colnames(element.cargo.count)[-c(1,2)],
                        clade=metadata$clade[match(sub("_.*", "", colnames(element.cargo.count)[-c(1,2)]),
                                                   metadata$new.strain)]),
        mapping=aes(color=clade),
        colors=c("Gh"="#F1F1F1",
                 "Ga"="#d8d8d8",
                 "GtA"="darkgrey",
                 "GtB"="#888888")
      ),
    set_sizes=
      upset_set_size(geom=geom_bar(aes(fill=category, x=group), width=0.7),
                     position="right") + 
      scale_y_continuous(limits=c(0, 150)) +
      scale_fill_manual(values=c("#e69f00", "#a58337")) +
      ylab("Total number\nof orthogroups") +
      theme(axis.title.x=element_text(size=7),
            axis.text.x=element_text(size=5),
            axis.ticks.x=element_line(),
            legend.key.size=unit(7, "pt"),
            legend.text=element_text(size=7),
            panel.grid=element_blank()),
    base_annotations=
      list(
        'Orthogroups in\nintersection'=
          intersection_size(mapping=aes(fill=category),
                            width=0.7,
                            text=list(size=1.5),
                            text_mapping=aes(colour='on_background',
                                             y=!!upset.text.position)) +
          scale_y_continuous(expand=expansion(mult=c(0, 0.2))) +
          scale_fill_manual(values=c("#e69f00", "#a58337"), guide="none")
      ),
    themes=
      upset_modify_themes(
        list(
          'Orthogroups in\nintersection'=upset_default_themes(
            axis.title.y=element_text(size=7),
            axis.text.y=element_blank(),
            panel.grid.major.y=element_blank(),
            panel.grid.minor.y=element_blank()
          ),
          intersections_matrix=theme(
            legend.position="none",
            axis.title.x=element_blank(),
            axis.text.y=element_text(size=5),
            panel.grid=element_blank()
          )
        )
      )
  ) +
  ggpreview(width=4, height=3.2)

#Make tree to plot alongside
tree.backbone.upset <- drop.tip(
  tree.backbone,
  c("Gt-19d1_s00091", "Gt-8d_s00067", "Gt-4e_s00054", "Gt-LH10_s00089")
)

#Plot base tree
gg.tree.upset <- ggtree(tree.backbone.upset, branch.length="none", lwd=0.2) +
  xlim(0, 8)

#Write to file
pdf(file=paste0("R:/GaeumannomycesStarships/shared_orthos-",
                Sys.Date(), ".pdf"),
    height=3.2, width=4.5)
((plot_spacer() / gg.tree.upset) + plot_layout(heights=c(1, 2.2)) |
  gg.upset.summary) + plot_layout(widths=c(0.1, 1))
dev.off()

#Proportion of specific genes
element.cargo.count %>% 
  pivot_longer(!c(orthogroup,category), names_to="element", values_to="count") %>% 
  filter(count > 0) %>%
  group_by(element, category) %>% 
  summarise(count=n()) %>% 
  pivot_wider(id_cols="element", names_from=category, values_from=count) %>% 
  mutate(prop=specific/(accessory+specific)) %>%
  arrange(desc(prop))


## Supplementary upset plot (accessory cargo) ##

#Arrange elements by clade
upset.clades <- data.frame(
  set=colnames(element.cargo.count)[-c(1,2)],
  clade=metadata$clade[match(sub("_.*", "", colnames(element.cargo.count)[-c(1,2)]),
                             metadata$new.strain)]) %>%
  arrange(match(clade, c('Gh', 'Ga', 'GtA', 'GtB')), clade)

#Filter cargo presence absence for accessory genes (i.e. in at least two elements)
element.cargo.accessory.count <- element.cargo.count %>%
  filter(category == "accessory") %>%
  dplyr::select(orthogroup, rev(tree.levels))

#Generate upset data for accessory genes
element.cargo.accessory.upset <- 
  upset_data(as.data.frame(element.cargo.accessory.count),
             colnames(element.cargo.accessory.count)[-1],
             sort_intersections_by=c("degree", "cardinality"),)

#Pull membership for each intersection
element.cargo.accessory.upset.members <- element.cargo.accessory.upset$presence %>% 
  group_by(intersection) %>% 
  dplyr::slice(1) %>%
  mutate(num=str_count(intersection, "-")+1) %>%
  arrange(num)
#Set order
element.cargo.accessory.upset.members$intersection <- 
  factor(element.cargo.accessory.upset.members$intersection,
         levels=rev(element.cargo.accessory.upset.members$intersection))
#Make labels for plot
element.cargo.accessory.members.labels <-  element.cargo.accessory.upset.members %>%
  group_by(num) %>%
  summarise(pos=mean(as.numeric(intersection)))

#Plot bar plot for number of elements per intersection
gg.upset.all.member.size <- ggplot(element.cargo.accessory.upset.members) +
  geom_bar(stat="identity",
           aes(x=intersection, y=num),
           width=0.7) +
  geom_text(data=element.cargo.accessory.members.labels,
            aes(x=pos, y=num, label=num),
            vjust=-0.5,
            size=1.5) +
  labs(y="Elements in\nintersection", x=NULL) +
  coord_cartesian(clip="off") +
  scale_y_continuous(expand=expansion(mult=c(0, 0.1))) +
  upset_default_themes() +
  theme(axis.title.y=element_text(size=7),
        axis.text=element_blank(),
        axis.ticks=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank()) +
  ggpreview(width=7.5, height=0.7)

#Plot upset plot for accessory cargo
gg.upset.accessory <-
  upset(
    as.data.frame(element.cargo.accessory.count),
    colnames(element.cargo.accessory.count)[-1],
    height_ratio=c(2.5),
    set_sizes=FALSE,
    sort_sets=FALSE,
    sort_intersections_by=c("degree", "cardinality"),
    matrix=
      (intersection_matrix(geom=geom_point(size=0.1),
                           segment=geom_segment(linewidth=0.2)) +
         scale_colour_manual(values=c("TRUE"='black', "FALSE"="white",
                                      'Gh'='#F1F1F1',
                                      'Ga'='#d8d8d8',
                                      'GtA'='darkgrey',
                                      'GtB'='#888888'))),
    stripes=
      upset_stripes(geom=geom_segment(linewidth=2.5),
                    data=upset.clades,
                    mapping=aes(color=clade),
                    colors=c('Gh'='#F1F1F1',
                             'Ga'='#d8d8d8',
                             'GtA'='darkgrey',
                             'GtB'='#888888')),
    base_annotations=list(
      'Orthogroups in\nintersection'=
        intersection_size(width=0.7,
                          text=list(size=1.5),
                          text_mapping=aes(colour='on_background',
                                           y=!!upset.text.position)) +
        scale_y_continuous(expand=expansion(mult=c(0, 0.1)))
    ),
    themes=upset_modify_themes(
      list(
        'Orthogroups in\nintersection'=upset_default_themes(
          axis.title.y=element_text(size=7),
          axis.text.y=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank()
        ),
        intersections_matrix=theme(
          legend.position="none",
          axis.title.x=element_blank(),
          axis.text.y=element_text(size=5),
          panel.grid=element_blank()
        )
      )
    )
  ) +
  ggpreview(width=7.5, height=3)

#Write to file
pdf(file=paste0("R:/GaeumannomycesStarships/shared_orthos_all-",
                Sys.Date(), ".pdf"),
    height=4, width=7.5)
(((plot_spacer() / gg.tree) + plot_layout(heights=c(2.2, 2.5)) |
    (gg.upset.all.member.size / gg.upset.accessory + plot_layout(heights=c(0.5, 2)))) + plot_layout(widths=c(0.1, 2)))
dev.off()


## Inset upset plot (species.level) ##

#Summarise presence-absence matrix by species
element.cargo.clade.count <- element.cargo.df %>%
  mutate(clade=metadata$clade[match(strain, metadata$strain)]) %>%
  filter(!is.na(orthogroup)) %>%
  dplyr::select(clade, orthogroup) %>%
  distinct() %>%
  mutate(value=as.numeric(1)) %>%
  pivot_wider(names_from="clade", values_fill=as.numeric(0)) %>%
  dplyr::select(orthogroup, c("GtA", "Ga", "GtB")) %>%
  dplyr::rename("<i>Gt</i>B"="GtB", "<i>Gt</i>A"="GtA", "<i>Ga</i>"="Ga")

#Plot upset plot at species-level
gg.upset.clades <-
  upset(
    as.data.frame(element.cargo.clade.count),
    colnames(element.cargo.clade.count)[-1],
    height_ratio=1,
    set_sizes=FALSE,
    sort_sets=FALSE,
    sort_intersections=FALSE,
    intersections=list(
      "<i>Gt</i>A", "<i>Gt</i>B", "<i>Ga</i>",
      c("<i>Gt</i>B", "<i>Ga</i>")
    ),
    matrix=
      (intersection_matrix(geom=geom_point(size=1)) +
         scale_colour_manual(
           values=c("TRUE"="black", "FALSE"="white",
                    "<i>Ga</i>"="#d8d8d8",
                    "<i>Gt</i>A"="darkgrey",
                    "<i>Gt</i>B"="#888888"),
           breaks=c("<i>Ga</i>", "<i>Gt</i>A", "<i>Gt</i>B"),
           labels=c(expression(italic("Ga")), expression(paste(italic("Gt"), "A")),
                    expression(paste(italic("Gt"), "B"))))
      ),
    stripes=
      upset_stripes(
        geom=geom_segment(linewidth=2.5),
        data=data.frame(set=c("<i>Ga</i>", "<i>Gt</i>A", "<i>Gt</i>B"),
                        clade=c("<i>Ga</i>", "<i>Gt</i>A", "<i>Gt</i>B")),
        mapping=aes(color=clade),
        colors=c('<i>Ga</i>'='#d8d8d8',
                 '<i>Gt</i>A'='darkgrey',
                 '<i>Gt</i>B'='#888888')
      ),
    base_annotations=list(
      'Intersect.\nratio'=
        intersection_ratio(
          text_mapping=aes(label=paste0(
            !!upset_text_percentage(digits=1),
            '\n(', !!get_size_mode('exclusive_intersection'), ')'
          )),
          bar_number_threshold=1,
          width=0.7,
          text=list(size=1.5)
        ) +
        coord_cartesian(clip="off") +
        scale_y_continuous(expand=expansion(mult=c(0, 0.5)))
    ),
    themes=upset_modify_themes(
      list(
        'Intersect.\nratio'=upset_default_themes(
          axis.title.y=element_text(size=7),
          axis.text.y=element_blank(),
          panel.grid.major.y=element_blank(),
          panel.grid.minor.y=element_blank()
        ),
        intersections_matrix=theme(
          legend.position=c(0.3, 1.5),
          legend.direction="horizontal",
          legend.key.size=unit(7, "pt"),
          legend.text=element_text(size=7, margin=margin(l=0.1)),
          axis.title.x=element_blank(),
          axis.text.y=ggtext::element_markdown(size=5),
          panel.grid=element_blank()
        )
      )
    )
  ) +
  ggpreview(width=1.5, height=1)

#Write to file
pdf(file=paste0("R:/GaeumannomycesStarships/shared_orthos_clades-",
                Sys.Date(), ".pdf"),
    height=1, width=1.5)
gg.upset.clades
dev.off()


###################################
## CARGO HIERARCHICAL CLUSTERING ##
###################################

#Normalise cargo gene presence-absence matrix
element.cargo.grid <- element.cargo.count %>%
  dplyr::select(-category) %>%
  column_to_rownames("orthogroup") %>%
  t() %>%
  scale()

#Make distance matrix
element.cargo.dist <- dist(element.cargo.grid, method="euclidean")

#Do hierarchical clustering
element.cargo.hclust <- hclust(element.cargo.dist, method="complete")

#Convert to trees
tree.clust <- as.phylo(element.cargo.hclust)

#Root element and hclust trees
outgroups <- c("Gt-19d1_s00091", "Gt-8d_s00067")
tree.backbone.rooted <- root(tree.backbone, 
                             outgroups,
                             resolve.root=TRUE,
                             edgelabel=TRUE)
tree.clust.rooted <- root(tree.clust,
                          outgroups,
                          resolve.root=TRUE,
                          edgelabel=TRUE)

#Make tanglegram
tanglegram <- cophylo(tree.backbone.rooted, tree.clust.rooted)

#Extract tanglegram trees
tree.backbone.untangled <- tanglegram$trees[[1]]
tree.clust.untangled <- tanglegram$trees[[2]]

#Plot trees
for (status in c("tree.backbone", "tree.clust")) {
  
  tree.untangled <- get(paste0(status, ".untangled"))
  other.untangled <- get(paste0(c("tree.backbone", "tree.clust")[which(!c("tree.backbone", "tree.clust") %in% status)], ".untangled"))
  
  #Plot base tree
  tree.tmp <- ggtree(tree.untangled, branch.length="none", lwd=0.2, ladderize=FALSE) +
    scale_y_continuous(expand=c(0, 0.5)) +
    theme(plot.title=element_text(hjust=0.5, size=7, face="bold"),
          plot.margin=margin(5.5, 0, 0, 0))
  
  #Flip hclust tree and add title
  if (status == "tree.clust") {
    
    tree.tmp2 <- tree.tmp +
      scale_x_reverse() +
      xlim(30, 0) +
      new_scale_colour() +
      geom_tiplab(geom="label",
                  fill="white",
                  size=1.8,
                  offset=-15,
                  linesize=0.2,
                  label.size=NA,
                  align=TRUE) +
      ggtitle("Cargo orthogroups hclust")
    
  } else {
    
    tree.tmp2 <- tree.tmp +
      xlim(0, 23) +
      new_scale_colour() +
      geom_tiplab(geom="label",
                  fill="white",
                  size=1.8,
                  offset=12,
                  hjust=1,
                  linesize=0.2,
                  label.size=NA,
                  align=TRUE) +
      ggtitle("Element kmer tree")
    
  }
  
  assign(paste0("gg.", status), tree.tmp2)
  
}

#Make dataframe for connecting lines
lines.df <- 
  data.frame(
    label1=gg.tree.backbone$data$label[gg.tree.backbone$data$isTip == TRUE],
    label2=gg.tree.clust$data$label[match(gg.tree.backbone$data$label[gg.tree.backbone$data$isTip == TRUE],
                                          gg.tree.clust$data$label[gg.tree.clust$data$isTip == TRUE])],
    y1=gg.tree.backbone$data$y[gg.tree.backbone$data$isTip == TRUE],
    y2=gg.tree.clust$data$y[match(gg.tree.backbone$data$label[gg.tree.backbone$data$isTip == TRUE],
                                  gg.tree.clust$data$label[gg.tree.clust$data$isTip == TRUE])]
  )

#Plot lines and labels
gg.lines <- ggplot(lines.df) +
  geom_segment(aes(x=0, y=y1, xend=1, yend=y2),
               lwd=0.2,
               lty="dashed") +
  scale_y_continuous(expand=c(0, 0.5)) +
  scale_x_continuous(expand=c(0, 0)) +
  theme_void() +
  theme(legend.position="none")

#Combine into tanglegram
gg.tanglegram <- plot_grid(gg.tree.backbone, gg.lines, gg.tree.clust,
                           ncol=3, rel_widths=c(1, 0.3, 1),
                           align="h", axis="bt")

ggpreview(gg.tanglegram, width=3, height=3)

#Write to file
pdf(file=paste0(dir.starship, "cargo_tanglegram-", Sys.Date(), ".pdf"),
    height=3, width=3)
gg.tanglegram
dev.off()

#Calculate topological distance of hierarchical clustering from element tree
RF.dist(tree.clust.rooted, tree.backbone.rooted, normalize=TRUE)
#Number of bipartitions that differ
print(paste0(RF.dist(unroot(tree.clust.rooted), unroot(tree.backbone.rooted), normalize=TRUE) *
               length(bitsplits(unroot(tree.backbone.rooted))[3]$freq), "/",
             length(bitsplits(unroot(tree.backbone.rooted))[3]$freq)))


###############################
## CARGO ACCUMULATION CURVES ##
###############################

#Function from https://github.com/SioStef/panplots
accumulation_curves <- function(data, curve="pan", iterations=100) {
  
  nr_rows <- nrow(data);
  nr_iterations <- iterations; #the number of iterations (100 by default) 
  #create empty matrix to store temp results
  temp <- matrix(data=NA,nrow=nr_rows,ncol=nr_iterations)
  
  if(curve == "core") {
    ## compute core_genome_accumulation_curve data for the number of iterations
    for(times in 1:nr_iterations){
      # random sampling of genomes
      for (i in 1: nr_rows){
        t=data[sample(nr_rows, i), ,drop=F]
        temp[i,times]=length(which(colSums(t) == i))
      }
    }
  } 
  
  if(curve == "pan") {
    ## compute gene_cluster_accumulation_curve data for the number of iterations
    for(times in 1:nr_iterations){
      # random sampling of genomes
      for (i in 1: nr_rows){
        t=data[sample(nr_rows, i), ,drop=F]
        temp[i,times]=length(which(colSums(t) > 0))
      }
    }
  } 
  
  if(curve == "uniq") {
    ## compute unique gene_cluster_accumulation_curve data for the number of iterations
    for(times in 1:nr_iterations){
      # random sampling of genomes
      for (i in 2: nr_rows){
        t=data[sample(nr_rows, i), ,drop=F]
        temp[i,times]=length(which(colSums(t) == 1))
      }
    }
  }
  
  # summerize permutation results using "matrixStats" library
  summary <- data.frame(genomes=c(1:nr_rows)) 
  summary$mean=rowMeans2(temp[,c(-1)])
  summary$sd=rowSds(temp[,c(-1)])
  summary$group=deparse(substitute(data))
  return(summary)
}

#Format abundance matrix
element.cargo.matrix <- t(element.cargo.count %>% 
                             dplyr::select(-category) %>% 
                             column_to_rownames("orthogroup"))

#Get data for accumulation curves
element.cargo.accumulation.df <- accumulation_curves(element.cargo.matrix, curve="pan")

pan.num <- element.cargo.accumulation.df$mean[which(element.cargo.accumulation.df$genomes == 28)]

#Make function for dynamic axis labels
addUnits <- function(n) {
  labels <- ifelse(n < 1000, n,
                   ifelse(n < 1e6, paste0(round(n/1e3), 'k')
                   )
  )
  return(labels)
}

#Plot accumulation curves
gg.accumulation <- ggplot(element.cargo.accumulation.df, aes(x=genomes, y=mean)) +
  geom_ribbon(aes(ymin=mean-sd, ymax=mean+sd),
              alpha=0.2, show.legend=FALSE) +
  geom_smooth(colour="black", se=FALSE,
              linewidth=0.5) +
  geom_point(size=0.5) +
  labs(x="Number of elements", y="Number of cargo genes") +
  scale_x_continuous(breaks=seq(0, max(element.cargo.accumulation.df$genomes), 2)) +
  scale_y_continuous(labels=addUnits) +
  theme_minimal() +
  theme(legend.position=c(0.15, 0.9),
        legend.title=element_blank(),
        panel.grid.minor.x=element_blank(),
        legend.key.size=unit(8, "pt"),
        legend.text=element_text(size=7, margin=margin(0, 5, 0, 0)),
        axis.title=element_text(size=7),
        axis.text=element_text(size=5)) +
  ggpreview(width=2, height=2)

#Write to file
pdf(paste0("R:/GaeumannomycesStarships/cargo_accumulation_curve-",
           Sys.Date(), ".pdf"),
    width=2, height=2)
gg.accumulation
dev.off()


###################
## GO ENRICHMENT ##
###################

#Read in AHRD annotations
ahrds.df <- do.call("rbind", lapply(
  c(Sys.glob(paste0("S:/CB-GENANNO-520_Mark_McMullan_Wheat_Take-all/Data_Package/*/ahrd_output.csv"))),
  function(fn) 
    data.frame(read.csv(fn, sep="\t", header=TRUE, quote="",
                        fill=TRUE, comment.char="#"))
)) %>%
  mutate(gene=sub("-", "", Protein.Accession))

#Add AHRD annotations to cargo and genes dataframes
element.cargo.go.df <- element.cargo.df %>% 
  left_join(ahrds.df, by="gene") %>%
  filter(Gene.Ontology.Term != "")
all.genes.go.df <- all.genes.df %>%
  mutate(Protein.Accession=sub("$", ".1", ID)) %>%
  left_join(ahrds.df, by="Protein.Accession") %>%
  filter(Gene.Ontology.Term != "")

#Test for GO enrichment per element
for (element in unique(element.cargo.go.df$elementID)) {
  
  universe <- all.genes.go.df[all.genes.go.df$new.strain == sub("_.*", "", element),]
  test.set <- element.cargo.go.df[element.cargo.go.df$elementID == element,]
  
  #Assign genes of interest amongst gene universe
  geneList <- factor(as.integer(universe$gene %in% test.set$gene))
  names(geneList) <- universe$gene
  
  #Format for topGO
  geneID2GO <- setNames(as.list(as.character(universe$Gene.Ontology.Term)), nm=universe$gene)
  
  #Run GO enrichment
  myGOdata <- new("topGOdata", description="Copy-number", ontology="BP",
                  allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
  
  #Fisher's exact test for significance
  resultFisher <- runTest(myGOdata, algorithm="weight01", statistic="fisher")
  
  #Summarise significant results
  sig.df <- GenTable(myGOdata, raw.p.value=resultFisher, topNodes=length(resultFisher@score)) %>%
    filter(raw.p.value < 0.05) %>%
    mutate(element=element)
  
  assign(paste0("sig.df.", element), sig.df)
  
}

#Combine results
bind_rows(mget(ls(pattern="sig.df.")))


################
## SCHEMATICS ##
################

#Combine TE and gene annotations
all.annotations.df <- bind_rows(all.genes.df, all.tes.df) %>%
  mutate(ID=paste0(sub("-", "", ID), ".1"))

#Get all element genes
element.genes <- elements.df %>%
  filter(category != "insert", category != "flank") %>%
  dplyr::rename(ID=gene) %>%
  left_join(all.annotations.df, by=c("ID", "strain", "new.strain")) %>%
  mutate(category=ifelse("transposable_element_gene" == biotype & !is.na(biotype),
                         "transposable_element_gene", category),
         category=sub("\\.", "gene", category),
         category=ifelse(is.na(orthogroups.stats$CAZyme[match(orthogroup, orthogroups.stats$orthogroup)]),
                         category, "CAZyme")) %>%
  dplyr::select(elementID, feat_id="ID", category, start, end, orthogroup) %>%
  mutate(seq_id=sub("Ga-3aA1_s00046,Ga-3aA1_s00047", "Ga-3aA1_s00046", elementID),
         seq_id=sub("^Ga-3aA1_s00047$", "Ga-3aA1_s00046", seq_id),
         seq_id=factor(seq_id, levels=tree.levels))


## FUNCTIONAL INFO ##
  
#Make dataframe of domains previously identified in element cargos/of interest
pfams <- data.frame(
  Pfam=c("PF12520", "PF01794", "PF08022", "PF08030",
         "PF01734", "PF17107", "PF05729", "PF00931", "PF12796", "PF00023", "PF06985"),
  domain.name=c("DUF3723", "Ferric_reduc", "FAD_binding_8", "NAD_binding_6",
                "Patatin", "SesA", "NACHT", "NB-ARC", "Ank_2", "Ank", "HET")
)

#Add domain names to AHRD dataframe
ahrds.df <- ahrds.df %>%
  mutate(Pfam=ifelse(grepl("Pfam", Human.Readable.Description),
                     sub("}.*", "", sub(".*Pfam:", "", Human.Readable.Description)), NA),
         domain.name=pfams$domain.name[match(Pfam, pfams$Pfam)],
         domain.name=ifelse(is.na(domain.name) & grepl("domain", Human.Readable.Description),
                            sub(".*=", "",
                                sub("domain-containing.*", "", Human.Readable.Description)),
                            domain.name))

#Add functional info to element genes dataframe
element.genes.func <- element.genes %>%
  mutate(domain=str_trim(ahrds.df$domain.name[match(feat_id, sub("-", "", ahrds.df$Protein.Accession))]))

#Most common genes/domains (in >=50% elements)
element.genes.domain.sum <- element.genes.func %>%
  filter(!is.na(domain)) %>%
  group_by(domain) %>%
  summarise(num.elements=n_distinct(elementID)) %>%
  arrange(desc(num.elements)) %>%
  mutate(prop=num.elements/length(unique(element.cargo.df$elementID))) %>%
  filter(prop >= 0.5)

#Plot bar graph of common domains in elements
gg.domain.sum <- ggplot(element.genes.domain.sum, aes(x=prop, y=fct_reorder(domain, prop))) +
  geom_bar(stat="identity") +
  scale_x_continuous(limits=c(0, 1),
                     expand=c(0, 0),
                     labels=label_percent()) +
  labs(x="Proportion of elements with domain", y=NULL, title="Domains in at least 50% of elements") +
  theme_minimal() +
  theme(legend.position="none",
        axis.title=element_text(size=6),
        axis.text=element_text(size=6),
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank(),
        plot.title=element_text(size=6, face="bold"),
        plot.margin=margin(5.5, 10, 5.5, 5.5)) +
  ggpreview(width=2.5, height=1.5)

#Most common domains in general
domains.df <- ahrds.df %>%
  mutate(strain=sub("_.*", "", Protein.Accession)) %>%
  filter(!is.na(domain.name)) %>%
  group_by(strain, domain.name) %>%
  summarise(num=n()) %>%
  arrange(strain, desc(num)) %>%
  slice_head(n=5)

#Plot bar graph of common domains across whole genome
gg.domain.all.sum <- ggplot(domains.df, aes(x=num, y=fct_reorder(domain.name, num))) +
  geom_boxplot(linewidth=0.3, outlier.size=0.3) +
  labs(x="Number of genes with domain", y=NULL, title="Most common domains in genome") +
  theme_minimal() +
  theme(legend.position="none",
        axis.title=element_text(size=6),
        axis.text=element_text(size=6),
        axis.ticks.y=element_blank(),
        panel.grid.major.y=element_blank(),
        plot.title=element_text(size=6, face="bold"),
        plot.margin=margin(5.5, 10, 5.5, 5.5)) +
  ggpreview(width=2.5, height=1.5)

#Write to file
pdf(paste0("R:/GaeumannomycesStarships/domains-",
           Sys.Date(), ".pdf"),
    width=5, height=1.5)
plot_grid(gg.domain.sum, gg.domain.all.sum, labels="auto")
dev.off()

#Look at function of common orthogroups (in >=50% elements)
common.genes <- element.cargo.accessory.upset.members$orthogroup[element.cargo.accessory.upset.members$num >= 10]
common.genes.df <- genes.df %>%
  filter(strain %in% c("GtLH10", "Gt-23d", "Gt-4e", "Gt-8d", "Gt-19d1", "Gt-CB1", "Gt-3aA1"),
         orthogroup %in% common.genes) %>%
  left_join(ahrds.df, by="gene") %>%
  arrange(orthogroup)
orthogroups.stats[orthogroups.stats$orthogroup %in% common.genes, c(1:11)]

#Read in BLAST results
blast.files <- Sys.glob(paste0(dir.starship, "blasts/*phi-base*.tsv"))
blast.files <- blast.files[file.size(blast.files) > 0]
blasts.df <- do.call("rbind", lapply(
  blast.files,
  function(fn)
    data.frame(read.csv(fn, sep="\t", header=FALSE))
)) %>%
  dplyr::rename(accession="V1", ID="V2", evalue="V3", bitscore="V4", pident="V5", length="V6") %>%
  filter(pident > 50) %>%
  arrange(ID, desc(bitscore)) %>%
  group_by(ID) %>%
  slice_head(n=1)

#Read in phi-base metadata
phibase.df <- read.csv(paste0(dir.starship, "blasts/phi-base_240801.csv"), skip=1)

#Match phi-base blast hits to orthogroups
genes.cargo.pb.df <- genes.df %>%
  mutate(biotype=orthogroups.stats$biotype[match(orthogroup, orthogroups.stats$orthogroup)]) %>%
  filter(is.na(biotype),
         gene %in% element.cargo.df$gene) %>%
  mutate(blast.hit=blasts.df$accession[match(gene, sub("-", "", blasts.df$ID))],
         pb.id=sub("#.*", "", sub(".*PHI:", "PHI:", blast.hit)),
         pb.gene=phibase.df$Gene[match(pb.id, phibase.df$PHI_MolConn_ID)],
         pb.function=phibase.df$Gene.Function[match(pb.id, phibase.df$PHI_MolConn_ID)],
         pb.species=phibase.df$Pathogen.species[match(pb.id, phibase.df$PHI_MolConn_ID)],
         pb.phen=phibase.df$Mutant.Phenotype[match(pb.id, phibase.df$PHI_MolConn_ID)]) %>%
  arrange(orthogroup)

#Find orthogroups with at least 50% of genes with a phi-base hit
orthogroups.pb <- genes.cargo.pb.df %>%
  group_by(orthogroup) %>%
  summarise(num=n(),
            num.pb=sum(!is.na(pb.gene))) %>%
  mutate(prop=num.pb/num) %>%
  filter(prop >= 0.5) %>%
  pull(orthogroup)
orthogroups.pb.df <- genes.cargo.pb.df %>%
  filter(orthogroup %in% orthogroups.pb,
         !is.na(pb.gene),
         !pb.gene %in% c("FvCpsA")) %>%
  mutate(pb.gene=sub(" \\(.*", "", pb.gene),
         pb.label=paste0(pb.gene, "-like")) %>%
  dplyr::select(orthogroup, pb.gene, pb.label) %>%
  unique() %>%
  print(n=30)

#Which CSEPs in elements
orthogroups.stats %>%
  filter(orthogroup %in% element.genes$orthogroup,
         !is.na(CSEP))
#Which CAZymes in elements
orthogroups.stats %>%
  filter(orthogroup %in% element.genes$orthogroup,
         !is.na(CAZyme))

#Make dataframe of presence of features previously identified in element cargos/of interest
func.sum.df <- data.frame(gene=c("DUF3723", "FRE", "PLP", "Spok (1-4)", "ToxA", "BGC", "NACHT domain-containing"),
                          presence=c(0, 0, 0, 0, 0, 0, 1))

#Plot table
gg.tab <- ggtexttable(func.sum.df, 
                      rows=NULL,
                      cols=c("Gene/feature", "Number"),
                      theme=ttheme(base_size=6,
                                   padding=unit(c(6, 2), "mm"))) %>%
  table_cell_bg(row=c(8), column=2,
              linewidth=2, fill="#CAE0AB", color="white") %>%
  table_cell_font(row=c(8), column=c(1, 2), size=6, face="bold")


## SCHEMATICS ##

#Add functions to element cargo and keep only most common/interesting categories for plotting
element.genes.func <- element.genes.func %>%
  mutate(category=ifelse(orthogroup %in% "N0.HOG0000073",
                         "GT2 CAZyme", category),
         category=ifelse(orthogroup %in% c("N0.HOG0000176", "N0.HOG0000260", "N0.HOG0010081"),
                         "FUG1", category),
         func=ifelse(grepl("NACHT", domain) & category == "gene",
                     "NACHT domain-containing",
                     NA),
         func=ifelse(is.na(orthogroups.stats$CSEP[match(orthogroup, orthogroups.stats$orthogroup)]),
                     func, "CSEP"),
         func=ifelse(category == "gene" & is.na(func),
                     orthogroups.pb.df$pb.label[match(orthogroup, orthogroups.pb.df$orthogroup)],
                     func))

#Get flanking repeats
element.flanks <- elements.df %>% 
  filter(category == "flank") %>%
  dplyr::rename(ID=gene) %>%
  dplyr::select(elementID, feat_id="ID", category, start, end) %>%
  mutate(category=str_match(feat_id, "\\|(.*?)\\|")[,2],
         seq_id=sub("Ga-3aA1_s00046,Ga-3aA1_s00047", "Ga-3aA1_s00046", elementID),
         seq_id=sub("^Ga-3aA1_s00047$", "Ga-3aA1_s00046", seq_id),
         seq_id=factor(seq_id, levels=tree.levels))

#Get nested element
element.nested <- elements.df %>%
  filter(category != "insert", elementID == "Ga-3aA1_s00046,Ga-3aA1_s00047" | elementID == "Ga-3aA1_s00047") %>%
  dplyr::select(elementID, category, start, end) %>%
  mutate(seq_id=sub("^Ga-3aA1_s00047$", "Ga-3aA1_s00046", elementID)) %>%
  filter(seq_id == "Ga-3aA1_s00046") %>%
  group_by(seq_id) %>%
  summarise(start=min(start), end=max(end)) %>%
  mutate(seq_id=factor(seq_id, levels=tree.levels))

#Combine genes and repeats
element.feats <- bind_rows(element.genes.func, element.flanks)

#Plot preliminary schematic and flip element so captains are always first
gg.element.schematic.tmp <- gggenomes(genes=element.feats, feats=element.nested) %>%
  flip_seqs(where(
    ~.x$seq_id %in%
      c("Gt-LH10_s00089", "Gt-LH10_s00088", "Gt-LH10_s00085", "Gt-LH10_s00079", "Gt-LH10_s00074",
        "Gt-4e_s00053", "Gt-4e_s00058", "Gt-4e_s00062",
        "Gt-23d_s00109", "Gt-23d_s00107", "Gt-23d_s00103", "Gt-23d_s00099",
        "Gt-19d1_s00091",
        "Ga-3aA1_s00044",
        "Ga-CB1_s00036",
        "Gh-1B17_s00001")
  ))

#Pull sequence and gene coordinates
element.seqs <- gg.element.schematic.tmp %>%
  get_seqs()
element.genes.flip <- gg.element.schematic.tmp %>%
  pull_genes()
#Format orthologous genes
element.orthos <- genes.df %>% 
  dplyr::select(cluster_id=orthogroup, feat_id=gene)

#Add all tracks to schematic
gg.element.schematic <- gggenomes(seqs=element.seqs, genes=element.genes.flip) |>
  add_clusters(element.orthos) +
  geom_link_line(colour="grey90", linewidth=0.3) +
  geom_feat(data=gg.element.schematic.tmp %>%
              pull_feats(),
            colour="#F7F056", linewidth=5) +
  geom_feat_text(data=gg.element.schematic.tmp %>%
                   pull_feats(),
                 label="Ga-3aA1_s00047", size=2, nudge_y=-0.7, nudge_x=100000) +
  geom_seq(linewidth=0.3) +
  geom_bin_label(data=gg.element.schematic.tmp %>%
                   pull_bins(),
                 size=2) +
  geom_gene(aes(fill=category),
            shape=0,
            colour=NA) +
  geom_label_repel(data=gg.element.schematic.tmp %>%
               pull_genes() %>%
               filter(!is.na(func)),
             aes(x=(x+xend)/2, y=y+0.15, label=func),
             nudge_y=0.3,
             size=1.5,
             fill="black",
             colour="white",
             fontface="bold",
             segment.colour="black",
             segment.size=0.3,
             min.segment.length=0,
             label.size=NA,
             label.padding=unit(1.5, "pt"),
             point.padding=NA,
             box.padding=0,
             direction="x") +
  scale_fill_manual(values=c("#E65518", "grey", "#AE76A3", "#7BAFDE",
                             "#F4A736", "dimgrey", "#FDDBC7",
                           "white", "white"),
                  breaks=c("cap", "gene", "GT2 CAZyme", "CAZyme",
                           "FUG1", "transposable_element_gene", "tyr"),
                  limits=c("cap", "gene", "GT2 CAZyme", "CAZyme",
                           "FUG1", "transposable_element_gene", "tyr"),
                  labels=c("captain", "gene", "GT2 CAZyme", "Other CAZyme",
                           "FUG1-like", "TE", "tyrosine recombinase")) +
  scale_x_continuous(position="top",
                     breaks=pretty_breaks(),
                     labels=label_number(scale=1e-3, suffix=" Kbp")) +
  coord_cartesian(clip="off")

#Get flanking repeats for custom shape plotting
element.repeats <- gg.element.schematic.tmp %>%
  get_seqs() %>%
  left_join(element.flanks, by="seq_id") %>%
  mutate(seq_id=factor(seq_id, levels=levels(element.genes$seq_id)),
         repeat.x=ifelse(xend-x > 0,
                         end.y-start.x,
                         end.x-start.y),
         repeat.xend=ifelse(xend-x > 0,
                            start.y-start.x,
                            end.x-end.y),
         repeat.pos=(repeat.x+repeat.xend)/2)

#Add repeats
gg.element.schematic.repeats <- gg.element.schematic +
  geom_point(data=element.repeats,
             aes(x=repeat.pos, y=y, shape="DR/TIR"),
             size=1, stroke=0.2) +
  scale_shape_manual(values=8) +
  guides(fill=guide_legend(order=1),
         shape=guide_legend(order=2,
                            override.aes=list(size=1.6, stroke=0.6))) +
  theme(legend.position=c(0.85, 0.75),
        legend.text=element_text(size=7, margin=margin(l=1, b=2)),
        legend.key.size=unit(8, "pt"),
        legend.title=element_blank(),
        legend.margin=margin(0,0,-10,0),
        strip.background=element_blank(),
        strip.text.y=element_text(size=8, face="bold", margin=margin(-5, 22, -5, 0)),
        plot.margin=margin(5.5, 5.5, 5.5, 5.5),
        panel.spacing=unit(0, "lines")) +
  # annotation_custom(ggplotGrob(gg.tab),
  #                   xmin=600000, xmax=700000,
  #                   ymin=0, ymax=6) +
  ggpreview(width=7, height=5)

#Write to file
pdf(paste0("R:/GaeumannomycesStarships/schematic-",
           Sys.Date(), ".pdf"),
    width=7, height=5)
(gg.tree | gg.element.schematic.repeats) + 
  plot_layout(widths=c(1, 12))
dev.off()


#############################################################
## Check closely related Gluck-Thaler et al. 2024 elements ##
#############################################################

GT2024.elements <- readxl::read_xlsx("R:/GaeumannomycesStarships/Gluck-Thaler2024/SupplementaryTables1-15.xlsx", skip=1, sheet="S7")

sordario.elements <- GT2024.elements %>%
  dplyr::rename(seqname=contigID, start=begin, ID=geneID, category=featureTag, elementID=starshipID) %>%
  filter(elementID %in% c("pyrory1_s09318", "pyrory1_s09315", "pyrory2_s09330",
                          "pyrory2_s09331", "pyrory2_s09332", "pyrory2_s09327",
                          "pyrory2_s09328", "spobra1_s09737"))

sordario.elements %>%
  group_by(elementID) %>%
  summarise(n())

#GtA gene: Gt8d_EIv1_0050330.1 Gt19d1_EIv1_0000150.1
#Check for Gt-19d1-spobra1 blast hit in elements
grep("spobra1_10718", sordario.elements$ID)
#GtA gene orthogroup
genes.df %>% filter(orthogroup == "N0.HOG0000251")

# write.csv(element.genes.func, "Papers/Gaeumannomyces_starships_paper/starfish_cargo.csv", row.names=FALSE)
