#Written in R v4.3.1
library(tidyverse)    #2.0.0
library(ape)          #5.7-1
library(ComplexUpset) #1.3.3
library(ggforce)      #0.4.2
library(gggenomes)    #1.0.0
library(ggtree)       #3.9.1
library(GO.db)        #3.17.0
library(matrixStats)  #1.3.0
library(patchwork)    #1.2.0
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
  labs(x="Genes only inside element",
       y="Genes with copies also outside element") +
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
  geom_hline(data=data.frame(y=c(2.5, 6.5, 9.5, 13.5)),
             aes(yintercept=y),
             colour="grey",
             linewidth=0.3) +
  geom_hline(data=data.frame(y=c(6.5)),
             aes(yintercept=y),
             linewidth=0.3) +
  scale_fill_manual(values=c("#4B854C", "#C3D6C3"),
                    labels=c("also outside", "only inside")) +
  coord_cartesian(clip="off") +
  labs(x="Genes", y=NULL, fill="Location of all\ngene copies") +
  theme(axis.title=element_blank(),
        axis.text.x=element_text(size=5),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position=c(0.5, 1.2),
        legend.key.size=unit(7, "pt"),
        legend.title=element_text(size=7, face="bold"),
        legend.text=element_text(size=7),
        panel.grid=element_blank(),
        panel.background=element_blank())

#Write to file
pdf(file=paste0("R:/GaeumannomycesStarships/shared_orthos_grid-",
                Sys.Date(), ".pdf"),
    height=3.2, width=4)
(gg.location.grid | gg.location) + 
  plot_layout(widths=c(1, 0.4)) +
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
  dplyr::select(orthogroup, category, everything())

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
                     position='right') + 
      scale_fill_manual(values=c("#e69f00", "#a58337")) +
      ylab("Total number\nof genes") +
      theme(axis.title.x=element_text(size=7),
            axis.text.x=element_text(size=5),
            axis.ticks.x=element_line(),
            legend.key.size=unit(7, "pt"),
            legend.text=element_text(size=7),
            panel.grid=element_blank()),
    base_annotations=
      list(
        'Genes in\nintersection'=
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
          'Genes in\nintersection'=upset_default_themes(
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

#Write to file
pdf(file=paste0("R:/GaeumannomycesStarships/shared_orthos-",
                Sys.Date(), ".pdf"),
    height=3.2, width=4)
gg.upset.summary
dev.off()

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
  dplyr::select(orthogroup, matches(upset.clades$set))

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
      'Genes in\nintersection'=
        intersection_size(width=0.7,
                          text=list(size=1.5),
                          text_mapping=aes(colour='on_background',
                                           y=!!upset.text.position)) +
        scale_y_continuous(expand=expansion(mult=c(0, 0.1)))
    ),
    themes=upset_modify_themes(
      list(
        'Genes in\nintersection'=upset_default_themes(
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
gg.upset.all.member.size / gg.upset.accessory + plot_layout(heights=c(0.5, 2))
dev.off()

#Most common genes (in > 10 elements)
element.genes %>%
  filter(orthogroup %in% element.cargo.accessory.upset.members$orthogroup[element.cargo.accessory.upset.members$num > 10]) %>%
  arrange(orthogroup)


## Inset upset plot (species.level) ##

#Summarise presence-absence matrix by species
element.cargo.clade.count <- element.cargo.df %>%
  mutate(clade=metadata$clade[match(strain, metadata$strain)]) %>%
  filter(!is.na(orthogroup)) %>%
  dplyr::select(clade, orthogroup) %>%
  distinct() %>%
  mutate(value=as.numeric(1)) %>%
  pivot_wider(names_from="clade", values_fill=as.numeric(0)) %>%
  dplyr::select(orthogroup, c("Ga", "GtA", "GtB")) %>%
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

pan.num <- accumulation.df$mean[which(accumulation.df$genomes == 28)]

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

#Read in tree
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

tree.levels <- gg.tree$data %>%
  filter(isTip) %>%
  arrange(desc(y)) %>%
  pull(label)

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
         category=ifelse(is.na(orthogroups.stats$CSEP[match(orthogroup, orthogroups.stats$orthogroup)]),
                         category, "CSEP"),
         category=ifelse(is.na(orthogroups.stats$CAZyme[match(orthogroup, orthogroups.stats$orthogroup)]),
                         category, "CAZyme")) %>%
  dplyr::select(elementID, feat_id="ID", category, start, end, orthogroup) %>%
  mutate(seq_id=sub("Ga-3aA1_s00046,Ga-3aA1_s00047", "Ga-3aA1_s00046", elementID),
         seq_id=sub("^Ga-3aA1_s00047$", "Ga-3aA1_s00046", seq_id),
         seq_id=factor(seq_id, levels=tree.levels))

#Read in BLAST results
blast.files <- Sys.glob(paste0(dir.starship, "blasts/*.tsv"))
blast.files <- blast.files[file.size(blast.files) > 0]
blasts.df <- do.call("rbind", lapply(
  blast.files,
  function(fn)
    data.frame(read.csv(fn, sep="\t", header=FALSE))
)) %>%
  dplyr::rename(accession="V1", ID="V2", evalue="V3", bitscore="V4", pident="V5", length="V6") %>%
  mutate(gene=recode(accession,
                     mp102_11352="FRE",
                     mp102_11357="PLP",
                     mp102_ofg737="NLR",
                     `JX560967.1`="Spok1",
                     `sp|B2AFA8.2|SPOK2_PODAN`="Spok2",
                     `MK521588.1`="Spok3",
                     `MK521589.1`="Spok4"))

#Make dataframe of domains previously identified in element cargos
pfams <- data.frame(
  Pfam=c("PF12520", "PF01794", "PF08022", "PF08030",
         "PF01734", "PF17107", "PF05729", "PF12796",
         "PF00023", "PF06985", "PF05729"),
  domain.name=c("DUF3723", "Ferric_reduc", "FAD_binding_8", "NAD_binding_6",
                "Patatin", "SesA", "NACHT", "Ank_2",
                "Ank", "HET", "NACHT")
)

#Add domain names to AHRD dataframe
ahrds.df <- ahrds.df %>%
  mutate(Pfam=ifelse(grepl("Pfam", Human.Readable.Description),
                     sub("}.*", "", sub(".*Pfam:", "", Human.Readable.Description)), NA),
         domain.name=pfams$domain.name[match(Pfam, pfams$Pfam)])

#Add functional info to element genes dataframe
element.genes <- element.genes %>%
  mutate(blast.hit=blasts.df$gene[match(feat_id, sub("-", "", blasts.df$ID))],
         pfam=ahrds.df$domain.name[match(feat_id, sub("-", "", ahrds.df$Protein.Accession))],
         category=ifelse(!is.na(blast.hit),
                         blast.hit, category),
         category=ifelse(!is.na(pfam) & category == "gene",
                         paste0(pfam, " domain"), category))

#Number of annotated features
table(element.genes$category)

#CSEPs in elements
orthogroups.stats %>%
  filter(orthogroup %in% element.genes$orthogroup,
         !is.na(CSEP))
#CAZymes in elements
orthogroups.stats %>%
  filter(orthogroup %in% element.genes$orthogroup,
         !is.na(CAZyme))

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
element.feats <- bind_rows(element.genes, element.flanks)

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

#Pull sequence coordinates
element.seqs <- gg.element.schematic.tmp %>% 
  get_seqs()

#Add all tracks to schematic
gg.element.schematic <- gggenomes(seqs=element.seqs) +
  geom_feat(data=gg.element.schematic.tmp %>%
               pull_feats(),
            colour="#F7F056", linewidth=5) +
  geom_feat_text(data=gg.element.schematic.tmp %>%
                   pull_feats(),
                 label="Ga-3aA1_s00047", size=2, nudge_y=-0.8, nudge_x=130000) +
  geom_seq(linewidth=0.3) +
  geom_bin_label(data=gg.element.schematic.tmp %>%
                   pull_bins(),
                 size=2) +
  geom_gene(data=gg.element.schematic.tmp %>%
              pull_genes(),
            aes(fill=category),
            shape=0,
            colour=NA) +
  scale_fill_manual(values=c("#E65518", "lightgrey", "#AE76A3", "#F4A736",
                             "#D1BBD7", "#4EB265", "dimgrey", "#7BAFDE",
                             "white", "white"),
                    breaks=c("cap", "gene", "CSEP", "CAZyme",
                             "NLR", "HET domain", "transposable_element_gene", "tyr"),
                    limits=c("cap", "gene", "CSEP", "CAZyme",
                             "NLR", "HET domain", "transposable_element_gene", "tyr"),
                    labels=c("cap", "gene", "CSEP", "CAZyme",
                             "NLR", "HET domain", "TE", "tyr")) +
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
        legend.text=element_text(size=7),
        legend.key.size=unit(8, "pt"),
        legend.title=element_blank(),
        legend.margin=margin(0,0,-10,0),
        strip.background=element_blank(),
        strip.text.y=element_text(size=8, face="bold", margin=margin(-5, 22, -5, 0)),
        panel.spacing=unit(0, "lines"))

#Write to file
pdf(paste0("R:/GaeumannomycesStarships/schematic-",
           Sys.Date(), ".pdf"),
    width=4.5, height=4)
(gg.tree | gg.element.schematic.repeats) +
  plot_layout(widths=c(1, 8))
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

# write.csv(cargo, "Papers/Gaeumannomyces_starships_paper/starfish_cargo.csv", row.names=FALSE)
