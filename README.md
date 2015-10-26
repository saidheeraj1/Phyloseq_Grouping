# Phyloseq_Grouping
This code is used while grouping the dataset
 library("phyloseq")
 
library("biom")

Grp1<-import_biom("C:/Users/SaiDheeraj/Desktop/Phyloseq/Project files/Mouse/QIIME.grp1.biom")

a=read_biom("C:/Users/SaiDheeraj/Desktop/Phyloseq/Project files/Mouse/QIIME.grp1.biom")

otumat = as(biom_data(a), "matrix")

taxmat = as.matrix(observation_metadata(a), rownames.force=TRUE)

OTU = otu_table(otumat, taxa_are_rows=TRUE)

TAX = tax_table(taxmat)

no1 = phyloseq(OTU, TAX)

library("ape")

tree = rtree(ntaxa(no1), rooted = TRUE, tip.label = taxa_names(no1))

new1 = phyloseq(OTU,TAX,tree)

sampleDF <- read.csv("C:/Users/SaiDheeraj/Desktop/Phyloseq/Project files/Mouse/Grp1.csv", row.names = "Metagenome_ID")

SD<-sample_data(sampleDF)

new1 = merge_phyloseq(new1, SD)

bacteria1<-subset_taxa(new1, taxonomy2=="Bacteria")

b1 = prune_taxa(taxa_sums(bacteria1) > 0, bacteria1)

OTU10 = names(sort(taxa_sums(b1), TRUE)[1:10])

b10=prune_taxa(OTU10, b1)

trial1 <- transform_sample_counts(b10, function(OTU) OTU/sum(OTU))

p<-plot_bar(trial1, "taxonomy3", fill = "taxonomy3", facet_grid=~Mouse.Type)

p + geom_bar(aes(color=taxonomy3, fill= taxonomy3), stat="identity", position="stack")


