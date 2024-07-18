library("data.table")
library("dplyr")
library("ggplot2")

# whole classification file 
sqantidir <- "/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/6_sqanti3/"
METADIR <- "/lustre/projects/Research_Project-MRC148213/Rosie/SFARIdevelopmentalgenomics/7_differential/3_DIU/metadata"
REFDIR <- "/lustre/projects/Research_Project-MRC148213/lsl693/reference/human"
class.files.names <- paste0(sqantidir, "WholeTargeted_cleaned_aligned_merged_collapsed_qced_RulesFilter_Whole_2reads2samples_classification_noMonoIntergenic.txt")

class.files <- fread(class.files.names)
knownGene <- class.files[class.files$novelGene == "Annotated Genes"] 
knownGene <- knownGene[,c("isoform", "chrom", "exons", "structural_category", "associated_gene", "associated_transcript", "novelGene")]
ensemblGene <- data.frame(fread(paste0(REFDIR, "/gencode.v40.annotation.geneannotation.txt")))
knownGene <- merge(knownGene, ensemblGene, by.x = "associated_gene", by.y = "GeneSymbol", all.x = T)

knownGene %>% select(associated_gene, Class) %>% distinct() %>% group_by(Class) %>% tally() %>% 
  ggplot(., aes(x = reorder(Class, -n), y = log10(n))) + geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
table(knownGene$Class)

proteinCodingGene <- knownGene[knownGene$Class == "protein_coding",]
lncRNAGene <- knownGene[knownGene$Class == "lncRNA",]

write.table(unique(knownGene$associated_gene), paste0(METADIR, "/WholeknownGenes.txt"),quote=F,row.names = F, col.names = F)
write.table(unique(proteinCodingGene$associated_gene), paste0(METADIR, "/WholeProteinCodingGenes.txt"),quote=F,row.names = F, col.names = F)
write.table(unique(lncRNAGene$associated_gene), paste0(METADIR, "/WholelncRNAGenes.txt"),quote=F,row.names = F, col.names = F)