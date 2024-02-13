library(ggplot2)
library(ggtranscript)
library(plyranges)
library(vroom)

y<-class.files$glob_targ_SQ
x1 <- read_gff ("/gpfs/mrc0/projects/Research_Project-MRC148213/Rosie/WholeTargeted/SQANTI/sqanti_with_ORF/WholeTargeted_with_ORF_corrected_renamed.gtf.cds.gff")
#x <- as.data.frame(x1 %>% filter(transcript_id %in% class.files$glob_targ_SQ[class.files$glob_targ_SQ$associated_gene == "CDH8","isoform"][1:10]))
dat2 <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/3_ficle_tc20bp/TRIO/Stats/TRIO_Exonskipping_tab.csv")
#dat2 <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/3_ficle_tc20bp/GRIK2/Stats/GRIK2_Exonskipping_tab.csv")
extensiveES <- dat2 %>% group_by(transcript_id) %>% tally() %>% filter(n > 40)
x <- as.data.frame(x1 %>% filter(transcript_id %in% 
                                   c(as.character(extensiveES$transcript_id[1:30]),
                                     y[y$associated_gene == "TRIO" & y$structural_category == "FSM","isoform"][1:10],
                                     c("ONT5_1471_4273","ONT5_1471_4815","ONT5_1471_4942","ONT5_1471_1636"))))

#head(y[y$associated_gene == "TRIO",] %>% arrange(-length))

#x <- x[x$transcript_id %in% c(unique(x$transcript_id)[1:20],"ONT16_2235_233","ONT16_2235_493"),]
#cpat<-vroom('/gpfs/mrc0/projects/Research_Project-MRC148213/sl693/RBFetal/2_cpat_tc20bp/WholeTargeted_fixed.ORF_prob.best.tsv'))
#x<-read_gff('~/isoform_analysis/data/examples/ANKRD12.gff')
#x<-as.data.frame(x)
#x<-x[which(x$transcript_id%in%y$isoform),]
#x<-x[which(x$transcript_id%in%cpat[which(cpat$Coding_prob>0.364),]$seq_ID),]
#x<-x[which(x$transcript_id%in%tmp),]

p2 <- plot_CDS(x)
plot_CDS <- function(x){
  x_exons <- x %>% dplyr::filter(type == "exon")
  head(x_exons)
  
  # extract cds
  x_cds <- x %>% dplyr::filter(type == "CDS")
  head(x_cds)
  
  x_cds_w_stop <- x_cds %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      end = ifelse(end == max(end), end + 3, end)
    ) %>%
    dplyr::ungroup()
  
  # add_utr() adds ranges that represent the UTRs
  x_cds_utr <- add_utr(
    x_exons,
    x_cds_w_stop,
    group_var = "transcript_id"
  )
  x_cds_utr_rescaled <-
    shorten_gaps(
      exons = x_cds_utr,
      introns = to_intron(x_cds_utr, "transcript_id"),
      group_var = "transcript_id"
    )
  
  z<-read_gff('/lustre/projects/Research_Project-MRC190311/longReadSeq/ONTRNA/SFARIdevelopmentalgenomics/0_output/gencode.v44.annotation.TRIO.gff3')
  z<-as.data.frame(z)
  z<-z[which(z$gene_name=='TRIO'),]
  
  z_exons <- z %>% dplyr::filter(type == "exon")
  z_exons %>% head()
  
  # extract cds
  z_cds <- z %>% dplyr::filter(type == "CDS")
  z_cds %>% head()
  
  z_cds_w_stop <- z_cds %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(
      end = ifelse(end == max(end), end + 3, end)
    ) %>%
    dplyr::ungroup()
  
  # add_utr() adds ranges that represent the UTRs
  z_cds_utr <- add_utr(
    z_exons,
    z_cds_w_stop,
    group_var = "transcript_id"
  )
  z_cds_utr_rescaled <-
    shorten_gaps(
      exons = z_cds_utr,
      introns = to_intron(z_cds_utr, "transcript_id"),
      group_var = "transcript_id"
    )
  
  
  x<-rbind(x, z[colnames(x)]) 
  
  p <- x_cds_utr_rescaled %>%
    dplyr::filter(type == "CDS") %>%
    ggplot(., aes(
      xstart = start,
      xend = end,
      y = transcript_id
    )) +
    geom_range(fill='blue') +
    geom_range(
      data = x_cds_utr_rescaled %>% dplyr::filter(type == "UTR"),
      height = 0.25,
      fill = "white"
    ) +
    geom_intron(
      data = to_intron(
        x_cds_utr_rescaled %>% dplyr::filter(type != "intron"),
        "transcript_id"
      ),
      arrow.min.intron.length = 110
    ) + 
    geom_range(data=z_cds_utr_rescaled %>% dplyr::filter(type == "CDS"),
               aes(
                 xstart = start,
                 xend = end,
                 y = transcript_id
               ),fill='blue') +
    geom_range(
      data = z_cds_utr_rescaled %>% dplyr::filter(type == "UTR"),
      height = 0.25,
      fill = "white"
    ) +
    geom_intron(
      data = to_intron(
        z_cds_utr_rescaled %>% dplyr::filter(type != "intron"),
        "transcript_id"
      ),
      arrow.min.intron.length = 110
    ) + theme_classic() + labs(x=NULL,y="Transcripts")
  
  
  return(p)
}
