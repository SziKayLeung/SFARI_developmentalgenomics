## ---------- Proteogenomics ----------

message("Number of RNA Transcripts to all genes:", nrow(class.files$glob_SQ))
message("Number of coding RNA isoforms:", length(unique(proteinInput$t2p.collapse.refined$pb_acc)))
message("Number of coding RNA isoforms:", length(unique(proteinInput$t2p.collapse.refined$corrected_acc)))
length(unique(class.files$ptarg_filtered$corrected_acc)) # -1 to not include "NA" 

setdiff(class.files$glob_SQ$isoform, c(as.character(proteinInput$no_orf$V1),as.character(proteinInput$cpat$seq_ID)))
# these are the isoforms that have been filtered from monoexonic isoforms from multiexonic genes
setdiff(c(as.character(proteinInput$no_orf$V1),as.character(proteinInput$cpat$seq_ID)),class.files$glob_SQ$isoform)

message("Number of transcripts with noORF: ", length(as.character(proteinInput$no_orf$V1)))
message("Number of transcripts with ORF: ", length(as.character(proteinInput$cpat$seq_ID)))
Rank1 <- proteinInput$mapped[proteinInput$mapped$orf_rank == "1",]
message("Number of transcripts with ORF ranked 1 and has stop codons: ", length(unique(Rank1[Rank1$has_stop_codon == "TRUE","transcript_id"])))
message("Number of transcripts with ORF ranked 1 and has stop codons, and coding potential > 0: ", length(unique(proteinInput$t2p.collapse.refined$pb_accs)))

RefinedORF <- proteinInput$cpat[proteinInput$cpat$seq_ID %in% proteinInput$t2p.collapse.refined$pb_accs,]
message("Number of transcripts with ORF ranked 1 and has stop codons, and coding potential > 0.364: ", nrow(RefinedORF[RefinedORF$Coding_prob > 0.364,]))
nrow(RefinedORF[RefinedORF$Coding_prob > 0.364,])/nrow(class.files$glob_SQ)

RefinedORFCoding <- RefinedORF[RefinedORF$Coding_prob > 0.364,]
message("Number of transcripts with ORF ranked 1 and has stop codons, and coding potential > 0.364, collapsed: ",
        length(unique(proteinInput$t2p.collapse.refined[proteinInput$t2p.collapse.refined$pb_accs %in% RefinedORFCoding$seq_ID,"corrected_acc"])))
590504/nrow(RefinedORF[RefinedORF$Coding_prob > 0.364,])

GSselectedcollapsedID <- proteinInput$t2p.collapse.refined[proteinInput$t2p.collapse.refined$pb_accs %in% RefinedORFCoding$seq_ID,"base_acc"]
CorrectedcollapsedID <- proteinInput$t2p.collapse.refined[proteinInput$t2p.collapse.refined$pb_accs %in% RefinedORFCoding$seq_ID,"corrected_acc"]

# note one transcript was filtered from SQANTI protein so off by 1 when calculating difference
message("Number of transcripts with ORF...collapsed: ", length(unique(setdiff(RefinedORFCoding$seq_ID, CorrectedcollapsedID))))
ORFCorrectedCollapsedID <- unique(setdiff(RefinedORFCoding$seq_ID, CorrectedcollapsedID))
class.files$targ_filtered %>% filter(isoform %in% ORFCorrectedCollapsedID) %>% group_by(structural_category) %>% tally()

class.files$protein_filtered_final <- class.files$protein_filtered[class.files$protein_filtered$pb %in% GSselectedcollapsedID,]
message("Number of protein products after filtering: ", nrow(class.files$protein_filtered_final))
