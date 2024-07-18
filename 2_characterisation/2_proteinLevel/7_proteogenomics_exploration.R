library("ggplot2")
library("cowplot")


message("Number of RNA transcripts: ", length(unique(input$t2p.collapse$pb_accs)))
message("Number of RNA isoforms: ", length(unique(input$t2p.collapse$base_acc)))

input$t2p.collapse %>% arrange(-numtxCollapsed) %>% select(base_acc, numtxCollapsed) %>% distinct() %>% 
  ggplot(., aes(x = numtxCollapsed)) + geom_histogram() +
  scale_x_log10() 

p2 <- input$t.class.files %>% group_by(associated_gene) %>% tally() %>% 
  ggplot(., aes(x = n)) + geom_histogram() + 
  scale_x_log10() + theme_classic() + 
  labs(x ="Number of transcripts", y = "Number of genes")

p3 <- input$t2p.collapse %>% arrange(-numtxCollapsed) %>% select(base_acc, numtxCollapsed) %>% distinct() %>% 
  left_join(input$t.class.files %>% select(isoform, associated_gene), by = c("base_acc" = "isoform")) %>%
  group_by(associated_gene) %>% tally() %>% 
  ggplot(., aes(x = n)) + geom_histogram() + 
  scale_x_log10() + theme_classic() + 
  labs(x ="Number of unique RNA isoforms", y = "Number of genes")

plot_grid(p2,p3)

# features of RNA isoforms that are uniquely (i.e. n = 1) associated to gene
uniqueIsoGene <- input$t2p.collapse %>% arrange(-numtxCollapsed) %>% select(base_acc, numtxCollapsed) %>% distinct() %>% 
  left_join(input$t.class.files %>% select(isoform, associated_gene), by = c("base_acc" = "isoform")) %>%
  group_by(associated_gene) %>% tally() %>% filter(n == 1) 

p4 <- input$t2p.collapse %>% arrange(-numtxCollapsed) %>% select(base_acc, numtxCollapsed) %>% distinct() %>% 
  left_join(input$t.class.files %>% select(isoform, structural_category, associated_gene), by = c("base_acc" = "isoform")) %>% 
  filter(associated_gene %in% uniqueIsoGene$associated_gene) %>% 
  group_by(structural_category) %>% tally() %>% 
  ggplot(., aes(x = structural_category, y = n)) + geom_bar(stat = "identity") + 
  theme_classic() + labs(x = "Structural Category", y = "Number of RNA isoforms")


totalT <- input$t.class.files %>% group_by(associated_gene) %>% tally() %>%
  `colnames<-`(c("associated_gene", "totalTranscripts"))

totalR <- input$t2p.collapse %>% arrange(-numtxCollapsed) %>% select(base_acc, numtxCollapsed) %>% distinct() %>% 
  left_join(input$t.class.files %>% select(isoform, associated_gene), by = c("base_acc" = "isoform")) %>%
  group_by(associated_gene) %>% tally() %>% 
  `colnames<-`(c("associated_gene", "totalRNAIsoforms"))

merge(totalR,totalT) %>% ggplot(., aes(x = totalRNAIsoforms, y = totalTranscripts)) + geom_point() +
  theme_classic() + labs(x = "Number of RNA isoforms in gene", y = "Number of transcripts in gene") + 
  ylim(0,4000) + xlim(0,4000)



dat <- input$t2p.collapse %>% arrange(-numtxCollapsed) %>% select(base_acc, numtxCollapsed) %>% distinct() %>% 
  left_join(input$t.class.files %>% select(isoform, structural_category, associated_gene), by = c("base_acc" = "isoform"))

dat %>% select(associated_gene, num) 

collapsedTranscripts <- setdiff(input$t2p.collapse$pb_accs,input$t2p.collapse$base_acc)
input$t.class.files[input$t.class.files$isoform %in% collapsedTranscripts, ] %>% group_by(structural_category) %>% tally() %>% 
  ggplot(., aes(x = structural_category, y = n)) + geom_bar(stat = "identity") +scale_y_continuous(labels = label_comma()) +
  theme_classic() + labs(x = "Structural category", y = "Number of transcripts collapsed")
