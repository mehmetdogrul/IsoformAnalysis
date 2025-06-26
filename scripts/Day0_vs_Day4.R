library(IsoformSwitchAnalyzeR)
library(readr)
library(ggplot2)
library(dplyr)
library(fgsea)
library(GSEABase)
library(tidyr)
library(biomaRt)

setwd("Documents/projects/rnaseq/scripts/Day0_Day4_analysis/")

### Loading the data
parentDir <- "../../stringtie_results/Day0_vs_Day4/"
conditionTable <- read.csv("../../stringtie_results/Day0_vs_Day4/sample_info.csv") # edit this to include ALL the comparasions, and you have a complete analysis!

stringtieQuant <- importIsoformExpression(
  parentDir, addIsofomIdAsColumn = TRUE, pattern = "Rep", readLength = 150
)

switchList <- importRdata(
  isoformCountMatrix = stringtieQuant$counts,
  isoformRepExpression = stringtieQuant$abundance,
  designMatrix = conditionTable,
  isoformExonAnnoation = "../../stringtie_results/Day0_vs_Day4/merged.gtf",
  showProgress = TRUE
)

### Filtering
switchListFiltered <- preFilter(
  switchAnalyzeRlist = switchList,
  geneExpressionCutoff = 1, # default
  isoformExpressionCutoff = 0, # default
  removeSingleIsoformGenes = TRUE # default
)
# Output: The filtering removed 53548 ( 56.03% of ) transcripts. There is now 42029 isoforms left

switchListFilteredStrict <- preFilter(
  switchAnalyzeRlist = switchList,
  geneExpressionCutoff = 10,
  isoformExpressionCutoff = 3,
  removeSingleIsoformGenes = TRUE # default
)
# Output: The filtering removed 89370 ( 93.51% of ) transcripts. There is now 6207 isoforms left

### Performing test for differential isoform usage with DEXSeq
switchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = switchListFiltered,
  reduceToSwitchingGenes = TRUE,
  showProgress = TRUE,
  alpha = 0.01
)

summary(switchListAnalyzed)

switchListAnalyzedStrict <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = switchListFilteredStrict,
  reduceToSwitchingGenes = TRUE,
  showProgress = TRUE
)

### Summarize switching features
extractSwitchSummary(switchListAnalyzed)
extractSwitchSummary(switchListAnalyzedStrict)
summaryResult <- summary(switchListAnalyzed)
# Question: Why is nrSwihces bigger than nrIsoforms?

### Predict alternative splicing
switchListAS <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = switchListAnalyzed,
  onlySwitchingGenes = FALSE # does not change anything due to previous restrictions
)

switchListASStrict <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = switchListAnalyzedStrict,
  onlySwitchingGenes = FALSE
)

table(switchListAS$AlternativeSplicingAnalysis$IR) 
table(switchListASStrict$AlternativeSplicingAnalysis$IR) 
switchListAS

### Genome-Wide analysis of alternative splicing
summary <- extractSplicingSummary(
  switchListAS,
  asFractionTotal = FALSE, # default
  plotGenes = FALSE # default
)

extractSplicingSummary(
  switchListASStrict,
  asFractionTotal = FALSE, # default
  plotGenes = FALSE # default
)

splicingEnrichment <- extractSplicingEnrichment(
  switchListAS,
  splicingToAnalyze = 'all',
  returnResult = FALSE,
  returnSummary = TRUE
)

splicingEnrichmentStrict <- extractSplicingEnrichment(
  switchListASStrict,
  splicingToAnalyze = 'all',
  returnResult = FALSE,
  returnSummary = TRUE
  
)

# check the parameters and save the plots

extractSplicingEnrichmentComparison(
  # do this after running having completed the whole pipeline
)

genomewide <- extractSplicingGenomeWide(
  switchListAS,
  featureToExtract = 'all',
  splicingToAnalyze = 'all',
  plot=TRUE,
  returnResult = FALSE,
  violinPlot = TRUE
)

### Volcano like plot:
volcano <- ggplot(data=switchListAnalyzed$isoformFeatures, aes(x=(dIF), y=-log10(isoform_switch_q_value))) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ), # default cutoff
    size=1
  ) +
  geom_hline(yintercept = -log10(0.01), linetype='dashed') + # default cutoff
  geom_vline(xintercept = c(-0.1, 0.1), linetype='dashed') + # default cutoff
  #facet_wrap( ~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='dIF', y='-Log10 ( Isoform Switch Q Value )') +
  theme_bw()

### Switch vs Gene changes:
hamburger <- ggplot(data=switchListAnalyzed$isoformFeatures, aes(x=gene_log2_fold_change, y=dIF)) +
  geom_point(
    aes( color=abs(dIF) > 0.1 & isoform_switch_q_value < 0.05), # default cutoff
    size=1
  ) + 
  #facet_wrap(~ condition_2) +
  #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
  geom_hline(yintercept = 0, linetype='dashed') +
  geom_vline(xintercept = 0, linetype='dashed') +
  scale_color_manual('Signficant\nIsoform Switch', values = c('black','red')) +
  labs(x='Gene log2 fold change', y='dIF') +
  theme_bw()

###   Save the results
dataDir <- "../../data_analysis/data/"
plotsDİr <- "../../data_analysis/plots/"

objectsToSave <- ls()
filePath <- file.path(dataDir, "Day0_vs_Day4.RData")
save(list = objectsToSave, file = filePath)
write.csv(summaryResult, "../../data_analysis/data/summary.csv", row.names = FALSE)


###   Save the plots
###   Create directories with sample_sorter.sh
ggsave("summary.png", summary, "png", "../../data_analysis/plots/Day0_vs_Day4/", create.dir = TRUE)
ggsave("enrich.png", splicingEnrichment, "png", "../../data_analysis/plots/Day0_vs_Day4/", create.dir = TRUE)
ggsave("genomewide.png", genomewide, "png", "../../data_analysis/plots/Day0_vs_Day4/", create.dir = TRUE)
ggsave("volcano.png", volcano, "png", "../../data_analysis/plots/Day0_vs_Day4/", create.dir = TRUE)
ggsave("hamburger.png", hamburger, "png", "../../data_analysis/plots/Day0_vs_Day4/", create.dir = TRUE)

###   pseudocode for gene-list
###   select conseqeunce
###   select significant
###   select unique
###
###   take genefold into account as:
###     repeat above
###     select the ones abs(x) < 0.1

###   form gene lists
###     for this, there are two problems: 1) I should form a distinct dataset for GSEA,
###     2) I have to pay attention to make sure that the samples are well aligned. They are disorganised between lists.
###     pseudocode
###       Pick ATTS.
###       Choose signifacnt ones.
###       Import whether up or down
###       Merge sets
###       Check out dataset
###       Filter, distinct
###       Check again
###       Divede into to sets as up and down
###       Along this process, match each variable well to avoid any confusion. Use of ID is crutial here. arrange() will slove this hopefully.

consequence <- switchListAS$AlternativeSplicingAnalysis %>% 
  select(isoform_id, ATTS) %>% 
  arrange(isoform_id)

significant <- switchListAS$isoformSwitchAnalysis %>% 
  select(isoform_id, padj, dIF) %>% 
  arrange(isoform_id)

geneNames <- switchListAS$isoformFeatures %>% 
  select(isoform_id,gene_ref,gene_id,gene_name, gene_log2_fold_change) %>% 
  arrange(isoform_id)

table(geneNames$isoform_id == significant$isoform_id) # TRUE
table(geneNames$isoform_id == consequence$isoform_id) # TRUE

significant_consequence_genes <- cbind(geneNames, consequence, significant)[,c(-6,-8)]

significant_consequence_genes.filtered <- significant_consequence_genes %>%
  filter(ATTS > 0,  padj< 5e-2) %>%
  mutate(
    regulation = case_when(
      dIF > 0.1 ~ "up",
      dIF < -0.1 ~ "down",
      TRUE ~ "not_significant")) %>% 
  select(isoform_id,gene_name,regulation, dIF, gene_log2_fold_change) %>% 
  distinct() %>% 
  na.omit()

atts_up <- significant_consequence_genes.filtered %>% 
  filter(regulation == "up") %>% 
  select(gene_name) %>% 
  distinct() %>% 
  arrange()

atts_down <- significant_consequence_genes.filtered %>% 
  filter(regulation == "down") %>% 
  select(gene_name) %>% 
  distinct() %>% 
  arrange()

write_csv(atts_up, "../../data_analysis/data/ATTS_up.csv", col_names = FALSE)
write_csv(atts_down, "../../data_analysis/data/ATTS_down.csv", col_names = FALSE)
###   organise the dvaid output

go_atts_up <- corrector(read.delim("../../data_analysis/data/ATTS_up_GO_annotation.txt"))
go_atts_down <- corrector(read.delim("../../data_analysis/data/ATTS_down_GO_annotation.txt"))

go_atts_up.filtered <- go_atts_up %>% 
  select(Term, Count, Genes, Benjamini) %>% 
  mutate(Term = str_replace(Term, ".*~", "")) %>% 
  filter(Benjamini < 5e-2)

go_atts_down.filtered <- go_atts_down %>% 
  select(Term, Count, Genes, Benjamini) %>% 
  mutate(Term = str_replace(Term, ".*~", "")) %>% 
  filter(Benjamini < 5e-2)

###   KEGG

kegg_atts_up <- corrector2(read.delim("../../data_analysis/data/ATTS_up_KEGG.txt"))

kegg_atts_down <- corrector2(read.delim("../../data_analysis/data/ATTS_down_KEGG.txt"))


###   GSEA

gmt_file <- "../../data_analysis/data/m5.go.v2024.1.Mm.symbols.gmt" 
gsea_gene_set <- getGmt(gmt_file)
pathways <- geneIds(gsea_gene_set)

atts_gsea_data <- significant_consequence_genes.filtered %>% 
  select(gene_name,gene_log2_fold_change) %>% 
  arrange(desc(gene_log2_fold_change)) %>% 
  distinct()

write_tsv(ATTS_gsea_data, "../../data_analysis/data/ATTS_gsea.tsv", col_names = FALSE)  
stats <- ATTS_gsea_data$gene_log2_fold_change
names(stats) <- ATTS_gsea_data$gene_name
fgsea_results <- fgsea(pathways, stats)
plotEnrichment(pathways[3557], stats) +
  labs(title = names(pathways)[3557])


###   DAVID with all isoforms

###   Form a gene list dataset
###   Pseudocode
###   Use switchAnalayzed
###   Select gene names, q values, dIF
###   Filter signifiacant ones
###   Check out if everything is OK
###   Divide into as down and up
###   Select unique ones
###   Check out if everything is OK

isoforms_david <- switchListAnalyzed$isoformFeatures %>% 
  select(gene_name, isoform_switch_q_value, dIF) %>%
  drop_na() %>% 
  filter(isoform_switch_q_value < 0.05) %>%
  mutate(
    regulation = case_when(
      dIF > 0.1 ~ "up",
      dIF < -0.1 ~ "down",
      TRUE ~ "not_significant")) %>% 
  arrange(gene_name)

isoforms_david_up <- isoforms_david %>% 
  filter(regulation == "up") %>% 
  select(gene_name) %>% 
  distinct() %>% 
  arrange(gene_name)

isoforms_david_down <- isoforms_david %>% 
  filter(regulation == "down") %>% 
  select(gene_name) %>% 
  distinct() %>% 
  arrange(gene_name)

write_csv(isoforms_david_up, "../../data_analysis/data/isoforms_david_up.csv", col_names = FALSE)
write_csv(isoforms_david_down, "../../data_analysis/data/isoforms_david_down.csv", col_names = FALSE)

go_isoforms_up <- corrector(read.delim("../../data_analysis/data/isoforms_up_go.txt"))
go_isoforms_down <- corrector(read.delim("../../data_analysis/data/isoforms_down_go.txt"))
# kegg_isoforms_up <- corrector2(read.delim("../../data_analysis/data/isoforms_up_kegg.txt"))
kegg_isoforms_down <- corrector2(read.delim("../../data_analysis/data/isoforms_down_kegg.txt"))

go_isoforms_up.filtered <- go_isoforms_up %>% 
  select(Term, Count, Genes, Benjamini) %>% 
  mutate(Term = str_replace(Term, ".*~", "")) %>% 
  filter(Benjamini < 5e-2)

go_isoforms_down.filtered <- go_isoforms_down %>% 
  select(Term, Count, Genes, Benjamini) %>% 
  mutate(Term = str_replace(Term, ".*~", "")) %>% 
  filter(Benjamini < 5e-2)



###   Manually formm gene lists

corrector <- function (data){
  df <- data %>% 
    select(Term, Count, Genes, Benjamini) %>% 
    mutate(Term = str_replace(Term, ".*~", ""))
  df
}

corrector2 <- function (data){
  df <- data %>% 
    select(Term, Count, Genes, Benjamini) %>% 
    mutate(Term = str_replace(Term, ".*:", ""))
  df
}



# genelister <- function (data, row_index){
#   df <- data[1:row_index,] %>% 
#     select(Term, Count, Genes, Benjamini) %>% 
#     mutate(Term = str_replace(Term, ".*~", ""))
#   df
# }


###   Output for David analysis: Signifiant genes
# write.csv(go_atts_down, "../../data_analysis/data/go_atts_down_final.csv")
# write.csv(go_atts_up, "../../data_analysis/data/go_atts_up_final.csv")
# write.csv(go_isoforms_down, "../../data_analysis/data/go_isoforms_down_final.csv")
# write.csv(go_isoforms_up, "../../data_analysis/data/go_isoforms_up_final.csv")
# write.csv(kegg_isoforms_down[1:5,], "../../data_analysis/data/kegg_isoforms_down_final.csv")



###   Creating a pathway vector so that we can colorize gene ontologies etc. in the hamburger plot
atts <- list(
  go_down = go_atts_down.filtered,
  go_up = go_atts_up.filtered
)

all_isoforms <- list(
  go_down = go_isoforms_down.filtered,
  go_up = go_isoforms_up.filtered,
  kegg_down = kegg_isoforms_down
)

david <- list(
  atts = atts,
  all_isoforms = all_isoforms
)

###   Convert human notation to mouse notation for david results using bioaRt

# human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://uswest.ensembl.org")
# mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host = "https://asia.ensembl.org")
# 
# human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://feb2024.archive.ensembl.org")
# mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://feb2024.archive.ensembl.org")
# 
# 
# gene_converter <- function (david_pointer, term_index){
#   genes <- david_pointer %>% 
#     filter(Term == Term[term_index])
#   genes <- genes$Genes
#   gene_vec <- strsplit(genes, ",\\s*")[[1]]
#   # print(strsplit(genes, ",\\s*")[[1]])
#   # print(genes)
#   getLDS(
#     attributes = c("hgnc_symbol"),
#     filters = "hgnc_symbol",
#     values = gene_vec,
#     mart = human_mart,
#     attributesL = c("mgi_symbol", "entrezgene_id"),
#     martL = mouse_mart,
#     uniqueRows = TRUE
#   )
# }
# gene_converter(david$all_isoforms$go_up, 1)


###   Conver gene names through string manuplation: uppercase to title

gene_converter <- function (david_pointer, term_index){
  genes <- david_pointer %>% 
    filter(Term == Term[term_index])
  genes <- genes$Genes
  genes_vector <- strsplit(genes, ",\\s*")[[1]]
  converted <- str_to_title(tolower(genes_vector))
  return(converted)
  
}
gene_converter(, 1)

hamburger_david_plotter <- function(david_pointer, regulation, term_index){
  
  df <- switchListAnalyzed$isoformFeatures
  genes <- gene_converter(david_pointer, term_index)
  
  if (regulation == "up"){
    df %>% mutate(
      pathway = ifelse(
        dIF > 0.1 & gene_name %in% genes, "pathway", 
        ifelse(
          abs(dIF) > 0.1 && isoform_switch_q_value < 0.05, "significant", "not significant")
        )
    )
  } 
  else if (regulation == "down"){
    df %>% mutate(
      pathway = ifelse(
        dIF > 0.1 & gene_name %in% genes, "pathway", 
        ifelse(
          abs(dIF) > 0.1 && isoform_switch_q_value < 0.05, "significant", "not significant")
        )
    )
  }
  else {
    stop("Invalid regulation argument")
  }
  
  ggplot(data=df, aes(x=gene_log2_fold_change, y=dIF)) +
    geom_point(
      aes(color=pathway, # default cutoff
      size=1)
      )+ 
    #facet_wrap(~ condition_2) +
    #facet_grid(condition_1 ~ condition_2) + # alternative to facet_wrap if you have overlapping conditions
    geom_hline(yintercept = 0, linetype='dashed') +
    geom_vline(xintercept = 0, linetype='dashed') +
    scale_color_manual(
        name = "Gene Category",
        values = c(
          "Not Significant" = "black",
          "Significant" = "red",
          "Pathway" = "blue"
          )
        ) +
    labs(x='Gene log2 fold change', y='dIF') +
    theme_bw()
}
hamburger_david_plotter(david$all_isoforms$go_down, "down", 1)

###   chatgpt 
hamburger_david_plotter <- function(david_pointer, regulation, term_index, gene_fold = FALSE) {
  df <- switchListAnalyzed$isoformFeatures
  genes <- gene_converter(david_pointer, term_index)
  
  # Flags
  is_pathway <- if (regulation == "up") {
    df$dIF > 0.1 & df$gene_name %in% genes
  } else if (regulation == "down") {
    df$dIF < -0.1 & df$gene_name %in% genes
  } else {
    stop("Invalid regulation input")
  }
  
  is_2d_sig <- df$gene_switch_q_value < 0.05 & abs(df$dIF) > 0.1 & df$isoform_switch_q_value < 0.05
  is_significant <- df$gene_switch_q_value > 0.05 & abs(df$dIF) > 0.1 & df$isoform_switch_q_value < 0.05
  
  # Assign category with correct priority
  df <- df %>%
    mutate(
      category = case_when(
        is_pathway ~ "Pathway",
        is_2d_sig ~ "2D Significant",
        is_significant ~ "Significant",
        TRUE ~ "Not Significant"
      )
    )
  
  # Separate for layered plotting
  df_pathway <- df %>% filter(category == "Pathway")
  df_other <- df %>% filter(category != "Pathway")
  
  ggplot() +
    # Base layer
    geom_point(data = df_other,
               aes(x = gene_log2_fold_change, y = dIF, color = category),
               size = 1.2, shape = 16) +
    
    # Pathway on top
    geom_point(data = df_pathway,
               aes(x = gene_log2_fold_change, y = dIF, color = category),
               size = 3, shape = 17) +
    
    # Threshold lines
    geom_hline(yintercept = c(-0.1, 0.1), linetype = 'dotdash', color = "gray40") +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = 'dotdash', color = "gray40") +
    
    # Colors
    scale_color_manual(
      name = "Gene Category",
      values = c(
        "Not Significant" = "#B0B0B0",
        "Significant" = "#D7263D",
        "2D Significant" = "#222222",
        "Pathway" = "#1E88E5"
      )
    ) +
    
    labs(
      x = 'Gene log2 fold change',
      y = 'dIF',
      title = david_pointer$Term[term_index]
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

  # Generate the plotrger_
hamburger_david_plotter(david$all_isoforms$go_up, "up", term_index = 5, gene_fold = FALSE)

###   remake of david hamburger


fn <- function(david_pointer, regulation, term_index) {
  genes <- gene_converter(david_pointer, term_index)
  df <- switchListAnalyzed$isoformFeatures
  
  df <- df %>%
    mutate(
      category = case_when(
        regulation == "up" & dIF > 0.1 & gene_name %in% genes & isoform_switch_q_value < 0.05 ~ "Pathway",
        regulation == "down" & dIF < -0.1 & gene_name %in% genes & isoform_switch_q_value < 0.05 ~ "Pathway",
        abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 & abs(gene_log2_fold_change) > 0.1 & gene_switch_q_value < 0.05 ~ "2D significant",
        abs(dIF) > 0.1 & isoform_switch_q_value < 0.05 ~ "Significant",
        TRUE ~ "Not significant"
      )
    )
  
  df_pathway <- df %>% filter(category == "Pathway")
  df_other <- df %>% filter(category != "Pathway")
  
  ggplot() +
    geom_point(data = df_other,
               aes(x = gene_log2_fold_change, y = dIF, color = category),
               size = 0.75, shape = 16) +
    geom_point(data = df_pathway,
               aes(x = gene_log2_fold_change, y = dIF, color = category),
               size = 2, shape = 17) +
    geom_hline(yintercept = c(-0.1, 0.1), linetype = 'dotdash', color = "gray40") +
    geom_vline(xintercept = c(-0.1, 0.1), linetype = 'dotdash', color = "gray40") +
    scale_color_manual(
      name = "Gene Category",
      values = c(
        "Not significant" = "#B0B0B0",
        "Significant" = "#D7263D",
        "2D significant" = "#414",
        "Pathway" = "#1E88E5"
      )
    ) +
    labs(
      x = 'Gene log2 fold change',
      y = 'dIF',
      title = david_pointer$Term[term_index]
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
  
}

fn(david$all_isoforms$go_down, "down", 9)

###   It looks like all the genes with abs(genefold) > 0.1 are 2D significant. Check if this is true
######    pseudocode
### create a new df
### select abs.genefold > 0.1         
### select abs.dIF > 0.1
### select gene_q_Value < 0.05
### select isoform_q_Value < 0.05
### create another df
### select abs.genefold > 0.1         
### select abs.dIF > 0.1
### select isoform_q_Value < 0.05

two_d_sig <- switchListAnalyzed$isoformFeatures %>% 
  filter(
    abs(dIF) > 0.1,
    isoform_switch_q_value < 0.05,
    abs(gene_log2_fold_change) > 0.1,
    gene_switch_q_value < 0.05
    )

two_d_sig.alt <- switchListAnalyzed$isoformFeatures %>% 
  filter(
    abs(dIF) > 0.1,
    isoform_switch_q_value < 0.05,
    abs(gene_log2_fold_change) > 0.1
  )


