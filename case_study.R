#################################################################################################################
######################### STATISTICAL META-ANALYSIS #############################################################
#################################################################################################################

##### Run ensemble meta-analysis #####

#### adj_pval_integrated_biopsy.txt is the output of a function adj_p_val_integration
integrated_adj_pvals<-read.table("case_study/adj_pval_integrated_biopsy.txt")

class<-rep(1, ncol(integrated_adj_pvals))

origin <- sub("_(adj_pval)$", "", colnames(integrated_adj_pvals))

ranked_gene_list<-run_ensembl_metanalysis(meta_dataframe = integrated_adj_pvals,class = class, origin = origin, 
                                          method = c("effect_size", "pvalue", "rank_product"), metric = "median")

write.table(ranked_gene_list, file="case_study/meta_analysis_ranked_gene_list.txt")

#####GSEA#####

gmt<-"case_study/c2.cp.reactome.v2023.2.Hs.symbols.gmt"

set.seed(42)

gsea<-compute_gsea(gene_list = ranked_gene_list, gmt_file = gmt)

set.seed(NULL)

gsea_ordered<-as.data.frame(gsea[order(gsea$pval, decreasing = F),])

gsea_ordered$leadingEdge <- sapply(gsea_ordered$leadingEdge, paste, collapse = ",")

write.table(gsea_ordered, file="case_study/gsea_results_ordered.txt", sep="\t")

##### Compute GSEA threshold #####

threshold<-compute_gsea_thresh(geneList = ranked_gene_list, fgsea_res = gsea, background = gmt)

write.table(threshold, file="case_study/gsea_threshold.txt", sep="\t")

#################################################################################################################
######################### NETWORK ANALYSIS #####################################################################
#################################################################################################################

##### Batch_adjustment_for_the_networks #####

##### Biopsy RNA-Seq disease #####

biopsy_rnaseq_disease_combined_expr_mat<-read.table("case_study/combined_symbol_expression_matrix_disease_biopsy_rnaseq_subset.txt", sep="\t")
biopsy_rnaseq_disease_batch<-read.table("case_study/batch_biopsy_rnaseq_disease.txt", sep="\t")
biopsy_rnaseq_disease_batch<-biopsy_rnaseq_disease_batch$x
biopsy_rnaseq_disease_samples<-read.table("case_study/samples_biopsy_rnaseq_disease.txt", sep="\t")
biopsy_rnaseq_disease_samples<-biopsy_rnaseq_disease_samples$x

adjusted_matrix_biopsy_rnaseq_disease<-multi_studies_adjust(expr_mat = biopsy_rnaseq_disease_combined_expr_mat, samples_label = biopsy_rnaseq_disease_samples, batch_labels = biopsy_rnaseq_disease_batch)
adjusted_matrix_biopsy_rnaseq_disease<-adjusted_matrix_biopsy_rnaseq_disease$x
write.table(adjusted_matrix_biopsy_rnaseq_disease, file = "case_study/adjusted_matrix_biopsy_rnaseq_disease.txt", sep="\t")

##### Biopsy RNA-Seq healthy #####

biopsy_rnaseq_healthy_combined_expr_mat<-read.table("case_study/combined_symbol_expression_matrix_healthy_biopsy_rnaseq_subset.txt", sep="\t")
biopsy_rnaseq_healthy_batch<-read.table("case_study/batch_biopsy_rnaseq_healthy.txt", sep="\t")
biopsy_rnaseq_healthy_batch<-biopsy_rnaseq_healthy_batch$x
biopsy_rnaseq_healthy_samples<-read.table("case_study/samples_biopsy_rnaseq_healthy.txt", sep="\t")
biopsy_rnaseq_healthy_samples<-biopsy_rnaseq_healthy_samples$x

adjusted_matrix_biopsy_rnaseq_healthy<-multi_studies_adjust(expr_mat = biopsy_rnaseq_healthy_combined_expr_mat, samples_label = biopsy_rnaseq_healthy_samples, batch_labels = biopsy_rnaseq_healthy_batch)
adjusted_matrix_biopsy_rnaseq_healthy<-adjusted_matrix_biopsy_rnaseq_healthy$x
write.table(adjusted_matrix_biopsy_rnaseq_healthy, file = "case_study/adjusted_matrix_biopsy_rnaseq_healthy.txt", sep="\t")

##### Biopsy microarray disease #####

biopsy_microarray_disease_combined_expr_mat<-read.table("case_study/combined_symbol_expression_matrix_disease_biopsy_micro_subset.txt", sep="\t")
biopsy_microarray_disease_batch<-read.table("case_study/biopsy_microarray_disease/batch_biopsy_micro_disease.txt", sep="\t")
biopsy_microarray_disease_batch<-biopsy_microarray_disease_batch$x
biopsy_microarray_disease_samples<-read.table("case_study/samples_biopsy_micro_disease.txt", sep="\t")
biopsy_microarray_disease_samples<-biopsy_microarray_disease_samples$x

adjusted_matrix_biopsy_microarray_disease<-multi_studies_adjust(expr_mat = biopsy_microarray_disease_combined_expr_mat, samples_label = biopsy_microarray_disease_samples, batch_labels = biopsy_microarray_disease_batch)
adjusted_matrix_biopsy_microarray_disease<-adjusted_matrix_biopsy_microarray_disease$x
write.table(adjusted_matrix_biopsy_microarray_disease, file = "case_study/adjusted_matrix_biopsy_microarray_disease.txt", sep="\t")

##### Biopsy microarray healthy #####

biopsy_microarray_healthy_combined_expr_mat<-read.table("case_study/combined_symbol_expression_matrix_healthy_biopsy_micro_subset.txt", sep="\t")
biopsy_microarray_healthy_batch<-read.table("case_study/batch_biopsy_micro_healthy.txt", sep="\t")
biopsy_microarray_healthy_batch<-biopsy_microarray_healthy_batch$x
biopsy_microarray_healthy_samples<-read.table("case_study/samples_biopsy_micro_healthy.txt", sep="\t")
biopsy_microarray_healthy_samples<-biopsy_microarray_healthy_samples$x

adjusted_matrix_biopsy_microarray_healthy<-multi_studies_adjust(expr_mat = biopsy_microarray_healthy_combined_expr_mat, samples_label = biopsy_microarray_healthy_samples, batch_labels = biopsy_microarray_healthy_batch)
adjusted_matrix_biopsy_microarray_healthy<-adjusted_matrix_biopsy_microarray_healthy$x
write.table(adjusted_matrix_biopsy_microarray_healthy, file = "case_study/adjusted_matrix_biopsy_microarray_healthy.txt", sep="\t")

##### Network_inference #####

matrices<-list.files(pattern="adjusted", include.dirs = FALSE)

for (i in 1:length(matrices)) {
  
  generatematrices<-get_ranked_consensus_matrix(gx_table = read.table(matrices[i], sep="\t"), iMethods = c("clr"),
                                                iEst = c("pearson"),
                                                iDisc=c("none"), ncores = 30, debug_output = TRUE, updateProgress = TRUE)
  
  
  #Parse ranked matrix and get bin_mat and edge_rank
  # Get edge rank list and binary inference matrix from edge rank matrix computed by get_ranked_consensus_matrix().
  # parse_edge_rank_matrix parses the edge rank matrix created by using the internal function get_ranked_consensus_matrix_matrix() to get a ranked edge list and a binary matrix.
  
  rankMat.parsed=parse_edge_rank_matrix(edge_rank_matrix = generatematrices, edge_selection_strategy = "default",
                                        mat_weights = "rank", topN = 10, debug_output = TRUE, updateProgress = TRUE)
  
  conGraph <- get_iGraph(rankMat.parsed$bin_mat)
  
  modified_name <- paste0("case_study/network_", gsub("adjusted_matrix_", "", matrices[i]))
  saveRDS(conGraph, file=sub(".txt$", ".rds", modified_name))
  
}


##### Compute modules #####

files <- list.files(pattern = "network")


list_of_modules_biopsy<-list()
for (i in 1:length(files)) {
  rds<-readRDS(files[i])
  #rds<-rds[[1]]
  list_of_modules_biopsy[[i]]<-get_modules(rds, method = "walktrap")
}
print("done")


module_list<-paste0("modules_", gsub(pattern = "network_", replacement = "", files))
module_list<-gsub(".rds", "", module_list)

names(list_of_modules_biopsy)<-module_list

save(list_of_modules_biopsy, file = "case_study/list_of_modules_biopsy.RData")


##### Enrichment analysis for the RNA-seq networks and visualisation #####


load("case_study/list_of_modules_biopsy.RData")

#names(list_of_modules_biopsy)
#[1] "modules_biopsy_microarray_disease" "modules_biopsy_microarray_healthy" "modules_biopsy_rnaseq_disease"     "modules_biopsy_rnaseq_healthy"    

network_rnaseq_disease<-readRDS("case_study/network_biopsy_rnaseq_disease.rds")
network_rnaseq_healthy<-readRDS("case_study/network_biopsy_rnaseq_healthy.rds")
modules_rnaseq_disease<-list_of_modules_biopsy[[3]]
modules_rnaseq_healthy<-list_of_modules_biopsy[[4]]

reactome_disease<-get_reactome_from_modules(modules = modules_rnaseq_disease, geneID="SYMBOL", pval_cutoff = 0.05, 
                                            outPath = paste0("case_study", "/pathway_results_disease"), layout = "overall")

reactome_healthy<-get_reactome_from_modules(modules = modules_rnaseq_healthy, geneID="SYMBOL", pval_cutoff = 0.05, 
                                            outPath = paste0("case_study", "/pathway_results_healthy"), layout = "overall")


bubble_disease<-get_bubbleplot_from_pathways(modules = modules_rnaseq_disease, geneID = "SYMBOL")
png("case_study/bubble_plot_disease_rnaseq.png")

bubble_healthy<-get_bubbleplot_from_pathways(modules = modules_rnaseq_healthy, geneID = "SYMBOL")
png("case_study/bubble_plot_healthy_rnaseq.png")

##### ESEA for RNA-seq disease and healthy #####

#install.packages("ESEA")
#library(ESEA)

#source("/nasdata/afederico/metanalysis_package_case_study/ESEA_function.R")
load("case_study/mypathwayEdge_KEGG.RData")

#build ESEA background

##### KEGG #####

#the edge background needs to be a matrix with two columns. Each row of the matrix is an edge

#pathwayEdge.db
#class is a character
#element i-th is "pathway name, \t,source,\t, edge (node1|node2)"

# source<-'kegg'
# mypathwayEdge.db<-c()
# for(i in 1:length(prior_data_kegg)) {
#   df<-read.csv(paste0(dir,prior_data_kegg[i]))
#   index_to_keep<-which(grepl(pattern = 'KEGG 2021' , x = df$r.source) | grepl(pattern = 'KEGG' , x = df$r.source))
#   df_cut<-df[index_to_keep,]
#   x<-as.vector(df_cut$g.Ensembl_ID)
#   pairs<-t(combn(x, 2))
#   pairs_right_noted<-c(paste0(pairs[,1],'|',pairs[,2]))
#   pathway_entry<-paste0(df_cut$p.name[1],'\t',source,'\t')
#   for(j in 1:length(pairs_right_noted)){
#     pathway_entry<-paste0(pathway_entry,pairs_right_noted[j],'\t')
#   }
#   mypathwayEdge.db[i]<-pathway_entry
# }

# load ESEA background 

#mypathwayEdge.db <- 

esea_background <- c()

for(i in 1:length(mypathwayEdge.db)){
  print(i)
  db.entry <- unlist(strsplit(mypathwayEdge.db[i], split = "\t"))
  df_temp <- data.frame(node1_ens=sapply(strsplit(db.entry[3:length(db.entry)], split = "|", fixed = TRUE), "[", 1), node2_ens=sapply(strsplit(db.entry[3:length(db.entry)], split = "|", fixed = TRUE), "[", 2), node1_sym = NA, node2_sym = NA)
  conv_table_n1 <- bitr(geneID = df_temp$node1_ens, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
  conv_table_n2 <- bitr(geneID = df_temp$node2_ens, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
  df_temp$node1_sym <- conv_table_n1$SYMBOL[match(df_temp$node1_ens, conv_table_n1$ENSEMBL)]
  df_temp$node2_sym <- conv_table_n2$SYMBOL[match(df_temp$node2_ens, conv_table_n2$ENSEMBL)]
  df_temp$edge <- paste(df_temp$node1_sym, df_temp$node2_sym, sep = "|")
  edges_recoded <- paste(df_temp$edge, collapse = "\t")
  pathway_entry <- paste(db.entry[1], db.entry[2], edges_recoded, sep = "\t")
  pathway_entry <- paste(pathway_entry, "\t", sep = "")
  esea_background[i] <- pathway_entry
  #Sys.sleep(5)
}


edge_list_biopsy_microarray_healthy <- readRDS("case_study/edge_rank_biopsy_microarray_healthy.rds")
edge_list_biopsy_rnaseq_healthy <- readRDS("case_study/edge_rank_biopsy_rnaseq_healthy.rds")
edge_list_biopsy_microarray_disease <- readRDS("case_study/edge_rank_biopsy_microarray_disease.rds")
edge_list_biopsy_rnaseq_disease <- readRDS("case_study/edge_rank_biopsy_rnaseq_disease.rds")

edge_list_biopsy_microarray_healthy <- gsub(edge_list_biopsy_microarray_healthy, pattern = ";", replacement = "|")
edge_list_biopsy_rnaseq_healthy <- gsub(edge_list_biopsy_rnaseq_healthy, pattern = ";", replacement = "|")
edge_list_biopsy_microarray_disease <- gsub(edge_list_biopsy_microarray_disease, pattern = ";", replacement = "|")
edge_list_biopsy_rnaseq_disease <- gsub(edge_list_biopsy_rnaseq_disease, pattern = ";", replacement = "|")

esea_vec_microarray_healthy <- c(length(edge_list_biopsy_microarray_healthy):1)
names(esea_vec_microarray_healthy) <- edge_list_biopsy_microarray_healthy
esea_vec_rnaseq_healthy <- c(length(edge_list_biopsy_rnaseq_healthy):1)
names(esea_vec_rnaseq_healthy) <- edge_list_biopsy_rnaseq_healthy
esea_vec_microarray_disease <- c(length(edge_list_biopsy_microarray_disease):1)
names(esea_vec_microarray_disease) <- edge_list_biopsy_microarray_disease
esea_vec_rnaseq_disease <- c(length(edge_list_biopsy_rnaseq_disease):1)
names(esea_vec_rnaseq_disease) <- edge_list_biopsy_rnaseq_disease

biopsy_microarray_healthy_esea <- compute_esea(EdgeCorScore = esea_vec_microarray_healthy, pathwayEdge.db = esea_background, pathway = "kegg", nperm = 100)
biopsy_rnaseq_healthy_esea <- compute_esea(EdgeCorScore = esea_vec_rnaseq_healthy, pathwayEdge.db = esea_background, pathway = "kegg", nperm = 100)
biopsy_microarray_disease_esea <- compute_esea(EdgeCorScore = esea_vec_microarray_disease, pathwayEdge.db = esea_background, pathway = "kegg", nperm = 100)
biopsy_rnaseq_disease_esea <- compute_esea(EdgeCorScore = esea_vec_rnaseq_disease, pathwayEdge.db = esea_background, pathway = "kegg", nperm = 100)


#################################################################################################################
######################### LATE INTEGRATION ######################################################################
#################################################################################################################

##### Disease #####

#names(list_of_modules_biopsy)
#[1] "modules_biopsy_microarray_disease" "modules_biopsy_microarray_healthy" "modules_biopsy_rnaseq_disease"     "modules_biopsy_rnaseq_healthy"    

modules_rnaseq_disease<-list_of_modules_biopsy[[3]]
modules_micro_disease<-list_of_modules_biopsy[[1]]

# Assuming rnaseq_disease is similar to micro_disease
# Transform micro_disease into the required format
micro_community_labels <- setNames(as.list(modules_micro_disease$membership), modules_micro_disease$names)
rnaseq_community_labels <- setNames(as.list(modules_rnaseq_disease$membership),modules_rnaseq_disease$names)

# Create a list of community labels for different views
community_labels_list <- list("micro" = micro_community_labels, "rnaseq" = rnaseq_community_labels)

n1<-max(modules_rnaseq_disease$membership)
n2<-max(modules_micro_disease$membership)

n_aggregate_communities_avg <- round((n1 + n2) / 2)

# Apply the aggregate_communities function
aggregated_results_disease <- aggregate_communities(community_labels_list, n_aggregate_communities = n_aggregate_communities_avg, max_iter = 50000)

# Print the results
print(aggregated_results_disease)

# Get the index of the maximum value in each row
most_probable_communities_disease <- get_most_probable_community(aggregated_results_disease)

##### Reactome from late integrated disease network #####

#####ORA FOR DISEASE#####

network_disease_rnaseq<-readRDS(file.path("case_study", "network_biopsy_rnaseq_disease.rds"))
network_disease_microarray<-readRDS(file.path("case_study", "network_biopsy_microarray_disease.rds"))

# Extract metrics from RNA-Seq network
rna_seq_metrics <- data.frame(
  name = V(network_disease_rnaseq)$name,
  degree = V(network_disease_rnaseq)$degree,
  betweenness = V(network_disease_rnaseq)$betweenness,
  closeness = V(network_disease_rnaseq)$closeness,
  eigenvector = V(network_disease_rnaseq)$eigenvector
)

# Extract metrics from microarray network
microarray_metrics <- data.frame(
  name = V(network_disease_microarray)$name,
  degree = V(network_disease_microarray)$degree,
  betweenness = V(network_disease_microarray)$betweenness,
  closeness = V(network_disease_microarray)$closeness,
  eigenvector = V(network_disease_microarray)$eigenvector
)

# Merge metrics from both networks
merged_metrics <- merge(rna_seq_metrics, microarray_metrics, by = "name", suffixes = c("_rnaseq", "_micro"))

# Add module information and compute weighted metrics
merged_metrics$module <- most_probable_communities_disease[merged_metrics$name]
module_contributions_disease <- aggregated_results_disease$view_contributions
  
merged_metrics <- merged_metrics %>%
  rowwise() %>%
  mutate(
    weighted_betweenness = betweenness_rnaseq * module_contributions_disease["rnaseq", module] +
      betweenness_micro * module_contributions_disease["micro", module],
    weighted_closeness = closeness_rnaseq * module_contributions_disease["rnaseq", module] +
      closeness_micro * module_contributions_disease["micro", module],
    weighted_degree = degree_rnaseq * module_contributions_disease["rnaseq", module] +
      degree_micro * module_contributions_disease["micro", module]
  )


# Filter top x% weighted betweenness within each module
filtered_genes <- merged_metrics %>%
  group_by(module) %>%
  filter(weighted_betweenness >= quantile(weighted_betweenness, 0.5)) %>%
  ungroup()

filtered_genes <- as.data.frame(filtered_genes)

# Load required libraries
library(clusterProfiler)
library(tidyverse)

# Assuming 'filtered_genes' is your dataframe
modules <- unique(filtered_genes$module) %>% 
  as.numeric() %>% 
  sort()
# Unique modules
ora_results <- list()

# Loop through each module
for (mod in modules) {
  # Filter genes by module
  genes <- filtered_genes %>% 
    filter(module == mod) %>% 
    pull(name)
  
  # Map gene symbols to Entrez IDs
  gene_symbols <- as.character(genes)
  
  # Map symbols to Entrez IDs
  entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns = "ENTREZID")
  
  # Remove rows with NA/invalid IDs
  entrez_ids <- entrez_ids %>% 
    filter(!is.na(ENTREZID)) %>% 
    pull(ENTREZID)
  
  # Perform ORA using Reactome
  reactome_ora <- enrichPathway(gene = entrez_ids, 
                                organism = 'human', 
                                pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.2)
  
  # Store the result
  ora_results[[paste0("Module_", mod)]] <- reactome_ora
}

# View the results for each module
ora_results


# Save the results to an Excel file
#install.packages("openxlsx")
library(openxlsx)
output_file <- "ORA_results_by_module_disease.xlsx"
wb <- createWorkbook()

for (sheet_name in names(ora_results)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, ora_results[[sheet_name]])
}

saveWorkbook(wb, file.path("case_study",output_file), overwrite = TRUE)

cat("Results saved to", output_file)

##### Reactome from late integrated disease network #####

merged_disease<-get_reactome_from_modules(modules = most_probable_communities_disease, geneID="SYMBOL", pval_cutoff = 0.05, 
                                          outPath = paste0("case_study", "/pathway_results_merged_disease"), layout = "overall")

bubble_merged_disease<-get_bubbleplot_from_pathways(modules = most_probable_communities_disease, geneID = "SYMBOL")

ggsave("case_study/bubble_plot_merged_disease.png")


##### Healthy #####

modules_rnaseq_healthy<-list_of_modules_biopsy[[4]]
modules_micro_healthy<-list_of_modules_biopsy[[2]]

# Assuming rnaseq_healthy is similar to micro_healthy
# Transform micro_healthy into the required format
micro_community_labels <- setNames(as.list(modules_micro_healthy$membership), modules_micro_healthy$names)
rnaseq_community_labels <- setNames(as.list(modules_rnaseq_healthy$membership),modules_rnaseq_healthy$names)

# Create a list of community labels for different views
community_labels_list <- list("micro" = micro_community_labels, "rnaseq" = rnaseq_community_labels)

n1<-max(modules_rnaseq_healthy$membership)
n2<-max(modules_micro_healthy$membership)

n_aggregate_communities_avg <- round((n1 + n2) / 2)

# Apply the aggregate_communities function
aggregated_results_healthy <- aggregate_communities(community_labels_list, n_aggregate_communities = n_aggregate_communities_avg, max_iter = 50000)

# Print the results
print(aggregated_results_healthy)

# Get the index of the maximum value in each row
most_probable_communities_healthy <- get_most_probable_community(aggregated_results_healthy)


##### ORA FOR healthy#####

network_healthy_rnaseq<-readRDS(file.path("case_study", "network_biopsy_rnaseq_healthy.rds"))
network_healthy_microarray<-readRDS(file.path("case_study", "network_biopsy_microarray_healthy.rds"))

# Extract metrics from RNA-Seq network
rna_seq_metrics <- data.frame(
  name = V(network_healthy_rnaseq)$name,
  degree = V(network_healthy_rnaseq)$degree,
  betweenness = V(network_healthy_rnaseq)$betweenness,
  closeness = V(network_healthy_rnaseq)$closeness,
  eigenvector = V(network_healthy_rnaseq)$eigenvector
)

# Extract metrics from microarray network
microarray_metrics <- data.frame(
  name = V(network_healthy_microarray)$name,
  degree = V(network_healthy_microarray)$degree,
  betweenness = V(network_healthy_microarray)$betweenness,
  closeness = V(network_healthy_microarray)$closeness,
  eigenvector = V(network_healthy_microarray)$eigenvector
)

# Merge metrics from both networks
merged_metrics <- merge(rna_seq_metrics, microarray_metrics, by = "name", suffixes = c("_rnaseq", "_micro"))

# Add module information and compute weighted metrics
merged_metrics$module <- most_probable_communities_healthy[merged_metrics$name]

module_contributions_healthy <- aggregated_results_healthy$view_contributions

merged_metrics <- merged_metrics %>%
  rowwise() %>%
  mutate(
    weighted_betweenness = betweenness_rnaseq * module_contributions_healthy["rnaseq", module] +
      betweenness_micro * module_contributions_healthy["micro", module],
    weighted_closeness = closeness_rnaseq * module_contributions_healthy["rnaseq", module] +
      closeness_micro * module_contributions_healthy["micro", module],
    weighted_degree = degree_rnaseq * module_contributions_healthy["rnaseq", module] +
      degree_micro * module_contributions_healthy["micro", module]
  )


# Filter top x% weighted betweenness within each module
filtered_genes <- merged_metrics %>%
  group_by(module) %>%
  filter(weighted_betweenness >= quantile(weighted_betweenness, 0.5)) %>%
  ungroup()

filtered_genes <- as.data.frame(filtered_genes)

# Load required libraries
#library(clusterProfiler)
#library(tidyverse)

# Assuming 'filtered_genes' is your dataframe
modules <- unique(filtered_genes$module) %>% 
  as.numeric() %>% 
  sort()
# Unique modules


ora_results <- list()

# Loop through each module
for (mod in modules) {
  # Filter genes by module
  genes <- filtered_genes %>% 
    filter(module == mod) %>% 
    pull(name)
  
  
  # Map gene symbols to Entrez IDs
  gene_symbols <- as.character(genes)
  
  # Map symbols to Entrez IDs
  entrez_ids <- AnnotationDbi::select(org.Hs.eg.db, keys = gene_symbols, keytype = "SYMBOL", columns = "ENTREZID")
  
  # Remove rows with NA/invalid IDs
  entrez_ids <- entrez_ids %>% 
    filter(!is.na(ENTREZID)) %>% 
    pull(ENTREZID)
  
  # Perform ORA using Reactome
  reactome_ora <- enrichPathway(gene = entrez_ids, 
                                organism = 'human', 
                                pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.2)
  
  # Store the result
  ora_results[[paste0("Module_", mod)]] <- reactome_ora
}

# View the results for each module
ora_results

# Save the results to an Excel file
#install.packages("openxlsx")
#library(openxlsx)
output_file <- "case_study/ORA_results_by_module_healthy.xlsx"
wb <- createWorkbook()

for (sheet_name in names(ora_results)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, ora_results[[sheet_name]])
}

saveWorkbook(wb, file.path("case_study",output_file), overwrite = TRUE)
cat("Results saved to", output_file)

##### Reactome from late integrated healthy network #####

merged_healthy<-get_reactome_from_modules(modules = most_probable_communities_healthy, geneID="SYMBOL", pval_cutoff = 0.05, 
                                          outPath = paste0("case_study", "/pathway_results_merged_healthy"), layout = "overall")

bubble_merged_healthy<-get_bubbleplot_from_pathways(modules = most_probable_communities_healthy, geneID = "SYMBOL")

ggsave("case_study/bubble_plot_merged_healthy.png")


#################################################################################################################
######################### SNF MULTI-OMICS INTEGRATION ###########################################################
#################################################################################################################


methylation_data_72h<-read.table("case_study/methylation_aggregated_exposure_72h.txt", sep="\t")

metadata_methylation<-xlsx::read.xlsx("case_study/metadata_methylation.xlsx", sheetIndex = 1)


##### Controls ######

# Extract rows where both conditions are met
metadata_methylation_72h_control <- metadata_methylation[metadata_methylation$Cytokines.exposure == "72h" & metadata_methylation$treatment == "C", ]

# Verify the result
head(metadata_methylation_72h_control)

# Step 1: Get the current column names
col_names <- colnames(methylation_data_72h)

# Step 2: Remove the 'X' prefix and '_median' suffix
new_col_names <- gsub("^X", "", col_names)          # Remove 'X' prefix
new_col_names <- gsub("_median$", "", new_col_names) # Remove '_median' suffix

# Step 3: Assign the new column names to the data frame
colnames(methylation_data_72h) <- new_col_names

# Verify the result
head(methylation_data_72h)

transcriptomics_data<-read.table("case_study/normalized_counts_matrix_deseq_SYMBOL.txt", sep="\t")

# Step 1: Get the current column names
col_names <- colnames(transcriptomics_data)

# Step 2: Remove the 'X' prefix and '_median' suffix
new_col_names <- gsub("^X", "", col_names)          # Remove 'X' prefix


# Step 3: Assign the new column names to the data frame
colnames(transcriptomics_data) <- new_col_names

metadata_transcriptomics<-xlsx::read.xlsx("case_study/metadata_transcriptomics.xlsx", sheetIndex = 1)

# Extract rows where both conditions are met
metadata_transcriptomics_72h_control <- metadata_transcriptomics[metadata_transcriptomics$Cytokines.exposure == "72h" & metadata_transcriptomics$treatment == "C", ]


methylation_data_72h_control<-methylation_data_72h[,colnames(methylation_data_72h)%in%metadata_methylation_72h_control$Sample.name]

transcriptomics_data_72h_control<-transcriptomics_data[,colnames(transcriptomics_data)%in%metadata_transcriptomics_72h_control$Sample.name]

data_list<-list(methylation_data_72h_control, transcriptomics_data_72h_control)

names(data_list)<-c("methylation_data_72h_control", "transcriptomics_data_72h_control")

snf_result_control_72h<-snf_based_integration(data_list)


##### IL4-IL13 #####

# Extract rows where both conditions are met
metadata_methylation_72h_IL4_IL13 <- metadata_methylation[metadata_methylation$Cytokines.exposure == "72h" & metadata_methylation$treatment == "IL4-IL13", ]

# Verify the result
head(metadata_methylation_72h_IL4_IL13)


# Verify the result
head(methylation_data_72h)


# Extract rows where both conditions are met
metadata_transcriptomics_72h_IL4_IL13 <- metadata_transcriptomics[metadata_transcriptomics$Cytokines.exposure == "72h" & metadata_transcriptomics$treatment == "IL4-IL13", ]

# Verify the result
head(metadata_transcriptomics_72h_IL4_IL13)

methylation_data_72h_IL4_IL13<-methylation_data_72h[,colnames(methylation_data_72h)%in%metadata_methylation_72h_IL4_IL13$Sample.name]

transcriptomics_data_72h_IL4_IL13<-transcriptomics_data[,colnames(transcriptomics_data)%in%metadata_transcriptomics_72h_IL4_IL13$Sample.name]

head(methylation_data_72h_IL4_IL13)

head(transcriptomics_data_72h_IL4_IL13)

data_list<-list(methylation_data_72h_IL4_IL13, transcriptomics_data_72h_IL4_IL13)

names(data_list)<-c("methylation_data_72h_IL4_IL13", "transcriptomics_data_72h_IL4_IL13")

snf_result_IL4_IL13_72h<-snf_based_integration(data_list)

##### LPS #####

# Extract rows where both conditions are met
metadata_methylation_72h_LPS <- metadata_methylation[metadata_methylation$Cytokines.exposure == "72h" & metadata_methylation$treatment == "LPS-INFgamma", ]

# Verify the result
head(metadata_methylation_72h_LPS)


# Verify the result
head(methylation_data_72h)

# Extract rows where both conditions are met
metadata_transcriptomics_72h_LPS <- metadata_transcriptomics[metadata_transcriptomics$Cytokines.exposure == "72h" & metadata_transcriptomics$treatment == "LPS-INFgamma", ]

# Verify the result
head(metadata_transcriptomics_72h_LPS)


methylation_data_72h_LPS<-methylation_data_72h[,colnames(methylation_data_72h)%in%metadata_methylation_72h_LPS$Sample.name]

transcriptomics_data_72h_LPS<-transcriptomics_data[,colnames(transcriptomics_data)%in%metadata_transcriptomics_72h_LPS$Sample.name]

head(methylation_data_72h_LPS)

head(transcriptomics_data_72h_LPS)

data_list<-list(methylation_data_72h_LPS, transcriptomics_data_72h_LPS)

names(data_list)<-c("methylation_data_72h_LPS", "transcriptomics_data_72h_LPS")

snf_result_LPS_72h<-snf_based_integration(data_list)

##### ESEA for multi-omics integrated networks #####

pruned_network_control_72h <- prune_snf_network(snf_result_control_72h, percentile_thr = 0.9)
pruned_network_IL4_IL13 <- prune_snf_network(snf_result_IL4_IL13_72h, percentile_thr = 0.9)
pruned_network_LPS <- prune_snf_network(snf_result_LPS_72h, percentile_thr = 0.9)

edge_list_control_72h <- get_ranked_edge_list(snf_result = pruned_network_control_72h)
edge_list_IL4_IL13 <- get_ranked_edge_list(snf_result = pruned_network_IL4_IL13)
edge_list_LPS <- get_ranked_edge_list(snf_result = pruned_network_LPS)

names(edge_list_control_72h) <- gsub(names(edge_list_control_72h), pattern = "_", replacement = "|")
names(edge_list_LPS) <- gsub(names(edge_list_LPS), pattern = "_", replacement = "|", fixed = TRUE)
names(edge_list_IL4_IL13) <- gsub(names(edge_list_IL4_IL13), pattern = "_", replacement = "|", fixed = TRUE)

control_72h_esea <- compute_esea(EdgeCorScore = edge_list_control_72h, pathwayEdge.db = esea_background, pathway = "kegg", nperm = 100)
IL4_IL13_72h_esea <- compute_esea(EdgeCorScore = edge_list_IL4_IL13, pathwayEdge.db = esea_background, pathway = "kegg", nperm = 100)
LPS_72h_esea <- compute_esea(EdgeCorScore = edge_list_LPS, pathwayEdge.db = esea_background, pathway = "kegg", nperm = 100)
