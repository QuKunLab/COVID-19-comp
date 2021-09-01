library(monocle3)
library(ggplot2)
library(dplyr)

file='cluster10_hvg_overlap2000_exp'
n_neighbors = 4
num_dim = 5
dir.create('trajectory')
dir.create(sprintf("trajectory/%s", file))
dir.create(sprintf("trajectory/%s/%s", file, n_neighbors))
dir.create(sprintf("trajectory/csv"))
tfidf <- read.csv(sprintf('%s.csv', file))[1:157,]
cluster <- read.csv('ann.csv')[1:157,]

expression_data <- t(tfidf[2:dim(tfidf)[2]])
cell_metadata <- data.frame(cluster)
gene_metadata <- data.frame(tfidf[2:dim(tfidf)[2],])
rownames(gene_metadata) <- colnames(tfidf)[2:dim(tfidf)[2]]
colnames(expression_data) <- rownames(cell_metadata)
cds <- new_cell_data_set(expression_data, cell_metadata = cell_metadata)

cds <- preprocess_cds(cds, method='PCA', num_dim = num_dim)
cds <- align_cds(cds, alignment_k=10, preprocess_method='PCA', alignment_group = "dataset")

cds <- reduce_dimension(cds, preprocess_method = 'Aligned', reduction_method='UMAP', umap.n_neighbors=n_neighbors)
cds <- cluster_cells(cds, reduction_method = 'UMAP', resolution=0.01)
cds <- learn_graph(cds)

pdf(file = sprintf("trajectory/%s/%s/trajectory_%s.pdf", file, n_neighbors, num_dim))
p<-plot_cells(cds, label_leaves=FALSE,label_groups_by_cluster=FALSE, label_cell_groups=FALSE, color_cells_by = "state", group_label_size=4,cell_size=1.2) + scale_color_manual(breaks = c("control", "moderate convalescent", "moderate onset", "severe convalescent", "severe onset"), values=c("#6699CC", "#1565C0", "#FF9900", '#6A1B9A','#990033')) + theme(legend.position = "right")
print(p)

p<-plot_cells(cds, label_leaves=FALSE,label_groups_by_cluster=FALSE, label_cell_groups=FALSE, color_cells_by = "cluster", group_label_size=4,cell_size=1.2)  # + scale_color_manual(breaks = c(1, 3, 2), values=c("#FF9900", "#6699CC", "#990033")) + theme(legend.position = "right")
print(p)

get_earliest_principal_node  <- function(cds, time_bin="control"){
  cell_ids <- which(colData(cds)[, "state"] == time_bin)
  
  
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
p<-plot_cells(cds, label_leaves=FALSE,label_groups_by_cluster=FALSE, label_cell_groups=FALSE, color_cells_by = "pseudotime", group_label_size=4,cell_size=1.2)

p<-plot_cells(cds, label_leaves=FALSE,label_groups_by_cluster=FALSE, label_cell_groups=FALSE, color_cells_by = "inflammatory_cytokine", group_label_size=4,cell_size=1.2)
print(p)

p<-plot_cells(cds, label_leaves=FALSE,label_groups_by_cluster=FALSE, label_cell_groups=FALSE, color_cells_by = "dataset", group_label_size=4,cell_size=1.2)

p<-plot_pc_variance_explained(cds)
print(p)

dev.off()

colData(cds)['cluster'] <- clusters(cds)
colData(cds)['pseudotime'] <- pseudotime(cds)
write.csv(colData(cds)[c('sampleID', 'state', 'cluster', 'inflammatory_cytokine', 'pseudotime')], sprintf("trajectory/csv/%s_%s_%s.csv", file, n_neighbors, num_dim), row.names = F)
