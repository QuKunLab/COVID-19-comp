setwd('~/workspace/projects/Infection_Virus/Codes.D05.Analysis02')
library(WGCNA)

options(stringsAsFactors = FALSE)
avg_exp = read.csv('./output/avg_exp.csv', row.names = 1)
avg_exp = t(avg_exp)
row.names(avg_exp)[1:5];colnames(avg_exp)[1:5]

powers = c(1:20)
sft = pickSoftThreshold(avg_exp, verbose = 5,powerVector=powers)

# Plot the results:
# sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = sft$powerEstimate;  # 10
adjacency = adjacency(avg_exp, power = softPower)

TOM = TOMsimilarity(adjacency,TOMType = "signed");
dissTOM = 1-TOM
geneTree = hclust(as.dist(dissTOM), method = "average");
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 25;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,minClusterSize = minModuleSize);
table(dynamicMods)


# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
# sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    main = "Gene dendrogram and module colors")

# Calculate eigengenes
MEList = moduleEigengenes(avg_exp, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")


# Call an automatic merging function
merge = mergeCloseModules(avg_exp, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

pdf(file = './figures/heatmap.WGCNA..dendro.pdf')
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

####################################################################################
new_gene_order = colnames(avg_exp)[geneTree$order]
avg_exp_new = avg_exp[, new_gene_order]
write.csv(avg_exp_new, file = './output/avg_exp.new.csv', quote = F)

write(mergedColors[geneTree$order], file = './output/WGCNA.colors.txt')
####################################################################################

moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;


save.image(file = './output/WGCNA.RData')

