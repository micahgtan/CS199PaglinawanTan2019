# Defining a function for finding undetected packages and checking if packages are loaded properly
pkgTest <- function(x)
{
	if (!require(x,character.only = TRUE))
	{
		install.packages(x,dep=TRUE)
		if(!require(x,character.only = TRUE)) stop("Package not found")
	}
}
# Usage: pkgTest("<undetected package>");

# Load the WGCNA package
library(WGCNA);
pkgTest("WGCNA");

# Allow multi-threading within WGCNA 
enableWGCNAThreads();

########## Yeast Data Set Preprocessing ##########

# Read in the yeast data set
combined = read.delim("combined.txt", stringsAsFactors=FALSE);

# Remove unneeded columns
combined = within(combined, rm(clb, alpha, cdc15, cdc28, elu));

# Transpose the data table
yeast = as.data.frame(t(combined));

# Update row and column names
names(yeast) = combined$X;
rownames(yeast) = names(combined);

# Remove the first row
yeast = yeast[-1,];

# Get the vector of the column names
yeastColumnNames = colnames(yeast);

# The data type of the elements in the data frame should be numeric or logical for the mean imputation function to work properly.
# One cannot safely convert factors directly to numeric.
# The as.character function must be applied first.
# Otherwise, the factors will be converted to their numeric storage values.

# Convert the elements of the data frame from integer type into character type
yeast[, c(yeastColumnNames)] = sapply(yeast[, c(yeastColumnNames)], as.character);

# Convert the elements of the data frame from character type into numeric type
yeast[, c(yeastColumnNames)] = sapply(yeast[, c(yeastColumnNames)], as.numeric);

# Mean imputation per column to fill up missing data
for(i in 1:ncol(yeast)) {
  yeast[ , i][is.na(yeast[ , i])] = mean(yeast[ , i], na.rm = TRUE)
}

# Cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
sampleTree = hclust(dist(yeast), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.size
GrWindow(12,9);
par(cex = 0.6);
par(mar = c(0,4,2,0));
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,cex.axis = 1.5, cex.main = 2);

# It appears there are 5 outliers
# Choose a height cut that will remove the offending samples, say 50
# Plot a line to show the cut
abline(h = 50, col = "red");

# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 50, minSize = 10);
table(clust);

# The yeast variable contains the samples we want to keep
keepSamples = (clust==1);
yeast = yeast[keepSamples, ];
nGenes = ncol(yeast);
nSamples = nrow(yeast);

########## Yeast Data Set Power Estimate ##########

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2));

# Call the network topology analysis function
sft = pickSoftThreshold(yeast, powerVector = powers, verbose = 5);

# Plot the results:
sizeGrWindow(9, 5);
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"));
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red");

# Check the computed soft thresholding power value
sft$powerEstimate;
 
# 5 was the value obtained
softPower = sft$powerEstimate;

########## Yeast Data Clustering and Module Detection ##########

##### Constructing the gene network and identifying modules 
##### This is the simplified function of WGCNA which uses hierarchical clustering
net = blockwiseModules(yeast, power = sft$powerEstimate, TOMType = "unsigned", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "yeastTOM", verbose = 3);

# Check the number of modules and module sizes 
table(net$colors);

# Open a graphics window
sizeGrWindow(12, 9);

# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors);

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05);

# Save the module assignment and module eigengene information 
moduleLabels = net$colors;
moduleColors = labels2colors(net$colors);
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
hcluster = moduleLabels;

# In order to proceed to clustering we need the dissimilarity matrix 
# Use the soft thresholding power value for adjacency matrix
adjacency = adjacency(yeast, power = softPower);

# Transform the adjacency matrix into TOM
TOM = TOMsimilarity(adjacency);

# Compute for the dissimilarity matrix
dissTOM = 1 - TOM;

# Since net function covered the hierarchical clustering, we now proceed to kmeans clustering
# We use the obtained number of clusters from hierarchical clustering which is 22
kclust = kmeans(dissTOM, 22);

# Save the clusters produced by kmeans clustering 
kcluster = kclust$cluster;

# We now proceed to DBSCAN

# Load the DBSCAN package
library(DBSCAN);
pkgTest("DBSCAN");

# Compute for the epsilon distance parameter
dbscan::kNNdistplot(dissTOM, k = 30);
# We get 1.4

# Use function DBSCAN, with minPts parameter 15 and 30
dclust_15 = dbscan(dissTOM, 1.4, 15);
dclust_30 = dbscan(dissTOM, 1.4, 30);

# Save the module assignment information
dcluster15 = dclust_15$cluster;
dcluster30 = dclust_30$cluster;

########## Yeast Data Simulation ##########

# Number of samples or microarrays in the data
no.obs = 72;

# Number of genes
nGenes1 = 6178;

# Getting the proportion of genes in the modules
# The last number in the list is the number of genes unassigned to a module
simulateProportions1 = c(929,871,801,408,383,354,163,161,153,152,132,132,125,107,94,82,80,76,74,59,52,39,751)/nGenes1;

# Sort the modules alphabetically
MEs = MEs[, c("ME0","ME1","ME2","ME3","ME4","ME5","ME6","ME7","ME8","ME9","ME10","ME11","ME12","ME13","ME14","ME15","ME16","ME17","ME18","ME19","ME20","ME21","ME22")];

# Remove "grey" module / genes not related to any of the modules
MEs = within(MEs, rm(ME0));

# Simulate data
simulatedYeastDataResult = simulateDatExpr(eigengenes=MEs, nGenes=nGenes1, modProportions=simulateProportions1, verbose=3); #signed?

# Get the simulated data set data frame
simulatedYeastData = as.data.frame(simulatedYeastDataResult$datExpr);

# Get the simulated data set module assignment
simulatedYeastDataModuleAssignment = simulatedYeastDataResult$setLabels;

########## Yeast Data Statistical Testing ##########

# Load the clv package for statistical testing tools
library(clv);
pkgTest("clv");

# Adding one to the integer vectors to avoid errors caused by "0"s
simulatedYeastDataModuleAssignment <- simulatedYeastDataModuleAssignment+1;
hcluster = hcluster+1;
kcluster = kcluster+1;
dcluster15 = dcluster15+1;
dcluster30 = dcluster30+1;

# Coerce the module assignment vectors into integers for the statistical testing tools to work properly.
simulatedYeastDataModuleAssignment <- as.integer(simulatedYeastDataModuleAssignment);
hcluster = as.integer(hcluster);
kcluster = as.integer(kcluster);
dcluster15 = as.integer(dcluster15);
dcluster30 = as.integer(dcluster30);

# Standard WGCNA
stdwgcna = std.ext(simulatedYeastDataModuleAssignment, hcluster);
# K means
stdkmeans = std.ext(simulatedYeastDataModuleAssignment, kcluster);
# DBSCAN
stddbscan15 = std.ext(simulatedYeastDataModuleAssignment, dcluster15);
stddbscan30 = std.ext(simulatedYeastDataModuleAssignment, dcluster30);

# Rand Index
randsimvswgcna = clv.Rand(stdwgcna);
randsimvskmeans = clv.Rand(stdkmeans);
randsimvsdbscan15 = clv.Rand(stddbscan15);
randsimvsdbscan30 = clv.Rand(stddbscan30);

# Jaccard Coefficient
jaccardsimvswgcna = clv.Jaccard(stdwgcna);
jaccardsimvskmeans = clv.Jaccard(stdkmeans);
jaccardsimvsdbscan15 = clv.Jaccard(stddbscan15);
jaccardsimvsdbscan30 = clv.Jaccard(stddbscan30);

# Similarity Index
similaritysimvswgcna = similarity.index(confusion.matrix(simulatedYeastDataModuleAssignment, hcluster));
similaritysimvskmeans = similarity.index(confusion.matrix(simulatedYeastDataModuleAssignment, kcluster));
similaritysimvsdbscan15 = similarity.index(confusion.matrix(simulatedYeastDataModuleAssignment, dcluster15));
similaritysimvsdbscan30 = similarity.index(confusion.matrix(simulatedYeastDataModuleAssignment, dcluster30));