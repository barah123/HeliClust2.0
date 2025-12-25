# HeliClust2.0  

![App Logo](./heliClust_logo.png)
---
HeliClust20 is an R package for streamlined hierarchical clustering of **metabolic biomarkers** and **microbiome/sequence (OTU/ASV) data**.  
It wraps common preprocessing steps (numeric conversion, scaling, distance calculation) and exposes simple, high-level functions for generating clustering-ready distance matrices and `hclust` objects.

HeliClust20 is designed to be:

- **Data-friendly** – accepts standard data frames or matrices with minimal formatting.
- **Method-flexible** – supports multiple distance and linkage methods through `vegan::vegdist()` and `stats::hclust()`.
- **Pipeline-ready** – returns base R objects (`matrix`, `dist`, `hclust`) that plug into existing workflows and visualizations.

---

## Installation

### From GitHub

- install.packages("remotes") # if needed
- remotes::install_github("barah123/HeliClust2.0", subdir = "HeliClust20")
- library(HeliClust20)
---

## Input data format

### Metabolite / biomarker data
- **Rows**: samples (e.g. subjects).
- **Columns**: numeric biomarkers (e.g. HOMA.IR, G30, G60, I30, …).
- Values should be numeric or coercible to numeric.
Example:
metabolites <- read.csv("lcs_v1.csv", row.names = 1, check.names = FALSE)

### OTU/ASV table for sequences

- **Rows**: samples.
- **Columns**: taxa / OTUs / ASVs.
- Entries: non-negative counts or relative abundances.

If your table has taxa as rows and samples as columns, transpose before using:

otu <- read.csv("otu_table_taxa_rows.csv", row.names = 1, check.names = FALSE)
otu <- t(otu) # rows = samples, columns = OTUs

---
## Core functions
### 1. `heliclust_dist()`
Compute a distance matrix suitable for hierarchical clustering:
dist_mat <- heliclust_dist(x,
dist_method = "euclidean",
scale_data = TRUE)
- `x`: data frame / matrix (rows = samples, columns = variables).
- `dist_method`: distance metric passed to `vegan::vegdist()` (e.g. `"euclidean"`, `"bray"`).
- `scale_data`: if `TRUE`, standardizes variables with `scale()`.
Returns a numeric matrix of pairwise distances.
---

### 2. `heliclust_metabolites()`

Hierarchical clustering for metabolic / biomarker data:
hc_met <- heliclust_metabolites(metabolites,
dist_method = "euclidean",
hclust_method = "ward.D2",
scale_data = TRUE)

plot(hc_met) # base R dendrogram

- Uses `heliclust_dist()` internally.
- Returns an `hclust` object created via `stats::hclust()`.

---

### 3. `heliclust_sequences()`
Hierarchical clustering for sequence / OTU data:
hc_seq <- heliclust_sequences(otu,
dist_method = "bray",
hclust_method = "average",
scale_data = FALSE)

plot(hc_seq)

Defaults are tuned for compositional sequence data (e.g. Bray–Curtis distances, average linkage), but can be customized.

---

## Example workflow

library(HeliClust20)

1. Load data
metabolites <- read.csv("lcs_v1.csv", row.names = 1, check.names = FALSE)

2. Distance matrix and clustering
dist_mat <- heliclust_dist(metabolites,
dist_method = "euclidean",
scale_data = TRUE)

hc_met <- heliclust_metabolites(metabolites,
dist_method = "euclidean",
hclust_method = "ward.D2",
scale_data = TRUE)

3. Visualizations
- plot(hc_met) # dendrogram
- heatmap(dist_mat) # simple distance heatmap

  ---
<div align="center">
  <img src="https://github.com/barah123/HeliClust2.0/blob/main/heat_1" width="45%" />
  <img src="https://github.com/barah123/HeliClust2.0/blob/main/heat_2" width="45%" />
  <img src="https://github.com/barah123/HeliClust2.0/blob/main/heat_3" width="45%" />
  <img src="https://github.com/barah123/HeliClust2.0/blob/main/heat_4" width="45%" />
</div>

  ---
4. Export distances (optional)
write.table(dist_mat, "metabolite_distances.txt",
sep = "\t", eol = "\n",
na = "", col.names = NA,
quote = FALSE, row.names = TRUE)




