fsom26 = readRDS("C:/Users/00ber/OneDrive/Desktop/VPC/human1/RDS/Edited/fsom_edited_26_meta.rds")
sample_file = read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/Info/MR242 Sample Information.csv")

comparisons = list(
  male_vs_female = list(Sex = list("male", "female")),
  male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female")),
  ctrl_vs_mibc = list(Disease = list("MIBC", "Ctrl")),
  nac_vs_no_nac = list(Disease = "MIBC", NAC = list("NAC", "No.NAC"))
)

sample_info = prepareSampleInfo("C:/Users/00ber/OneDrive/Desktop/VPC/human1/Info/MR242 Sample Information.csv",
                                "Sample.Name", "File.Name", comparisons)

disease_groups = list("Control" = c("female_Ctrl_X", "male_Ctrl_X"),
                      "MIBC" = c("female_MIBC_No.NAC", "female_MIBC_NAC", "male_MIBC_No.NAC", "male_MIBC_NAC"))
umap = plotGroupUMAPs(fsom, sample_info, disease_groups, num_cells = 1000, seed = 41)

colors <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4",
            "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", "#008080", "#E6BEFF",
            "#9A6324", "#FF4500", "#800000", "#AAFFC3", "#808000", "#FFD8B1",
            "#000075", "#808080", "#0000FF", "#000000", "#A9A9A9", "#DC143C",
            "#7FFF00", "#8B4513")

labels <- c("Naive CD8 T cells", "CD8 Temra", "CD8 Tscm", "CD8 Tcm", "CD8 Tem",
            "CD8 T cell other", "Naive CD4 T cells", "CD4 Temra", "CD4 Tscm", "CD4 Tcm",
            "CD4 Tem", "CD4 T cell other", "Treg CD4 T cells", "IgD+ B cells", "IgD- B cells",
            "Plasmablasts", "CD56hi NK cells", "CD56lo NK cells", "NK T cells", "gdT cells", "cMo",
            "intMo", "pMo", "cDC", "pDC", "Undefined")

umap_full = plotUMAP(fsom, num_cells = 1000, labels = labels, colors = colors, seed = 43)

umaps = plotGroupUMAPs(fsom, sample_info, disease_groups, umap = umap_full,
                       num_cells = 30000, seed = 43)

design = makeDesignMatrix(sample_info)
contrasts = makeContrastsMatrix(sample_info, comparisons)
counts = makeCountMatrix(fsom, sample_info)
s = doDAAnalysis(design, counts, contrasts, sample_info, norm_method = "none", dir_tables = "unnorm")

###
file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
file2 <- system.file("extdata", "aggregate.fcs", package = "flowFun")
agg <- flowCore::read.FCS(file)
agg2 <- flowCore::read.FCS(file2)

agg_dt2 <- getTableFromFCS(list(file,file2))

attr(agg_dt, "clustered") <- colnames(agg_dt)[c(14,15,16)]
attributes(agg_dt)

cols_to_use <- c(10, 12:14, 15, 17:22, 25:34)+2

fsom2 <- flowSOMWrapper(agg_dt2,
                      cols_to_cluster = cols_to_use,
                      num_clus = 10,
                      seed = 42,
                      fsom_file = NULL)

plotUMAPNew(fsom, num_cells = 500)
plotMetaclusterMFIsNew(fsom)

merge <- editTableMetaclusters(fsom, new_labels = c("3" = "test"), cluster_assignments = c("96" = "test"))

merge %>%
 tidytable::select(c(Cluster, Metacluster, Meta_original, cell_id))

merge2 <- editTableMetaclusters(merge, new_labels = c("test" = "9", "7" = "8"), cluster_assignments = c("56" = "3"))

merge2 %>%
 tidytable::select(c(Cluster, Metacluster, Meta_original, cell_id))

filter <- createFilteredAggregateTable(fsom, num_cells = 1000, metaclusters = c(1,2,3,4,5))

pca <- doPCA(filter, cols_to_use)

plotPCAScree(pca)

fsom_pca <- clusterTableWithPCA(filter, pca, 5, 10)

id_cols <- agg_dt %>%
 tidytable::select(.id, cell_id) %>%
 data.frame(check.names = FALSE)


input <- flowCore::read.FCS("C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw/Ab_PHA_Ctrl AWB2.fcs")
file <- "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw/Ab_PHA_Ctrl AWB2.fcs"
mat <- read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/auto_final.csv", check.names = FALSE)[-1]
colnames(mat) = sub(" :: .*", "", colnames(mat))
table <- doPreprocessing(file, "BUV496-A", mat)

### DE
gagd <- file.path("C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw",
                  list("Ab_PHA_Ctrl AWB4.fcs", "Ab_PHA_Ctrl AWB9.fcs",
                       "Ab_PHA_MIBC MR93.fcs", "Ab_PHA_MIBC MR66.fcs",
                       "Ab_PHA_MIBC MR110.fcs", "Ab_PHA_MIBC MR75.fcs"))
tab <- getTableFromFCS(gagd, num_cells = 10000)
file <- system.file("extdata", "sample_aggregate.fcs", package = "flowFun")
file2 <- system.file("extdata", "aggregate.fcs", package = "flowFun")

p = doPreprocessing(ex, ld_channel = "BUV496-A", compensation = mat,
                     transformation_type = "logicle")

ex = "C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw/Ab_PHA_Ctrl MR201.fcs"
ex = flowCore::read.FCS(ex, truncate_max_range = TRUE)
ex <- PeacoQC::RemoveMargins(ex, c("FSC-A", "SSC-A"))
ex <- PeacoQC::RemoveDoublets(ex)

gate = flowDensity::flowDensity(obj = ex,
                                channels = c("FSC-A", "SSC-A"),
                                position = c(TRUE, NA),
                                percentile = c(0.45, NA))
filter = slot(gate, "filter")
ex_g <- flowDensity::getflowFrame(gate)

t = flowDensity::flowDensity(agg, channels = c("FSC-A", "SSC-A"),
                             position = c(FALSE, NA),
                             filter = filter)

flowDensity::plotDens(ex, c("FSC-A", "SSC-A"))
flowDensity::plotDens(gate, c("FSC-A", "SSC-A"))

ex_c <- flowCore::compensate(ex_g, spillover = mat)
transformation <- flowCore::estimateLogicle(ex_c, channels = "BUV496-A")
ex_c <- flowCore::transform(ex_c, transformation)

ex_l <- flowDensity::flowDensity(obj = ex_c,
                                 channels = c("BUV496-A", "FSC-A"),
                                 position = c(FALSE, NA))

flowDensity::plotDens(ex_c, c("BUV496-A", "FSC-A"))
flowDensity::plotDens(ex_l, c("BUV496-A", "FSC-A"))

###

# Stratified sampling of raw data
# dt <- getTableFromFCS("C:/Users/00ber/OneDrive/Desktop/VPC/human1/Data/Raw")
table <- system.file("extdata", "raw_samples.rds", package = "flowFun")
table <- readRDS(table)

# Get all file names
file_names <- dt %>%
  tidytable::pull(.id) %>%
  unique()

# Read in compensation matrix
mat <- read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/auto_final.csv", check.names = FALSE)[-1]
colnames(mat) = sub(" :: .*", "", colnames(mat))

# Test out preprocessing parameters
live <- flowCore::rectangleGate(.gate = list("BUV496-A"= c(-1, 1.75),
                                             "FSC-A" = c(0, 280000)))

flowDensity::deGate(flowCore::read.FCS(gagd[2]), "SSC-A")

previewPreprocessing(gagd[2],
                     ld_channel = "BUV496-A",
                     compensation = mat,
                     transformation_type = "logicle",
                     nmad = 3.5,
                     live_gate = live)


# Preprocess all samples
preprdt <- doPreprocessing(tab,
                           ld_channel = "BUV496-A",
                           compensation = mat,
                           transformation_type = "logicle",
                           nmad = 3.5)

# Clustering
cols_to_use <- c(10, 12:14, 15, 17:22, 25:34)+2
fsom_dt <- flowSOMWrapper(preprdt,
                          cols_to_cluster = cols_to_use,
                          num_clus = 20,
                          seed = 42,
                          fsom_file = NULL)

comparisons = list(
 male_vs_female = list(Sex = list("male", "female")),
 male_vs_female_mibc = list(Disease = "MIBC", Sex = list("male", "female")),
 ctrl_vs_mibc = list(Disease = list("MIBC", "Ctrl")),
 nac_vs_no_nac = list(Disease = "MIBC", NAC = list("NAC", "No.NAC"))
)

sample_file = read.csv("C:/Users/00ber/OneDrive/Desktop/VPC/human1/Info/MR242 Sample Information.csv")

sample_info = prepareSampleInfo("C:/Users/00ber/OneDrive/Desktop/VPC/human1/Info/MR242 Sample Information.csv",
                               "Sample.Name", "File.Name", comparisons)

design = makeDesignMatrix(sample_info)
contrasts = makeContrastsMatrix(sample_info, comparisons)
counts = makeCountMatrix(fsom_dt)

da = doDAAnalysis(design, counts, contrasts, sample_info, norm_method = "none")
de = doDEAnalysisTable(fsom_dt, sample_info, design, contrasts, counts, "PHA-L")
