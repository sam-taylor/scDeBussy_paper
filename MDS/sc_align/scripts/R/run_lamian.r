# LAMIAN MDS ANALYSIS
# CellAlignDTW Project | March 2025

library(Lamian)
library(anndata)
library(Matrix)

suppressMessages(library(Lamian))
# vignette("Lamian")
load("MDS/sc_align/data/processed/man_tree_data.rda")
data(man_tree_data)
set.seed(12345)
res <- infer_tree_structure(
    pca = man_tree_data[["pca"]],
    expression = man_tree_data[["expression"]],
    cellanno = man_tree_data[["cellanno"]],
    origin.marker = c("CD34"),
    number.cluster = 5,
    xlab = "Principal component 1",
    ylab = "Principal component 2"
)



Lamian::plotmclust(res,
    cell_point_size = 0.5,
    x.lab = "Principal component 1",
    y.lab = "Principal component 2"
)

tryCatch(
    {
        # ---- Data Loading ----
        adata <- read_h5ad("young_old_mds.h5ad")
        message("Loaded ", ncol(adata), " cells with ", nrow(adata), " genes")

        # ---- Format Conversion ----
        expr_mat <- t(adata$layers$counts) # Genes x Cells
        cell_anno <- data.frame(
            cell = colnames(adata),
            pseudotime = adata$obs$palantir_pseudotime,
            group = adata$obs$age_group
        )

        # ---- Run Lamian ----
        set.seed(123)
        result <- lamian.test(
            expr = expr_mat,
            cellanno = cell_anno[, c("cell", "pseudotime")],
            pseudotime = "pseudotime",
            design = model.matrix(~group, data = cell_anno),
            test.type = "variable"
        )

        # ---- Save Outputs ----
        saveRDS(result, "lamian_result.rds")
        write.csv(result$stat, "lamian_pseudotime_stats.csv")
        message("Analysis completed. Outputs saved.")
    },
    error = function(e) {
        message("ERROR: ", conditionMessage(e))
        dump.frames("lamian_dump", to.file = TRUE)
    }
)
