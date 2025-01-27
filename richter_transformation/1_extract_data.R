library(Seurat)
library(Matrix)

base_dir = "../data/Nadeu2022_NatMed_scRNAseq_data/"
setwd(base_dir)
in_path = "seurat_objects"

# files with tumor and normal cells 
batch_1 = readRDS(file.path(in_path, "patient_63/2.seurat_object_filtered.rds"))
batch_2 = readRDS(file.path(in_path, "3.seurat_filtered.rds"))

cat("Merging data for tumor and normal cells...\n")
batch_1$donor_id = batch_1$donor_id.x
batch_1$donor_id.x = NULL
batch_1$donor_id.y = NULL
seurat_obj = merge(batch_1, batch_2)


cat("Writing output files for tumor and normal cells...\n")
outpath = "filtered"
if (!file.exists(outpath)) {
    dir.create(outpath)
}
counts <- GetAssayData(seurat_obj, layer="counts")
writeMM(counts, file.path(outpath, "matrix.mtx"))
write.table(rownames(counts), file = file.path(outpath, "genes.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)
write.table(colnames(counts), file = file.path(outpath, "barcodes.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)
write.table(seurat_obj@meta.data, file = file.path(outpath, "metadata.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)


# files with tumor cells only with annotations
cat("Merging annotated tumor cells...\n")
rds_files = file.path(in_path, c("6.seurat_annotated_12.rds", 
                                "6.seurat_annotated_19.rds", 
                                "6.seurat_annotated_3299.rds", 
                                "6.seurat_annotated_365.rds",
                                "patient_63/3.seurat_annotated.rds"))
seurat_list = list()
for (rds_file in rds_files) {
    seurat_obj = readRDS(rds_file)
    if ("donor_id.x" %in% colnames(seurat_obj@meta.data)) {
        donor_id = seurat_obj$donor_id.x
        seurat_obj$donor_id = donor_id
        seurat_obj$donor_id.x = NULL
        seurat_obj$donor_id.y = NULL
    }
    print(unique(seurat_obj$donor_id))
    seurat_list = c(seurat_list, seurat_obj)
}
seurat_obj = merge(seurat_list[[1]], seurat_list[-1])


cat("Writing combined output files for annotated tumor cells...\n")
outpath = "combined"
if (!file.exists(outpath)) {
    dir.create(outpath)
}
counts <- GetAssayData(seurat_obj, layer="counts")
writeMM(counts, file.path(outpath, "matrix.mtx"))
write.table(rownames(counts), file = file.path(outpath, "genes.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)
write.table(colnames(counts), file = file.path(outpath, "barcodes.tsv"), sep = "\t", quote = FALSE, col.names = FALSE)
write.table(seurat_obj@meta.data, file = file.path(outpath, "metadata.tsv"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = TRUE)
