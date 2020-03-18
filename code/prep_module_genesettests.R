# packages
library("data.table")
library("here")


# constants
path_wgcna_geneMod <- "/projects/jonatan/pub-perslab/18-mousebrain/18-mousebrain_7/tables/ClusterName_prior_191127a_cell_cluster_module_genes.csv.gz"
path_df_ortho <- here("data","gene_annotation.hsapiens_mmusculus_unique_orthologs.GRCh37.ens_v91.txt.gz")

# load files
dt_geneMod <- fread(path_wgcna_geneMod)
dt_ortho <- fread(path_df_ortho)

# map genes from ENSMUS to ENSG
dt_mapping <- data.table(
  ENSMUS = dt_geneMod$ensembl
)

idx_ENSG <- match(dt_geneMod$ensembl,dt_ortho$mmusculus_homolog_ensembl_gene)

dt_mapping <- dt_mapping[!is.na(idx_ENSG)]

idx_ENSG <- idx_ENSG[!is.na(idx_ENSG)]

dt_mapping$ENSG <- dt_ortho$ensembl_gene_id[idx_ENSG]

dt_geneMod$ENSG <- dt_mapping$ENSG[match(dt_geneMod$ensembl, dt_mapping$ENSMUS)]

# pull out module genes as a list of vectors
list_modGenes_ENSG <- lapply(unique(dt_geneMod$module), function(module) {
  vec_mod_ENSG <- dt_geneMod$ENSG[dt_geneMod$module==module]
  vec_mod_ENSG <- vec_mod_ENSG[!is.na(vec_mod_ENSG)]
  vec_mod_ENSG
  })

names(list_modGenes_ENSG) <- unique(dt_geneMod$module)

# save to file
saveRDS(object = list_modGenes_ENSG, file = here("data", "list_vec_WGCNA_modgenes_ENSG.RDS"))
