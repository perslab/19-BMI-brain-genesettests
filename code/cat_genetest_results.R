# contatenate es geneset enrichment results
library("here")
library("data.table")
library("openxlsx")

vec_files = dir(here("output"), pattern = "23_genes_wilcoxon_genetestOuts_200613.csv.gz", full.names = T)

dt_out = Reduce(x=lapply(vec_files, fread), f=rbind)

file.out = here("output", "all_BMI_23_genes_wilcoxon_genetestOuts_200613.csv")
fwrite(dt_out, file.out)
system2(command = "gzip", args = file.out)

write.xlsx(dt_out, here("output", "all_BMI_23_genes_wilcoxon_genetestOuts_200613.xlsx"))
