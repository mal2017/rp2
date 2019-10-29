library(DESeq2)
library(tximport)

# get info
samps <- snakemake@params[["names"]]
quant_sfs <- unlist(snakemake@input[["qsf"]])
names(quant_sfs) <- samps
output <- snakemake@output[[1]]
tx2gene <- readRDS(snakemake@input[["tx2g"]])
conds <- snakemake@params[["conds"]]
ctrl <- snakemake@params[["ctrl"]]

# make colData
cd <- as.data.frame(t(as.data.frame(conds)))[samps,, drop=F]

# import from Salmon
txi <- tximport(quant_sfs, type = "salmon", tx2gene = tx2gene)

# generate formula
formula <- as.formula(paste0("~",paste(colnames(cd),collapse=" +")))

# make dds
dds <- DESeqDataSetFromTximport(txi, cd, formula)

# save dds
saveRDS(dds, output)

