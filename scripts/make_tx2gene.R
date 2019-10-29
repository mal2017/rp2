library(GenomicFeatures)
library(readr)

gtf <- snakemake@input[["gtf"]]
output <- snakemake@output[[1]]

# make txdb
txdb <- makeTxDbFromGFF(gtf)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")

tx2gene$GENEID <- vapply(strsplit(tx2gene$GENEID,"\\."),`[[`,1,FUN.VALUE=c("a"))

saveRDS(tx2gene,output)
