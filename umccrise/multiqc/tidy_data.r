#!/usr/bin/env Rscript

library("TidyMultiqc")
library("arrow")
#library("paws")

if (length(args)==0) {
  stop("At least one argument for multiqc_data.json input must be specified", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "multiqc_data.parquet"
}

df = load_multiqc(args[1])
write_parquet(df, args[2])
