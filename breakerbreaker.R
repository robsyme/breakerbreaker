#!/usr/bin/env Rscript

library(argparser, quietly=TRUE)


p <- arg_parser("Identify Patterns of Anomalous Reads")
p <- add_argument(p,
    "bam",
    help = "BAM file to parse",
    type = "character")

argv <- parse_args(p)


