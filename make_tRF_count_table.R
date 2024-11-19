#!/usr/bin/env Rscript
## Checks if the necessary packages are installed and installs them if they are not
if (!requireNamespace("argparse", quietly = TRUE)) {
  install.packages("argparse")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  install.packages("readxl")
}

library("argparse")
library("readr")
library("readxl")

## Parse arguments
parser <- ArgumentParser()
parser$add_argument("-m", "--mintmap", type="character",
                    help="The path of the folder containing MINTmap output files",
                    required=TRUE)
parser$add_argument("-b", "--multicov", type="character",
                    help="The path of the folder containing BEDtools multicov output files",
                    required=TRUE)
parser$add_argument("-f", "--fasta", type="character",
                    help="The path of the FASTA file containing the sequences of the 5' leader and 3' trailer tRFs (Human_hg38lift_tsRNA_and_5leader_final_edited.fasta)",
                    required=TRUE)
parser$add_argument("-n", "--nci60", type="character",
                    help="The path of the Excel file containing the NCI-60 metadata (NCI60_metadata.xlsx)",
                    required=TRUE)
parser$add_argument("-o", "--output", type="character",
                    help="The path of the output directory",
                    required=TRUE)

args <- parser$parse_args()

mintmap_out_path <- args$mintmap
if (!dir.exists(mintmap_out_path)) {
  stop("The path of the folder containing MINTmap output files does not exist.")
}

bedtools_out_path <- args$multicov
if (!dir.exists(bedtools_out_path)) {
  stop("The path of the folder containing BEDtools multicov output files does not exist.")
}

leader_trailer_tRF_seq_path <- args$fasta
file <- file.info(leader_trailer_tRF_seq_path)
if (file$isdir) {
  leader_trailer_tRF_seq_path <- file.path(leader_trailer_tRF_seq_path, "Human_hg38lift_tsRNA_and_5leader_final_edited.fasta")
}
if (!file.exists(leader_trailer_tRF_seq_path)) {
  stop("The FASTA file containing the sequences of the 5' leader and 3' trailer tRFs does not exist.")
}

NCI60_metadata_path <- args$nci60
file <- file.info(NCI60_metadata_path)
if (file$isdir) {
  NCI60_metadata_path <- file.path(NCI60_metadata_path, "NCI60_metadata.xlsx")
}
if (!file.exists(NCI60_metadata_path)) {
  stop("The Excel file containing the NCI-60 metadata does not exist.")
}

output_path <- args$output
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
}
if (!dir.exists(output_path)) {
  stop("The output directory does not exist and could not be created.")
}



## Import MINTmap output - tRFs

mintmap_out_files <- list.files(mintmap_out_path, pattern = "-MINTmap_v2.0-alpha-exclusive-tRFs.expression.txt")
mintmap_out_matrix <- read_delim(file.path(mintmap_out_path, mintmap_out_files[1]), 
                                 delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 6)

mintmap_out_matrix <- mintmap_out_matrix[,c(1:4)]
colnames(mintmap_out_matrix)[4] <- gsub("-MINTmap_v2.0-alpha-exclusive-tRFs.expression.txt","",mintmap_out_files[1])

i <- 2

while(i<=length(mintmap_out_files)){
  tmp <- read_delim(file.path(mintmap_out_path, mintmap_out_files[i]), 
                    delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 6)
  tmp <- tmp[,c(1:4)]
  colnames(tmp)[4] <- gsub("-MINTmap_v2.0-alpha-exclusive-tRFs.expression.txt","",mintmap_out_files[i])
  mintmap_out_matrix <- merge(x = mintmap_out_matrix, y = tmp, by = c("License Plate", "tRF sequence",  "tRF type(s)"), all = TRUE)
  i <- i + 1
  rm(tmp)
}

rm(i)

mintmap_out_matrix[is.na(mintmap_out_matrix)] <- 0

colnames(mintmap_out_matrix)[1:3] <- c("id", "sequence", "tRNA_derived_type")

## Import BEDtools multicov output - tsRNAs and 5' leader sequences
bedtools_out_files <- list.files(bedtools_out_path, pattern = ".txt")

bedtools_out_matrix <- read_delim(file.path(bedtools_out_path, bedtools_out_files[1]), 
                                  delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

bedtools_out_matrix <- bedtools_out_matrix[,c(1,4)]
colnames(bedtools_out_matrix)[2] <- gsub(".txt","",bedtools_out_files[1])

i <- 2

while(i<=length(bedtools_out_files)){
  tmp <- read_delim(file.path(bedtools_out_path, bedtools_out_files[i]), 
                    delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  tmp <- tmp[,c(1,4)]
  colnames(tmp)[2] <- gsub(".txt","",bedtools_out_files[i])
  bedtools_out_matrix <- merge(x = bedtools_out_matrix, y = tmp, by = "X1")
  i <- i + 1
  rm(tmp)
}

rm(i)

colnames(bedtools_out_matrix)[1] <- "id"
bedtools_out_matrix$tRNA_derived_type <- 0

r <- 1

while(r <= nrow(bedtools_out_matrix)){
  if(substring(bedtools_out_matrix[r,1],1,2)=="5P"){
    bedtools_out_matrix[r,ncol(bedtools_out_matrix)] <-"5_leader_tRF"
  } else {bedtools_out_matrix[r,ncol(bedtools_out_matrix)] <- "3_trailer_tRF"}
  r <- r + 1
}

rm(r)

leader_trailer_seq <- read_delim(leader_trailer_tRF_seq_path, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

seq <-  leader_trailer_seq[grep(pattern = ">", leader_trailer_seq$X1) + 1,1]
colnames(seq) <- "sequence"

colnames(leader_trailer_seq)[1] <- "id"
leader_trailer_seq <- leader_trailer_seq[grep(pattern = ">", leader_trailer_seq$id),]

leader_trailer_seq <- cbind(leader_trailer_seq, seq)
leader_trailer_seq$id <- gsub(">", "", leader_trailer_seq$id)

bedtools_out_matrix <- merge(x = bedtools_out_matrix, y = leader_trailer_seq, by = "id")

bedtools_out_matrix <- bedtools_out_matrix[,c(1,ncol(bedtools_out_matrix),ncol(bedtools_out_matrix)-1,2:c(ncol(bedtools_out_matrix)-2))]

## Combine MINTmap and BEDtools outputs and extract

final_table <- rbind(mintmap_out_matrix, bedtools_out_matrix)

NCI60_metadata <- read_excel(NCI60_metadata_path)
NCI60_metadata <- cbind(NCI60_metadata, col_names = colnames(final_table[,c(4:ncol(final_table))]))

colnames(final_table)[4:ncol(final_table)] <- NCI60_metadata$`Sample name`

write.table(final_table,
            file = file.path(output_path, "tRFs_harmonized_raw_count_table.txt"),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


