## Import MINTmap output - tRFs
library(readr)

mintmap_out_path <- "/mintmap_output_NCI60/" # <- Insert path folder containing MINTmap output files

mintmap_out_files <- list.files(mintmap_out_path, pattern = "-MINTmap_v2.0-alpha-exclusive-tRFs.expression.txt")
mintmap_out_matrix <- read_delim(paste(mintmap_out_path, mintmap_out_files[1], sep = ""), 
                                 delim = "\t", escape_double = FALSE, trim_ws = TRUE, skip = 6)

mintmap_out_matrix <- mintmap_out_matrix[,c(1:4)]
colnames(mintmap_out_matrix)[4] <- gsub("-MINTmap_v2.0-alpha-exclusive-tRFs.expression.txt","",mintmap_out_files[1])

i <- 2

while(i<=length(mintmap_out_files)){
  tmp <- read_delim(paste(mintmap_out_path, mintmap_out_files[i], sep = ""), 
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
bedtools_out_path <- "/bedtools_multicov_output_NCI60/" # <- Insert path folder containing BEDtools output files
bedtools_out_files <- list.files(bedtools_out_path, pattern = ".txt")

bedtools_out_matrix <- read_delim(paste(bedtools_out_path, bedtools_out_files[1], sep = ""), 
                                  delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)

bedtools_out_matrix <- bedtools_out_matrix[,c(1,4)]
colnames(bedtools_out_matrix)[2] <- gsub(".txt","",bedtools_out_files[1])

i <- 2

while(i<=length(bedtools_out_files)){
  tmp <- read_delim(paste(bedtools_out_path, bedtools_out_files[i], sep = ""), 
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

leader_trailer_tRF_seq_path <- "/Human_hg38lift_tsRNA_and_5leader_final_edited.fasta" # <- Add path FASTA file reporting the sequences of the 5' leader and 3' trailer tRFs

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
library(readxl)

output_path <- "/NCI60_data/" # <- Add path output folder

final_table <- rbind(mintmap_out_matrix, bedtools_out_matrix)

NCI60_metadata <- read_excel("/NCI60_metadata.xlsx") # <- Add path excel file containing NCI-60 metadata
NCI60_metadata <- cbind(NCI60_metadata, col_names = colnames(final_table[,c(4:ncol(final_table))]))

colnames(final_table)[4:ncol(final_table)] <- NCI60_metadata$`Sample name`

write.table(final_table,
            file = paste(output_path,"NCI60_tRFs_harmonized_raw_count_table.txt", sep = ""),
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


