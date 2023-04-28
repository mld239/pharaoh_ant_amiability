#Packages
library(phylotools)
library(Biostrings)
library(deepredeff)
library(purrr)
library(stringr)
library(tidyverse)


### Making code into function ###

#Making function
nucleotide_to_AA_translation <- function(input_file) {
  
  #Storing Assembling
  assembly <- read.fasta(file = input_file)
  
  #Making into compatible DNAstringset
  sequence_vector <- assembly$seq.text
  names(sequence_vector) <- assembly$seq.name
  seq_DNAStringset <- DNAStringSet(x=sequence_vector, start=NA, end=NA, width=NA, use.names=TRUE)
  
  #Translating
  seq_translation <- translate(x=seq_DNAStringset, genetic.code=GENETIC_CODE, if.fuzzy.codon="solve")
  
  #Converting back to dataframe
  seq_translation_dataframe <- aasset_to_df(aas=seq_translation)
  colnames(seq_translation_dataframe) <- c("seq.name", "seq.text")
  
  #Converting dataframe to fasta and saving to specific name
  dir.create("amino_acid_sequences")
  species_code <- str_split_i(string = input_file, pattern = "/", i = 3) %>%
    str_split_i(pattern = "_", i = 1)
  output_filename <- paste("./amino_acid_sequences/", species_code, "_genomeAssembly.fna.transdecoder.pep", sep = "")
  dat2fasta(dat=seq_translation_dataframe, outfile = output_filename)
}


#Making list of genomes
genome_assembly_list <- list.files(path = "./genomes_cleaned", pattern = "*.fna", all.files = FALSE,
           full.names = TRUE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)


#Making list of genomes that Transdecoder worked on
Transdecoder_worked <- list.files(path = "./nucleotide_sequences/pep_files", pattern = "*.pep", all.files = FALSE,
                                   full.names = TRUE, recursive = FALSE,
                                   ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)
#Isolating species codes
Transdecoder_worked_species <- str_split_i(Transdecoder_worked, pattern = "/", i = 4) %>% 
  str_split_i(pattern = "_", i = 1)

genome_list_species <- str_split_i(genome_assembly_list, pattern = "/", i = 3) %>%
  str_split_i(pattern = "_", i = 1)


# Only worked w/ first element in Transdecoder_worked #
#desired_files <- genome_list_species[!grepl(Transdecoder_worked_species, genome_list_species)]

#Getting species code on non-worked files
desired_files <- grep(paste(Transdecoder_worked_species,collapse="|"), 
                        genome_list_species, value=TRUE, invert=TRUE)
desired_files <- paste("./genomes_cleaned/",desired_files, "_genomeAssembly.fna",sep="")

possibly_nucleotide_to_AA_translation <- possibly(nucleotide_to_AA_translation, otherwise = "error", quiet = TRUE )
map(desired_files, possibly_nucleotide_to_AA_translation)

#Moving Transdecoder output to AA translations directory as well
file_name <- str_split_i(Transdecoder_worked, pattern = "/", i = 4)
file.copy(from = paste0(Transdecoder_worked),
          to = paste0("./amino_acid_sequences/", file_name))

