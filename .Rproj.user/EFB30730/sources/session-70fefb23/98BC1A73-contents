#Packages
library(phylotools)
library(Biostrings)
library(deepredeff)

install.packages("deepredeff")

#Making into compatible DNAstringset
sequence_vector <- sinv_assembly$seq.text
names(sequence_vector) <- sinv_assembly$seq.name
seq_DNAStringset <- DNAStringSet(x=sequence_vector, start=NA, end=NA, width=NA, use.names=TRUE)

#Translating
seq_translation <- translate(x=seq_DNAStringset, genetic.code=GENETIC_CODE, if.fuzzy.codon="error")

#Converting back to dataframe
seq_translation_dataframe <- aasset_to_df(aas=seq_translation)

#Converting dataframe to fasta and saving to specific name
dir.create("amino_acid_sequences")
dat2fasta(dat=seq_translation_dataframe, outfile = "./amino_acid_sequences/sinv_AA_translation.fasta")

#Storing Assembling
sinv_assembly <- read.fasta(file = "./genomes_cleaned/sinv_genomeAssembly.fna")


