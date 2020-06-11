#Libraries use for this activity
library("seqinr")
library("rBLAST")
library("R.utils")

#Downloading the E.coli CDS sequence from the Ensembl FTP page
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz",
              destfile = "ecoli.fa.gz")

#Uncompressing the file
gunzip("ecoli.fa.gz")

#Creating the blast database 
makeblastdb("ecoli.fa",dbtype="nucl", "-parse_seqids")

#Download the sample file
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa",
              destfile = "sample.fa")

#Read the sample file into R
d <- read.fasta("sample.fa")
mygene <- d[[3]]
mygene

#Length of my gene sequence in bp
str(mygene)
length(mygene)

#Proportion of GC content
seqinr::GC(mygene)

#Functions to create blast
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/mutblast_functions.R",
              destfile = "mutblast.R")
source("mutblast.R")

#Blast search for E. coli genes that matches best with my gene sequence
res <- myblastn_tab(myseq = mygene, db = "ecoli.fa")

#View my Blast results
res
View(res)

str(res)
head(res)

#Making point mutations to my_gene sequence
mygene_mutation <- mutator(mygene,20)
res <- myblastn_tab(myseq = mygene_mutation, db = "ecoli.fa")
res

#Writing the blast index
write.fasta(mygene,names = "my gene",file.out = "mygene.fa")
makeblastdb(file="mygene.fa",dbtype = "nucl")

#Mutation testing with 100 mismatches
mygene_mutation <- mutator(myseq=mygene,100)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#Mutation testing with 150 mismatches
mygene_mutation <- mutator(myseq=mygene,150)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#Mutation testing with 200 mismatches
mygene_mutation <- mutator(myseq=mygene,200)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#Mutation testing with 250 mismatches
mygene_mutation <- mutator(myseq=mygene,250)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#Mutation testing with 300 mismatches
mygene_mutation <- mutator(myseq=mygene,450)
res <- myblastn_tab(myseq = mygene_mutation, db = "mygene.fa")
res

cointoss <- function(mygene_mutation) {
  sample(c(0,1),1,replace = TRUE)
}

mean(replicate(100,cointoss(mygene_mutation)))

#Plot for bit score
plot(res$bitscore)
