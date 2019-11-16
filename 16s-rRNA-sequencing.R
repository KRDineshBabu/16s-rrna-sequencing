install.packages("filesstrings")
install.packages("seqinr")
install.packages("seqRFLP")

# Bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("sangerseqR")

install.packages("devtools")
library(devtools)
install_github("roblanf/sangeranalyseR")

library(filesstrings)
library(sangerseqR)
library(sangeranalyseR)
library(seqinr)
library(seqRFLP)

seq_path = "C:\\Users\\asckrdb\\OneDrive - Temasek Polytechnic\\Soil benefical bacteria\\Result\\Sequencing\\sanger sequencing\\H and I"
output = "C:\\Users\\asckrdb\\OneDrive - Temasek Polytechnic\\Soil benefical bacteria\\Result\\Sequencing\\sanger sequencing\\H and I\\output"

seq_files = list.files(path = seq_path, pattern = ".ab1")
file.rename(seq_files, paste0(gsub("F_.*","F.ab1",seq_files)))
file.rename(seq_files, paste0(gsub("R_.*","R.ab1",seq_files)))
seq_files_fwd = list.files(path = seq_path, pattern = "F.ab1")


for (file in seq_files_fwd){
  setwd(seq_path)
  seq_read = readsangerseq(file)
  file_name = basename(file)
  filepdf = give_ext(file_name,"pdf", replace = TRUE)
  basecall =  makeBaseCalls(seq_read, ratio = 0.33)
  setwd(output)
  chromato =  chromatogram(basecall, width = 100, height = 2, trim5 = 250, trim3 = 300, showcalls = "both", filename = filepdf)
  rmext = sub(pattern = ".ab1", '',file_name)
  primeryseq = primarySeq(basecall)
  
  if(nchar(primeryseq)>810){
    #trim = sread(primeryseq)
    # write.fasta(sequences = trim, names = give_ext(rmext, ".prim"), file.out = give_ext(rmext, ".trim.fasta") )
    
    trim_primery = substring(primeryseq, 250, 800)
    write.fasta(sequences = trim_primery, names = give_ext(rmext, ".prim"), file.out = give_ext(rmext, ".prim.fasta") )
    secondary = secondarySeq(basecall)
    trim_sec = substring(secondary, 250, 800)
    write.fasta(sequences = trim_sec, names = give_ext(rmext, ".sec"), file.out = give_ext(rmext, "sec.fasta") )
  }
}

setwd(seq_path)   
seq_files_rev = list.files(path = seq_path, pattern = "R.ab1")

for (file in seq_files_rev){
  setwd(seq_path)
  
  seq_read = readsangerseq(file)
  file_name = basename(file)
  file = give_ext(file_name,"pdf", replace = TRUE)
  basecall =  makeBaseCalls(seq_read, ratio = 0.33)
  setwd(output)
  chromato =  chromatogram(basecall, width = 100, height = 2, trim5 = 250, trim3 = 300, showcalls = "both", filename = file)
  rmext = sub(pattern = ".ab1", '',file_name)
  primeryseq = primarySeq(basecall)
  rev_comp = reverseComplement(primeryseq)
  if(nchar(rev_comp)>810){
    trim_primery = substring(rev_comp, 150, 800)
    write.fasta(sequences = trim_primery, names = give_ext(rmext, ".prim"), file.out = give_ext(rmext, ".prim.fasta") )
    secondary = secondarySeq(basecall)
    trim_sec = substring(secondary, 150, 800)
    write.fasta(sequences = trim_sec, names = give_ext(rmext, ".sec"), file.out = give_ext(rmext, "sec.fasta") )
  }
}

file_list = list.files(path = output, pattern = "prim.fasta")
setwd(output)
if(!exists("dataset")){
  
  
  for (file in file_list){
    
    fastaFile <- readDNAStringSet(file)
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    df <- data.frame(seq_name, sequence)
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- df
      
    }
    else{
      # if the merged dataset does exist, append to it
      temp_dataset <-df
      dataset<-rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
  }
  
} else {
  
  rm(dataset)
  
  for (file in file_list){
    
    fastaFile <- readDNAStringSet(file)
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    df <- data.frame(seq_name, sequence)
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- df
      
    }
    
    # if the merged dataset does exist, append to it
    else {
      temp_dataset <-df
      dataset<-rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
  }
}

dataset_fasta =  dataframe2fas(dataset)
detach(package : seqinr)
subDir <- "output_combine"

if (file.exists(subDir)){
  setwd(file.path(output, subDir))
} else {
  dir.create(file.path(output, subDir))
  setwd(file.path(output, subDir))
}
write.fasta(dataset_fasta, file = "dataset_fasta_primery.fasta")

setwd(output)
file_list = list.files(path = output, pattern = "sec.fasta")


if(!exists("dataset")){
  
  for (file in file_list){
    
    fastaFile <- readDNAStringSet(file)
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    df <- data.frame(seq_name, sequence)
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- df
      
    }
    
    # if the merged dataset does exist, append to it
    else{
      temp_dataset <-df
      dataset<-rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
  }
  
} else {
  
  rm(dataset)
  
  for (file in file_list){
    
    fastaFile <- readDNAStringSet(file)
    seq_name = names(fastaFile)
    sequence = paste(fastaFile)
    df <- data.frame(seq_name, sequence)
    
    # if the merged dataset doesn't exist, create it
    if (!exists("dataset")){
      dataset <- df
      
    }
    
    # if the merged dataset does exist, append to it
    else{
      temp_dataset <-df
      dataset<-rbind(dataset, temp_dataset)
      rm(temp_dataset)
    }
  }
}

dataset_fasta =  dataframe2fas(dataset)

if (file.exists(subDir)){
  setwd(file.path(output, subDir))
} else {
  dir.create(file.path(output, subDir))
  setwd(file.path(output, subDir))
}
write.fasta(dataset_fasta, file = "dataset_fasta_secondary.fasta")
