#Description: This is a script to compare experimentally determined deconvoluted masses to a computationally digested list of transcripts.
#First, digests a fasta file of transcript sequences and calculates the resultant masses.
#Next, compares and scores each transcript based upon a similarity score to observed deconvoluted masses.
#Author: Eric J. Wolf, New England Biolabs
#Date: January 24, 2022

#Load Dependencies
library(seqinr)
library(stringr)
library(hash)

Nucleotides <- hash(keys = c("C","U","G","A"), 
                    values = c(305.0413, 306.0253, 345.0474, 329.0525)) #Hash of individual nucleoside masses


seq1 <- seqinr::read.fasta("refMrna.fa", as.string = TRUE, forceDNAtolower = FALSE) #Input the fasta file for searching
seq1 <- lapply(seq1, function(x) toupper(x)) #Make all sequences in the fasta file upper case
seq1 <- lapply(seq1, function(x) str_replace_all(x, "T", "U")) #Convert any DNA 'T's to RNA 'U's

CleavageProducts <- lapply(seq1, function(x) strsplit(x, "(?<=[U])(?=[A|G])", perl = TRUE)) #Completely digest the fasta file with the specificity of hRNase4 
rm(seq1) #Get rid of the original sequence file to save memory

#A function to give each cleavage product a name and filer the list by length >=4
AnnotateT1Products <- function(CleavageProducts) {
  length.v  <- nchar(CleavageProducts[[1]])
  cumsum.v <- cumsum(length.v)
  inclusion.bool <- length.v >= 4
  seqlist <- list(CleavageProducts[[1]][inclusion.bool])
  names(seqlist[[1]]) <- as.character(cumsum.v[inclusion.bool])
  return(seqlist)
}

RNase4 <- lapply(CleavageProducts, function(x) AnnotateT1Products(x)) #A list of annotated cleavage products produced by the AnnotatedT1Products fucntion
rm(CleavageProducts) #Get rid of the original list to save memory 

#A function to detemine the mass of each cleavage product in the RNase4 list
CountMass <- function(myHash, mySeq) {
  count <- 18.0106 - 79.9663 #plus water - phosphate
  for (key in keys(myHash)) { #add the mass of the number of occurances of each nucleotide in each cleavage product
    count <- count + str_count(mySeq, key)*myHash[[key]] #sum the massof the each nucleotide in each cleavage product
  }
  return(count)
}

RNase4mass <- lapply(RNase4, function(x) CountMass(Nucleotides, unlist(x))) #A list of all masses of cleavage products
rm(RNase4) #Get rid of the original list to save memory 

#A function to compare the experimentally derived masses to theoretical intensities.
#Takes an experimental mass v intensity dataframe as an input and the RNase4mass list.
MatchMasses <- function(Experimental, Theoretical, cutoff, outputType) {
  i <- 0
  Intensity <- 0
  all <- length(unique(Theoretical))
  TheoMatches <- c()
  for (Mass in Experimental$Mass) { #loop through experimentally derived masses
    i <- i + 1
    min <- which.min(abs(Mass - Theoretical)) #find the theoretical mass which best matches the experimental mass
    PPM <- abs(10^6*((Theoretical[min] - Mass)/Mass)) #calculate the ppm mass error
    if (PPM < cutoff) { #if the theoretical mass is within a PPM tolerance include it as a match
      Intensity <- Intensity + Experimental[i,"Intensity"]
      TheoMatches <- append(TheoMatches, Theoretical[min])
    }
  }
  if (outputType == "Score") { #generate a score of the product of the fraction of Intensity explained and the fraction of all theoretical oligos identified
    FractionIntensity <- (Intensity/sum(Experimental$Intensity)) * (length(unique(TheoMatches))/all)
  } 
  #Other metrics which may be useful to assess the identity of a deconvoluted mass spectrum
  else if (outputType == "NumTotal") {
    FractionIntensity <- all
  }
  else if (outputType == "FractionFound") {
    FractionIntensity <- (length(unique(TheoMatches))/all)
  }
  else if (outputType == "FractionIntensity") {
    FractionIntensity <- (Intensity/sum(Experimental$Intensity))
  }
  else if (outputType == "Intensity") {
    FractionIntensity <- (Intensity)
  }
  else {
    stop("You have entered an invalid output type")
  }
  return(FractionIntensity)
}

#Run the functions against a mass v intensity dataset
Score <- lapply(RNase4mass, function(x) MatchMasses(MatchedMasses, x, 5, "Score"))

#Make a dataframe of the transcript scores
Score <- data.frame(Score = unlist(Score), Transcripts = names(Score))

#Output the scores for processing outside of R
dir.create("Tables")

write.table(Score, "Tables/TranscriptScores.txt", sep = "\t", quote = FALSE, row.names = FALSE)
