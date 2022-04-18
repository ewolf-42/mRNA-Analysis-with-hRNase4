#Script to generate cleavage products from a list of 1000 randomly selected human transcripts with endoribonucleases previously utilized for RNA LC-MS/MS
#Author: Eric J. Wolf, New England Biolabs
#Date: January 24, 2022

#Load in libraries
library(seqinr)
library(stringr)
library(hash)
library(ggplot2)
library(cowplot)
library(data.table)

#Load in Fasta File Function 
InputFasta <- function(fasta) {
  seq1 <- seqinr::read.fasta(fasta, as.string = TRUE, forceDNAtolower = FALSE)
  seq1 <- lapply(seq1, function(x) toupper(x))
  seq1 <- lapply(seq1, function(x) str_replace_all(x, "T", "U"))
  return(seq1)
}

#Cleavage motifs with enzyme name
ribonucleases <- list(c("(?<=[G])"), c("(?<=[U|C])"), c("(?<=[A|U|C])(?=[U])"), c("(?=[U])"), c("(?<=[U|A|C])(?=[U])","(?<=[C])(?=[U|G])"), c("(?<=[C])(?=[U|G|A])"), c("(?<=[G])(?=[U])"), c("(?<=[U])(?=[A|G])"))
names(ribonucleases) <- c("RNaseT1", "RNaseA", "MC1-2021", "MC1-2015", "Cusativin-2021", "Cusativin-2017", "Colicin-E5", "RNase4")

#Function to Annotate the position of each product in the context of each sequence
AnnotateT1Products <- function(CleavageProducts,IncludeLengths = TRUE) {
  seqlist <- list()
  for (Enzyme in names(CleavageProducts)) {
    if (IncludeLengths == TRUE) {
      length.v  <- nchar(CleavageProducts[[Enzyme]][[1]])
      inclusion.bool <- length.v >= 4  & length.v < 40
      seqlist[[Enzyme]] <- CleavageProducts[[Enzyme]][[1]][inclusion.bool]
    }
    else {
      seqlist[[Enzyme]] <- list(CleavageProducts)
    }
  }
  return(seqlist)
}

#Function to cleave products
Cleavageproducts <- function (ribonucleases, FastaFile){
  count <- 0
  CleavageProducts <- list()
  for (ribo in ribonucleases) {
    count <- count + 1
    Split <- list()
    for (i in 1:length(ribo)) {
      print(i)
      if (i == 1) {
        Split <- str_split(unlist(FastaFile), pattern = ribo[i])
      }
      else {
        Split <- str_split(unlist(Split), pattern = ribo[i])
        Split<-list(unlist(Split))
      }
    }
    CleavageProducts[[names(ribonucleases)[count]]] <- Split
  }
  CleavageProducts [["NA"]] <- NULL
  return(CleavageProducts)
}

#Function to make data frames 
MakeDataframes <- function(list1) {
  Outputlist <- list()
  i <- 0
  for (element in list1) {
    i <- i + 1
    df <- data.frame(sequences = element, stringsAsFactors = FALSE)
    Outputlist[[names(list1)[i]]] <- df
  }
  return(Outputlist)
}

#Function to count all nucleotides in Fasta File with Duplicate check 
GetNucRedundancy <- function(seqs, Removeduplicates = FALSE) {
  nucs <- c("A","C", "G", "U")
  for (Enzyme in names(seqs)) {
    seqs[[Enzyme]]$NucCompose <- ""
      for (nuc in nucs) {
        nuccount <- str_count(seqs[[Enzyme]][,"sequences"], pattern = nuc)
        nucstring <- paste(nuc, as.character(nuccount), sep = "")
        seqs[[Enzyme]][,"NucCompose"] <- paste(seqs[[Enzyme]][,"NucCompose"],nucstring, sep = ":")
      }
      seqs[[Enzyme]]$Duplicate <- duplicated(seqs[[Enzyme]]$NucCompose,fromLast = TRUE)
      seqs[[Enzyme]]$Duplicate2 <- duplicated(seqs[[Enzyme]]$NucCompose)
      seqs[[Enzyme]]$DuplicatedCheck <- FALSE
      seqs[[Enzyme]][seqs[[Enzyme]]$Duplicate == TRUE | seqs[[Enzyme]]$Duplicate2 == TRUE, "DuplicatedCheck"] <- TRUE
      seqs[[Enzyme]]$Duplicate <- duplicated(seqs[[Enzyme]]$sequences)
      seqs[[Enzyme]]$Duplicate2 <- duplicated(seqs[[Enzyme]]$sequences, fromLast = TRUE)
      seqs[[Enzyme]] <- seqs[[Enzyme]][seqs[[Enzyme]]$Duplicate == FALSE & seqs[[Enzyme]]$Duplicate2 == FALSE,]
      seqs[[Enzyme]]$Duplicate <- NULL 
      seqs[[Enzyme]]$Duplicate2 <- NULL 
      
  }
  return(seqs)
}

#Compute Coverages For Each Transcript With Each Enzyme
ComputeCoverage <- function(Products, lengths, Isomeric) {
  output <- list()
  for (Product in names(Products)) {
    for (Enzyme in names(Products[[Product]])) {
      if (Isomeric == "Total") {
        output[[Enzyme]][[Product]] <- sum(nchar(Products[[Product]][[Enzyme]][, "sequences"]))/lengths[Product]
      }
      if (Isomeric == "Unique") {
        output[[Enzyme]][[Product]] <- sum(nchar(Products[[Product]][[Enzyme]][Products[[Product]][[Enzyme]]$DuplicatedCheck == FALSE, "sequences"]))/lengths[Product]
      }
    }
  }
  return(output)
}

#Function to combine lists into a dataframe for graphing
CombineLists <- function(list1) {
  out <- data.frame(Enzyme = "", Coverage = 0, stringsAsFactors = FALSE)
  for(Enzyme in names(list1)) {
    add <- data.frame(Enzyme = Enzyme, Coverage = unlist(list1[[Enzyme]]))
    out <- merge(out, add, all = TRUE)
  }
  out <- out[out$Enzyme != "",]
  return(out)
}

#Function To Plot Coverage Distribution
PlotCoverageBox <- function(df) {
  plot <- ggplot(data = df, aes(x = Enzyme, y=Coverage, fill = Enzyme, colour = Enzyme)) +
    geom_boxplot(fill = NA, outlier.size = 0.2) +
    facet_grid(cols = vars(Unique)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          axis.title = element_text(size = 8,family = "sans", colour = "Black"),
          legend.key = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(size = 8,family = "sans", colour = "Black"),
          axis.text.x = element_text(size = 8,family = "sans", colour = "Black", vjust = 0.5, hjust = 1, angle = 90),
          title = element_text(size = 8, family = "sans", colour = "Black")) +
    xlab("Enzyme") +
    ylab("Theoretical Sequence Coverage")
  return(plot)
}

####Run the Functions######
#Fasta input
FastaFile <- InputFasta("refMrna.fa")
FastaFile <- FastaFile[names(FastaFile) %like% "NM"]
FastaFile <- FastaFile[nchar(FastaFile) < 5000]
set.seed(250)
select <- FastaFile[sample(1:length(FastaFile), 1000, replace =F)]
lengths <- nchar(select)

#Cleave products
CleavageProducts <- lapply(select, function(x) Cleavageproducts(ribonucleases, x))

#Apply Annotation function 
EnzymeCollection <- lapply(CleavageProducts, function(x) AnnotateT1Products(x,TRUE)) 

#Create enzyme collection data frame 
ProductList <- lapply(EnzymeCollection, function(x) MakeDataframes(x))
Products <- lapply(ProductList, function(x) GetNucRedundancy(x, TRUE))

#Calculate Coverages
CoverageOut <- ComputeCoverage(Products, lengths, "Total")
CoverageOut.unique <- ComputeCoverage(Products, lengths, "Unique")
CoverageOut.df <- CombineLists(CoverageOut)
CoverageOut.unique.df <- CombineLists(CoverageOut.unique)
CoverageOut.unique.df$Unique <- "Unique"
CoverageOut.df$Unique <- "Total"
CoverageOut.df <- merge(CoverageOut.df, CoverageOut.unique.df, all = TRUE)


CoverageBox <- PlotCoverageBox(CoverageOut.df)

dir.create("Graphs")

pdf("Graphs/CoverageBox.2.pdf")  
ggdraw() +
  draw_plot(CoverageBox, x = 0, y = 0, width = 1/2, height = 1/3)
dev.off()

