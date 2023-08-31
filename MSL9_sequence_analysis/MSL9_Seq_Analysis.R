#This is the code used to perform the following steps of analyzing sequence variation of the M-SL9 epitope:
#   A. Import FASTA files obtained from BV-BRC containing the M1-coding nucleotide sequence
#   B. Translate the nucleotide sequence to amino acids and crop at just the M-SL9 variant seq plus preceding codon.
#   C. Make sequence logos for each time period and flu subtype.

###  !!!! READ ME !!!!
###  BEFORE RUNNING THIS PROGRAM, you must do these 2 things:
###  1. Set the working directory to a folder that contains all of your target fastas to turn into seqlogos in one batch.
###  2. Edit the code in STEP 3 so that the number of seqlogo_N lines shown equals the number of your target fasta files. 

library(ggseqlogo)
library(tools)
library(seqinr)
library(bioseq)
library(ggplot2)
library(ggeasy)
library(stringr)
library(cowplot)

#STEP 1: Pull the file names for all FASTA files and store them as a character vector
fluDNAFastaVector <- list.files()


#STEP 2: Use a for loop to read in FASTA files; test for sequences >300 nt starting with Met (M1 protein);
#pull out 10 codons starting 89 nt from start (M-SL9 opening reading frame); and convert a.a. into seqLogos:
Nmax <- length(fluDNAFastaVector)
N=1

for (eachFileName in fluDNAFastaVector) {
  if (N<=Nmax)
DNAvectorN <- assign(file_path_sans_ext(fluDNAFastaVector[N]), read_fasta(fluDNAFastaVector[N], type="DNA"))
M1TransN <- seq_translate(DNAvectorN, codon_frame = 1)
M1MetTestN <- ifelse(startsWith(M1TransN,"M"), TRUE, FALSE)
DNAvectorMetClean <- DNAvectorN[M1MetTestN]
DNAvectorMetClean_asVec <- as.vector(DNAvectorMetClean)
M1lengthtestN <- ifelse(nchar(DNAvectorMetClean_asVec) > 300, TRUE, FALSE)
DNAvectorM1Clean <- DNAvectorMetClean[M1lengthtestN]
AAvectorN <- assign(paste0(file_path_sans_ext(fluDNAFastaVector[N]),"_aa"), seq_translate(DNAvectorM1Clean, codon_frame = 89))
AAvectorCropN <- assign(paste0(file_path_sans_ext(fluDNAFastaVector[N]),"_aa_crop"), seq_crop_position(AAvectorN, position_in = 1, position_out = 10))
AAvectorCropN <- assign(paste0(file_path_sans_ext(fluDNAFastaVector[N]),"_aa_vector"), as.vector(AAvectorCropN))
BooleanVectorN <- ifelse(nchar(AAvectorCropN) == 10, TRUE, FALSE)
AAvectorCropCleanN <- assign(paste0(file_path_sans_ext(fluDNAFastaVector[N]),"_aa_vector_clean"), AAvectorCropN[BooleanVectorN])
seqlogoN <- assign(paste0(file_path_sans_ext(fluDNAFastaVector[N]),"_seqlogo"),
                 ggplot()
                  + geom_logo(AAvectorCropCleanN, method="p", font = "helvetica_bold")
                  + theme_logo()
                  + xlab("")
                  + ylab("")
                  + theme(legend.position = "none")
                  + scale_x_continuous(breaks=1:10, labels=c(-1,1,2,3,4,5,6,7,8,9))
                  + theme(axis.text = element_text(size = 11))
                  + theme(axis.text.x = element_text(vjust = 2))
                  + ggtitle(str_split_i(file_path_sans_ext(fluDNAFastaVector[N]),"_M1",1))
                  + theme(plot.title = element_text(hjust = 0.4)) 
                  + theme(plot.title=element_text(size=12))
                   )
seqlogo_Numbered <- assign(paste0("seqlogo_",N), seqlogoN)
N <- N+1
}


#STEP 3: graph all the seqlogos in a grid
#Commment/uncomment lines as needed to make the number of seqlogo lines equal to the number of fasta files
cowplot::plot_grid(
  seqlogo_1
   ,seqlogo_2
   ,seqlogo_3
   ,seqlogo_4
   ,seqlogo_5
   ,seqlogo_6
  # ,seqlogo_7
  # ,seqlogo_8
  # ,seqlogo_9
  # ,seqlogo_10
  
  #Choose how many columns you want to display your seqlogo plots in:
  ,ncol = 1
                  )
