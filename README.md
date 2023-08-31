# Data and code to support Hogan et al. Nature Immunology 2023

Full reference: 

Michael J. Hogan, Nikita Maheshwari, Bridget E. Begg, Annalisa Nicastri, Emma J. Hedgepeth, Hiromi Muramatsu, Norbert Pardi, Michael A. Miller, Shanelle P. Reilly, Laurent Brossay, Kristen W. Lynch, Nicola Ternette, and Laurence C. Eisenlohr. "Cryptic MHC-E epitope from influenza elicits a potent cytolytic T cell response." Nature Immunology (2023)

## TCR-beta sequence data and code to analyze in R

Explanation: six C57Bl/6 mice were infected with influenza A virus PR8 and lungs were collected on day 9 post-infection. Antigen-specific or naive T cells were distinguished by staining with NP366/Db and M-SL9/Qa-1 tetramers and FACS-sorted, and then genomic DNA was isolated and shipped to Adaptive Biotechnologies for mouse TCR-beta sequencing. Sequence data is available as a .tsv file for each combination of mouse and T cell specificity (18 total). The accompanying R code shows how to generate the PCA plot (Extended Data Fig. 8c), Pearson correlation coefficients (Extended Data Fig. 8d), and "Spectratype" graphs (Fig. 3e) plotting the average frequency of TCRBV gene usage by CDR3 length.

## M-SL9 sequence analysis

Explanation: This R code (single file) shows how to analyze fasta files containing matrix protein 1-coding nucleotide sequences obtained from the BV-BRC database to pull out the amino acids corresponding to M-SL9 (and its variants) from influenza A virus strains of interest and assemble seqlogo diagrams. This made use of the ggseqlogo package referenced below:

https://omarwagih.github.io/ggseqlogo/
Wagih, Omar. ggseqlogo: a versatile R package for drawing sequence logos. Bioinformatics (2017).
https://doi.org/10.1093/bioinformatics/btx469

