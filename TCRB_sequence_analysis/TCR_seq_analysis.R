library(immunarch)

#FIRST set working directory to directory where ImmunoSEQ_TCRB_sequences folder is saved.

#load all data into R in Immunarch format using Immunarch function
#FYI: all sequences are already in frame and must have been filtered by ImmunoSEQ
immdata <- repLoad("ImmunoSEQ_TCRB_sequences")

#Analyze TCRB gene usage and statistics, first for all T cells and then focusing only on flu-specific T cells:
imm_VgeneUsage_ALL <- geneUsage(immdata$data, "musmus.trbv", .norm = T, .type = "segment")
imm_VgeneUsage_flu <- geneUsage(immdata$data[c(7:18)], "musmus.trbv", .norm = T, .type = "segment")

#Generate heatmap of Pearson correlation coefficients:
imm_geneUsage_flu_cor <- geneUsageAnalysis(imm_VgeneUsage_flu, .method = "cor",.cor = "pearson", .verbose = F)
vis(imm_geneUsage_flu_cor, .title = "Gene usage correlation", .leg.title = "Pearson coefficient", .text.size = 1.5)

#Generate principal component analysis plot:
imm_geneUsage_PCA_ALL <- geneUsageAnalysis(imm_VgeneUsage_ALL, .method = "pca", .verbose = F)
PCAplot_withLabels <- vis(imm_geneUsage_PCA_ALL, .title = "Gene usage PCA", .leg.title = "PCA", .text.size = 1.5, .text=TRUE)
PCAplot_noLabels <- vis(imm_geneUsage_PCA_ALL, .title = "Gene usage PCA", .leg.title = "PCA", .text.size = 1.5, .text=FALSE)
PCAplot_withLabels
PCAplot_noLabels


#Spectratyping average gene usage frequencies for NAIVE T CELLS, samples 1 to 6 (files beginning with "A"):
n=1
for (n in 1:6) {
  n_countTable = assign(paste0("countTable",n), spectratype(immdata$data[[n]], .quant = "count", .col = "aa+v"))
  n_freqTable <- assign(paste0("freqTable",n), mutate(n_countTable, Val=(Val/sum(n_countTable$Val))/6)) #converts sequence counts into frequencies (out of a total of 1/6 since I will add these together to make a total of 1)
  n_freqTable <- assign(paste0("freqTable",n), mutate(n_freqTable, Length_gene = paste0(Length,"_",Gene)))
}

  joined1to6 <- bind_rows(    #combining frequency tables for samples 1-6; frequencies should now add up to 1
    freqTable1,
    freqTable2,
    freqTable3,
    freqTable4,
    freqTable5,
    freqTable6
    ) 
  joined1to6 <- arrange(joined1to6, Length_gene, desc(Val)) #arranged in order of ascending a.a. length and descending frequency
  
r=1
for (r in 1:(nrow(joined1to6)-1)) {
  if(joined1to6$Length_gene[r] == joined1to6$Length_gene[r+1]) {
      joined1to6$Val[r+1] <- joined1to6$Val[r] + joined1to6$Val[r+1]  #sums together the frequencies of each unique combination of Gene Usage and CDR3 a.a. length and stores in the last row of each unique combination
    }
}
  joined1to6_cum_desc = arrange(joined1to6, Length_gene, desc(Val))  #arranged table of building cumulative values in order of ascending a.a. length and descending frequency
  joined1to6_cum_desc
  
  joined1to6_unique = joined1to6_cum_desc[!duplicated(joined1to6_cum_desc$Length_gene),]  #de-duplicated table so only the sum is displayed for each unique combination of Gene Usage and CDR3 a.a. length
  
  joined1to6_unique_byVal = arrange(joined1to6_unique, desc(Val))  #arranges in order of descending frequency

  naiveSpectra <- vis(joined1to6_unique_byVal) + xlim(10,20)
  naiveSpectra #makes and prints an averaged spectratype graph for naive T cell samples
  


#Spectratyping average for NP366-specific T CELLS, samples 7 to 12 (files beginning with "B")
  n=1
  for (n in 7:12) {
    n_countTable = assign(paste0("countTable",n), spectratype(immdata$data[[n]], .quant = "count", .col = "aa+v"))
    n_freqTable <- assign(paste0("freqTable",n), mutate(n_countTable, Val=(Val/sum(n_countTable$Val))/6)) #converts sequence counts into frequencies (out of a total of 1/6 since I will add these together to make a total of 1)
    n_freqTable <- assign(paste0("freqTable",n), mutate(n_freqTable, Length_gene = paste0(Length,"_",Gene)))
  }
  
  joined7to12 <- bind_rows(    #combining frequency tables for samples 7-12; frequencies should now add up to 1
    freqTable7,
    freqTable8,
    freqTable9,
    freqTable10,
    freqTable11,
    freqTable12
  ) 
  joined7to12 <- arrange(joined7to12, Length_gene, desc(Val)) #arranged in order of ascending a.a. length and descending frequency
  
  r=1
  for (r in 1:(nrow(joined7to12)-1)) {
    if(joined7to12$Length_gene[r] == joined7to12$Length_gene[r+1]) {
      joined7to12$Val[r+1] <- joined7to12$Val[r] + joined7to12$Val[r+1]  #sums together the frequencies of each unique combination of Gene Usage and CDR3 a.a. length and stores in the last row of each unique combination
    }
  }
  joined7to12_cum_desc = arrange(joined7to12, Length_gene, desc(Val))  #arranged table of building cumulative values in order of ascending a.a. length and descending frequency
  joined7to12_cum_desc
  
  joined7to12_unique = joined7to12_cum_desc[!duplicated(joined7to12_cum_desc$Length_gene),]  #de-duplicated table so only the sum is displayed for each unique combination of Gene Usage and CDR3 a.a. length
  
  joined7to12_unique_byVal = arrange(joined7to12_unique, desc(Val))  #arranges in order of descending frequency
  
  NP366Spectra <- vis(joined7to12_unique_byVal) + xlim(10,20)
  NP366Spectra #makes and prints an averaged spectratype graph for NP366 T cell samples
  
  
#Spectratyping average for M-SL9-specific T CELLS, samples 13-18 (files beginning with "C")
  n=1
  for (n in 13:18) {
    n_countTable = assign(paste0("countTable",n), spectratype(immdata$data[[n]], .quant = "count", .col = "aa+v"))
    n_freqTable <- assign(paste0("freqTable",n), mutate(n_countTable, Val=(Val/sum(n_countTable$Val))/6)) #converts sequence counts into frequencies (out of a total of 1/6 since I will add these together to make a total of 1)
    n_freqTable <- assign(paste0("freqTable",n), mutate(n_freqTable, Length_gene = paste0(Length,"_",Gene)))
  }
  
  joined13to18 <- bind_rows(    #combining frequency tables for samples 13-18; frequencies should now add up to 1
    freqTable13,
    freqTable14,
    freqTable15,
    freqTable16,
    freqTable17,
    freqTable18
  ) 
  joined13to18 <- arrange(joined13to18, Length_gene, desc(Val)) #arranged in order of ascending a.a. length and descending frequency
  
  r=1
  for (r in 1:(nrow(joined13to18)-1)) {
    if(joined13to18$Length_gene[r] == joined13to18$Length_gene[r+1]) {
      joined13to18$Val[r+1] <- joined13to18$Val[r] + joined13to18$Val[r+1]  #sums together the frequencies of each unique combination of Gene Usage and CDR3 a.a. length and stores in the last row of each unique combination
    }
  }
  joined13to18_cum_desc = arrange(joined13to18, Length_gene, desc(Val))  #arranged table of building cumulative values in order of ascending a.a. length and descending frequency
  joined13to18_cum_desc
  
  joined13to18_unique = joined13to18_cum_desc[!duplicated(joined13to18_cum_desc$Length_gene),]  #de-duplicated table so only the sum is displayed for each unique combination of Gene Usage and CDR3 a.a. length
  
  joined13to18_unique_byVal = arrange(joined13to18_unique, desc(Val))  #arranges in order of descending frequency
  
  MSL9Spectra <- vis(joined13to18_unique_byVal) + xlim(10,20)
  MSL9Spectra #makes and prints an averaged spectratype graph for M-SL9 T cell samples
  