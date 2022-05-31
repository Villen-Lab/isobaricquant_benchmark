#define input and output file name
file_path_in <- "data/a18182_MS2_1.mokapot.psms.txt"
file_path_out <- "data/a18182_MS2_1_hits_oldComet.csv"



#load function to recreate old Comet peptide sequence output
peptide_gsub_Comet_old <- function(string){
  if(!is.character(string)) stop("'input' needs to be a character vector")
  
  string <- gsub("[15.9949]", "*", string, fixed=T)
  string <- gsub("n[42.0106]", "n#", string, fixed=T)
  
  return(string)
}



#prepare data, by filtering Comet hits with Percolator output
data <- read.table(file = file_path_in, sep = '\t', header = TRUE)

#filter comet output based on q-value
data <- data[data$mokapot.q.value<=0.01,]

#sort by ScanNr
data <- data[order(data$ScanNr),]

#create row number column
data$PeptideId <- 1:nrow(data)

#extract charge state
data$Charge <- as.integer(gsub("^.*_.*_.*_.*_(.*)_.*$", "\\1", data$SpecId))

#export IsobaricQuant hits input file
data_exp <- data.frame(SearchID = 1,
                       PeptideID = data$PeptideId,
                       ModifiedPeptideSeq = peptide_gsub_Comet_old(data$Peptide),
                       Reference = gsub("^([^\\\t]*).*$", "\\1", data$Proteins, perl=T),
                       Charge = data$Charge,
                       ScanNumber = data$ScanNr,
                       mz = (data$ExpMass+data$Charge*1.0073-1.0073)/data$Charge)

write.table(data_exp, file = file_path_out, col.names = F, row.names = F, sep = ",")
