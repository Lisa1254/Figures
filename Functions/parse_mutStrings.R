#' Functions for parsing mutations strings in different styles
#' In the scales plot script, input data features a string as a code indicating
#' an initial amino acid, (optional chain character), position, and final aa.
#' Different software for in silico analysis can provide these strings in
#' different formats. 
#' 
#' Here are two self-drafted functions that can help parse the strings
#' to their component parts.
#' For both functions, the output will be a modification of an initial dataframe
#' to add new columns representing each part, where
#' ir = initial residue, pos = position, and fr = final residue
#' 
#' Also included is a function takes as input a dataframe with the 
#' initial amino acid and final after mutation, and outputs 
#' a vector with the row index of all scores that represent mutations 
#' to the same amino acid as initially expressed. 
#' The function recognizes "=" as final residue indicating mutation to self. 
#' AA can be expressed as 1 or 3 letters.

## Using digits as the position marker to split ----

#This was drafted using the expected format according to HGVS nomenclature
# http://varnomen.hgvs.org/recommendations/general/
#3-letter codes are indicated as preferred, but often 1-letter codes get used instead.
#As protein, codes will have "p." as prefix.
#This can be easily modified to parse DNA or other types of sequence variants
#This assumes that the only digits present are the position, and are all sequential without any other characters in between
#inputs are df = dataframe with mutation information to be modified and col = column index of mutation

parse_p_hgvs <- function(df, col) {
  ir <- unlist(lapply(df[,col], function(x){gsub("p.", "", strsplit(x, "[[:digit:]]+")[[1]][1])}))
  fr <- unlist(lapply(df[,col], function(x){strsplit(x, "[[:digit:]]+")[[1]][2]}))
  pos <- unlist(lapply(df[,col], function(x){stringr::str_extract(x, "[[:digit:]]+")}))
  df <- data.frame(df, ir, pos=as.numeric(pos), fr)
  return(df)
}

#

## Using number of characters in each position ----

#When the ddG scores were calculated using the Position Scan command of foldX software, the mutation code was expressed in the output file by three letter amino acid code for initial residue, chain character, position, and one letter amino acid code for the mutation residue
#The software would sometimes use different protonation states of histidine, describing it as H1S or H2S instead of HIS. The digits prevent the use of the above function, so instead this one makes use of the consistent number of characters for each part
#Variabless are:
#df = dataframe with mutation code,
#col = column number that has mutation information
#nchar1 = number of characters for each initial residue (usually will be 1 or 3)
#nchar2 = number of characters used to express final residue (usually will be 1 or 3)
#chain = FALSE is the default, but can change to true if the mutation string has the chain included (at the position directly after the initial amino acid)
#e.g. SERB123P has nchar1=3, nchar2=1, and chain=TRUE
split_mutation_string <- function(df, col, nchar1, nchar2, chain=FALSE) {
  ir <- substr(df[,col], 1, nchar1)
  fr <- substr(df[,col], nchar(df[,col])- nchar2 + 1, nchar(df[,col]))
  if (chain == TRUE) {
    pos <- as.numeric(substr(df[,col], nchar1+2, nchar(df[,col])-nchar2))
  } else {
    pos <- as.numeric(substr(df[,col], nchar1+1, nchar(df[,col])-nchar2))
  }
  df <- cbind(df, ir, pos, fr)
  return(df)
}

#

## Mutation to Self ----

#Function to identify rows with mutation to self in dataset to delete or filter if desired. Will recognize "=" character as self mutation, or when ir and fr are the same. Input for either can be either 3 letter or 1 letter amino acid codes

mut2self.pro <- function(df, ir_col, fr_col) {
  AA1 = c('G', 'A', 'L', 'V', 'I', 'P', 'R', 'T', 'S', 'C', 'M', 'K', 'E', 'Q', 'D', 'N', 'W', 'Y', 'F', 'H')
  AA3 = c('GLY', 'ALA', 'LEU', 'VAL', 'ILE', 'PRO', 'ARG', 'THR', 'SER', 'CYS', 'MET', 'LYS', 'GLU', 'GLN', 'ASP', 'ASN', 'TRP', 'TYR', 'PHE', 'HIS')
  rows.self <- which(df[,fr_col] == "=")
  df[,ir_col] <- toupper(df[,ir_col])
  df[,fr_col] <- toupper(df[,fr_col])
  for (i in 1:nrow(df)) {
    if (df[i,ir_col] == df[i,fr_col]) {
      rows.self <- c(rows.self, i)
    } else if ((nchar(df[i,ir_col]) == 1) & (nchar(df[i,fr_col]) == 3)) {
      if (plyr::mapvalues(df[i,ir_col], AA1, AA3, warn_missing = F) == df[i,fr_col]) {
        rows.self <- c(rows.self, i)
      }
    } else if ((nchar(df[i,ir_col]) == 3) & (nchar(df[i,fr_col]) == 1)) {
      if (plyr::mapvalues(df[i,ir_col], AA3, AA1, warn_missing = F) == df[i,fr_col]) {
        rows.self <- c(rows.self, i)
      }
    }
  }
  return(rows.self)
}



