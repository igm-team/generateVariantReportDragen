Filter.by.Flag <- function(data, type, strict = FALSE) {
  #check for column name
  stopifnot(length(which(colnames(data) == "Flag"))>0)
  if (strict) {
    result <- data[data["Flag"] == type,]
  } else {
    result <- data[grep(type, data$"Flag") ,] 
  }
  # put NA rows back if any
  result = rbind(result,data[data["Flag"] == "NA",])
}

normalized.name <- function(names) {
  names.new <- gsub("[( )#&-]",".",names)
  names.new
}

Filter.for.compound.heterozygote <- function(data) {
  #check for correct columns, just in case the format has been change
  columns <- c("Var Ctrl Freq #1 & #2 (co-occurance)","Genotype (mother) (#1)", 
               "Genotype (father) (#1)","Genotype (mother) (#2)","Genotype (father) (#2)",
                "Function (#1)","Function (#2)",
               "Ctrl MAF (#1)", "Evs Eur Maf (#1)", "Evs Afr Maf (#1)",
               "ExAC global maf (#1)", "ExAC afr maf (#1)", "ExAC amr maf (#1)",
               "ExAC eas maf (#1)","ExAC sas maf (#1)","ExAC fin maf (#1)",
               "ExAC nfe maf (#1)","ExAC oth maf (#1)",
               "Ctrl MAF (#2)", "Evs Eur Maf (#2)", "Evs Afr Maf (#2)",
               "ExAC global maf (#2)", "ExAC afr maf (#2)", "ExAC amr maf (#2)",
               "ExAC eas maf (#2)","ExAC sas maf (#2)","ExAC fin maf (#2)",
               "ExAC nfe maf (#2)","ExAC oth maf (#2)",
               "Minor Hom Ctrl (#1)","Evs All Genotype Count (#1)","ExAC global gts (#1)",
               "Evs Filter Status (#1)","ExAC global gts (#1)",
               "Minor Hom Ctrl (#2)","Evs All Genotype Count (#2)","ExAC global gts (#2)",
               "Evs Filter Status (#2)","ExAC global gts (#2)" 
              )
  
  #make sure all columns are present
  stopifnot(length(setdiff(normalized.name(columns),colnames(data))) ==0)
  
  if (dim(data)[1] ==0) { return(data)}
  
  #step 1: 
  data   <- data[is.na(data[normalized.name("Var Ctrl Freq #1 & #2 (co-occurance)")])
                 | data[normalized.name("Var Ctrl Freq #1 & #2 (co-occurance)")] == 0,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 2: 
  data <- data[(is.na(data[normalized.name("Genotype (mother) (#1)")])
                | data[normalized.name("Genotype (mother) (#1)")] != "hom")          &
                (is.na(data[normalized.name("Genotype (father) (#1)")])
                | data[normalized.name("Genotype (father) (#1)")] != "hom")          &
                (is.na(data[normalized.name("Genotype (mother) (#2)")])
                | data[normalized.name("Genotype (mother) (#2)")] != "hom")          &
                (is.na(data[normalized.name("Genotype (father) (#2)")])
                | data[normalized.name("Genotype (father) (#2)")] != "hom") 
               ,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 3:
  Index <- grep("^SYNONYMOUS",data[normalized.name("Function (#1)")][,1])
  if (length(Index) >0) {
    data <- data[-Index,]
  }
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- grep("^INTRON_EXON",data[normalized.name("Function (#1)")][,1])
  if (length(Index) >0) {
    data <- data[-Index,]
  }
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- grep("^SYNONYMOUS",data[normalized.name("Function (#2)")][,1])
  if (length(Index) >0) {
    data <- data[-Index,]
  }
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- grep("^INTRON_EXON",data[normalized.name("Function (#2)")][,1])
  if (length(Index) >0) {
    data <- data[-Index,]
  }
  if (dim(data)[1] ==0) { return(data)}
  
  
  #step 4:
  data <- data[(is.na(data[normalized.name("Ctrl MAF (#1)")])
                | data[normalized.name("Ctrl MAF (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("Evs Eur Maf (#1)")])
                  | data[normalized.name("Evs Eur Maf (#1)")] < 0.005)        & 
                 (is.na(data[normalized.name("Evs Afr Maf (#1)")])
                  | data[normalized.name("Evs Afr Maf (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC global maf (#1)")])
                  | data[normalized.name("ExAC global maf (#1)")] < 0.005)    &
                 (is.na(data[normalized.name("ExAC afr maf (#1)")])
                  | data[normalized.name("ExAC afr maf (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC amr maf (#1)")])
                  | data[normalized.name("ExAC amr maf (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC eas maf (#1)")])
                  | data[normalized.name("ExAC eas maf (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC sas maf (#1)")])
                  | data[normalized.name("ExAC sas maf (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC fin maf (#1)")])
                  | data[normalized.name("ExAC fin maf (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC nfe maf (#1)")])
                  | data[normalized.name("ExAC nfe maf (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC oth maf (#1)")])
                  | data[normalized.name("ExAC oth maf (#1)")] < 0.005)        &
                 (is.na(data[normalized.name("Ctrl MAF (#2)")])
                  | data[normalized.name("Ctrl MAF (#2)")] < 0.005)        &
                 (is.na(data[normalized.name("Evs Eur Maf (#2)")])
                  | data[normalized.name("Evs Eur Maf (#2)")] < 0.005)        & 
                 (is.na(data[normalized.name("Evs Afr Maf (#2)")])
                  | data[normalized.name("Evs Afr Maf (#2)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC global maf (#2)")])
                  | data[normalized.name("ExAC global maf (#2)")] < 0.005)    &
                 (is.na(data[normalized.name("ExAC afr maf (#2)")])
                  | data[normalized.name("ExAC afr maf (#2)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC amr maf (#2)")])
                  | data[normalized.name("ExAC amr maf (#2)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC eas maf (#2)")])
                  | data[normalized.name("ExAC eas maf (#2)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC sas maf (#2)")])
                  | data[normalized.name("ExAC sas maf (#2)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC fin maf (#2)")])
                  | data[normalized.name("ExAC fin maf (#2)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC nfe maf (#2)")])
                  | data[normalized.name("ExAC nfe maf (#2)")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC oth maf (#2)")])
                  | data[normalized.name("ExAC oth maf (#2)")] < 0.005)
               ,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 5:
  data <- data[(is.na(data[normalized.name("Minor Hom Ctrl (#1)")])
                | data[normalized.name("Minor Hom Ctrl (#1)")] == 0),]
  if (dim(data)[1] ==0) {return(data)}
  
  Index <- is.na(data[normalized.name("Evs All Genotype Count (#1)")]) 
  for (i in 1:length(Index)) {
    if (!Index[i]) { #not is NA
      str <-as.character(data[i,normalized.name("Evs All Genotype Count (#1)")])
      not.done <- (length(grep("A1A1=[[123456789]]",str))==0) & (length(grep("A2A2=[[123456789]]",str))==0) & 
        (length(grep("A3A3=[[123456789]]",str))==0) 
      
      if (not.done) { #not found 
        genotypes <- strsplit(str,"/")[[1]]
        genotype.count <- strsplit(genotypes,"=")[[1]][2]
        not.done <- genotype.count == "0"
      }
      if (not.done) { #the last condition has been checked
        Index[i] <- TRUE
      }
    }
  }
  data <- data[Index,]
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- is.na(data[normalized.name("ExAC global gts (#1)")]) 
  for (i in 1:length(Index)) {
    if (!Index[i]) { #not is NA
      str <- as.character(data[i,normalized.name("ExAC global gts (#1)")])
      not.done <- !is.na(str)
      if (not.done) {
        not.done <- (length(grep("/0'$",str))>0)
      }
      
      if (not.done) {
        #print(str)
        genotype.count <- strsplit(str,"/")[[1]]
        if ((length(genotype.count) ==4) & (length(grep("NA",genotype.count[1])) >0)) {
          not.done <- (genotype.count[3] == "0") & (length(grep("^0",genotype.count[4])) >0)
        }
      }
      if (not.done) { #the last condition has been checnke
        Index[i] <- TRUE
      }
    }
  }
  data <- data[Index,]
  if (dim(data)[1] ==0) { return(data)}
  
  data <- data[((is.na(data[normalized.name("Evs Filter Status (#1)")])
                  | data[normalized.name("Evs Filter Status (#1)")] != "FAIL"))
               ,]
  if (dim(data)[1] ==0) { return(data)}
  
  data <- data[(is.na(data[normalized.name("Minor Hom Ctrl (#2)")])
                | data[normalized.name("Minor Hom Ctrl (#2)")] == 0),]
  if (dim(data)[1] ==0) {return(data)}
  
  Index <- is.na(data[normalized.name("Evs All Genotype Count (#2)")]) 
  for (i in 1:length(Index)) {
    if (!Index[i]) { #not is NA
      str <-as.character(data[i,normalized.name("Evs All Genotype Count (#2)")])
      not.done <- (length(grep("A1A1=[[123456789]]",str))==0) & (length(grep("A2A2=[[123456789]]",str))==0) & 
        (length(grep("A3A3=[[123456789]]",str))==0) 
      
      if (not.done) { #not found 
        genotypes <- strsplit(str,"/")[[1]]
        genotype.count <- strsplit(genotypes,"=")[[1]][2]
        not.done <- genotype.count == "0"
      }
      if (not.done) { #the last condition has been checked
        Index[i] <- TRUE
      }
    }
  }
  data <- data[Index,]
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- is.na(data[normalized.name("ExAC global gts (#2)")]) 
  for (i in 1:length(Index)) {
    if (!Index[i]) { #not is NA
      str <- as.character(data[i,normalized.name("ExAC global gts (#2)")])
      not.done <- !is.na(str)
      if (not.done) {
        not.done <- (length(grep("/0'$",str))>0)
      }
      
      if (not.done) {
        #print(str)
        genotype.count <- strsplit(str,"/")[[1]]
        if ((length(genotype.count) ==4) & (length(grep("NA",genotype.count[1])) >0)) {
          not.done <- (genotype.count[3] == "0") & (length(grep("^0",genotype.count[4])) >0)
        }
      }
      #browser("chet")
      if (not.done) { #the last condition has been checnke
           Index[i] <- TRUE
      }
    }
  }
  data <- data[Index,]
  if (dim(data)[1] ==0) { return(data)}
  
  data <- data[((is.na(data[normalized.name("Evs Filter Status (#2)")])
                 | data[normalized.name("Evs Filter Status (#2)")] != "FAIL"))
                ,]
  data
}

Filter.for.homozygous <- function(data) {
  #check for correct columns
  columns <- c("Percent Read Alt (child)","Genotype Qual GQ (child)", 
               "Rms Map Qual MQ (child)","Evs Filter Status",
               "Samtools Raw Coverage (mother)", "Genotype (mother)", "Genotype (father)",
               "Function", "Ctrl MAF", "Evs Eur Maf", "Evs Afr Maf","ExAC global maf",
               "ExAC afr maf","ExAC amr maf","ExAC eas maf","ExAC sas maf",
               "ExAC fin maf","ExAC nfe maf","ExAC oth maf","Minor Hom Ctrl",
               "Evs All Genotype Count", "ExAC global gts")
  
  #make sure all columns are present
  stopifnot(length(setdiff(normalized.name(columns),colnames(data))) ==0)
  
  if (dim(data)[1] ==0) { return(data)}
  #step 1: 
  data   <- data[is.na(data[normalized.name("Percent Read Alt (child)")])
                 | data[normalized.name("Percent Read Alt (child)")] > 0.799,]
  if (dim(data)[1] ==0) { return(data)}
  #step 2: 
  data <- data[(is.na(data[normalized.name("Genotype Qual GQ (child)")])
                | data[normalized.name("Genotype Qual GQ (child)")] > 20)        &
                 (is.na(data[normalized.name("Rms Map Qual MQ (child)")])
                  | data[normalized.name("Rms Map Qual MQ (child)")] > 40)       &
                 (is.na(data[normalized.name("Evs Filter Status")])
                  | data[normalized.name("Evs Filter Status")] != "FAIL") 
               ,]
  if (dim(data)[1] ==0) { return(data)}
  #step 3:
  data <- data[(is.na(data[normalized.name("Genotype (mother)")])
                  | data[normalized.name("Genotype (mother)")] == "het")          &
                 (is.na(data[normalized.name("Genotype (father)")])
                  | data[normalized.name("Genotype (father)")] == "het")  
               ,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 4:
  Index <- grep("^SYNONYMOUS",data$"Function")
  if (length(Index) >0) {
    data <- data[-Index,]
  }
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- grep("^INTRON_EXON",data$"Function")
  if (length(Index) >0) {
    data <- data[-Index,]
  }
  if (dim(data)[1] ==0) { return(data)}
  
  #step 5:
  data <- data[(is.na(data[normalized.name("Ctrl MAF")])
                | data[normalized.name("Ctrl MAF")] < 0.005)        &
                 (is.na(data[normalized.name("Evs Eur Maf")])
                  | data[normalized.name("Evs Eur Maf")] < 0.005)        & 
                 (is.na(data[normalized.name("Evs Afr Maf")])
                  | data[normalized.name("Evs Afr Maf")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC global maf")])
                  | data[normalized.name("ExAC global maf")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC afr maf")])
                  | data[normalized.name("ExAC afr maf")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC amr maf")])
                  | data[normalized.name("ExAC amr maf")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC eas maf")])
                  | data[normalized.name("ExAC eas maf")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC sas maf")])
                  | data[normalized.name("ExAC sas maf")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC fin maf")])
                  | data[normalized.name("ExAC fin maf")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC nfe maf")])
                  | data[normalized.name("ExAC nfe maf")] < 0.005)        &
                 (is.na(data[normalized.name("ExAC oth maf")])
                  | data[normalized.name("ExAC oth maf")] < 0.005)  
               ,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 6:
  data <- data[(is.na(data[normalized.name("Minor Hom Ctrl")])
                | data[normalized.name("Minor Hom Ctrl")] == 0),]
  if (dim(data)[1] ==0) {return(data)}
  
  Index <- is.na(data[normalized.name("Evs All Genotype Count")]) 
  for (i in 1:length(Index)) {
    if (!Index[i]) { #not is NA
      str <-as.character(data[i,normalized.name("Evs All Genotype Count")])
      not.done <- (length(grep("A1A1=[[123456789]]",str))==0) & (length(grep("A2A2=[[123456789]]",str))==0) & 
        (length(grep("A3A3=[[123456789]]",str))==0)
      
      if (not.done) { #not found 
        genotypes <- strsplit(str,"/")[[1]]
        genotype.count <- strsplit(genotypes,"=")[[1]][2]
        not.done <- genotype.count == "0"
      }
      if (not.done) { #the last condition has been checked
        Index[i] <- TRUE
      }
    }
  }
  data <- data[Index,]
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- is.na(data[normalized.name("ExAC global gts")]) 
  for (i in 1:length(Index)) {
    if (!Index[i]) { #not is NA
      str <- as.character(data[i,normalized.name("ExAC global gts")])
      not.done <- !is.na(str)
      if (not.done) {
        not.done <- (length(grep("/0'$",str))>0)
      }
      
      if (not.done) {
        #print(str)
        genotype.count <- strsplit(str,"/")[[1]]
        if ((length(genotype.count) ==4) & (length(grep("NA",genotype.count[1])) >0)) {
          not.done <- (genotype.count[3] == "0") & (length(grep("^0",genotype.count[4])) >0)
        }
      }
    #browser("hom")
        if (not.done) { #the last condition has been checnke
            Index[i] <- TRUE
        }
    }
  }
  data <- data[Index,]
  data
}

Filter.for.hemizygous <- function(data) {
  #check for correct columns
  columns <- c("Percent Read Alt (child)","Genotype Qual GQ (child)", 
               "Rms Map Qual MQ (child)","Evs Filter Status",
               "Samtools Raw Coverage (mother)", "Genotype (mother)", "Genotype (father)",
               "Function", "Ctrl MAF", "Evs Eur Maf", "Evs Afr Maf","ExAC global maf",
               "ExAC afr maf","ExAC amr maf","ExAC eas maf","ExAC sas maf",
               "ExAC fin maf","ExAC nfe maf","ExAC oth maf","Minor Hom Ctrl",
               "Evs All Genotype Count", "ExAC global gts")
  
  #make sure all columns are present
  stopifnot(length(setdiff(normalized.name(columns),colnames(data))) ==0)
  
  if (dim(data)[1] ==0) { return(data)}
  #step 1: 
  data   <- data[is.na(data[normalized.name("Percent Read Alt (child)")])
                 | data[normalized.name("Percent Read Alt (child)")] > 0.7999,]
  if (dim(data)[1] ==0) { return(data)}
  #step 2: 
  data <- data[(is.na(data[normalized.name("Genotype Qual GQ (child)")])
                | data[normalized.name("Genotype Qual GQ (child)")] > 20)        &
                 (is.na(data[normalized.name("Rms Map Qual MQ (child)")])
                  | data[normalized.name("Rms Map Qual MQ (child)")] > 40)       &
                 (is.na(data[normalized.name("Evs Filter Status")])
                  | data[normalized.name("Evs Filter Status")] != "FAIL") 
               ,]
  if (dim(data)[1] ==0) { return(data)}
  #step 3:
  data <- data[(is.na(data[normalized.name("Samtools Raw Coverage (mother)")])
                  | data[normalized.name("Samtools Raw Coverage (mother)")] > 9) &
                 (is.na(data[normalized.name("Genotype (mother)")])
                  | data[normalized.name("Genotype (mother)")] == "het")          &
                 (is.na(data[normalized.name("Genotype (father)")])
                  | data[normalized.name("Genotype (father)")] == "hom ref")  
               ,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 4:
  Index <- grep("^SYNONYMOUS",data$"Function")
  if (length(Index) >0) {
    data <- data[-Index,]
  }
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- grep("^INTRON_EXON",data$"Function")
  if (length(Index) >0) {
    data <- data[-Index,]
  }
  if (dim(data)[1] ==0) { return(data)}
  
  #step 5:
  data <- data[(is.na(data[normalized.name("Ctrl MAF")])
                | data[normalized.name("Ctrl MAF")] < 0.005)        &
                 (is.na(data[normalized.name("Evs Eur Maf")])
                | data[normalized.name("Evs Eur Maf")] < 0.005)        & 
                 (is.na(data[normalized.name("Evs Afr Maf")])
                | data[normalized.name("Evs Afr Maf")] < 0.005)        &
                (is.na(data[normalized.name("ExAC global maf")])
                | data[normalized.name("ExAC global maf")] < 0.005)        &
                (is.na(data[normalized.name("ExAC afr maf")])
                | data[normalized.name("ExAC afr maf")] < 0.005)        &
                (is.na(data[normalized.name("ExAC amr maf")])
                | data[normalized.name("ExAC amr maf")] < 0.005)        &
                (is.na(data[normalized.name("ExAC eas maf")])
                | data[normalized.name("ExAC eas maf")] < 0.005)        &
                (is.na(data[normalized.name("ExAC sas maf")])
                | data[normalized.name("ExAC sas maf")] < 0.005)        &
                (is.na(data[normalized.name("ExAC fin maf")])
                | data[normalized.name("ExAC fin maf")] < 0.005)        &
                (is.na(data[normalized.name("ExAC nfe maf")])
                | data[normalized.name("ExAC nfe maf")] < 0.005)        &
                (is.na(data[normalized.name("ExAC oth maf")])
                | data[normalized.name("ExAC oth maf")] < 0.005)  
               ,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 6:
  data <- data[(is.na(data[normalized.name("Minor Hom Ctrl")])
                | data[normalized.name("Minor Hom Ctrl")] == 0),]
  if (dim(data)[1] ==0) {return(data)}
  
  Index <- is.na(data[normalized.name("Evs All Genotype Count")]) 
  for (i in 1:length(Index)) {
    if (!Index[i]) { #not is NA
      str <-as.character(data[i,normalized.name("Evs All Genotype Count")])
      not.done <- (length(grep("A1A1=0",str))==0) & (length(grep("A2A2=0",str))==0) & 
        (length(grep("A3A3=0",str))==0) & (length(grep("/A=\\d/",str))==0) &
        (length(grep("/C=\\d/",str))==0) & (length(grep("/G=\\d/",str))==0)    & 
        (length(grep("/T=\\d/",str))==0)  
      
      if (not.done) { #not found 
        genotypes <- strsplit(str,"/")[[1]]
        genotype.count <- strsplit(genotypes,"=")[[1]][2]
        not.done <- genotype.count == "0"
      }
      if (not.done) { #the last condition has been checked
        Index[i] <- TRUE
      }
    }
  }
  data <- data[Index,]
  if (dim(data)[1] ==0) { return(data)}
  
  Index <- is.na(data[normalized.name("ExAC global gts")]) 
  for (i in 1:length(Index)) {
    if (!Index[i]) { #not is NA
      str <- as.character(data[i,normalized.name("ExAC global gts")])
      not.done <- !is.na(str)
      if (not.done) {
        not.done <- (length(grep("/0'$",str))>0)
      }
          
      if (not.done) {
        #print(str)
        genotype.count <- strsplit(str,"/")[[1]]
        if ((length(genotype.count) ==4) & (length(grep("NA",genotype.count[1])) >0)) {
          not.done <- (genotype.count[3] == "0") & (length(grep("^0",genotype.count[4])) >0)
        }
      }
    #browser("hemi")
        if (not.done) { #the last condition has been checnke
            Index[i] <- TRUE
        }
    }
  }
  data <- data[Index,]
  data
}

Filter.for.denovo <- function(data) {
  #check for correct columns
  columns <- c("Percent Read Alt (child)","Genotype Qual GQ (child)", "Qual By Depth QD (child)",
               "Rms Map Qual MQ (child)","Genotype Qual GQ (child)","Qual (child)",
               "Samtools Raw Coverage (child)","Evs Filter Status", "Samtools Raw Coverage (mother)",
               "Samtools Raw Coverage (father)","QC Fail Ctrl","QC Fail Ctrl","Ctrl MAF",
               "Evs All Genotype Count","ExAC global gts")
  
  #make sure all columns are present
  stopifnot(length(setdiff(normalized.name(columns),colnames(data))) ==0)
  
  #steps can be merged or splitted
  #but merged steps can usually perform fasted
  
  #step 4:
  data <- data[(is.na(data[normalized.name("QC Fail Ctrl")])
                | data[normalized.name("QC Fail Ctrl")] == 0)                   &
                (is.na(data[normalized.name("Ctrl MAF")])
                  | data[normalized.name("Ctrl MAF")]  == 0)                    &
                is.na(data[normalized.name("Evs All Genotype Count")])          &
                is.na(data[normalized.name("ExAC global gts")])
               ,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 1: 
  data   <- data[is.na(data[normalized.name("Percent Read Alt (child)")])
                 | data[normalized.name("Percent Read Alt (child)")] > 0.1999,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 2: 
  data <- data[(is.na(data[normalized.name("Genotype Qual GQ (child)")])
                  | data[normalized.name("Genotype Qual GQ (child)")] > 20)      &
               (is.na(data[normalized.name("Qual By Depth QD (child)")]) 
                  | data[normalized.name("Qual By Depth QD (child)")] > 2)       &
               (is.na(data[normalized.name("Rms Map Qual MQ (child)")])
                  | data[normalized.name("Rms Map Qual MQ (child)")] > 40)       &
               (is.na(data[normalized.name("Genotype Qual GQ (child)")]) 
                  | data[normalized.name("Genotype Qual GQ (child)")] > 20)      &
               (is.na(data[normalized.name("Qual (child)")]) 
                  | data[normalized.name("Qual (child)")] > 50)                  &
               (is.na(data[normalized.name("Samtools Raw Coverage (child)")])
                  | data[normalized.name("Samtools Raw Coverage (child)")] > 9) &
               (is.na(data[normalized.name("Evs Filter Status")])
                  | data[normalized.name("Evs Filter Status")] != "FAIL") 
              ,]
  if (dim(data)[1] ==0) { return(data)}
  
  #step 3: //de novo flag is not filtered here. Which can be filtered separately
  data <- data[(!is.na(data[normalized.name("Samtools Raw Coverage (mother)")])
                & data[normalized.name("Samtools Raw Coverage (mother)")] > 9)      &
                (!is.na(data[normalized.name("Samtools Raw Coverage (father)")])
                  & data[normalized.name("Samtools Raw Coverage (father)")] > 9)
               ,]
  
  data
}

parse.hemihomo.count.exac <- function(s) {
  if (is.na(s)) {
    return(0)
  } else {
    has.4.genotypes <- (length(grep("na",s))>0)
    fields <- strsplit(gsub("na/","",gsub("'","",s)),"/")[[1]]
    fields <- as.numeric(fields)  
    if (has.4.genotypes) { #na/heterozygous/hemizygous/homozygous
      return (fields[2] + fields[3]) 
    } else { #hom-ref/het/hom
      return(fields[3])
    }
  }
}

parse.allele.count.exac <- function(s) {
  if (is.na(s)) {
    return(0)
  } else {
    has.4.genotypes <- (length(grep("na",s))>0)
    fields <- strsplit(gsub("na/","",gsub("'","",s)),"/")[[1]]
    fields <- as.numeric(fields)  
    if (has.4.genotypes) { #na/heterozygous/hemizygous/homozygous
       return (fields[1] + fields[2] + 2* fields[3]) 
    } else { #hom-ref/het/hom
      return (min(2*fields[1]+fields[2], fields[2]+ 2* fields[3]))  
    }
  }
}

parse.hemihomo.count.evs <- function(s) {
  if (is.na(s)) {
    return(0)
  } else {
    fields <- strsplit(s,"/")[[1]]
    if (length(grep('r',s)) == 0) { #snv, change the string so it has the same format as indels
      alt.allele <- substr(s,1,1)
      s <- gsub(alt.allele,"_",s) # make it a non alpha first
      s <- gsub("[[:alpha:]]","r",s) # change the other one to ref
      s <- gsub("_","a",s) # change the alt to "a".
      fields <- strsplit(s,"/")[[1]] #split again since s has been normlized
    }
    hemi.hom.count <- 0
    for (i in 1: length(fields)) {
      current.record <- strsplit(fields[i],"=")[[1]]
      current.type = gsub("[[:digit:]]","",current.record[1]) #works for indels only, no effect on snv
      current.count <- as.numeric(current.record[2])
      if (length(grep("r",current.type)) == 0) { #counting only aa or a
        hemi.hom.count <- hemi.hom.count +  current.count
      } 
    }
    return (hemi.hom.count)
  }
}

parse.allele.count.evs <- function(s) {
  if (is.na(s)) {
    return(0)
  } else {
    fields <- strsplit(s,"/")[[1]]
    if (length(grep('r',s)) == 0) { #snv, change the string so it has the same format as indels
      alt.allele <- substr(s,1,1)
      s <- gsub(alt.allele,"_",s) # make it a non alpha first
      s <- gsub("[[:alpha:]]","r",s) # change the other one to ref
      s <- gsub("_","a",s) # change the alt to "a".
      fields <- strsplit(s,"/")[[1]] #split again since s has been normlized
    }
      
    hom.allele.count <- 0
    het.allele.count <- 0
    hom.ref.allele.count <- 0
    for (i in 1: length(fields)) {
        current.record <- strsplit(fields[i],"=")[[1]]
        current.type = gsub("[[:digit:]]","",current.record[1]) #works for indels only, no effect on snv
        current.count <- as.numeric(current.record[2])
        if (length(grep("rr",current.type)) > 0) {
          hom.ref.allele.count <- hom.ref.allele.count + 2 * current.count
        } else if (length(grep("aa",current.type)) > 0) {
          hom.allele.count <- hom.allele.count + 2 * current.count
        } else if (length(grep("ar",current.type)) > 0){
          het.allele.count <- het.allele.count + current.count
        } else if (length(grep("r",current.type)) > 0) {
          hom.ref.allele.count <- hom.ref.allele.count + current.count 
        } else { # a
          hom.allele.count <- hom.allele.count + current.count
        }
    }
    return (min(het.allele.count + hom.allele.count, het.allele.count +  hom.ref.allele.count ))
  }
}

#temp <- data[data["variant.type"] == "indel" & !is.na(data["evs.all.genotype.count"]),c("ref.allele","alt.allele","evs.all.genotype.count")]

Filter.by.function <- function(x) {
  return ((x=="stop_gained") | (x =="stop_lost") | (x =="start_lost")
          | (x=="frame_shift") | (x== "splice_site_donor") | (x=="splice_site_acceptor") | (x=="codon_deletion") | (x=="codon_insertion") | (x=="codon_change_plus_codon_deletion") | (x=="codon_change_plus_codon_insertion") | (x=="exon_deleted"))
}

Filter.by.hemihomo.count <- function(data,threshold,is.comphet = false) {
  #get igm genotype count
  if (is.comphet) {
    stopifnot(length(which(colnames(data) == normalized.name("major hom ctrl (#1)")))>0)
    stopifnot(length(which(colnames(data) ==  normalized.name("minor hom ctrl (#1)")))>0)
    genotypecount.1 <- cbind(data[normalized.name("major hom ctrl (#1)")] , data[normalized.name("minor hom ctrl (#1)")] )
    genotypecount.1 <- apply(genotypecount.1,1,function(x) min(x))  
    
    stopifnot(length(which(colnames(data) == normalized.name("major hom ctrl (#2)")))>0)
    stopifnot(length(which(colnames(data) ==  normalized.name("minor hom ctrl (#2)")))>0)
    genotypecount.2 <- cbind(data[normalized.name("major hom ctrl (#2)")] , data[normalized.name("minor hom ctrl (#2)")] )
    genotypecount.2 <- apply(genotypecount.2,1,function(x) min(x))  
    
    #get evs genotype count
    column.name = normalized.name("evs.all.genotype.count (#1)")
    stopifnot(length(which(colnames(data) == column.name))>0)
    tmp <- sapply(data[column.name], as.character)
    evs.genotypecount.1<- sapply(tmp, parse.hemihomo.count.evs)
    
    column.name = normalized.name("evs.all.genotype.count (#2)")
    stopifnot(length(which(colnames(data) == column.name))>0)
    tmp <- sapply(data[column.name], as.character)
    evs.genotypecount.2<- sapply(tmp, parse.hemihomo.count.evs)
    
    #get exac genotype count
    column.name = normalized.name("exac.global.gts (#1)")
    tmp <- sapply(data[column.name], as.character)
    exac.genotypecount.1 <- sapply(tmp, parse.hemihomo.count.exac)
    
    column.name = normalized.name("exac.global.gts (#2)")
    tmp <- sapply(data[column.name], as.character)
    exac.genotypecount.2 <- sapply(tmp, parse.hemihomo.count.exac)
    
    all <- ((genotypecount.1 + evs.genotypecount.1 + exac.genotypecount.1) <= threshold) 
    all <- all & ((genotypecount.2 + evs.genotypecount.2 + exac.genotypecount.2) <= threshold)
    
  } else {
    stopifnot(length(which(colnames(data) == "major.hom.ctrl"))>0)
    stopifnot(length(which(colnames(data) == "minor.hom.ctrl"))>0)
    genotypecount <- cbind(data["major.hom.ctrl"] , data["minor.hom.ctrl"] )
    genotypecount <- apply(genotypecount,1,function(x) min(x))  
    
    #get evs genotype count
    stopifnot(length(which(colnames(data) == "evs.all.genotype.count"))>0)
    tmp <- sapply(data["evs.all.genotype.count"], as.character)
    evs.genotypecount<- sapply(tmp, parse.hemihomo.count.evs)
    
    #get exac genotype count
    stopifnot(length(which(colnames(data) == "exac.global.gts"))>0)
    tmp <- sapply(data["exac.global.gts"], as.character)
    exac.genotypecount <- sapply(tmp, parse.hemihomo.count.exac)
    
    all <- ((genotypecount + evs.genotypecount + exac.genotypecount) <= threshold)  
  }
  
  result <- data[all,]
}

Filter.by.Allele.Count <- function(data, threshold) {  
  #get igm allele count
    browser()
  stopifnot(length(which(colnames(data) == "major.hom.ctrl"))>0)
  stopifnot(length(which(colnames(data) == "het.ctrl"))>0)
  stopifnot(length(which(colnames(data) == "minor.hom.ctrl"))>0)
  allelecount <- cbind(2*data["major.hom.ctrl"] + data["het.ctrl"], data["het.ctrl"] + 2*data["minor.hom.ctrl"] )
  allelecount <- apply(allelecount,1,function(x) min(x))
  
  #get exac allele count
  stopifnot(length(which(colnames(data) == "evs.all.genotype.count"))>0)
  
  tmp <- sapply(data["evs.all.genotype.count"], as.character)
  evs.allelecount <- sapply(tmp, parse.allele.count.evs)
  
  tmp <- sapply(data["exac.global.gts"], as.character)
  exac.allelecount <- sapply(tmp, parse.allele.count.exac)
  
  qc.fail.ctrl <- sapply(data["qc.fail.ctrl"], as.numeric)
  all <- ((allelecount + evs.allelecount + exac.allelecount + qc.fail.ctrl) <= threshold)
  result <- data[all,]
}


Filter.for.tier2 <- function(data, is.comphet = false) {
  if (is.comphet) {
    #check for correct columns
    columns <- c("hgmd variant class (#1)","hgmd variant class (#2)","hgmd flanking count (#1)", "hgmd flanking count (#2)","clinvar flanking count (#1)","clinvar clinical significance (#1)","function (#1)","function (#2)","clingen haploinsufficiencydesc (#1)","clinvar pathogenic indel count (#1)", "clinvar pathogenic cnv count (#1)", "clinvar pathogenic snv splice count (#1)", "clinvar pathogenic snv nonsense count (#1)")
   
    #make sure all columns are present
    stopifnot(length(setdiff(normalized.name(columns),colnames(data))) ==0)
    
    #inclusion rule 1: 
    r1.1 <- sapply(data[normalized.name("hgmd variant class (#1)")], as.character)
    r1.1 <- sapply(r1.1, function(x) length(grep("dm",x)) > 0)
    r1.2 <- sapply(data[normalized.name("hgmd variant class (#2)")], as.character)
    r1.2 <- sapply(r1.2, function(x) length(grep("dm",x)) > 0) 
    r1 <- r1.1 | r1.2

    #inclusion rule 2:
    suppresswarnings(temp1 <- sapply(data[normalized.name("hgmd flanking count (#1)")], as.numeric))
    temp1[is.na(temp1)] <- 0
    suppresswarnings(temp2 <- sapply(data[normalized.name("clinvar flanking count (#1)")], as.numeric))
    temp2[is.na(temp2)] <- 0
    suppresswarnings(temp3 <- sapply(data[normalized.name("hgmd flanking count (#2)")], as.numeric))
    temp3[is.na(temp3)] <- 0
    suppresswarnings(temp4 <- sapply(data[normalized.name("clinvar flanking count (#2)")], as.numeric))
    temp4[is.na(temp4)] <- 0
    r2 <- temp1 > 0 | temp2 > 0 | temp3 > 0 | temp4 > 0

    #inclusion rule 3:
    temp <- sapply(data[normalized.name("clinvar clinical significance (#1)")], as.character)
    r3 = (temp == "pathogenic") | (temp == "likely pathogenic") | (temp == "likely pathogenic;pathogenic")
    
    functional.1 <- Filter.by.function(sapply(data[normalized.name("function (#1)")], as.character))
    functional.2 <- Filter.by.function(sapply(data[normalized.name("function (#2)")], as.character))

    #inclusion rule 4:
    r4 <- (functional.1 | functional.2) & (!eval(parse(text=paste0("data$",normalized.name("clingen haploinsufficiencydesc (#1)")))) %in% c("na","no evidence","little evidence"))
    
    #inclusion rule 5:
    suppresswarnings(temp <- sapply(data[normalized.name("clinvar pathogenic indel count (#1)")], as.numeric))
    temp[is.na(temp)] <- 0
    r5 <- (functional.1 | functional.2) & (temp > 0)
    
    #inclusion rule 6:
    suppresswarnings(temp <- sapply(data[normalized.name("clinvar pathogenic snv splice count (#1)")], as.numeric))
    temp[is.na(temp)] <- 0
    r6 <- (functional.1 | functional.2) & (temp > 0)
    
    #inclusion rule 7:
    suppresswarnings(temp <- sapply(data[normalized.name("clinvar pathogenic snv nonsense count (#1)")], as.numeric))
    temp[is.na(temp)] <- 0
    r7 <- (functional.1 | functional.2) & (temp > 0)

    #inclusion rule 8:
    suppresswarnings(temp <- sapply(data[normalized.name("clinvar pathogenic cnv count (#1)")], as.numeric))
    temp[is.na(temp)] <- 0
    r8 <- (functional.1 | functional.2) & (temp > 0)
    
  } else {
    #check for correct columns
    columns <- c("hgmd variant class","hgmd flanking count","clinvar flanking count","clinvar clinical significance","function","clingen haploinsufficiencydesc","clinvar pathogenic indel count", "clinvar pathogenic cnv count", "clinvar pathogenic snv splice count", "clinvar pathogenic snv nonsense count")
    
    #make sure all columns are present
    stopifnot(length(setdiff(normalized.name(columns),colnames(data))) ==0)
   
    #inclusion rule 1: 
    r1 <- sapply(data[normalized.name("hgmd variant class")], as.character)
    r1 <- sapply(r1, function(x) length(grep("dm",x)) > 0)
    
    #inclusion rule 2:
    suppresswarnings(temp1 <- sapply(data[normalized.name("hgmd flanking count")], as.numeric))
    temp1[is.na(temp1)] <- 0
    
    suppresswarnings(temp2 <- sapply(data[normalized.name("clinvar flanking count")], as.numeric))
    temp2[is.na(temp2)] <- 0
    r2 <- temp1 > 0 | temp2 > 0
    
    #inclusion rule 3:
    temp <- sapply(data[normalized.name("clinvar clinical significance")], as.character)
    r3 = (temp == "pathogenic") | (temp == "likely pathogenic") | (temp == "likely pathogenic;pathogenic")
    
    functional <- Filter.by.function(sapply(data[normalized.name("function")], as.character))
  
    #inclusion rule 4:
    r4 <- functional & (!eval(parse(text=paste0("data$",normalized.name("clingen haploinsufficiencydesc")))) %in% c("na","no evidence","little evidence"))
    
    #inclusion rule 5:
    suppresswarnings(temp <- sapply(data[normalized.name("clinvar pathogenic indel count")], as.numeric))
    temp[is.na(temp)] <- 0
    r5 <- functional & (temp > 0)
    
    #inclusion rule 6:
    suppresswarnings(temp <- sapply(data[normalized.name("clinvar pathogenic snv splice count")], as.numeric))
    temp[is.na(temp)] <- 0
    r6 <- functional & (temp > 0)
    
    #inclusion rule 7:
    suppresswarnings(temp <- sapply(data[normalized.name("clinvar pathogenic snv nonsense count")], as.numeric))
    temp[is.na(temp)] <- 0
    r7 <- functional & (temp > 0)  

    #inclusion rule 8:
    suppresswarnings(temp <- sapply(data[normalized.name("clinvar pathogenic cnv count")], as.numeric))
    temp[is.na(temp)] <- 0
    r8 <- functional & (temp > 0)  

  }
  
  r.all <- r1 | r2 | r3 | r4 | r5 | r6 | r7 | r8
  data <- data[r.all,]
  data
}
