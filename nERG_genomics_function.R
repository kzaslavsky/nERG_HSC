#genomics_functions for nERG project
#provided as is
#helpful for parsing HGVS and non-HGVS variant strings
#updates variant positions to latest reference transcripts
#maps to hg38
#creates VCF file
#Author: Kirill Zaslavsky
#version 0.1


## a useful function: rev() for strings
strReverse <- function(x)
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse="")

compareVarType <- function(var_chosen, cdna_change_oldvar)
{
  if(grepl(pattern = "^NM", var_chosen$variation_set$cdna_change))
  {
    #print("inside main if")
    #transcript.id.clinvar <- sub("^(NM.*?)(\\().*$", "\\1", var_chosen$variation_set$variation_name)
    cdna_change_altvar <- sub("^[^:]*:([^ ]+).*", "\\1", var_chosen$variation_set$cdna_change)
    
    print(cdna_change_altvar)
    
    var_type_alt <- parseHGVS(cdna_change_altvar)$type
    print(var_type_alt)
    
    var_type_old <- parseHGVS(cdna_change_oldvar)$type
    print(var_type_old)
    
    if(var_type_alt == var_type_old)
    {print("true")
      return(TRUE)}
    
    
  }else
  {
    print("no")
    return(FALSE)
  }
  return(FALSE)
}

#works better than VarType due to how entrez database is structured
compareVarTypeTitle <- function(var_chosen, cdna_change_oldvar, ensembl = ensembl)
{
  if(grepl(pattern = "^NM", var_chosen$title))
  {
    #print("inside main if") "NM_001256789.3(CACNA1F):c.3133dup (p.Leu1045fs)"
    transcript.id.clinvar <- sub("^(NM.*?)(\\().*$", "\\1", var_chosen$title)
    gene_name <- sub("^[^\\(]*\\(([^\\)]+)\\).*", "\\1", var_chosen$title)
    cdna_change_altvar <- sub("^[^:]*:([^ ]+).*", "\\1", var_chosen$title)
    
    print(paste0("altVar ", cdna_change_altvar))
    print(paste0("oldVar ", cdna_change_oldvar))
    
    var_alt <- parseHGVS(cdna_change_altvar)
    print(var_alt)
    
    var_old <- parseHGVS(cdna_change_oldvar)
    print(var_old)
    
    if(var_alt$type == var_old$type)
    {print("true")
      return(TRUE)}
    
    else if(var_old$type == "insertion" & var_alt$type == "singleDuplication" )
    {
      srch <- entrez_search(db="nucleotide", term = transcript.id.clinvar)
      transcript.fasta <- entrez_fetch(db="nucleotide", id = srch$ids, rettype = "fasta")
      tx.fasta <- read.fasta(file = textConnection(transcript.fasta))
      var_alt_fromfasta <- toupper(tx.fasta[[1]][var_alt$start+1])
      
      strand <- getStrandInfo(gene_name, ensembl)$strand
      var_alt_fromfasta_sense <- checkSense(var_alt_fromfasta, strand)
      
      if(var_alt_fromfasta_sense == var_old$variant)
      {
        print("clinvar suggested duplication is same as query insertion variant")
        return(TRUE)
      }
      
    }
    
  }else
  {
    print("no")
    return(FALSE)
  }
  return(FALSE)
}

#check DNA strand and get complement if necessary
checkSense <- function(dna, strand){
  #check if strand is 1 or -1
  if (strand == 1)
  {
    return(dna)
  }
  if (strand == -1)
  {
    antisense <- chartr("ACGT", "TGCA", dna)
    if(nchar(dna) > 1)
    {
      print(antisense)
      antisense <- strReverse(antisense)
    }
    return(antisense)
  }
}



#populate clinvar table to merge
processVariantClinvar <- function(var_chosen)
{
  #var_chosen <- entrez_summary(db="clinvar", id = var_id) 
  
  #print(paste("var chosen", var_id, var_chosen, sep = ", "))
  clinvar_to_paste <- init_clinvar_paste()
  
  clinvar_to_paste[1,"clinvar_uid"] <- var_chosen$uid
  clinvar_to_paste[1,"clinvar_title"] <- var_chosen$title
  
  if(grepl(pattern = "^NM", var_chosen$variation_set$cdna_change))
  {
    transcript.id.clinvar <- sub("^(NM.*?)(\\().*$", "\\1", var_chosen$variation_set$variation_name)
    cdna_change <- sub(".*:(.*?) .*", "\\1", var_chosen$variation_set$variation_name)
  }else
  {
    transcript.id.clinvar <- sub("^(NM.*?)(\\().*$", "\\1", var_chosen$variation_set$variation_name)
    cdna_change <- sub(".*:(.*?) .*", "\\1", var_chosen$variation_set$cdna_change)
  }
  
  clinvar_to_paste[1, "clinvar_transcript_id"] <- transcript.id.clinvar
  clinvar_to_paste[1,"clinvar_cdna_change"] <- cdna_change
  clinvar_to_paste[1,"clinvar_assembly_name"] <- var_chosen[6]$variation_set$variation_loc[[1]]$assembly_name[1]
  clinvar_to_paste[1,"clinvar_chr"] <- var_chosen[6]$variation_set$variation_loc[[1]]$chr[1]
  clinvar_to_paste[1,"clinvar_start"] <- var_chosen[6]$variation_set$variation_loc[[1]]$start[1]
  clinvar_to_paste[1,"clinvar_stop"] <- var_chosen[6]$variation_set$variation_loc[[1]]$stop[1]
  clinvar_to_paste[1,"clinvar_protchange"] <- var_chosen$protein_change
  clinvar_to_paste[1,"clinvar_clinsigreview"] <- var_chosen$clinical_significance$review_status
  clinvar_to_paste[1,"clinvar_clinsig"] <- var_chosen$clinical_significance$description
  
  return(clinvar_to_paste)
}

#initialize dummy variable to populate with clinvar data
init_clinvar_paste <- function()
{
  clinvar_to_paste <- data.frame(clinvar_uid = as.character(),
                                 clinvar_title = as.character(),
                                 clinvar_transcript_id = as.character(),
                                 clinvar_cdna_change = as.character(),
                                 clinvar_assembly_name = as.character(),
                                 clinvar_chr = as.character(),
                                 clinvar_start = as.character(),
                                 clinvar_stop = as.character(),
                                 clinvar_protchange = as.character(),
                                 clinvar_clinsigreview = as.character(),
                                 clinvar_clinsig = as.character()
  )
  
}

#split the gene variant into gene name, cDNA variant, PTN variant; PTN variant will frequently return NA
splitGeneVariant <- function(qgene) {
  qgene2 <- str_split(qgene, pattern = " ")[[1]]
  qgene2 <- qgene2[qgene2 != ""]
  return(qgene2)
}

#get info on whether transcript is on sense or anti-sense strand - important for calculating position
getStrandInfo <- function(gene_name,ensembl = ensembl)
{
  # Retrieve information about RS1 and its exons
  strand <- getBM(attributes = c('ensembl_gene_id', 
                                 'strand'
                                 #'hgnc_symbol'
                                 #'ensembl_transcript_id'
                                 #'transcript_mane_select',
                                 #'entrezgene_trans_name',
                                 #'transcript_is_canonical',
  ),
  filters = 'hgnc_symbol',
  values = gene_name,
  mart = ensembl)
  
}

#prevents code from stopping if biomart connection drops
getMart <- function(biomart, species) {
  e <- simpleError("")
  while(is(e, "simpleError")) {
    e <- tryCatch(useEnsembl(biomart, species), error = function(x) x)
  }
  e
}

#find clinvar match
findMatchClinvar <- function(qgene, mdf2_dump, ensembl = ensembl)
{
  # split variant into gene name and cDNA change
  qgene2 <- splitGeneVariant(qgene)
  if (qgene2[2] == "NA") {print ("cDNA variant is NA, exiting")
    return(list(var_info = NULL, mdf2_dump = qgene))
  }
  
  #search for variation IDs
  srch <- entrez_search(db="clinvar", term = qgene) 
  
  #exit if cannot find variant and just return thousands of ids related to gene name only
  pattern <- "^[A-Z0-9]+\\[All Fields\\](?!$)"
  if(!grepl(pattern, srch$QueryTranslation, perl = TRUE))
  {
    print(paste0(qgene2[2], " not found, exiting"))
    return(list(var_info = NULL, mdf2_dump = mdf2_dump))
  }
  
  #find the clinvar query where the cdna change is an exact match
  #sometimes clinvar query returns similar cdna changes e.g., 598 and 589
  #only do this if the clinvar query returned an answer
  if(!is_empty(srch$ids))
  { 
    for (id in srch$ids)
    {
      #load var_info
      var_info <- entrez_summary(db="clinvar", id = id) 
      
      #exit if cannot find in clinvar
      if(is_empty(var_info$variation_set$cdna_change)) 
      { print(paste0("No cDNA change recorded for variant uid ", id))
        next} 
      
      
      #store cdna change in 'var_change' variable
      var_change <- var_info$variation_set$cdna_change
      # print(id)
      # print(var_change)
      # print(qgene2)
      # 
      # find matching variation id to to final clinvar query
      if (length(var_change) > 1 & any(var_change %in% qgene2[2]))
      {
        var_id <- id
        print(paste0("exact match & more than 2 variants returned for ", var_id, ", will dump"))
        mdf2_dump <- rbind(mdf2_dump, data_frame(Variant = var_info$title))
        #print(mdf2.dump)
        next
      }else if (var_change == qgene2[2]) 
      {
        var_id <- id
        #transcript.id.clinvar <- sub("^(NM.*?)(\\().*$", "\\1", var_info$variation_set$variation_name)
        print(paste0("exact match 2 ", var_id))
        return(list(var_info = var_info, mdf2_dump = mdf2_dump))
        break
      }else if (grepl("^NM", var_change))
      {
        print(paste0("Clinvar Return Alternative Variant ", var_change, " instead of ", qgene2[2]))
        
        #check to ensure same variant type
        if (compareVarTypeTitle(var_info, qgene2[2], ensembl))
        {
          print("Old Var Type matches Alt Var Type")
          return(list(var_info = var_info, mdf2_dump = mdf2_dump))
          break
        }else
        {
          print("Old Var Type DOES NOT MATCH Alt Var Type")
          #return(list(var_info = NULL, mdf2_dump = mdf2_dump))
          next
        }
      }
      else if (grepl("^NM", var_info$title) & var_change != qgene2[2])
      {
        print(paste0("Clinvar Return Alternative Variant ", var_info$title, " instead of ", qgene2[2]))
        
        #check to ensure same variant type
        if (compareVarTypeTitle(var_info, qgene2[2], ensembl))
        {
          print("Old Var Type matches Alt Var Type")
          return(list(var_info = var_info, mdf2_dump = mdf2_dump))
          break
        }else
        {
          print("Old Var Type DOES NOT MATCH Alt Var Type")
          #return(list(var_info = NULL, mdf2_dump = mdf2_dump))
          next
        }
      }
      else if (grepl("^NC", var_change))
      {
        print(paste0("Clinvar Return Alternative Variant ", var_change, " instead of ", qgene2[2]))
        
        
        return(list(var_info = var_info, mdf2_dump = mdf2_dump))
        break
      }
      
      else 
      {
        print("incorrect variant")
        var_id <- NULL
        #no match / incorrect variant
      }
      
    } 
    return(list(var_info = NULL, mdf2_dump = mdf2_dump))
    
  }else { print("Not found in clinvar")
    mdf2_dump <- rbind(mdf2_dump, data_frame(Variant = qgene))
    return(list(var_info = NULL, mdf2_dump = mdf2_dump))
  }
}

#parse HGVS strings - adapted from Jochen Weile's function
parseHGVS <- function(strings,aacode=c(NA,1,3)) {
  
  #Check that parameters are valid
  if (!is.character(strings)) {
    stop("Input for 'parse' function must be a character vector! Found '",class(strings),"' instead.")
  }
  
  aacode <- aacode[[1]]
  if (!is.na(aacode) && !(aacode %in% c(1,3))) {
    warning("Invalid aacode parameter, defaulting to NA!")
    aacode <- NA
  }
  
  #Helper function: turns a list of lists (lol) in to a dataframe
  to.df <- function(lol) {
    colnames <- unique(do.call(c,lapply(lol,names)))
    columns <- lapply(colnames,function(cn) sapply(lol,function(row) {
      if (cn %in% names(row)) row[[cn]] else NA
    }))
    names(columns) <- colnames
    empty <- which(sapply(columns,function(xs)all(is.na(xs))))
    columns[empty] <- NULL
    do.call(data.frame,columns)
  }
  
  # ###
  # # Binds matrices of same size together to a 3D matrix, analogously
  # # to cbind and rbind.
  # #
  # zbind <- function(...) {
  # 	x <- list(...)
  # 	y <- array(0,dim=c(nrow(x[[1]]),ncol(x[[1]]),length(x)),dimnames=dimnames(x[[1]]))
  # 	for (i in 1:length(x)) y[,,i] <- x[[i]]
  # 	y
  # }
  
  
  ###
  # Function to *locally* excise regex groups from string vectors.
  # I.e. only extract the first occurrence of each group within each string.
  # x = string vector
  # re = regular expression with groups
  #
  extract.groups <- function(x, re) {
    matches <- regexpr(re,x,perl=TRUE)
    start <- attr(matches,"capture.start")
    end <- start + attr(matches,"capture.length") - 1
    do.call(cbind,lapply(1:ncol(start), function(i) {
      sapply(1:nrow(start),function(j){
        if (start[j,i] > -1) substr(x[[j]],start[j,i],end[j,i]) else NA
      })
    }))
  }
  
  # ###
  # # Function to *globally* excise regex groups from string vectors.
  # # x = string vector
  # # re = regular expression with groups
  # #
  # global.extract.groups <- function(x,re) {
  #     all.matches <- gregexpr(re,x,perl=TRUE)
  #     mapply(function(matches,x) {
  #         start <- attr(matches,"capture.start")
  #         end <- start + attr(matches,"capture.length") - 1
  #         apply(zbind(start,end),c(1,2),function(pos) substr(x,pos[[1]],pos[[2]]) )
  #     },matches=all.matches,x=x,SIMPLIFY=FALSE)
  # }
  
  ###
  # Helper function to split multi-mutant bodies into their individual 
  # elements. Returns a vector of strings containing these elements. 
  # An attribute "multi" is attached to the vector, detailing the type 
  # of multi-mutant
  # 
  splitMulti <- function(body) {
    if (regexpr("\\[.+\\];\\[.+\\]",body) > 0) {
      out <- strsplit(substr(body,2,nchar(body)-1),"\\];\\[")[[1]]
      attr(out,"multi") <- "trans"
    } else if (regexpr("\\[.+\\(;\\).+\\]",body) > 0) {
      out <- strsplit(substr(body,2,nchar(body)-1),"\\(;\\)")[[1]]
      attr(out,"multi") <- "unknown"
    } else if (regexpr("\\[.+;.+\\]",body) > 0) {
      out <- strsplit(substr(body,2,nchar(body)-1),";")[[1]]
      attr(out,"multi") <- "cis"
    } else {
      out <- body
      attr(out,"multi") <- "single"
    }
    return(out)
  }
  
  ###
  # Helper function:
  # Given an HGVS body and a list of regexes corresponding to types,
  # find the (first) matching type.
  findType <- function(body,types) {
    i <- 0
    found <- done <- FALSE
    while (!found && !done) {
      found <- regexpr(types[[i <- i+1]],body) > 0
      done <- i >= length(types)
    }
    if (found) {
      return(names(types)[[i]])
    } else {
      return("invalid")
    }
  }
  
  out <- lapply(strings,function(s) {
    
    if (regexpr("^[gcnmrp]\\.",s) < 1) {
      return(list(list(hgvs=s,subject="invalid",type="invalid")))
    }
    
    body <- substr(s,3,nchar(s))
    
    subbodies <- splitMulti(body)
    
    subjects <- c(
      g="genomic",c="coding",n="noncoding",
      m="mitochondrial",r="rna",p="protein"
    )
    subject <- subjects[[substr(s,1,1)]]
    
    if (subject=="genomic") {
      
      types <- c(
        substitution="\\d+[ACGT]>[ACGT]", singledeletion="^\\d+del$",
        deletion="\\d+_\\d+del$",inversion="\\d+_\\d+inv",
        duplication="\\d+_\\d+dup",insertion="\\d+_\\d+ins[ATCG]+",
        conversion="\\d+_\\d+con\\d+_\\d+",delins="\\d+_\\d+delins[ATCG]+",
        amplification="\\d+_\\d+\\[\\d+\\]"
      )
      
      phasing <- attr(subbodies,"multi")
      isMulti <- length(subbodies) > 1
      
      lapply(1:length(subbodies), function(i.multi) {
        body <- subbodies[[i.multi]]
        
        type <- findType(body,types)
        
        if (type == "substitution") {
          groups <- extract.groups(body,"(\\d+)([ACGT])>([ACGT])")[1,]
          position <- as.integer(groups[[1]])
          ancestral <- groups[[2]]
          variant <- groups[[3]]
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=position,ancestral=ancestral,variant=variant))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=position,
                        ancestral=ancestral,variant=variant))
          }
          
        } else if (type == "singledeletion") {
          groups <- extract.groups(body,"(\\d+)del")[1,]
          position <- as.integer(groups[[1]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=position))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=position))
          }
          
        } else if (type == "deletion") {
          groups <- extract.groups(body,"(\\d+)_(\\d+)del")
          start <- as.integer(groups[[1]])
          end <- as.integer(groups[[2]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,end=end))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,end=end))
          }
          
        } else if (type == "inversion") {
          groups <- extract.groups(body,"(\\d+)_(\\d+)inv")
          start <- as.integer(groups[[1]])
          end <- as.integer(groups[[2]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,end=end))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,end=end))
          }
          
        } else if (type == "duplication") {
          groups <- extract.groups(body,"(\\d+)_(\\d+)dup")
          start <- as.integer(groups[[1]])
          end <- as.integer(groups[[2]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,end=end))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,end=end))
          }
          
        } else if (type == "insertion") {
          groups <- extract.groups(body,"(\\d+)_(\\d+)ins([ATCG]+)")
          start <- as.integer(groups[[1]])
          end <- as.integer(groups[[2]])
          if (abs(end-start)!=1) {
            warning("Invalid insertion definition: 
							Start and end positions must be adjacent!")
          }
          variant <- groups[[3]]
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,end=end,variant=variant))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,end=end,
                        variant=variant))
          }
          
        } else if (type == "conversion") {
          groups <- extract.groups(body,"(\\d+)_(\\d+)con(\\d+)_(\\d+)")
          start <- as.integer(groups[[1]])
          end <- as.integer(groups[[2]])
          tStart <- as.integer(groups[[3]])
          tEnd <- as.integer(groups[[4]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,end=end,
                        templateStart=tStart,templateEnd=tEnd))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,end=end,
                        templateStart=tStart,templateEnd=tEnd))
          }
          
        } else if (type == "delins") {
          groups <- extract.groups(body,"(\\d+)_(\\d+)delins([ATCG]+)")
          start <- as.integer(groups[[1]])
          end <- as.integer(groups[[2]])
          variant <- groups[[3]]
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,end=end,variant=variant))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,end=end,
                        variant=variant))
          }
          
        } else if (type == "amplification") {
          groups <- extract.groups(body,"(\\d+)_(\\d+)\\[(\\d+)\\]")
          start <- as.integer(groups[[1]])
          end <- as.integer(groups[[2]])
          copies <- as.integer(groups[[3]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,end=end,copies=copies))
          } else {
            return(list(hgvs=s,subject=subject,type=type,
                        start=start,end=end,copies=copies))
          }
          
        } else {
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,type="invalid"))
          } else {
            return(list(hgvs=s,subject=subject,type="invalid"))
          }
        }
        
      })
      
    } else if (subject=="coding") {
      #coding needs to be handled separately from genomic, as the syntax may differ
      #e.g. it allows for offset descriptions relative to exon-intron borders
      
      types <- c(
        substitution="\\d+([+-]\\d+)?[ACGT]>[ACGT]", 
        singledeletion="^\\d+([+-]\\d+)?del[A-Z]*$",
        singleDuplication="^\\d+([+-]\\d+)?dup[A-Z]*$", # New type for single nucleotide duplication
        deletion="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?del$",
        inversion="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?inv",
        duplication="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?dup",
        insertion="^(\\d+)([+-]\\d+)?(_\\d+([+-]\\d+)?)?ins([A-Z]*)$",
        conversion="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?con\\d+([+-]\\d+)?_\\d+([+-]\\d+)?",
        delins="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?delins[ATCG]+",
        amplification="\\d+([+-]\\d+)?_\\d+([+-]\\d+)?\\[\\d+\\]"
      )
      
      phasing <- attr(subbodies,"multi")
      isMulti <- length(subbodies) > 1
      
      lapply(1:length(subbodies), function(i.multi) {
        body <- subbodies[[i.multi]]
        
        type <- findType(body,types)
        
        if (type == "substitution") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?([ACGT])>([ACGT])")[1,]
          position <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          ancestral <- groups[[3]]
          variant <- groups[[4]]
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=position,startIntron=intronOffset,ancestral=ancestral,
                        variant=variant))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=position,
                        startIntron=intronOffset,ancestral=ancestral,variant=variant))
          }
          
        } else if (type == "singledeletion") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?del")[1,]
          position <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=position,startIntron=intronOffset))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=position,
                        startIntron=intronOffset))
          }
        } else if (type == "singleDuplication") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?dup")[1,]
          position <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=position,startIntron=intronOffset))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=position,
                        startIntron=intronOffset))
          }
        } else if (type == "deletion") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?del")
          start <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          end <- as.integer(groups[[3]])
          intronOffset2 <- as.integer(groups[[4]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,startIntron=intronOffset,
                        end=end,endIntron=intronOffset2))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,
                        startIntron=intronOffset,end=end,endIntron=intronOffset2))
          }
        } else if (type == "inversion") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?inv")
          start <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          end <- as.integer(groups[[3]])
          intronOffset2 <- as.integer(groups[[4]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,startIntron=intronOffset,
                        end=end,endIntron=intronOffset2))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,
                        startIntron=intronOffset,end=end,endIntron=intronOffset2))
          }
        } else if (type == "duplication") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?dup")
          start <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          end <- as.integer(groups[[3]])
          intronOffset2 <- as.integer(groups[[4]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,startIntron=intronOffset,
                        end=end,endIntron=intronOffset2))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,
                        startIntron=intronOffset,end=end,endIntron=intronOffset2))
          }
        } else if (type == "insertion") {
          # groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?ins([ATCG]+)")
          groups <- extract.groups(body, "^(\\d+)([+-]\\d+)?(_\\d+([+-]\\d+)?)?ins([A-Z]*)$")
          start <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          end <- as.integer(groups[[3]])
          intronOffset2 <- as.integer(groups[[4]])
          variant <- groups[[5]]
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,startIntron=intronOffset,
                        end=end,endIntron=intronOffset2,variant=variant))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,
                        startIntron=intronOffset,end=end,endIntron=intronOffset2,variant=variant))
          }
        } else if (type == "conversion") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?con(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?")
          start <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          end <- as.integer(groups[[3]])
          intronOffset2 <- as.integer(groups[[4]])
          tStart <- as.integer(groups[[5]])
          intronOffset3 <- as.integer(groups[[6]])
          tEnd <- as.integer(groups[[7]])
          intronOffset4 <- as.integer(groups[[8]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,startIntron=intronOffset,
                        end=end,endIntron=intronOffset2,templateStart=tStart,
                        templateStartIntron=intronOffset3,templateEnd=tEnd,
                        templateEndIntron=intronOffset4))
          } else {
            return(list(hgvs=s,subject=subject,type=type,
                        start=start,startIntron=intronOffset,
                        end=end,endIntron=intronOffset2,templateStart=tStart,
                        templateStartIntron=intronOffset3,templateEnd=tEnd,
                        templateEndIntron=intronOffset4))
          }
        } else if (type == "delins") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?delins([ATCG]+)")
          start <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          end <- as.integer(groups[[3]])
          intronOffset2 <- as.integer(groups[[4]])
          variant <- groups[[5]]
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,startIntron=intronOffset,
                        end=end,endIntron=intronOffset2,variant=variant))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,
                        startIntron=intronOffset,end=end,endIntron=intronOffset2,variant=variant))
          }
        } else if (type == "amplification") {
          groups <- extract.groups(body,"(\\d+)([+-]\\d+)?_(\\d+)([+-]\\d+)?\\[(\\d+)\\]")
          start <- as.integer(groups[[1]])
          intronOffset <- as.integer(groups[[2]])
          end <- as.integer(groups[[3]])
          intronOffset2 <- as.integer(groups[[4]])
          copies <- as.integer(groups[[3]])
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=start,startIntron=intronOffset,
                        end=end,endIntron=intronOffset2,copies=copies))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=start,
                        startIntron=intronOffset,end=end,endIntron=intronOffset2,copies=copies))
          }
        } else {
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type="invalid"))
          } else {
            return(list(hgvs=s,subject=subject,type="invalid"))
          }
        }
        
      })
      
    } else if (subject=="protein") {
      
      one2three <- c(A="Ala",C="Cys",D="Asp",E="Glu",F="Phe",G="Gly",H="His",
                     I="Ile",K="Lys",L="Leu",M="Met",N="Asn",P="Pro",Q="Gln",R="Arg",
                     S="Ser",T="Thr",V="Val",W="Trp",Y="Tyr",`*`="Ter")
      three2one <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",
                     Gly="G",His="H",Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",
                     Ser="S",Thr="T",Trp="W",Tyr="Y",Val="V",Ter="*")
      codes <- paste(c(one2three,three2one[-21],"\\*"),collapse="|")
      
      types <- list(
        synonymous1="^=$",
        synonymous2=paste0("^(",codes,")(\\d+)=$"),
        substitution=paste0("^(",codes,")(\\d+)(",codes,")$"),
        singledeletion=paste0("^(",codes,")(\\d+)del$"),
        deletion=paste0("^(",codes,")(\\d+)_(",codes,")(\\d+)del$"),
        duplication=paste0("^(",codes,")(\\d+)_(",codes,")(\\d+)dup$"),
        insertion=paste0("^(",codes,")(\\d+)_(",codes,")(\\d+)ins((",codes,")+)$"),
        delins=paste0("^(",codes,")(\\d+)_(",codes,")(\\d+)delins((",codes,")+)$"),
        frameshift1=paste0("^(",codes,")(\\d+)fs$"),
        frameshift2=paste0("^(",codes,")(\\d+)(",codes,")fs(Ter|\\*)(\\d+)$")
      )
      
      phasing <- attr(subbodies,"multi")
      isMulti <- length(subbodies) > 1
      
      lapply(1:length(subbodies), function(i.multi) {
        body <- subbodies[[i.multi]]
        
        type <- findType(body,types)
        
        if (type == "synonymous1") {
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,type=type))
          } else {
            return(list(hgvs=s,subject=subject,type="synonymous"))
          }
        } else if (type == "synonymous2") {
          groups <- extract.groups(body,types$synonymous2)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          if (aa1 %in% c(one2three,three2one)) {
            if (is.na(aacode)) {
              #do nothing
            } else if (aacode == 1) {
              if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
            } else if (aacode ==3) {
              if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
            } else {
              #this should never happen, as it's supposed to be detected at start of function
              stop("Invalid aacode. If you see this, report this as a bug!")
            }
            if (isMulti) {
              return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                          type="synonymous",start=pos,ancestral=aa1))
            } else {
              return(list(hgvs=s,subject=subject,type="synonymous",start=pos,
                          ancestral=aa1))
            }
          } else {#not valid amino acid
            if (isMulti) {
              return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                          type="invalid"))
            } else {
              return(list(hgvs=s,subject=subject,type="invalid"))
            }
          }
        } else if (type == "substitution") {
          groups <- extract.groups(body,types$substitution)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          aa2 <- groups[[3]]
          if (aa1 %in% c(one2three,three2one) && aa2 %in% c(one2three,three2one)) {
            if (is.na(aacode)) {
              #do nothing
            } else if (aacode == 1) {
              if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
              if (nchar(aa2) == 3) aa2 <- three2one[[aa2]]
            } else if (aacode ==3) {
              if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
              if (nchar(aa2) == 1) aa2 <- one2three[[aa2]]
            } else {
              #this should never happen, as it's supposed to be detected at start of function
              stop("Invalid aacode. If you see this, report this as a bug!")
            }
            if (isMulti) {
              return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                          type=type,start=pos,ancestral=aa1,variant=aa2))
            } else {
              return(list(hgvs=s,subject=subject,type=type,start=pos,
                          ancestral=aa1,variant=aa2))
            }
          } else {#not valid amino acid
            if (isMulti) {
              return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                          type="invalid"))
            } else {
              return(list(hgvs=s,subject=subject,type="invalid"))
            }
          }
          
        } else if (type == "singledeletion") {
          groups <- extract.groups(body,types$singledeletion)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          if (is.na(aacode)) {
            #do nothing
          } else if (aacode == 1) {
            if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
          } else if (aacode == 3) {
            if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
          } else {
            #this should never happen, as it's supposed to be detected at start of function
            stop("Invalid aacode. If you see this, report this as a bug!")
          }
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=pos,ancestral=aa1))
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=pos,ancestral=aa1))
          }
        } else if (type == "deletion") {
          groups <- extract.groups(body,types$deletion)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          aa2 <- groups[[3]]
          pos2 <- as.integer(groups[[4]])
          if (is.na(aacode)) {
            #do nothing
          } else if (aacode == 1) {
            if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
            if (nchar(aa2) == 3) aa2 <- three2one[[aa2]]
          } else if (aacode == 3) {
            if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
            if (nchar(aa2) == 1) aa2 <- one2three[[aa2]]
          } else {
            #this should never happen, as it's supposed to be detected at start of function
            stop("Invalid aacode. If you see this, report this as a bug!")
          }
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=pos,ancestral=aa1,end=pos2,ancestral2=aa2))	
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=pos,
                        ancestral=aa1,end=pos2,ancestral2=aa2))	
          }
          
        } else if (type == "duplication") {
          groups <- extract.groups(body,types$duplication)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          aa2 <- groups[[3]]
          pos2 <- as.integer(groups[[4]])
          if (is.na(aacode)) {
            #do nothing
          } else if (aacode == 1) {
            if (nchar(aa1) == 3) aa1 <- three2one[[aa1]]
            if (nchar(aa2) == 3) aa2 <- three2one[[aa2]]
          } else if (aacode == 3) {
            if (nchar(aa1) == 1) aa1 <- one2three[[aa1]]
            if (nchar(aa2) == 1) aa2 <- one2three[[aa2]]
          } else {
            #this should never happen, as it's supposed to be detected at start of function
            stop("Invalid aacode. If you see this, report this as a bug!")
          }
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=pos,ancestral=aa1,end=pos2,ancestral2=aa2))	
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=pos,
                        ancestral=aa1,end=pos2,ancestral2=aa2))	
          }
          
        } else if (type == "insertion") {
          groups <- extract.groups(body,types$insertion)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          aa2 <- groups[[3]]
          pos2 <- as.integer(groups[[4]])
          insert <- groups[[5]]
          #TODO: Implement code conversion 
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=pos,ancestral=aa1,end=pos2,ancestral2=aa2,variant=insert))	
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=pos,
                        ancestral=aa1,end=pos2,ancestral2=aa2,variant=insert))	
          }
          
        } else if (type == "delins") {
          groups <- extract.groups(body,types$delins)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          aa2 <- groups[[3]]
          pos2 <- as.integer(groups[[4]])
          insert <- groups[[5]]
          #TODO: Implement code conversion 
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type=type,start=pos,ancestral=aa1,end=pos2,ancestral2=aa2,variant=insert))	
          } else {
            return(list(hgvs=s,subject=subject,type=type,start=pos,
                        ancestral=aa1,end=pos2,ancestral2=aa2,variant=insert))	
          }
          
        } else if (type == "frameshift1") {
          groups <- extract.groups(body,types$frameshift1)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          #TODO: Implement code conversion 
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type="frameshift",start=pos,ancestral=aa1))
          } else {
            return(list(hgvs=s,subject=subject,type="frameshift",start=pos,
                        ancestral=aa1))
          }
          
        } else if (type == "frameshift2") {
          groups <- extract.groups(body,types$frameshift2)
          aa1 <- groups[[1]]
          pos <- as.integer(groups[[2]])
          aa2 <- groups[[3]]
          term <- as.integer(groups[[5]])
          #TODO: Implement code conversion 
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type="frameshift",start=pos,ancestral=aa1,variant=aa2,end=term))
          } else {
            return(list(hgvs=s,subject=subject,type="frameshift",start=pos,
                        ancestral=aa1,variant=aa2,end=term))
          }
          
        } else {#unmatched type
          if (isMulti) {
            return(list(hgvs=s,subject=subject,phasing=phasing,multiPart=i.multi,
                        type="invalid"))
          } else {
            return(list(hgvs=s,subject=subject,type="invalid"))
          }
        }
        
      })
      
    } else if (subject=="noncoding") {
      #FIXME: These need to be list of lists to match postprocessing
      return(list(list(hgvs=s,subject="not_implemented",type="not_implemented")))
    } else if (subject=="mitochondrial") {
      return(list(list(hgvs=s,subject="not_implemented",type="not_implemented")))
    } else if (subject=="rna") {
      return(list(list(hgvs=s,subject="not_implemented",type="not_implemented")))
    } else {#unmatched subject, shouldn't happen
      stop("Unmatched subject! If you see this, report it as a bug!")
    }
  })
  
  #demote multimutants to enforce simple list structure
  multiLengths <- sapply(out,length)
  ids <- do.call(c,lapply(1:length(multiLengths),
                          function(i) if (multiLengths[[i]]==1) as.character(i) else paste0(i,".",1:multiLengths[[i]])
  ))
  out2 <- do.call(c,out)
  names(out2) <- ids
  
  
  
  to.df(out2)
}


# function to manually calculate genomic position from a selected canonical transcipt
calcGenomicPosition <- function(var_loc, exons.ens.list, transcripts.ens, canon_transcript,ensembl)
{
  # Initialize variables to keep track of the cumulative exon length and the genomic position
  cumulative_exon_length <- 0
  genomic_position <- 0
  #print(canon_transcript)
  
  # find/define canonical transcript
  # start with labeled canonical transcript
  # use canonical transcript to find cDNA position, then derive genomic position
  # if there is no canonical transcript defined, or it differs from the one used from variant notation, figure out using ccds
  # print(exons.ens.list)
  for (exon_list in exons.ens.list)
  {
    if (unique(exon_list$transcript_is_canonical == 1))
    {
      #print(exon_list)
      #print(canon_transcript)
      canon_transcript <- exon_list
      break
    }
  }
  #print(canon_transcript)
  #print(nrow(canon_transcript))
  
  # Define the cDNA position - for example, 589
  # in c.XXXC>T notation, the position is relative to the start codon of coding sequence, which starts at a specific position
  # therefore, cDNA position is actually - translation start + XXX
  # first attempt is to try and define it on the canonical transcript
  # if canonical transcript is linked to NCBI refseq transcript via MANE project, the 'transcipt_mane_select' will have a value
  # if there is no MANE transcript, define transcript using ccds - to be added later
  transcr.term <- subset(transcripts.ens[grepl("^NM", transcripts.ens$transcript_mane_select),], 
                         transcript_is_canonical == 1)$transcript_mane_select
  
  cds_interval <- getCDS(transcr.term)
  
  srch <- entrez_search(db="nucleotide", term = transcr.term) # Find refseq NCBI entry for MANE transcript
  
  # pull transcript genbank from NCBI to find CDS coordinates
  transcript.xml <- entrez_fetch(db="nucleotide", id = srch$ids, rettype = "gb", retmode = "xml",
                                 parsed = TRUE)
  transcript.xml2 <- xmlParse(transcript.xml)
  
  # pull out coding sequence coordinates
  cds_features <- xpathSApply(transcript.xml2, "//GBFeature[GBFeature_key='CDS']")
  
  # Extracting GBinterval_from values
  
  cds_interval <- lapply(cds_features, function(feature) {
    val_from <- as.numeric(xpathSApply(feature, "./GBFeature_intervals/GBInterval/GBInterval_from", xmlValue))
    val_to <- as.numeric(xpathSApply(feature, "./GBFeature_intervals/GBInterval/GBInterval_to", xmlValue))
    acc_num <- xpathSApply(feature, "./GBFeature_intervals/GBInterval/GBInterval_accession", xmlValue)
    # If there are multiple intervals, this will extract all of them
    df.ret <- (data.frame(accnum = acc_num, from = as.numeric(val_from), to = as.numeric(val_to)))
    return(df.ret)
  })[[1]]
  
  # defines the cdna position
  cdna_position <- var_loc + cds_interval[1,"from"] - 1
  
  # Loop through the canonical transcript exons to find where the cDNA position falls
  
  for(i in 1:nrow(canon_transcript)) {
    exon_length <- canon_transcript$exon_chrom_end[i] - canon_transcript$exon_chrom_start[i] + 1
    #print(paste0("exon length = " ,exon_length))
    if (cumulative_exon_length + exon_length >= cdna_position) {
      # Calculate the genomic position
      if (canon_transcript$strand[i] == 1) {  # Forward strand
        genomic_position <- canon_transcript$exon_chrom_start[i] + (cdna_position - cumulative_exon_length) - 1
        #print(paste0("strand is 1 for variant ", qgene))
        
      } else {  # Reverse strand
        genomic_position <- canon_transcript$exon_chrom_end[i] - (cdna_position - cumulative_exon_length) + 1
        #print(canon_transcript[[1]]$exon_chrom_end[i])
        #print(paste0("strand is -1 for variant ", qgene))
      }
      break
    }
    cumulative_exon_length <- cumulative_exon_length + exon_length
    #print(cumulative_exon_length)
  }
  
  return(genomic_position)
}

#getCDS interval (start, end) for given accession number
getCDS <- function(transcr.term)
{
  srch <- entrez_search(db="nucleotide", term = transcr.term) # Find refseq NCBI entry for MANE transcript
  
  # pull transcript genbank from NCBI to find CDS coordinates
  transcript.xml <- entrez_fetch(db="nucleotide", id = srch$ids, rettype = "gb", retmode = "xml",
                                 parsed = TRUE)
  transcript.xml2 <- xmlParse(transcript.xml)
  
  # pull out coding sequence coordinates
  cds_features <- xpathSApply(transcript.xml2, "//GBFeature[GBFeature_key='CDS']")
  
  # Extracting GBinterval_from values
  
  cds_interval <- lapply(cds_features, function(feature) {
    val_from <- as.numeric(xpathSApply(feature, "./GBFeature_intervals/GBInterval/GBInterval_from", xmlValue))
    val_to <- as.numeric(xpathSApply(feature, "./GBFeature_intervals/GBInterval/GBInterval_to", xmlValue))
    acc_num <- xpathSApply(feature, "./GBFeature_intervals/GBInterval/GBInterval_accession", xmlValue)
    # If there are multiple intervals, this will extract all of them
    df.ret <- (data.frame(accnum = acc_num, from = as.numeric(val_from), to = as.numeric(val_to)))
    return(df.ret)
  })[[1]]
  
  return(cds_interval)
}

#function to make VCF line
constructVCF_line <- function(CHROM, POS, REF, ALT, FORMAT = "GT:DP", NA00001 = "1:150", INFO = ".")
{data_frame(CHROM = CHROM,
            POS = POS,
            ID = ".",
            REF = REF,
            ALT = ALT,
            QUAL = ".",
            FILTER = ".",
            INFO = INFO,
            FORMAT = "GT:DP",
            NA00001 = "1:150")
}

addPtID_toVCF <- function(vcf_table)
{
  for (i in 1:dim(vcf_table)[1])
  {
    if(vcf_table[i, INFO] == "." | is.na(vcf_table[i, INFO]))
    {
      infopaste <- paste0("PtID=", vcf_table[i,PtID])
    }else
    {
      infopaste = paste0("PtID=",vcf_table[i,PtID], ";",vcf_table[i, INFO])
    }
    vcf_table[i, INFO := infopaste]
    #print(i)
    #print(vcf_table[i, INFO])
  }
  return(vcf_table)
}

removePtID_fromVCF <- function(vcf_table) {
  # Loop through each row of the VCF table
  for (i in 1:nrow(vcf_table)) {
    # Check if the INFO field starts with "PtID="
    if (startsWith(vcf_table$INFO[i], "PtID=")) {
      # Remove the PtID= part including its value and the following semicolon
      vcf_table$INFO[i] <- sub("^PtID=[^;]+;?", "", vcf_table$INFO[i])
    }
  }
  return(vcf_table)
}



# get exons for all transcripts with refseq matching IDs via MANE or CCDS
getExonsList <- function(exons.ens, transcripts.ens, transcripts.ens.list)
{
  exons.ens.list <- lapply(transcripts.ens.list, function(tx)
  {
    for (transcr in tx)
    {
      #print(exons)
      exons.ens.list <- exons.ens[exons.ens$"ensembl_transcript_id" == transcr,]
      
      return(exons.ens.list)
    }
  })
  
  #####
  # Order exons in ascending order, and label canoncial transcript if present
  # functions to product appropriate VCF file output do check for presence of canonical transcript
  for (i in 1:length(exons.ens.list))
  {
    exons.ens.list[[i]] <- exons.ens.list[[i]][order(exons.ens.list[[i]]$rank),]
    exons.ens.list[[i]]$transcript_is_canonical <- 0
    if(unique(exons.ens.list[[i]]$ensembl_transcript_id) == transcripts.ens[!is.na(transcripts.ens[transcripts.ens$transcript_is_canonical==1, "ensembl_transcript_id"]),"ensembl_transcript_id"])
    {
      exons.ens.list[[i]]$transcript_is_canonical <- 1
    }
  }
  
  return(exons.ens.list)
}

#make a list of transcripts that have matching refseq NM accession numbers via MANE or ccds match
ExtractMatchedRefseqEnsTranscripts <- function(transcripts.ens)
{
  # if(any(grepl("^NM", transcripts.ens$transcript_mane_select))) {
  #   transcripts.ens.list <- as.list(transcripts.ens[grepl("^NM", transcripts.ens$transcript_mane_select) ,
  #                                                   "ensembl_transcript_id"])
  # } else 
  # {
  transcripts.ens.list <- as.list(transcripts.ens[grepl("^NM", transcripts.ens$transcript_ccds_match) ,
                                                  "ensembl_transcript_id"])
  # }
  return(transcripts.ens.list)
}

#function to retrieve exons & transcripts associated with each gene, as well as match to ncbi accession numbers by ccds
retrieveGeneInfo <- function(qgene2, ensembl) {
  exons <- getBM(attributes = c('ensembl_gene_id',  'transcription_start_site','exon_chrom_start', 'exon_chrom_end', 'strand', 'ensembl_transcript_id', 'rank', 'cds_start', 'cds_end', 'cds_length'),
                 filters = 'hgnc_symbol',
                 values = qgene2[[1]],
                 mart = ensembl)
  
  transcripts.ens <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id', 'ccds', 'transcript_mane_select', 'entrezgene_trans_name', 'ensembl_peptide_id', 'transcript_is_canonical', 'ensembl_peptide_id', 'chromosome_name'),
                           filters = 'hgnc_symbol',
                           values = qgene2[[1]],
                           mart = ensembl)
  # remove transcripts with inappropriately labeled chromosome names
  transcripts.ens <- transcripts.ens[!grepl(pattern = "^H",transcripts.ens$chromosome_name), ]
  
  
  
  # if there is no refeq ID linked, search using CCDS
  # first get ccds associated with each transcript - working
  # find the refseq transcirpt for those CCDS - working
  transcripts.ens$transcript_ccds_match <- rep("", dim(transcripts.ens)[1])
  transcripts.ens$transcript_ccds_cap <- rep("", dim(transcripts.ens)[1])
  ccds.ids <- unique(transcripts.ens[transcripts.ens$ccds != "",]$ccds)
  
  print(paste0("CCDS IDS ", ccds.ids))
  if(any(!is.na(ccds.ids)))
  {
    
    for (ccds.id in ccds.ids)
    {
      srch.ccds <- entrez_search(db = "nucleotide", term = ccds.id)  
      for (srch.ccds.id in srch.ccds$ids)
      {
        #print(srch.ccds.id)
        tx.id.search <- entrez_summary(db = "nuccore", id = srch.ccds.id)
        if (grepl(pattern = "^NM", tx.id.search$caption))
        {
          #NM.id <- data_frame(CCDS.ID = ccds.id, NM.ID = tx.id.search$caption)
          transcripts.ens[transcripts.ens$ccds == ccds.id, ]$transcript_ccds_match = tx.id.search$accessionversion
          transcripts.ens[transcripts.ens$ccds == ccds.id, ]$transcript_ccds_cap = tx.id.search$caption
          print(paste0("found match for ccds id ", ccds.id))
          break
        }
      }
    }
  }
  
  #failsafe to ensure MANE transcript is not overwritten if there are multiple CCDS - ref accnum matches
  transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), ]$transcript_ccds_match <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), ]$transcript_mane_select
  accession <- sub("\\..*$", "", transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), ]$transcript_ccds_match)
  transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), ]$transcript_ccds_cap <- accession
  
  
  return(list(exons = exons, transcripts = transcripts.ens))
}


#get Fasta seeuqnce for specific accesion number
getFasta <- function(tx)
{
  tx.srch <- entrez_search(db="nucleotide", term = tx)
  tx.fasta <- entrez_fetch(db="nucleotide", id = tx.srch$ids, rettype = "fasta")
  tx.fasta <- read.fasta(file = textConnection(tx.fasta))
  return(tx.fasta)
}

#find differences in transcript sequences - useful to troubleshooting if transcripts differ between accession numbers slightly
compareTranscr <- function(vec1, vec2){
  # Determining the length of the shorter vector
  min_length <- min(length(vec1), length(vec2))
  
  # Finding positions where the vectors differ
  diff_positions <- which(vec1[1:min_length] != vec2[1:min_length])
  
  return_diff <- data_frame(pos_ens = numeric(),
                            ens = character(),
                            var = character())
  
  # Displaying the differences
  for (i in 1:length(diff_positions)) {
    return_diff[i, "pos_ens"] <- diff_positions[i]
    return_diff[i, "ens"] <- vec1[diff_positions[i]]
    return_diff[i, "var"] <- vec2[diff_positions[i]]
    
    # cat(sprintf("Position %d: %s (vec1) vs %s (vec2)\n", pos, vec1[pos], vec2[pos]))
  }
  
  # # If the vectors are of different lengths, handle the extra elements
  # if (length(vec1) != length(vec2)) {
  #   longer_vec <- ifelse(length(vec1) > length(vec2), vec1, vec2)
  #   extra_positions <- (min_length+1):length(longer_vec)
  #   for (pos in extra_positions) {
  #     cat(sprintf("Extra element in longer vector at position %d: %s\n", pos, longer_vec[pos]))
  #   }
  # }
  
  return(return_diff)
}

#calculate genomic position using ensembl functions
calcGenomicPositionEns <- function(var_loc, canon_transcript)
{
  ## Define the position within the CDS of the transcript
  var_cds <- IRanges(start = c(var_loc), 
                     names = rep(canon_transcript, 1))
  var_transcr <- cdsToTranscript(var_cds, ahEdb)
  var_gen <- transcriptToGenome(var_transcr, ahEdb)
  return(var_gen)
}

#calculate genomic position for large deletions using ensembl functions
calcGenomicPositionEns.del <- function(var_start, var_end, canon_transcript)
{
  ## Define the position within the CDS of the transcript
  var_cds <- IRanges(start = c(var_start), end = c(var_end),
                     names = rep(canon_transcript, 1))
  var_transcr <- cdsToTranscript(var_cds, ahEdb)
  var_gen <- transcriptToGenome(var_transcr, ahEdb)
  return(var_gen)
}


checkVarVsMismatch <- function(var_parsed, mmt)
{
  #check var_parsed$start is not in range of patternStart and PatternEnd
  #returns TRUE if variant position is NOT in the mismatch region
  #return FALSE if it is - therefore, variant is likely incorrect now
  for(i in 1:dim(mmt)[1])
  {
    if(mmt[i,"PatternStart"] <= var_parsed$start & var_parsed$start <= mmt[i,"PatternEnd"])
    {
      return(FALSE)
    }
  }
  return(TRUE)
}


getRefSeq <- function(server, var_start, var_end, chr)
{
  
  #pulls ref seq from specified coordinates (hg38)
  
  #construct query to pull ref seq
  ens.query <- paste0("/sequence/region/human/", 
                      chr,
                      ":",
                      var_start,
                      "..",
                      var_end
                      #":",``
                      #"?"
  )
  
  seq <- GET(paste(server, ens.query, sep = ""), content_type("text/plain"))
  
  stop_for_status(seq)
  
  
  REF.seq <- print(content(seq))
  return(REF.seq)
}

#parse intronic formats written as IVS (old standard)
parseIVS <- function(variant_string)
{
  pattern <- "^IVS(\\d+)([-+]?\\d+)([A-Z])>([A-Z])$"
  
  # Use str_match from the stringr package to extract parts
  parsed_parts <- stringr::str_match(variant_string, pattern)
  
  # Create the data frame
  variant_data <- data.frame(
    variant_name = variant_string,
    type = "intervening sequence",
    IVS_number = as.numeric(parsed_parts[,2]),
    IntronOffset = as.numeric(parsed_parts[,3]),
    ancestral = parsed_parts[,4],
    variant = parsed_parts[,5],
    stringsAsFactors = FALSE # Avoid factors to make it easier to work with the data
  )
  
  return(variant_data)
}

#gets genomic position for simple SNV, deletion, insertion, duplication(single) i.e. only start coord needed
getGenomicPosition.simple <- function(var_tx_rs_accnum, var_loc, transcripts.ens, exons.ens.list, ensembl)
{
  if(any(grepl(var_tx_rs_accnum,transcripts.ens$transcript_ccds_match)))
  {
    #tx ref accnum matches ens version by ccds
    # ensembl way
    canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == var_tx_rs_accnum), "ensembl_transcript_id"]
    
    #possible that ccds_match will return multiple ens transcripts. take the one that matches MANE or otherwise the first one
    if(length(canon_transcript) > 1 & 
       transcripts.ens[which(transcripts.ens$transcript_is_canonical ==1), "transcript_mane_select"] == var_tx_rs_accnum)
    {
      canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_mane_select == var_tx_rs_accnum), "ensembl_transcript_id"]
    }
    else
    {
      canon_transcript <- canon_transcript[1]
    }
    
    genomic_position_ens_result <- calcGenomicPositionEns(var_loc, canon_transcript)
    print("supplied var ref transcr matches ref transcr in ensembl table; proceed with esnembl way of deriving position")
    
  }
  else if(var_tx_rs_accnum != "NA")
  {
    #tx ref accnum does not match ens canon_tx
    #manual way
    #pull NCBI accession number linked by CCDS to accession of supplied reference transcript in esnembl db (ie latest version)
    accession <- sub("\\..*$", "", var_tx_rs_accnum)
    pattern <- paste0("^", accession)
    ens_vers_transcr <- unique(transcripts.ens[grep(pattern, transcripts.ens$transcript_ccds_cap), "transcript_ccds_match"])
    use_canon_tx<-0
    
    #if there is no match, assign canoncial to ens_ver_transcr
    if(is_empty(ens_vers_transcr))
    {
      ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_mane_select"]
      use_canon_tx<-1
      print("no match in ensembl, using ensembl canonical transcript")
    }
    
    #pull var_ref_transcr fasta and cds_interval
    var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
    var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
    var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
    var_tx_rs_accnum.cds <- paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = "")
    
    #pull ens_ref_transr fasta and cds interval
    ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
    ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
    ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
    ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
    
    # #pairwise align full tx - unnecessary
    # ens_vs_var_transcr <- pairwiseAlignment(paste0(ens_vers_transcr.fasta[[1]], collapse = ""),
    #                                         paste0(var_tx_rs_accnum.fasta[[1]], collapse = ""))
    
    
    #pairwaise aling CDS
    #pattern is ens
    #subject is var
    
    ens_vs_var_transcr.cds <- pairwiseAlignment(
      ens_vers_transcr.cds,
      var_tx_rs_accnum.cds
    )
    
    # check if CDS length equal
    
    if (var_tx_rs_accnum.cds_len == ens_vers_transcr.cds_len)
    {
      
      #cDNA position on the ens tx
      #cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
      
      #check if percent identity is 100% i.e. CDS are identical between versions
      if (pid(ens_vs_var_transcr.cds) == 100 )
      {
        #standard ensembl way to calc genomic position
        canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
        if(use_canon_tx == 1)
        {
          canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
        }
        genomic_position_ens_result <- calcGenomicPositionEns(var_loc, canon_transcript)
        
      }
      else
      {
        #check if var_parsed$start is not in mismatch table
        mmt <- mismatchTable(ens_vs_var_transcr.cds)
        if (checkVarVsMismatch(var_parsed, mmt))
        {
          #standard ensembl way to calc genomic position
          canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
          genomic_position_ens_result <- calcGenomicPositionEns(var_loc, canon_transcript)
          
        }
        else
        {
          #error out
          print("ref seqeunce changed, variant is likely incorrect")
          
          #dump qgene to some sort of log
          return(NULL)
          
        }
      }
    }
    else 
    {
      #CDS length is not equal
      
      canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
      var_start_adjusted <- adjustVarLoc2(var_loc, ens_vs_var_transcr.cds)
      
      genomic_position_ens_result <- calcGenomicPositionEns(var_start_adjusted, canon_transcript)
      
      #new cDNA position on the ens tx
      # cdna_position <- var_loc_adjusted + ens_vers_transcr.cds_int[1,"from"] - 1
      
      
      print("ref CDS sequence length changed, position adjusted for indels between transcripts")  
      
    }
    
  }
  else if(var_tx_rs_accnum == "NA")
  {
    #no tx ref accnum supplied, use ens canon_tx
    canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
    #ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_ccds_match"]
    
    
    genomic_position_ens_result <- calcGenomicPositionEns(var_loc, canon_transcript)
    
    # genomic_position_man <- calcGenomicPosition(var_loc =  var_parsed$start, 
    #                                             exons.ens.list =  exons.ens.list, 
    #                                             transcripts.ens =  transcripts.ens, 
    #                                             canon_transcript = canon_transcript, 
    #                                             ensembl = ensembl)
    # 
    # if(genomic_position_man == start(genomic_position_ens_result[[1]]) ) #matches; likely correct transcript chosen
    # {
    #   print("No var trascript ref supplied, used canonical transcript - both ens and manual ways align")
    # }
    # else
    # {
    #   print("no match between ens and man")
    #   return(NULL)
    # }
  }
  #return(genomic_position_ens_result)
  return(list(gen_pos = genomic_position_ens_result, canon_tx = canon_transcript))
  
}

#gets genomic position for deletions, need both start and end coords
getGenomicPosition.del <- function(var_tx_rs_accnum, var_loc_start, var_loc_end, transcripts.ens, exons.ens.list, ensembl)
{
  if(any(grepl(var_tx_rs_accnum,transcripts.ens$transcript_ccds_match)))
  {
    #tx ref accnum matches ens version by ccds
    # ensembl way
    canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == var_tx_rs_accnum), "ensembl_transcript_id"]
    
    #possible that ccds_match will return multiple ens transcripts. take the one that matches MANE or otherwise the first one
    if(length(canon_transcript) > 1 & 
       transcripts.ens[which(transcripts.ens$transcript_is_canonical ==1), "transcript_mane_select"] == var_tx_rs_accnum)
    {
      canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_mane_select == var_tx_rs_accnum), "ensembl_transcript_id"]
    }
    else
    {
      canon_transcript <- canon_transcript[1]
    }
    
    genomic_position_ens_result <- calcGenomicPositionEns.del(var_loc_start, var_loc_end, canon_transcript)
    print("supplied var ref transcr matches ref transcr in ensembl table; proceed with esnembl way of deriving position")
    
  }
  else if(var_tx_rs_accnum != "NA")
  {
    #tx ref accnum does not match ens canon_tx
    #manual way
    #pull NCBI accession number linked by CCDS to accession of supplied reference transcript in esnembl db (ie latest version)
    accession <- sub("\\..*$", "", var_tx_rs_accnum)
    pattern <- paste0("^", accession)
    ens_vers_transcr <- unique(transcripts.ens[grep(pattern, transcripts.ens$transcript_ccds_cap), "transcript_ccds_match"])
    use_canon_tx<-0
    
    #if there is no match, assign canoncial to ens_ver_transcr
    if(is_empty(ens_vers_transcr))
    {
      ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_mane_select"]
      use_canon_tx<-1
      print("no match in ensembl, using ensembl canonical transcript")
    }
    
    #pull var_ref_transcr fasta and cds_interval
    var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
    var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
    var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
    var_tx_rs_accnum.cds <- paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = "")
    
    #pull ens_ref_transr fasta and cds interval
    ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
    ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
    ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
    ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
    
    # #pairwise align full tx - unnecessary
    # ens_vs_var_transcr <- pairwiseAlignment(paste0(ens_vers_transcr.fasta[[1]], collapse = ""),
    #                                         paste0(var_tx_rs_accnum.fasta[[1]], collapse = ""))
    
    
    #pairwaise aling CDS
    #pattern is ens
    #subject is var
    
    ens_vs_var_transcr.cds <- pairwiseAlignment(
      ens_vers_transcr.cds,
      var_tx_rs_accnum.cds
    )
    
    # check if CDS length equal
    
    if (var_tx_rs_accnum.cds_len == ens_vers_transcr.cds_len)
    {
      
      #cDNA position on the ens tx
      #cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
      
      #check if percent identity is 100% i.e. CDS are identical between versions
      if (pid(ens_vs_var_transcr.cds) == 100 )
      {
        #standard ensembl way to calc genomic position
        canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
        if(use_canon_tx == 1)
        {
          canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
        }
        genomic_position_ens_result <- calcGenomicPositionEns.del(var_loc_start, var_loc_end, canon_transcript)
        
      }
      else
      {
        #check if var_parsed$start is not in mismatch table
        mmt <- mismatchTable(ens_vs_var_transcr.cds)
        if (checkVarVsMismatch(var_parsed, mmt))
        {
          #standard ensembl way to calc genomic position
          canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
          genomic_position_ens_result <- calcGenomicPositionEns.del(var_loc_start, var_loc_end, canon_transcript)
          
        }
        else
        {
          #error out
          print("ref seqeunce changed, variant is likely incorrect")
          
          #dump qgene to some sort of log
          return(NULL)
          
        }
      }
    }
    else 
    {
      #CDS length is not equal
      
      canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
      var_start_adjusted <- adjustVarLoc2(var_loc_start, ens_vs_var_transcr.cds)
      var_end_adjusted <- adjustVarLoc2(var_loc_end, ens_vs_var_transcr.cds)
      
      
      genomic_position_ens_result <- calcGenomicPositionEns.del(var_start_adjusted, var_end_adjusted, canon_transcript)
      
      #new cDNA position on the ens tx
      # cdna_position <- var_loc_adjusted + ens_vers_transcr.cds_int[1,"from"] - 1
      
      
      print("ref CDS sequence length changed, position adjusted for indels between transcripts")  
      
    }
    
  }
  else if(var_tx_rs_accnum == "NA")
  {
    #no tx ref accnum supplied, use ens canon_tx
    canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
    #ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_ccds_match"]
    
    
    genomic_position_ens_result <- calcGenomicPositionEns.del(var_loc_start, var_loc_end, canon_transcript)
    
  }
  #return(genomic_position_ens_result)
  return(list(gen_pos = genomic_position_ens_result, canon_tx = canon_transcript))
  
}

adjustVarLoc <- function(var_loc, ens_vs_var_transcr.cds)
{
  #adjust location based on insertions or deletions in supplied var tx relative to ens tx
  
  #this function works for simple cases, but fails when there are lots of indels
  
  #ens is pattern, var is subj
  
  
  ins <- insertion(indel(ens_vs_var_transcr.cds))
  dels <- deletion(indel(ens_vs_var_transcr.cds))
  
  #insertions in ens tx suggest ens tx is longer than var tx
  #for every insertion before variant location, add with of insertion from variant location
  #this should give corrected variant location on ens tx
  for (i in 1:length(ins))
  {
    # print(i)
    if(length(ins[[i]]) > 0)
    {
      for (k in 1:length(ins[[i]]))
      {
        #  print(k)
        if(!is.na(start(ins[[i]])[k]))
        {
          #print(c(i,k))
          if (var_loc > start(ins[[i]][k]))
          {
            var_loc <- var_loc + width(ins[[i]][k])
            #       print(var_loc)
          }
        }
      }
    }
  }
  
  #deletions in ens tx suggest ens tx is shorter than var tx
  #for every deletion before variant location, sub width of deletion to variant location
  #this should give corrected variant location on ens tx
  for (i in 1:length(dels))
  {
    # print(i)
    if(length(dels[[i]]) > 0)
    {
      for (k in 1:length(dels[[i]]))
      {
        #   print(k)
        if(!is.na(start(dels[[i]])[k]))
        {
          #     print("a")
          if (var_loc > start(dels[[i]])[k])
          {
            var_loc <- var_loc - width(dels[[i]][k])
            #         print(var_loc)
          }
        }
      }
    }
  }
  #ir1
  #v01
  #  print(paste0("var loc is ", var_loc))
  return(var_loc)
  
}



adjustVarLoc2 <- function(var_loc, ens_vs_var_transcr_cds) {
  # Extract indels
  ins <- insertion(indel(ens_vs_var_transcr_cds))
  dels <- deletion(indel(ens_vs_var_transcr_cds))
  
  # Initialize a data frame to store indels with type
  indel_df <- data.frame(
    start = integer(),
    width = integer(),
    type = character()
  )
  
  # Add insertions to the data frame
  if (length(ins) > 0 && length(width(ins[[1]])) > 0) {
    indel_df <- rbind(indel_df, data.frame(
      start = start(ins[[1]]),
      width = width(ins[[1]]),
      type = 'insertion'
    ))
  }
  
  # Add deletions to the data frame
  if (length(dels) > 0 && length(width(dels[[1]])) > 0) {
    indel_df <- rbind(indel_df, data.frame(
      start = start(dels[[1]]),
      width = width(dels[[1]]),
      type = 'deletion'
    ))
  }
  
  # Sort indel_df by start positions
  indel_df <- indel_df[order(indel_df$start), ]
  
  # Adjust var_loc based on sorted indels
  for (i in seq_len(nrow(indel_df))) {
    if (var_loc > indel_df$start[i]) {
      if (indel_df$type[i] == "insertion") {
        var_loc <- var_loc + indel_df$width[i]
      } else if (indel_df$type[i] == "deletion") {
        var_loc <- var_loc - indel_df$width[i]
      }
    }
  }
  
  return(var_loc)
}


infer_cDNA_change <- function(qgene2)
{
  #parse ptn variant
  var_parsed <- parseHGVS(qgene2[[3]])
  var_tx_rs_accnum <- qgene2[[4]]
  
  if(grepl("\\?", qgene2[3])){
    print("false")
    return("none")}
  
  #get cds seq for translation
  var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
  var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
  var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
  var_tx_rs_accnum.cds <- DNAString(paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = ""))
  var_tx_rs_accnum.cds.translation <- translate(var_tx_rs_accnum.cds)
  
  
  #confirm position of variant matches translation
  #convert 3-letter code to 1-letter code
  refAA <- as.character(convert_aa_3to1(var_parsed$ancestral))
  varAA <- as.character(convert_aa_3to1(var_parsed$variant))
  
  # Check if refAA matches translated CDS at variant position
  if (as.character(var_tx_rs_accnum.cds.translation[var_parsed$start]) != refAA) {
    print("Ref AA does not match AA in translated CDS at variant position.")
    return(NULL)
  }
  
  # Extract codon at position
  codon_start <- 1 + 3 * (var_parsed$start - 1)
  ref_codon <- as.character(substr(var_tx_rs_accnum.cds, codon_start, codon_start + 2))
  
  # Use seqinr::translate() to check all possible single nucleotide mutations
  possible_mutations <- c()
  for (pos in 1:3) {
    for (base in c("A", "C", "G", "T")) {
      if (substr(ref_codon, pos, pos) != base) {
        mutated_codon <- paste0(substr(ref_codon, 1, pos - 1), base, substr(ref_codon, pos + 1, 3))
         #print(mutated_codon)
         #print(as.character(translate(DNAString(mutated_codon),no.init.codon = TRUE )))
        if (as.character(translate(DNAString(mutated_codon), no.init.codon= TRUE)) == varAA) {
          cDNA_position <- codon_start + pos - 1
          possible_mutations <- c(possible_mutations, paste0("c.", cDNA_position, substr(ref_codon, pos, pos), ">", base))
          #print(possible_mutations)
        }
      }
    }
  }
  
  # Output proposed cDNA changes
  if (length(possible_mutations) == 0) {
    print("No possible cDNA change could be found that results in the specified AA change.")
  } else {
    print(paste0("Possible cDNA changes: ", toString(possible_mutations)))
  }
  
  return(possible_mutations)
  
}

# Create a mapping from three-letter to one-letter amino acid codes
aa_mapping <- c(
  Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
  Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
  Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
  Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
  Ter = "*", "*" = "*"
)

# Function to convert a vector of three-letter codes to one-letter codes
convert_aa_3to1 <- function(three_letter_vector) {
  one_letter_vector <- aa_mapping[three_letter_vector]
  return(one_letter_vector)
}


#annotate PM3 and PP1 based on ClinVar submission comments
Annotate_CV_PM3_PP1 <- function(merged_dbNERG_clinvar_test, clinvar_ids_to_test =NULL, existing_tbl = NULL, output_CVtable_PM3_PP1_path = "CV_PM3_PM1.txt")
{
  print("test")
  
  if (!is.null(existing_tbl)) {
    merged_dbNERG_clinvar_test <- existing_tbl
    print("Using supplied table")
  }
  
  if (is.null(clinvar_ids_to_test)) {
    clinvar_ids_to_test <- unique(merged_dbNERG_clinvar_test[!is.na(INFO_dbNSFP_clinvar_id), INFO_dbNSFP_clinvar_id])
    print("Extracting all ClinVar IDs")
  }
  
  if (length(clinvar_ids_to_test) == 0) {
    print("No ClinVar IDs to process")
    return(NULL)
  }
  
  len <- length(clinvar_ids_to_test)
  i <- 1
  for (cid in clinvar_ids_to_test)
  {
    print(paste("processing ", i, " of ", len))
    i <- i+1
    
    if(!(merged_dbNERG_clinvar_test[INFO_dbNSFP_clinvar_id == 2095068, clinvar_comments] == ""))
    {
      print(paste0("Clinvar ID ", cid, "already processed"))
      next
    }
    print(paste0("working on clinvar ID ", cid))
    clinvar_search <- entrez_search("clinvar",  cid, retmax = 186)
    clinvar_summary <- entrez_summary(db = "clinvar", id = clinvar_search$ids[1])
    variant_accession <- clinvar_summary[3]
    
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&id=", variant_accession)
    
    
    response <- GET(url)
    xml_data <- content(response, "text")
    xml_file <- read_xml(xml_data)
    
    comment_nodes <- xml_find_all(xml_file, ".//Comment")
    comments <- xml_text(comment_nodes)
    
    for (cmt in comments)
    {
      #print(cmt)
      merged_dbNERG_clinvar_test[INFO_dbNSFP_clinvar_id == cid,clinvar_comments := paste0(merged_dbNERG_clinvar_test[INFO_dbNSFP_clinvar_id == cid,clinvar_comments],";" ,cmt)]
      if(grepl("the data is consistent with being in trans", cmt))
      {
        #assign PM3 clinvar moderate at minimum
        merged_dbNERG_clinvar_test[INFO_dbNSFP_clinvar_id == cid,clinvar_pm3_transseg := "Moderate"]
        #print("Moderate")
        #check if "at least one" is present
        if(grepl("In at least one individual the data is consistent with being in trans", cmt))
        {
          #assign PM3 clinvar strong if 2
          merged_dbNERG_clinvar_test[INFO_dbNSFP_clinvar_id == cid,clinvar_pm3_transseg := "Strong"]
          
          #print("strong")
          #assign PM3 clinvar very strong if 4
          
        }
        
      }
      if(grepl("segregate with disease", cmt))
      {
        merged_dbNERG_clinvar_test[INFO_dbNSFP_clinvar_id == cid,clinvar_pp1_diseaseseg := "Supporting"]
      }
      
    }
    merged_dbNERG_clinvar_test[, clinvar_comments := sub(";", "", clinvar_comments, fixed = TRUE)]
    
    # write.table(merged_dbNERG_clinvar_test, file = output_vcf_path, row.names = FALSE, sep = "\t")
    write.table(merged_dbNERG_clinvar_test, file = output_CVtable_PM3_PP1_path, row.names = FALSE, sep = "\t")
    
    
  }
  return(merged_dbNERG_clinvar_test)
}





## determine if variant if P/LP in clinVar
findClinVar.PLP <- function(gene_query)
{
  srch <- entrez_search(db="clinvar", term = gene_query) 
  if(!is_empty(srch$ids))
  {
    for (srchid in srch$ids)
    {
      var_info <- entrez_summary(db="clinvar", id = srchid) 
      
      #check for exact match to gene query as clinvar search is too permissive
      substring_gene_query <- substring(gene_query, 1, regexpr(":", gene_query) - 1)
      
      if (substring(var_info[5], 1, nchar(substring_gene_query)) == substring_gene_query) {
        print(var_info$germline_classification$description)
        if(tolower(var_info$germline_classification$description) %in% tolower(c("Pathogenic", "Likely Pathogenic", "Pathogenic/Likely Pathogenic", "Pathogenic/Likely pathogenic")))
        {
          print("match")
          return(list(clinvar_id = var_info$uid, clinvar_class = var_info$germline_classification$description, variation = test_inferred[i]))
        }
      }
      
      
    }
    
    
  }
}



## annotate PS1 criterion with Strong if match found in clinVar
testPS1 <- function(db)
{
  PTN.changes <- as.data.table(db[,c("Gene_Name", "HGVS.c", "HGVS.p", "Feature_ID", "INFO_PtID")])
  PTN.changes[,SnpEff_Variant := paste(Gene_Name, HGVS.c, HGVS.p, Feature_ID, sep = " ")]
  PTN.changes[,`:=` (PS1 = character(),
                     PS1.clinVarID = character())]
  genevars <- PTN.changes[!is_empty(HGVS.p) & !(is.na(HGVS.p)) & HGVS.p != "NA" & HGVS.p != "" & !grepl("del", HGVS.p) & !grepl("fs", HGVS.p) &!grepl("dup", HGVS.c), SnpEff_Variant]
  
  for (genevar in genevars)
  {
    test_genevar <- splitGeneVariant(genevar)
    print(paste0("testing ", genevar))
    print(test_genevar)
    test_inferred <- infer_cDNA_change(test_genevar)
    
    if(length(test_inferred) > 1)
    {
      for (i in 1:length(test_inferred))
      {
        if(test_inferred[i] != test_genevar[[2]])
        {
          PS.test <- findClinVar.PLP(paste0(test_genevar[[4]], ":", test_inferred[i]))
          print(PS.test[[1]])
          if(!is.null(PS.test))
          {PTN.changes[genevar == SnpEff_Variant, `:=` (PS1 = "Strong", 
                                                       PS1.clinVarID = PS.test[[1]]) ]}
        }
      }
    }
    
  }
  return(PTN.changes)
}

#test PM5 criterion
testPM5 <- function(db, old_table = NULL, output_PTN_table_path = "PTN_PM5.txt")
{
  
  #create PTN.changes table to populate unless it exists already
  if(!is_empty(old_table))
  {PTN.changes <- old_table
  print("using old table")
  }
  else {
    print("creating new table")
    PTN.changes <- as.data.table(db[,c("Gene_Name", "HGVS.c", "HGVS.p", "Feature_ID", "INFO_PtID")])
    PTN.changes[,SnpEff_Variant := paste(Gene_Name, HGVS.c, HGVS.p, Feature_ID, sep = " ")]
    PTN.changes[,`:=` (PM5 = character(),
                       PM5.clinVarID = character())]
    
  }
  
  
  genevars <- PTN.changes[!is_empty(HGVS.p) & !(is.na(HGVS.p)) & HGVS.p != "NA" & HGVS.p != "" & !grepl("del", HGVS.p) & !grepl("fs", HGVS.p) &!grepl("dup", HGVS.c), SnpEff_Variant]
  
  for (genevar in genevars)
  {
    test_genevar <- splitGeneVariant(genevar)
    print(paste0("testing ", genevar))
    #print(test_genevar)
    
    
    
    PM5.test <- infer_cDNA_change.PM5(test_genevar)
    print(PM5.test)
    
    if(!is.null(PM5.test)){
      print(PM5.test[[1]])
      PTN.changes[genevar == SnpEff_Variant, `:=` (PM5 = "Moderate",
                                                   PM5.clinVarID = PM5.test[[1]]) ]
    }
    
    
    
  }
  
  write.table(PTN.changes, file = output_PTN_table_path, row.names = FALSE, sep = "\t")
  
  return(PTN.changes)
}


#infer all possible cDNA changes that can result in diff amino acids and test them
infer_cDNA_change.PM5 <- function(qgene2)
{
  #parse ptn variant
  var_parsed <- parseHGVS(qgene2[[3]])
  var_tx_rs_accnum <- qgene2[[4]]
  
  if(grepl("\\?", qgene2[3])){
    print("false")
    return("none")}
  
  #get cds seq for translation
  var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
  var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
  var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
  var_tx_rs_accnum.cds <- DNAString(paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = ""))
  var_tx_rs_accnum.cds.translation <- translate(var_tx_rs_accnum.cds)
  
  
  #confirm position of variant matches translation
  #convert 3-letter code to 1-letter code
  refAA <- as.character(convert_aa_3to1(var_parsed$ancestral))
  varAA <- as.character(convert_aa_3to1(var_parsed$variant))
  
  # Check if refAA matches translated CDS at variant position
  if (as.character(var_tx_rs_accnum.cds.translation[var_parsed$start]) != refAA) {
    print("Ref AA does not match AA in translated CDS at variant position.")
    return(NULL)
  }
  
  # Extract codon at position
  codon_start <- 1 + 3 * (var_parsed$start - 1)
  ref_codon <- as.character(substr(var_tx_rs_accnum.cds, codon_start, codon_start + 2))
  
  # Use seqinr::translate() to check all possible single nucleotide mutations
  possible_mutations <- c()
  for (pos in 1:3) {
    for (base in c("A", "C", "G", "T")) {
      if (substr(ref_codon, pos, pos) != base) {
        mutated_codon <- paste0(substr(ref_codon, 1, pos - 1), base, substr(ref_codon, pos + 1, 3))
        #print(mutated_codon)
        #print(as.character(translate(DNAString(mutated_codon),no.init.codon = TRUE )))
        if (as.character(translate(DNAString(mutated_codon), no.init.codon= TRUE)) != varAA) {
          cDNA_position <- codon_start + pos - 1
          var.hgvs <- paste0(var_tx_rs_accnum, ":", paste0("c.", cDNA_position, substr(ref_codon, pos, pos), ">", base))
          print(var.hgvs)
          PM5.test <- findClinVar.PLP(var.hgvs)
          if(!is.null(PM5.test))
          {
            return(PM5.test)
          }
          
          
          #possible_mutations <- c(possible_mutations, paste0("c.", cDNA_position, substr(ref_codon, pos, pos), ">", base))
          #print(possible_mutations)
        }
      }
    }
  }
  return(NULL)
}

#getClinVar comments and other metadata for clinvar entries
getClinVarComments <- function(query_term)
{
  clinvar_search <- entrez_search("clinvar",  query_term, retmax = 186)
  
  clinvar_ids <- character()
  clinvar_interp <- character()
  clinvar_comments <- character()
  clinvar_varloc_assemblyName <- character()
  clinvar_varloc_chr <- character()
  clinvar_varloc_start <- character()
  clinvar_varloc_stop <- character()
  
  i = 1
  for (srchid in clinvar_search$ids)
  {
    #clinvar_summary <- entrez_summary(db = "clinvar", id = srchid)  
    var_info <- entrez_summary(db="clinvar", id = srchid) 
    
    
    variant_accession <- var_info[3]
    url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&id=", variant_accession)
    
    response <- GET(url)
    xml_data <- content(response, "text")
    xml_file <- read_xml(xml_data)
    
    comment_nodes <- xml_find_all(xml_file, ".//Comment")
    
    clinvar_ids[i] <- srchid
    clinvar_interp[i] <- var_info$germline_classification$description
    if(!is_empty( xml_text(comment_nodes)))
    {clinvar_comments[i] <- xml_text(comment_nodes)}
    else{clinvar_comments[i] <- ""}
    clinvar_varloc_assemblyName[i] <- var_info$variation_set$variation_loc[[1]]$assembly_name
    clinvar_varloc_chr[i] <- var_info$variation_set$variation_loc[[1]]$chr
    clinvar_varloc_start[i] <- var_info$variation_set$variation_loc[[1]]$start
    clinvar_varloc_stop[i] <- var_info$variation_set$variation_loc[[1]]$stop
    
    print(i)
    i<-i+1
  }
  
  return(list(clinvar_ids, clinvar_interp, clinvar_comments, clinvar_varloc_assemblyName, clinvar_varloc_chr, clinvar_varloc_start, clinvar_varloc_stop))
}


#extract exon numbers
extract_exon_numbers <- function(input_string) {
  # Split the string on the dash character
  parts <- strsplit(input_string, "-")[[1]]
  
  # Initialize variables for first and last number
  start_exon <- as.numeric(parts[1])
  end_exon <- NA
  
  # Check if there is a second part and assign it to last_number if it exists
  if (length(parts) > 1) {
    end_exon <- as.numeric(parts[2])
  } 
  
  # Return a named vector with the first and last number
  return(c(Start_exon = start_exon, End_exon = end_exon))
}


#parse complex structural variants with or without uncertain breakpoints
parse_hgvs_universal <- function(hgvs_string) {
  # Define the regex pattern with proper groups for each part
  pattern <- "c\\.\\((-?\\d*|\\?)\\+?(-?\\d*?)_(-?\\d*|\\?)\\-?(-?\\d*?)\\)_\\((-?\\d*|\\*1|\\?)\\+?(-?\\d*?)_(-?\\d*|\\?)\\-?(-?\\d*?)\\)(del|dup|ins)"
  
  # Check if the pattern matches the string
  matches <- regmatches(hgvs_string, regexec(pattern, hgvs_string))
  match <- matches[[1]]
  
  if(is_empty(match))
  {
    print(paste0("Cannot parse variant to cDNA: ", hgvs_string))
    return(NULL)
  }
  
  # Verify a match occurred and capture groups are correct
  if (length(match) == 10) {
    # Helper function to interpret components
    interpret_component <- function(x) {
      if (x == "?" || x == "") {
        "uncertain"
      } else if (x == "*1") {
        "3UTR_start"
      } else {
        as.integer(x)  # Direct conversion to integer
      }
    }
    
    # Extract and interpret each component
    downstream_start_exon <- interpret_component(match[2])
    downstream_start_exon_offset <- interpret_component(match[3])
    upstream_start_exon <- interpret_component(match[4])
    upstream_start_exon_offset <- interpret_component(match[5])
    downstream_stop_exon <- interpret_component(match[6])
    downstream_stop_exon_offset <- interpret_component(match[7])
    upstream_stop_exon <- interpret_component(match[8])
    upstream_stop_exon_offset <- interpret_component(match[9])
    variant_type <- tolower(match[10])  # Convert to lowercase
    
    # Compile all extracted data into a list
    result <- list(
      Type = 'cDNA',
      Downstream_Start_Exon = downstream_start_exon,
      Downstream_Start_Exon_Offset = downstream_start_exon_offset,
      Upstream_Start_Exon = upstream_start_exon,
      Upstream_Start_Exon_Offset = upstream_start_exon_offset,
      Downstream_Stop_Exon = downstream_stop_exon,
      Downstream_Stop_Exon_Offset = downstream_stop_exon_offset,
      Upstream_Stop_Exon = upstream_stop_exon,
      Upstream_Stop_Exon_Offset = upstream_stop_exon_offset,
      Variant_Type = variant_type
    )
    
    return(result)
  } else {
    return("Invalid HGVS format or the pattern does not match.")
  }
}

#get structural variant coordinates from a structural variant  no in clinvar
#returns chr start end
getSVcoords <- function(genevar)
{
  genevar.split <- splitGeneVariant(genevar)
  
  #get gene info
  geneInfo <- retrieveGeneInfo(genevar.split,ensembl)
  exons.ens <- as.data.table(geneInfo$exons)
  transcripts.ens <- as.data.table(geneInfo$transcripts)
  
  exon_SV_start_end <- extract_exon_numbers(genevar.split[3])
  if(is.atomic(exon_SV_start_end))
  {
    SV_start_exon <- exon_SV_start_end[1]
    SV_end_exon <- exon_SV_start_end[1]
  }
  else
  {
    SV_start_exon <- exon_SV_start_end$Start_exon
    SV_end_exon <- exon_SV_start_end$End_exon
  }
  

  #use canonical transcript (for now) - may need to handle exceptions
  transcript_ref <- transcripts.ens[transcript_is_canonical == 1, ensembl_transcript_id]
  exons_ref <- exons.ens[ensembl_transcript_id == transcript_ref,]
  chr <- transcripts.ens[ensembl_transcript_id == transcript_ref, chromosome_name]

  
  #estimate chrom start and end based on transcript data
  if(!is.na(exon_SV_start_end[1]))
  {exon_start_coord <- exons_ref[rank == exon_SV_start_end[1], exon_chrom_start]}
  
  if(!is.na(exon_SV_start_end[2]))
  {exon_end_coord <- exons_ref[rank == exon_SV_start_end[2], exon_chrom_end]}
  else{
    exon_end_coord <- exons_ref[rank == exon_SV_start_end[1], exon_chrom_end]
    SV_end_exon <- SV_start_exon
  }
  
  result <- list(
    chr = chr,
    SV_start_coord = exon_start_coord,
    SV_exon_start = SV_start_exon,
    SV_end_coord = exon_end_coord,
    SV_exon_end = SV_end_exon,
    Var_Type = genevar.split[4]
  )
  
  return(result)
}

#given exon list of a transcript, derive translation start site
findTranslationStartLoc <- function(exons_ref){
  for (i in 1:dim(exons_ref)[1])
  {
    if (!is.na(exons_ref[rank == i, cds_start]))
    {
      if(exons_ref[1,strand] == -1){
        TSS_loc <- exons_ref[i,exon_chrom_start] + exons_ref[i,cds_end] 
        return(c(TSS_loc, i))
      }
      #double-check this
      else if(exons_ref[1,strand] == 1){
        TSS_loc <- exons_ref[i,exon_chrom_end] - exons_ref[i,cds_end] 
        return(c(TSS_loc, i))
      }
    }
  }
}

#functions for parsing annotated VCF after SnpEff, SnpSift, and SpliceAI
# Function to parse VCF header and extract annotations
extract_info_annotations <- function(file_path) {
  lines <- read_lines(file_path)
  info_lines <- lines[str_detect(lines, "^##INFO=")]
  annotations <- str_extract(info_lines, "(?<=ID=)[^,]+")
  return(annotations)
}

# Function to parse the INFO column into separate annotation columns
parse_info_column <- function(vcf_data, annotations) {
  # Expand INFO column into multiple rows per annotation
  # Add a row identifier to keep track of original rows
  vcf_data <- vcf_data %>%
    mutate(row_id = row_number())
  
  # Expand INFO column into multiple rows per annotation
  info_expanded <- vcf_data %>%
    mutate(INFO = str_split(INFO, ";")) %>%
    unnest_longer(INFO) %>%
    separate(INFO, into = c("Annotation", "Value"), sep = "=", fill = "right", extra = "merge") %>%
    pivot_wider(names_from = Annotation, values_from = Value, values_fn = list(Value = dplyr::first), names_prefix = "INFO_", id_cols = "row_id")
  
  # Ensure all expected annotations have columns, even if they are missing
  missing_annotations <- setdiff(paste0("INFO_", annotations), names(info_expanded))
  for (col in missing_annotations) {
    info_expanded[[col]] <- NA
  }
  
  # Join back with original data to ensure no loss of rows
  # Join back with original data to ensure no loss of rows and all columns are retained
  vcf_data <- left_join(dplyr::select(vcf_data, -INFO), info_expanded, by = "row_id") %>%
    dplyr::select(-row_id)
  
  return(vcf_data)
}

# Main function to parse the VCF file
parse_vcf <- function(vcf_file_path) {
  # Extract annotations from the VCF file header
  annotations <- extract_info_annotations(vcf_file_path)
  
  # Read the VCF data, excluding the metadata header
  vcf_data <- read_tsv(vcf_file_path, comment = "##", col_names = TRUE, col_types = cols(.default = "c"))
  if ("#CHROM" %in% names(vcf_data)) {
    names(vcf_data)[which(names(vcf_data) == "#CHROM")] <- "CHROM"
  }
  
  # Parse the INFO column into separate annotation columns
  vcf_parsed <- parse_info_column(vcf_data, annotations)
  
  return(vcf_parsed)
}

# Adjusted function to extract sub-annotation titles from the VCF header
extract_ann_titles <- function(vcf_file_path) {
  lines <- read_lines(vcf_file_path)
  ann_line <- lines[str_detect(lines, "^##INFO=<ID=ANN")]
  # Extract the part of the line after "Functional annotations:"
  ann_titles_str <- str_extract(ann_line, "(?<=Functional annotations: ).*")
  # Split the titles into a vector and trim spaces and remove quotes
  ann_titles <- str_split(ann_titles_str, "\\|", simplify = FALSE) %>%
    unlist() %>%
    map_chr(~str_trim(.)) %>%
    map_chr(~str_remove_all(., "'"))
  return(ann_titles)
}



# Function to parse the INFO_ANN column into individual sub-annotation columns
# Adjusted function to parse the INFO_ANN column with modifications
# Function to parse INFO_ANN column
parse_info_ann <- function(vcf_data, info_ann_titles) {
  # Ensure that info_ann_titles has exactly 16 elements
  if (length(info_ann_titles) != 16) {
    stop("info_ann_titles must contain exactly 16 titles for sub-annotations.")
  }
  
  vcf_data %>%
    # Generate a row identifier to keep track of original rows
    mutate(row_id = row_number()) %>%
    # Expand INFO_ANN into rows
    separate_rows(INFO_ANN, sep = ",") %>%
    # Separate the sub-annotations into individual columns
    separate(INFO_ANN, into = info_ann_titles, sep = "\\|", fill = "right") %>%
    # Optionally, handle conversion of empty strings to NA or other preprocessing
    mutate(across(.cols = starts_with("ANN"), .fns = ~na_if(.x, ""))) %>%
    # Remove rows where all annotations are NA (if needed)
    # filter_if(.predicate = starts_with("ANN"), .f = any_vars(!is.na(.)))
    dplyr::select(-row_id) # Optionally remove the row identifier if it's no longer needed
}

# Adjusted function to extract sub-annotation titles from the VCF header
extract_SpliceAI_titles <- function(vcf_file_path) {
  lines <- read_lines(vcf_file_path)
  metadata <- lines[str_detect(lines, "^##INFO=<ID=SpliceAI")]
  # Extract the part of the line after "Functional annotations:"
  titles <- sub('.*Format: (.+?)">', '\\1', metadata)
  #titles <- unlist(strsplit(titles, "\\|"))
  # Split the titles into a vector and trim spaces and remove quotes
  titles <- str_split(titles, "\\|", simplify = FALSE) %>%
    unlist() %>%
    map_chr(~str_trim(.)) %>%
    map_chr(~str_remove_all(., "'"))
  
  titles <- paste("SpliceAI", titles, sep = "_")
  
  return(titles)
}

# Function to parse INFO_SpliceAI column
parse_info_spliceAI <- function(vcf_data, info_spliceAI_titles) {
  if (length(info_spliceAI_titles) != 10) {  # Ensure there are exactly 10 titles for sub-annotations
    stop("info_spliceAI_titles must contain exactly 10 titles for sub-annotations.")
  }
  
  vcf_data %>%
    # Generate a row identifier to keep track of original rows
    mutate(row_id = row_number()) %>%
    # Expand INFO_SpliceAI into rows by separating on commas
    separate_rows(INFO_SpliceAI, sep = ",") %>%
    # Separate the sub-annotations into individual columns
    separate(INFO_SpliceAI, into = info_spliceAI_titles, sep = "\\|", fill = "right") %>%
    # Optionally, handle conversion of empty strings to NA or other preprocessing
    mutate(across(starts_with("SpliceAI"), ~na_if(.x, ""))) %>%
    # Remove the row identifier if it's no longer needed
    dplyr::select(-row_id)
}

#function to interpret all leftover variants - usually for deletions or very copmlex cdna hgvs strings
Annotate_Unannotated <- function(saved_vcf_file){
  
  ##create DT to export
  melted_df_anno <- as.data.table(saved_vcf_file)
  setDT(melted_df_anno)[, POS := as.character(POS)]
  melted_df_anno[, `:=`(
    GeneName = character(),
    clinvar_ids = character(),
    clinvar_interp = character(),
    clinvar_comments = character(),
    clinvar_varloc_assemblyName = character(),
    clinvar_varloc_chr = character(),
    clinvar_varloc_start = character(),
    clinvar_varloc_stop = character(),
    infer_SVloc_chr = character(),
    infer_SVloc_start = character(),
    infer_SVloc_end = character(),
    infer_SVloc_exon_start = character(),
    infer_SVloc_exon_end = character(),
    infer_SVloc_type = character()
  )]
  
  melted_df_anno <- melted_df_anno[is.na(CHROM),]
  
  #VCF_to_paste <- init_VCF_paste()
  #merged_df_annoVCF <- cbind(melted_df_anno, VCF_to_paste)
  
  genevars <- melted_df_anno[,Gene_Variant]
  
  i <- 1
  for (genevar in genevars)
  {
    print(i)
    i <- i+1
    genevar.split <- splitGeneVariant(genevar)
    melted_df_anno[Gene_Variant == genevar, ]
    
    #first search clinvar
    query_term <- paste(genevar.split[1:4], sep = " ", collapse = " ")
    clinvar_query <- getClinVarComments(query_term)
    
    #clinvar returns something
    if(!is_empty(clinvar_query[[1]]))  {
      melted_df_anno[Gene_Variant == genevar, `:=`(clinvar_ids = paste(clinvar_query[[1]], sep = ", ", collapse = ", "),
                                                   clinvar_interp = paste(clinvar_query[[2]], sep = ", ", collapse = ", "),
                                                   clinvar_comments = paste(clinvar_query[[3]], sep = ", ", collapse = ", "),
                                                   clinvar_varloc_assemblyName = paste(clinvar_query[[4]], sep = ", ", collapse = ", "),
                                                   clinvar_varloc_chr = paste(clinvar_query[[5]], sep = ", ", collapse = ", "),
                                                   clinvar_varloc_start = paste(clinvar_query[[6]], sep = ", ", collapse = ", "),
                                                   clinvar_varloc_stop = paste(clinvar_query[[7]], sep = ", ", collapse = ", "),
                                                   CHROM = clinvar_query[[5]][1],
                                                   POS = clinvar_query[[6]][1])]
    }
    
    #does not exist in clinvar and labeled as deletion or duplication
    else if(is_empty(clinvar_query[[1]]) & genevar.split[4] %in% c("deletion", "duplication")){
      print(paste0("manually inferring coordinates for ", genevar))
      SV_coords <- getSVcoords(genevar)
      melted_df_anno[Gene_Variant == genevar, `:=`(infer_SVloc_chr = SV_coords$chr,
                                                   infer_SVloc_start = SV_coords$SV_start_coord,
                                                   infer_SVloc_end = SV_coords$SV_end_coord,
                                                   infer_SVloc_exon_start = SV_coords$SV_exon_start,
                                                   infer_SVloc_exon_end = SV_coords$SV_exon_end,
                                                   infer_SVloc_type = SV_coords$Var_Type,
                                                   CHROM = SV_coords$chr,
                                                   POS = SV_coords$SV_start_coord)]
      
    }
    
    #does not return clinvar and is not labeled as deletion with possible HGVS exon boundaries
    #i.e. deletion given in c.(A+B_C-D)_(E+F_G_H)del|dup|ins
    else if(is_empty(clinvar_query[[1]]) & (genevar.split[2] != "NA") & (!genevar.split[4] %in% c("deletion", "duplication")) ){
      print(paste0("manually inferring coordinates for SV in HGVS format for ", genevar))
      SV_coords <- getSVcoords.hgvs(genevar)
      #if variant cannot be parsed, insert dummy variables (to be removed later)
      if(is.null(SV_coords))
      {
        melted_df_anno[Gene_Variant == genevar, `:=`(
          CHROM = paste0("DUMMYCHR", i),
          POS = paste0("DUMMYPOS", i))]
        next
      }
      melted_df_anno[Gene_Variant == genevar, `:=`(infer_SVloc_chr = SV_coords$chr,
                                                   infer_SVloc_start = SV_coords$SV_start_coord,
                                                   infer_SVloc_end = SV_coords$SV_end_coord,
                                                   infer_SVloc_exon_start = SV_coords$SV_exon_start,
                                                   infer_SVloc_exon_end = SV_coords$SV_exon_end,
                                                   infer_SVloc_type = SV_coords$Var_Type,
                                                   CHROM = SV_coords$chr,
                                                   POS = SV_coords$SV_start_coord)]
    }
    
    #other variants do not apply
    else{
      melted_df_anno[Gene_Variant == genevar, `:=`(
        CHROM = paste0("DUMMYCHR", i),
        POS = paste0("DUMMYPOS", i))]
      
      next
      
    }
    
    
    #if no clinvar
    #if(is.na(merged_df_annoVCF[Gene_Variant == genevar, clinvar_ids]))
    
    
    
    #print(var_info$variation_set$variation_loc)
    #print(var_info$germline_classification$description)
  }
  
  return(melted_df_anno)
  
}

#infer SV coordinates from HGVS string
getSVcoords.hgvs <- function(genevar)
{
  #parse variant
  genevar.split <- splitGeneVariant(genevar)
  genevar.parsed <- parse_hgvs_universal(genevar.split[2])
  if(is.null(genevar.parsed))
  {return(NULL)}
  
  #get gene info
  geneInfo <- retrieveGeneInfo(genevar.split,ensembl)
  exons.ens <- as.data.table(geneInfo$exons)
  transcripts.ens <- as.data.table(geneInfo$transcripts)
  
  #try to use the correct transcript acccession, get chr of transcript
  if(genevar.split[4] != "NA"){
    accession <- sub("\\..*$", "", genevar.split[4])
    transcript_ref <- transcripts.ens[transcript_ccds_cap == accession, ensembl_transcript_id]
    exons_ref <- exons.ens[ensembl_transcript_id == transcript_ref,]
    chr <- transcripts.ens[ensembl_transcript_id == transcript_ref, chromosome_name]
  }
  #use canonical transcript if no transcript accession supplied
  else if (genevar.split[4] == "NA"){
    transcript_ref <- transcripts.ens[transcript_is_canonical == 1, ensembl_transcript_id]
    exons_ref <- exons.ens[ensembl_transcript_id == transcript_ref,]
    chr <- transcripts.ens[ensembl_transcript_id == transcript_ref, chromosome_name]
  }
  
  #define cds start based on cds_start and cds_end; returns TSS and exon rank
  translation_start_loc <- findTranslationStartLoc(exons_ref)
  
  #derive deletion start
  SV_start <- findSVStart(exons_ref, translation_start_loc, genevar.parsed)
  
  #derive deletion end
  SV_end <- findSVEnd(exons_ref, translation_start_loc, genevar.parsed)
  
  #parse variant type
  if(genevar.parsed$Variant_Type == "del"){variation_type <- "deletion"}
  if(genevar.parsed$Variant_Type == "dup"){variation_type <- "duplication"}
  if(genevar.parsed$Variant_Type == "ins"){variation_type <- "insertion"}
  
  result <- list(
    chr = chr,
    SV_start_coord = SV_start$SV_Start,
    SV_exon_start = SV_start$Exon_start,
    SV_end_coord = SV_end$SV_End,
    SV_exon_end = SV_end$Exon_End,
    Var_Type = variation_type
  )
  
  
  return(result)
}

#Derive SV end position and exon
findSVEnd <- function(exons_ref, translation_start_loc, genevar.parsed)
{
  if(genevar.parsed$Downstream_Stop_Exon != "3UTR_start"){
    for (i in translation_start_loc[2]:dim(exons_ref)[1])
    {
      #print(i)
      if (genevar.parsed$Downstream_Stop_Exon == exons_ref[rank == i, cds_end])
      {
        #print(exons_ref[rank == i, cds_end])
        return(list(SV_End = exons_ref[rank == i, exon_chrom_end], Exon_End = i))
      }
    }
  }
  else if(genevar.parsed$Downstream_Stop_Exon == "3UTR_start"){
    return(list(SV_End = exons_ref[rank == dim(exons_ref)[1], exon_chrom_end], Exon_End = exons_ref[rank == dim(exons_ref)[1], rank]))
  }
}


#Derive SV start position and exon
findSVStart <- function(exons_ref, translation_start_loc, genevar.parsed)
{
  #if start is before TSS
  if(exons_ref[1,strand] == -1 & genevar.parsed$Upstream_Start_Exon < 0){
    #
    SV_start <- translation_start_loc[1] - genevar.parsed$Upstream_Start_Exon - 1
    if(SV_start == exons_ref[rank == translation_start_loc[2], exon_chrom_end])
    {
      return(list(SV_Start = SV_start, Exon_start = translation_start_loc[2]))
    }
  }
  # if start exon is after TSS
  else if(exons_ref[1,strand] == -1 & genevar.parsed$Upstream_Start_Exon > 0){
    SV_start <- exons_ref[cds_start == genevar.parsed$Upstream_Start_Exon, exon_chrom_start]
    exon_rank <- exons_ref[cds_start == genevar.parsed$Upstream_Start_Exon, rank]
    return(list(SV_Start = SV_start, Exon_start = translation_start_loc[2]))
  }
  #same for antisense strand
  #if start is before TSS
  else if(exons_ref[1,strand] == 1 & genevar.parsed$Upstream_Start_Exon < 0){
    #
    SV_start <- translation_start_loc[1] + genevar.parsed$Upstream_Start_Exon - 1
    if(SV_start == exons_ref[rank == translation_start_loc[2], exon_chrom_end])
    {
      return(list(SV_Start = SV_start, Exon_start = translation_start_loc[2]))
    }
  }
  #doublecheck this
  else if(exons_ref[1,strand] == 1 & genevar.parsed$Upstream_Start_Exon > 0){
    #SV_start <- translation_start_loc[1] + genevar.parsed$Upstream_Start_Exon + 1
    SV_start <- exons_ref[cds_start == genevar.parsed$Upstream_Start_Exon, exon_chrom_start]
    exon_rank <- exons_ref[cds_start == genevar.parsed$Upstream_Start_Exon, rank]
    
    return(list(SV_Start = SV_start, Exon_start = exon_rank))
    
  }
}


# Main function
processClinVar <- function(genevars, ensembl = ensembl, meltedDF = melted_df, mdf2_dump, output_file, old_table = NULL) {
  # initiate clinvar_to_paste
  clinvar_to_paste <- init_clinvar_paste()
  
  # create melted_df to populate unless it exists already
  if(!is_empty(old_table))
  {melted_df_merge2 <- old_table
  print("using old table")
  }
  else {
    print("creating new table")
    melted_df_merge2 <- cbind(meltedDF, clinvar_to_paste)
  }
  
  # main processing loop
  for (qgene in genevars) {
    print(paste0("Processing Variant ", qgene))
    
    if(!is.na(any(melted_df_merge2[Gene_Variant == qgene, clinvar_uid])))
    {
      print("variant already processed")
      next
    }
    
    # query clinVar - return summary object with gene
    result <- findMatchClinvar(qgene, mdf2_dump, ensembl)
    clinvar_summary <- result$var_info
    mdf2.dump <- result$mdf2_dump
    
    # populate info for pasting
    if(!is_empty(clinvar_summary)) {
      clinvar_to_paste <- processVariantClinvar(clinvar_summary)
    } else {
      next
    }
    
    # merge
    if(!is_empty(clinvar_to_paste)) {
      melted_df_merge2[Gene_Variant == qgene, names(clinvar_to_paste) := clinvar_to_paste]
    } else {
      next
    }
    
    # Save the current state of melted_df_merge2 to a file
    write.table(melted_df_merge2, file = output_file, row.names = FALSE, sep = "\t")
    # or use saveRDS(melted_df_merge2, file = output_file) for an R data format
  }
  
  return(list(melted_df_merge2 = melted_df_merge2, mdf2_dump = mdf2.dump))
}


# Main function - create VCF entries in metled table that can then be extracted for analysis as VCF file
CreateVCF <- function(genevars, ensembl = ensembl, meltedDF = melted_df, mdf2_dump, output_file, old_table = NULL) {
  
  # initiate clinvar_to_paste
  VCF_to_paste <- init_VCF_paste()
  
  #create melted_df to populate unless it exists already
  if(!is_empty(old_table))
  {melted_df_merge2 <- old_table
  print("using old table")
  }
  else {
    print("creating new table")
    melted_df_merge2 <- cbind(melted_df, VCF_to_paste)
  }
  
  
  
  
  #main processing loop
  for (qgene in genevars) {
    print(paste0("Processing Variant ", qgene))
    VCF_to_paste <- init_VCF_paste() # reinitialize empty VCF to paste to prevent accidental copying of old VCF data
    
    
    if(any(!is.na(melted_df_merge2[Gene_Variant == qgene, CHROM])))
    {
      print("variant already processed")
      next
    }
    
    # build VCF entry - one per variant
    result <- Build_VCF_entry(qgene, mdf2_dump, ensembl)
    VCF_to_paste <- result
    #VCF_summary <- result$var_info when get mdf2 dump working
    #mdf2.dump <- result$mdf2_dump
    
    # populate info for pasting
    # if(!is_empty(VCF_summary)) {
    #   VCF_to_paste <- processVariantVCF(VCF_summary)
    # } else {
    #   next
    # }
    
    # merge
    if(!is_empty(VCF_to_paste)) {
      print(VCF_to_paste)
      melted_df_merge2[Gene_Variant == qgene, names(VCF_to_paste) := VCF_to_paste]
    } else {
      print("cannot merge")
      next
    }
    
    # Save the current state of melted_df_merge2 to a file
    write.table(melted_df_merge2, file = output_vcf_path, row.names = FALSE, sep = "\t")
    # or use saveRDS(melted_df_merge2, file = output_file) for an R data format
  }
  
  return(list(melted_df_merge2 = melted_df_merge2, mdf2_dump = mdf2.dump))
}


# make VCF entry to add to file
# parses cDNA variant type
# calls different processing functions depending on variant type
Build_VCF_entry <- function(qgene, mdf2_dump, ensembl)
{
  print(paste0("building ", qgene))
  #####
  
  qgene2 <- splitGeneVariant(qgene) # 1 is hgnc, 2 is cdna, 3 is ptn
  if(qgene2[2] == "NA" & qgene2[3] == "NA") {
    print("no cDNA or Protein var")
    return(NULL)} #stop if no cDNA variant present
  
  if(grepl(pattern = "^IVS", qgene2[[2]]))
  {
    var_parsed <- parseIVS(qgene2[[2]])
  }
  
  if(qgene2[2] == "NA" & qgene2[3] != "NA") {
    print("Protein var only, imputing cDNA change")
    inferred_cDNA <- infer_cDNA_change(qgene2)
    var_parsed <- parseHGVS(inferred_cDNA[1])} #stop if no cDNA v
  
  else
  {var_parsed <- parseHGVS(qgene2[[2]])} # parse variant
  
  
  # 
  
  
  # stop if cannot parse variant
  if (var_parsed$type == "invalid" & qgene2[[4]] != "deletion")
  {
    print("cannot parse variant - invalid")
    return(NULL)
  }
  
  if (var_parsed$type == "not_implemented" & qgene2[[4]] != "deletion")
  {
    print("cannot parse variant - not implemented")
    return(NULL)
  }
  
  
  
  # Purpose is to construct vcf entries for all types of variables
  # Retrieve information about gene and all of its exons
  #####
  
  geneInfo <- retrieveGeneInfo(qgene2,ensembl)
  exons.ens <- geneInfo$exons
  transcripts.ens <- geneInfo$transcripts
  
  # get exons for all transcripts with refseq id - first try MANE / canonical
  # if no MANE transcript, then use ccds match
  
  # print("geneInfo")
  # print(geneInfo)
  # print("transcripts.ens")
  # print(transcripts.ens)
  # print("exons.ens")
  # print(exons.ens)
  
  
  transcripts.ens.list <- ExtractMatchedRefseqEnsTranscripts(transcripts.ens)
  
  # Get exons for transcripts that have a refseq ID (from mane or from ccds) as a list
  
  exons.ens.list <- getExonsList(exons.ens, transcripts.ens, transcripts.ens.list)
  
  # print("transcripts.ens.list")
  # print(transcripts.ens.list)
  # 
  # print("exons.ens.list")
  # print(exons.ens.list)
  
  
  # enter VCF for point variant
  if (var_parsed$type  == "substitution" & is.null(var_parsed$startIntron))
  {
    var_vcf_toReturn <- build_VCF_entry_SNV(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  if (var_parsed$type  == "substitution" & !is.null(var_parsed$startIntron))
  {
    var_vcf_toReturn <- build_VCF_entry_SNV_intronic(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  if (var_parsed$type  == "intervening sequence")
  {
    var_vcf_toReturn <- build_VCF_entry_SNV_intronic_IVS(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  else if (var_parsed$type == "singleDuplication")
  {
    var_vcf_toReturn <- build_VCF_entry_singleDup(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  else if (var_parsed$type == "singledeletion")
  {
    var_vcf_toReturn <- build_VCF_entry_singleDel(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  else if (var_parsed$type == "delins")
  {
    var_vcf_toReturn <- build_VCF_entry_delins(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  else if (var_parsed$type == "duplication")
  {
    var_vcf_toReturn <- build_VCF_entry_duplication(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  else if (var_parsed$type == "deletion")
  {
    var_vcf_toReturn <- build_VCF_entry_deletion(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  else if (var_parsed$type == "insertion")
  {
    var_vcf_toReturn <- build_VCF_entry_insertion(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  else if (var_parsed$type == "invalid" & qgene2[[4]] == "deletion")
  {
    var_vcf_toReturn <- build_VCF_entry_SV(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
  }
  
  return(var_info = var_vcf_toReturn)
  
}


# build SNV intronic VCF entry
build_VCF_entry_SNV_intronic <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
{
  print(paste0("inside intronic function for ", qgene2[[2]]))
  
  #if reference transcript is known
  var_tx_rs_accnum <- qgene2[[4]]
  
  #note that for intronic variants, if on antisense / neg strand, must do opposite
  if(exons.ens.list[[1]]$strand[1] == 1)
  {
    var_intron_shift <- var_parsed$startIntron
  }
  else if (exons.ens.list[[1]]$strand[1] == -1)
  {
    var_intron_shift <- -(var_parsed$startIntron)
  }
  
  #####
  # if(any(grepl(var_tx_rs_accnum,transcripts.ens$transcript_ccds_match)))
  # {
  #   #tx ref accnum matches ens version by ccds
  #   # ensembl way
  #   canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == var_tx_rs_accnum), "ensembl_transcript_id"]
  #   
  #   #possible that ccds_match will return multiple ens transcripts. take the one that matches MANE or otherwise the first one
  #   if(length(canon_transcript) > 1 & 
  #      transcripts.ens[which(transcripts.ens$transcript_is_canonical ==1), "transcript_mane_select"] == var_tx_rs_accnum)
  #   {
  #     canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_mane_select == var_tx_rs_accnum), "ensembl_transcript_id"]
  #   }
  #   else
  #   {
  #     canon_transcript <- canon_transcript[1]
  #   }
  #   
  #   genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
  #   print("supplied var ref transcr matches ref transcr in ensembl table; proceed with esnembl way of deriving position")
  #   
  #   
  # }
  # else if(var_tx_rs_accnum != "NA")
  # {
  #   #tx ref accnum does not match ens canon_tx
  #   accession <- sub("\\..*$", "", var_tx_rs_accnum)
  #   pattern <- paste0("^", accession)
  #   ens_vers_transcr <- unique(transcripts.ens[grep(pattern, transcripts.ens$transcript_ccds_cap), "transcript_ccds_match"])
  #   use_canon_tx<-0
  #   
  #   #if there is no match, assign canoncial to ens_ver_transcr
  #   if(is_empty(ens_vers_transcr))
  #   {
  #     ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_mane_select"]
  #     use_canon_tx<-1
  #     print("no match in ensembl, using ensembl canonical transcript")
  #   }
  #   
  #   #pull var_ref_transcr fasta and cds_interval
  #   var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
  #   var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
  #   var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
  #   var_tx_rs_accnum.cds <- paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = "")
  #   
  #   #pull ens_ref_transr fasta and cds interval
  #   ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
  #   ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
  #   ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
  #   ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
  #   
  #   ens_vs_var_transcr.cds <- pairwiseAlignment(
  #     ens_vers_transcr.cds,
  #     var_tx_rs_accnum.cds
  #   )
  #   # check if CDS length equal
  #   
  #   if (var_tx_rs_accnum.cds_len == ens_vers_transcr.cds_len)
  #   {
  #     
  #     #cDNA position on the ens tx
  #     #cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
  #     
  #     #check if percent identity is 100% i.e. CDS are identical between versions
  #     if (pid(ens_vs_var_transcr.cds) == 100 )
  #     {
  #       #standard ensembl way to calc genomic position
  #       canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
  #       if(use_canon_tx == 1)
  #       {
  #         canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
  #       }
  #       genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
  #       
  #     }
  #     else
  #     {
  #       #check if var_parsed$start is not in mismatch table
  #       mmt <- mismatchTable(ens_vs_var_transcr.cds)
  #       if (checkVarVsMismatch(var_parsed, mmt))
  #       {
  #         #standard ensembl way to calc genomic position
  #         canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
  #         genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
  #         
  #       }
  #       else
  #       {
  #         #error out
  #         print("ref seqeunce changed, variant is likely incorrect")
  #         
  #         #dump qgene to some sort of log
  #         return(NULL)
  #         
  #       }
  #     }
  #   }
  #   else 
  #   {
  #     #CDS length is not equal
  #     
  #     canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
  #     var_start_adjusted <- adjustVarLoc(var_parsed$start, ens_vs_var_transcr.cds)
  #     
  #     
  #     genomic_position_ens_result <- calcGenomicPositionEns(var_start_adjusted, canon_transcript)
  #     
  #     #new cDNA position on the ens tx
  #     # cdna_position <- var_loc_adjusted + ens_vers_transcr.cds_int[1,"from"] - 1
  #     
  #     
  #     print("ref CDS sequence length changed, position adjusted for indels between transcripts")  
  #     
  #   }
  #   
  # }
  # else if(var_tx_rs_accnum == "NA")
  # {
  #   #no tx ref accnum supplied, use ens canon_tx
  #   canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
  #   #ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_ccds_match"]
  #   
  #   
  #   genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
  #   
  #   #get tx seq
  #   # ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
  #   # ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
  #   # ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
  #   # ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
  #   # 
  #   
  #   #cDNA position on the ens tx
  #   #cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
  #   
  #   
  #   #manually extract genomic position
  #   genomic_position_man <- calcGenomicPosition(var_loc =  var_parsed$start, 
  #                                                     exons.ens.list =  exons.ens.list, 
  #                                                     transcripts.ens =  transcripts.ens, 
  #                                                     canon_transcript = canon_transcript, 
  #                                                     ensembl = ensembl)
  #   
  #   #check if ensembl and manual extractions agree, then check if ref seq from ens matches transcript
  #   
  #   
  #   if(genomic_position_man == genomic_position_ens && as.data.frame(checkGenRef)$value == var_parsed$ancestral) #matches; likely correct transcript chosen
  #   {
  #     print("No var trascript ref supplied, used canonical transcript - both ens and manual ways align")
  #   }
  #   else
  #   {
  #     print("no match between ens and man")
  #     return(NULL)
  #   }
  # }
  
  # Main function - create VCF entries in metled table that can then be extracted for analysis as VCF file
  CreateVCF <- function(genevars, ensembl = ensembl, meltedDF = melted_df, mdf2_dump, output_file, old_table = NULL) {
    
    # initiate clinvar_to_paste
    VCF_to_paste <- init_VCF_paste()
    
    #create melted_df to populate unless it exists already
    if(!is_empty(old_table))
    {melted_df_merge2 <- old_table
    print("using old table")
    }
    else {
      print("creating new table")
      melted_df_merge2 <- cbind(melted_df, VCF_to_paste)
    }
    
    
    
    
    #main processing loop
    for (qgene in genevars) {
      print(paste0("Processing Variant ", qgene))
      VCF_to_paste <- init_VCF_paste() # reinitialize empty VCF to paste to prevent accidental copying of old VCF data
      
      
      if(any(!is.na(melted_df_merge2[Gene_Variant == qgene, CHROM])))
      {
        print("variant already processed")
        next
      }
      
      # build VCF entry - one per variant
      result <- Build_VCF_entry(qgene, mdf2_dump, ensembl)
      VCF_to_paste <- result
      #VCF_summary <- result$var_info when get mdf2 dump working
      #mdf2.dump <- result$mdf2_dump
      
      # populate info for pasting
      # if(!is_empty(VCF_summary)) {
      #   VCF_to_paste <- processVariantVCF(VCF_summary)
      # } else {
      #   next
      # }
      
      # merge
      if(!is_empty(VCF_to_paste)) {
        print(VCF_to_paste)
        melted_df_merge2[Gene_Variant == qgene, names(VCF_to_paste) := VCF_to_paste]
      } else {
        print("cannot merge")
        next
      }
      
      # Save the current state of melted_df_merge2 to a file
      write.table(melted_df_merge2, file = output_vcf_path, row.names = FALSE, sep = "\t")
      # or use saveRDS(melted_df_merge2, file = output_file) for an R data format
    }
    
    return(list(melted_df_merge2 = melted_df_merge2, mdf2_dump = mdf2.dump))
  }
  
  
  # make VCF entry to add to file
  # parses cDNA variant type
  # calls different processing functions depending on variant type
  Build_VCF_entry <- function(qgene, mdf2_dump, ensembl)
  {
    print(paste0("building ", qgene))
    #####
    
    qgene2 <- splitGeneVariant(qgene) # 1 is hgnc, 2 is cdna, 3 is ptn
    if(qgene2[2] == "NA" & qgene2[3] == "NA") {
      print("no cDNA or Protein var")
      return(NULL)} #stop if no cDNA variant present
    
    if(grepl(pattern = "^IVS", qgene2[[2]]))
    {
      var_parsed <- parseIVS(qgene2[[2]])
    }
    
    if(qgene2[2] == "NA" & qgene2[3] != "NA") {
      print("Protein var only, imputing cDNA change")
      inferred_cDNA <- infer_cDNA_change(qgene2)
      var_parsed <- parseHGVS(inferred_cDNA[1])} #stop if no cDNA v
    
    else
    {var_parsed <- parseHGVS(qgene2[[2]])} # parse variant
    
    
    # 
    
    
    # stop if cannot parse variant
    if (var_parsed$type == "invalid" & qgene2[[4]] != "deletion")
    {
      print("cannot parse variant - invalid")
      return(NULL)
    }
    
    if (var_parsed$type == "not_implemented" & qgene2[[4]] != "deletion")
    {
      print("cannot parse variant - not implemented")
      return(NULL)
    }
    
    
    
    # Purpose is to construct vcf entries for all types of variables
    # Retrieve information about gene and all of its exons
    #####
    
    geneInfo <- retrieveGeneInfo(qgene2,ensembl)
    exons.ens <- geneInfo$exons
    transcripts.ens <- geneInfo$transcripts
    
    # get exons for all transcripts with refseq id - first try MANE / canonical
    # if no MANE transcript, then use ccds match
    
    # print("geneInfo")
    # print(geneInfo)
    # print("transcripts.ens")
    # print(transcripts.ens)
    # print("exons.ens")
    # print(exons.ens)
    
    
    transcripts.ens.list <- ExtractMatchedRefseqEnsTranscripts(transcripts.ens)
    
    # Get exons for transcripts that have a refseq ID (from mane or from ccds) as a list
    
    exons.ens.list <- getExonsList(exons.ens, transcripts.ens, transcripts.ens.list)
    
    # print("transcripts.ens.list")
    # print(transcripts.ens.list)
    # 
    # print("exons.ens.list")
    # print(exons.ens.list)
    
    
    # enter VCF for point variant
    if (var_parsed$type  == "substitution" & is.null(var_parsed$startIntron))
    {
      var_vcf_toReturn <- build_VCF_entry_SNV(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    if (var_parsed$type  == "substitution" & !is.null(var_parsed$startIntron))
    {
      var_vcf_toReturn <- build_VCF_entry_SNV_intronic(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    if (var_parsed$type  == "intervening sequence")
    {
      var_vcf_toReturn <- build_VCF_entry_SNV_intronic_IVS(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    else if (var_parsed$type == "singleDuplication")
    {
      var_vcf_toReturn <- build_VCF_entry_singleDup(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    else if (var_parsed$type == "singledeletion")
    {
      var_vcf_toReturn <- build_VCF_entry_singleDel(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    else if (var_parsed$type == "delins")
    {
      var_vcf_toReturn <- build_VCF_entry_delins(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    else if (var_parsed$type == "duplication")
    {
      var_vcf_toReturn <- build_VCF_entry_duplication(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    else if (var_parsed$type == "deletion")
    {
      var_vcf_toReturn <- build_VCF_entry_deletion(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    else if (var_parsed$type == "insertion")
    {
      var_vcf_toReturn <- build_VCF_entry_insertion(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    else if (var_parsed$type == "invalid" & qgene2[[4]] == "deletion")
    {
      var_vcf_toReturn <- build_VCF_entry_SV(qgene2, var_parsed, transcripts.ens, exons.ens.list, ensembl)
    }
    
    return(var_info = var_vcf_toReturn)
    
  }
  
  
  # build SNV intronic VCF entry
  build_VCF_entry_SNV_intronic <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside intronic function for ", qgene2[[2]]))
    
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    
    #note that for intronic variants, if on antisense / neg strand, must do opposite
    if(exons.ens.list[[1]]$strand[1] == 1)
    {
      var_intron_shift <- var_parsed$startIntron
    }
    else if (exons.ens.list[[1]]$strand[1] == -1)
    {
      var_intron_shift <- -(var_parsed$startIntron)
    }
    
    #####
    # if(any(grepl(var_tx_rs_accnum,transcripts.ens$transcript_ccds_match)))
    # {
    #   #tx ref accnum matches ens version by ccds
    #   # ensembl way
    #   canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == var_tx_rs_accnum), "ensembl_transcript_id"]
    #   
    #   #possible that ccds_match will return multiple ens transcripts. take the one that matches MANE or otherwise the first one
    #   if(length(canon_transcript) > 1 & 
    #      transcripts.ens[which(transcripts.ens$transcript_is_canonical ==1), "transcript_mane_select"] == var_tx_rs_accnum)
    #   {
    #     canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_mane_select == var_tx_rs_accnum), "ensembl_transcript_id"]
    #   }
    #   else
    #   {
    #     canon_transcript <- canon_transcript[1]
    #   }
    #   
    #   genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #   print("supplied var ref transcr matches ref transcr in ensembl table; proceed with esnembl way of deriving position")
    #   
    #   
    # }
    # else if(var_tx_rs_accnum != "NA")
    # {
    #   #tx ref accnum does not match ens canon_tx
    #   accession <- sub("\\..*$", "", var_tx_rs_accnum)
    #   pattern <- paste0("^", accession)
    #   ens_vers_transcr <- unique(transcripts.ens[grep(pattern, transcripts.ens$transcript_ccds_cap), "transcript_ccds_match"])
    #   use_canon_tx<-0
    #   
    #   #if there is no match, assign canoncial to ens_ver_transcr
    #   if(is_empty(ens_vers_transcr))
    #   {
    #     ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_mane_select"]
    #     use_canon_tx<-1
    #     print("no match in ensembl, using ensembl canonical transcript")
    #   }
    #   
    #   #pull var_ref_transcr fasta and cds_interval
    #   var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
    #   var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
    #   var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
    #   var_tx_rs_accnum.cds <- paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = "")
    #   
    #   #pull ens_ref_transr fasta and cds interval
    #   ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
    #   ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
    #   ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
    #   ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
    #   
    #   ens_vs_var_transcr.cds <- pairwiseAlignment(
    #     ens_vers_transcr.cds,
    #     var_tx_rs_accnum.cds
    #   )
    #   # check if CDS length equal
    #   
    #   if (var_tx_rs_accnum.cds_len == ens_vers_transcr.cds_len)
    #   {
    #     
    #     #cDNA position on the ens tx
    #     #cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
    #     
    #     #check if percent identity is 100% i.e. CDS are identical between versions
    #     if (pid(ens_vs_var_transcr.cds) == 100 )
    #     {
    #       #standard ensembl way to calc genomic position
    #       canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #       if(use_canon_tx == 1)
    #       {
    #         canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
    #       }
    #       genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #       
    #     }
    #     else
    #     {
    #       #check if var_parsed$start is not in mismatch table
    #       mmt <- mismatchTable(ens_vs_var_transcr.cds)
    #       if (checkVarVsMismatch(var_parsed, mmt))
    #       {
    #         #standard ensembl way to calc genomic position
    #         canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #         genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #         
    #       }
    #       else
    #       {
    #         #error out
    #         print("ref seqeunce changed, variant is likely incorrect")
    #         
    #         #dump qgene to some sort of log
    #         return(NULL)
    #         
    #       }
    #     }
    #   }
    #   else 
    #   {
    #     #CDS length is not equal
    #     
    #     canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #     var_start_adjusted <- adjustVarLoc(var_parsed$start, ens_vs_var_transcr.cds)
    #     
    #     
    #     genomic_position_ens_result <- calcGenomicPositionEns(var_start_adjusted, canon_transcript)
    #     
    #     #new cDNA position on the ens tx
    #     # cdna_position <- var_loc_adjusted + ens_vers_transcr.cds_int[1,"from"] - 1
    #     
    #     
    #     print("ref CDS sequence length changed, position adjusted for indels between transcripts")  
    #     
    #   }
    #   
    # }
    # else if(var_tx_rs_accnum == "NA")
    # {
    #   #no tx ref accnum supplied, use ens canon_tx
    #   canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
    #   #ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_ccds_match"]
    #   
    #   
    #   genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #   
    #   #get tx seq
    #   # ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
    #   # ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
    #   # ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
    #   # ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
    #   # 
    #   
    #   #cDNA position on the ens tx
    #   #cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
    #   
    #   
    #   #manually extract genomic position
    #   genomic_position_man <- calcGenomicPosition(var_loc =  var_parsed$start, 
    #                                                     exons.ens.list =  exons.ens.list, 
    #                                                     transcripts.ens =  transcripts.ens, 
    #                                                     canon_transcript = canon_transcript, 
    #                                                     ensembl = ensembl)
    #   
    #   #check if ensembl and manual extractions agree, then check if ref seq from ens matches transcript
    #   
    #   
    #   if(genomic_position_man == genomic_position_ens && as.data.frame(checkGenRef)$value == var_parsed$ancestral) #matches; likely correct transcript chosen
    #   {
    #     print("No var trascript ref supplied, used canonical transcript - both ens and manual ways align")
    #   }
    #   else
    #   {
    #     print("no match between ens and man")
    #     return(NULL)
    #   }
    # }
    
    
    #####
    #construct output
    
    genomic_position_ens_result <- getGenomicPosition.simple(var_tx_rs_accnum,var_loc = var_parsed$start,
                                                             exons.ens.list = exons.ens.list, transcripts.ens = transcripts.ens, ensembl)
    
    
    var_intron_loc_adjust <- start(genomic_position_ens_result$gen_pos[[1]]) + var_intron_shift
    
    #construct query to pull ref seq
    ens.query <- paste0("/sequence/region/human/", 
                        as.character(runValue(seqnames(genomic_position_ens_result$gen_pos[[1]]))),
                        ":",
                        var_intron_loc_adjust,
                        "..",
                        var_intron_loc_adjust
                        #":",``
                        #"?"
    )
    
    seq <- GET(paste(server, ens.query, sep = ""), content_type("text/plain"))
    
    stop_for_status(seq)
    
    
    REF.seq <- print(content(seq))
    
    #build VCF file entries
    CHROM = as.character(runValue(seqnames(genomic_position_ens_result$gen_pos[[1]])))
    POS = var_intron_loc_adjust
    REF = strsplit(REF.seq, NULL)[[1]][1]
    ALT = checkSense(var_parsed$variant, exons.ens.list[[1]]$strand[1])
    # INFO.VCF <- paste0(
    #   "SVTYPE=DEL;END=",
    #   end(genomic_position_ens_result[[1]])+1
    # )
    
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
    #print(paste(CHROM, POS, REF, ALT, sep = " "))
    #print(var_toAdd)
    return(var_toAdd)
  }
  
  # build SNV intronic VCF entry
  build_VCF_entry_SNV_intronic_IVS <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside intronic IVS function for ", qgene2[[2]]))
    
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    
    #note that for intronic variants, if on antisense / neg strand, must do opposite
    if(exons.ens.list[[1]]$strand[1] == 1)
    {
      var_intron_shift <- var_parsed$IntronOffset
    }
    else if (exons.ens.list[[1]]$strand[1] == -1)
    {
      var_intron_shift <- -(var_parsed$IntronOffset)
    }
    
    #note that for intronic variants, if on antisense / neg strand, must do opposite
    if(var_parsed$IntronOffset < 0) # position relative to start of 2nd exon
    {
      exon_rank <- var_parsed$IVS_number + 1
      exon_ref_loc <- as.data.table(exons.ens.list[[1]])[rank == exon_rank, cds_start]
    }
    else if (var_parsed$IntronOffset > 0) #position relative to end of 1st exon
    {
      exon_rank <- var_parsed$IVS_number
      exon_ref_loc <- as.data.table(exons.ens.list[[1]])[rank == exon_rank, cds_end]
    }
    
    exon_ref_genome_loc <- getGenomicPosition.simple(var_tx_rs_accnum = var_tx_rs_accnum,
                                                     var_loc = exon_ref_loc,
                                                     transcripts.ens = transcripts.ens, 
                                                     exons.ens.list = exons.ens.list,
                                                     ensembl = ensembl)
    
    
    
    
    #####
    #construct output
    
    var_intron_loc_adjust <- start(exon_ref_genome_loc$gen_pos[[1]]) + var_intron_shift
    
    #construct query to pull ref seq
    ens.query <- paste0("/sequence/region/human/", 
                        as.character(runValue(seqnames(exon_ref_genome_loc$gen_pos[[1]]))),
                        ":",
                        var_intron_loc_adjust,
                        "..",
                        var_intron_loc_adjust
                        #":",``
                        #"?"
    )
    
    seq <- GET(paste(server, ens.query, sep = ""), content_type("text/plain"))
    
    stop_for_status(seq)
    
    
    REF.seq <- print(content(seq))
    
    #build VCF file entries
    CHROM = as.character(runValue(seqnames(exon_ref_genome_loc$gen_pos[[1]])))
    POS = var_intron_loc_adjust
    REF = strsplit(REF.seq, NULL)[[1]][1]
    ALT = checkSense(var_parsed$variant, exons.ens.list[[1]]$strand[1])
    # INFO.VCF <- paste0(
    #   "SVTYPE=DEL;END=",
    #   end(genomic_position_ens_result[[1]])+1
    # )
    
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
    #print(paste(CHROM, POS, REF, ALT, sep = " "))
    #print(var_toAdd)
    return(var_toAdd)
  }
  
  # build singleDup VCF entry
  build_VCF_entry_singleDup <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside singleDup function for ", qgene2[[2]]))
    
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    
    
    
    genomic_position_ens_result <- getGenomicPosition.simple(var_tx_rs_accnum,var_loc = var_parsed$start,
                                                             exons.ens.list = exons.ens.list, transcripts.ens = transcripts.ens, ensembl)
    
    genomic_position_ens <- start(genomic_position_ens_result$gen_pos[[1]]) #exact pos of duplicated base
    canon_transcript <- genomic_position_ens_result$canon_tx
    #####
    
    # if(any(grepl(var_tx_rs_accnum,transcripts.ens$transcript_ccds_match)))
    # {
    #   #tx ref accnum matches ens version by ccds
    #   # ensembl way
    #   canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == var_tx_rs_accnum), "ensembl_transcript_id"]
    #   
    #   #possible that ccds_match will return multiple ens transcripts. take the one that matches MANE or otherwise the first one
    #   if(length(canon_transcript) > 1 & 
    #      transcripts.ens[which(transcripts.ens$transcript_is_canonical ==1), "transcript_mane_select"] == var_tx_rs_accnum)
    #   {
    #     canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_mane_select == var_tx_rs_accnum), "ensembl_transcript_id"]
    #   }
    #   else
    #   {
    #     canon_transcript <- canon_transcript[1]
    #   }
    #   
    #   genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start+1, canon_transcript)
    #   genomic_position_ens <- start(genomic_position_ens_result[[1]])
    #   
    #   #cDNA position on the ens tx
    #   cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
    #   
    #   print("supplied var ref transcr matches ref transcr in ensembl table; proceed with esnembl way of deriving position")
    #   
    # }
    # else if(var_tx_rs_accnum != "NA")
    # {
    #   #tx ref accnum does not match ens canon_tx
    #   #manual way
    #   #pull NCBI accession number linked by CCDS to accession of supplied reference transcript in esnembl db (ie latest version)
    #   accession <- sub("\\..*$", "", var_tx_rs_accnum)
    #   pattern <- paste0("^", accession)
    #   ens_vers_transcr <- unique(transcripts.ens[grep(pattern, transcripts.ens$transcript_ccds_cap), "transcript_ccds_match"])
    #   use_canon_tx<-0
    #   
    #   #if there is no match, assign canoncial to ens_ver_transcr
    #   if(is_empty(ens_vers_transcr))
    #   {
    #     ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_ccds_match"]
    #     use_canon_tx<-1
    #     print("no match in ensembl, using ensembl canonical transcript")
    #   }
    #   
    #   #pull var_ref_transcr fasta and cds_interval
    #   var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
    #   var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
    #   var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
    #   var_tx_rs_accnum.cds <- paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = "")
    #   
    #   # var_fun <- "NM_001256789.3"
    #   # 
    #   # var_fun.fasta <- getFasta(var_fun)
    #   # var_fun.cds_int <- getCDS(var_tx_rs_accnum)
    #   # var_fun.cds_len <- (var_fun.cds_int[,3] - var_fun.cds_int[,2])
    #   # 
    #   
    #   #pull ens_ref_transr fasta and cds interval
    #   ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
    #   ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
    #   ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
    #   ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
    #   
    #   # #pairwise align full tx
    #   # ens_vs_var_transcr <- pairwiseAlignment(paste0(ens_vers_transcr.fasta[[1]], collapse = ""),
    #   #                                         paste0(var_tx_rs_accnum.fasta[[1]], collapse = ""))
    #   
    #   #pattern is ens
    #   #subject is var
    #   
    #   ens_vs_var_transcr.cds <- pairwiseAlignment(
    #     ens_vers_transcr.cds,
    #     var_tx_rs_accnum.cds
    #   )
    #   
    #   
    #   # check if CDS length equal
    #   
    #   if (var_tx_rs_accnum.cds_len == ens_vers_transcr.cds_len)
    #   {
    #     
    #     #cDNA position on the ens tx
    #     cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
    #     
    #     #check if percent identity is 100% i.e. CDS are identical between versions
    #     if (pid(ens_vs_var_transcr.cds) == 100 )
    #     {
    #       #standard ensembl way to calc genomic position
    #       canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #       if(use_canon_tx == 1)
    #       {
    #         canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
    #       }
    #       genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start+1, canon_transcript)
    #       genomic_position_ens <- start(genomic_position_ens_result[[1]])
    #     }
    #     else
    #     {
    #       #check if var_parsed$start is not in mismatch table
    #       mmt <- mismatchTable(ens_vs_var_transcr.cds)
    #       if (checkVarVsMismatch(var_parsed, mmt))
    #       {
    #         #standard ensembl way to calc genomic position
    #         canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #         genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start+1, canon_transcript)
    #         genomic_position_ens <- start(genomic_position_ens_result[[1]])
    #       }
    #       else
    #       {
    #         #error out
    #         print("ref seqeunce changed, variant is likely incorrect")
    #         
    #         #dump qgene to some sort of log
    #         return(NULL)
    #         
    #       }
    #     }
    #   }
    #   else 
    #   {
    #     #CDS length is not equal
    #     
    #     canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #     var_loc_adjusted <- adjustVarLoc(var_parsed$start, ens_vs_var_transcr.cds)
    #     
    #     genomic_position_ens_result <- calcGenomicPositionEns(var_loc_adjusted, canon_transcript)
    #     genomic_position_ens <- start(genomic_position_ens_result[[1]])
    #     
    #     #new cDNA position on the ens tx
    #     cdna_position <- var_loc_adjusted + ens_vers_transcr.cds_int[1,"from"] - 1
    #     
    #     
    #     print("ref CDS sequence length changed, position adjusted for indels between transcripts")  
    #     
    #   }
    # }
    # else if(var_tx_rs_accnum == "NA")
    # {
    #   #no tx ref accnum supplied, use ens canon_tx
    #   
    #   canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
    #   
    #   genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #   genomic_position_ens <- start(genomic_position_ens_result[[1]])
    #   
    #   #extract ref sequence at extracted genomic location - to check if matches variant
    #   #checkGenRef <- getSeq(BSgenome.Hsapiens.NCBI.GRCh38, genomic_position_ens_result)
    #   
    #   #new cDNA position on the ens tx
    #   cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
    #   
    #   
    #   #manually extract genomic position
    #   genomic_position_man <- calcGenomicPosition(var_loc =  var_parsed$start, 
    #                                               exons.ens.list =  exons.ens.list, 
    #                                               transcripts.ens =  transcripts.ens, 
    #                                               canon_transcript = canon_transcript, 
    #                                               ensembl = ensembl)
    #   
    #   if(genomic_position_man == genomic_position_ens) #matches; likely correct transcript chosen
    #   {
    #     print("No var trascript ref supplied, used canonical transcript - both ens and manual ways align")
    #   }
    #   else
    #   {
    #     print("no match between ens and man")
    #     return(NULL)
    #   }
    # }
    
    #####
    
    
    #extract ref sequence at extracted genomic location - to check if matches variant
    
    #note direction of dup will change based on strand
    if(exons.ens.list[[1]]$strand[1] == 1)
    {
      var_shift <- 1
      REF.fullseq <- getRefSeq(server, genomic_position_ens-var_shift, genomic_position_ens,   
                               as.character(runValue(seqnames(genomic_position_ens_result[[1]]))))
      POS = genomic_position_ens
    }
    else if (exons.ens.list[[1]]$strand[1] == -1)
    {
      var_shift <- -1
      REF.fullseq <- getRefSeq(server, genomic_position_ens+var_shift, genomic_position_ens,  
                               as.character(runValue(seqnames(genomic_position_ens_result[[1]]))))
      POS <- genomic_position_ens+var_shift
    }
    
    
    REF.seq <- substr(REF.fullseq, 1,1)
    ALT.seq <- REF.fullseq
    
    #check if checkGenRef matches REF.seq
    # if(as.character(unlist(checkGenRef)) == REF.seq)
    # {print("variant mapped correctly")}
    # 
    CHROM = as.character(runValue(seqnames(genomic_position_ens_result$gen_pos[[1]])))
    #POS = genomic_position_ens
    REF = REF.seq
    ALT =  ALT.seq
    #print(paste(CHROM, POS, REF, ALT, sep = " "))
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
    
    #print(var_toAdd)
    return(var_toAdd)
  }
  
  # build singleDel VCF entry
  build_VCF_entry_singleDel <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside singelDel function for ", qgene2[[2]]))
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    
    gen_pos_result <- getGenomicPosition.simple(var_tx_rs_accnum, var_parsed$start, transcripts.ens, exons.ens.list ,ensembl)
    genomic_position_ens <- start(gen_pos_result$gen_pos[[1]])
    canon_transcript <- gen_pos_result$canon_tx
    
    #adjust based on strand orientation
    #note that for deletions, if on antisense / neg strand, must calc Gen Position from "end"
    # if(exons.ens.list[[1]]$strand[1] == 1)
    # {
    #   genomic_position_ens.adj <- genomic_position_ens-1
    # }
    # else if (exons.ens.list[[1]]$strand[1] == -1)
    # {
    #   genomic_position_ens.adj <- genomic_position_ens+1
    # }
    
    genomic_position_ens.adj <- genomic_position_ens-1
    
    
    #extract ref sequence at extracted genomic location - to check if matches variant
    REF.fullseq <- getRefSeq(server, genomic_position_ens.adj, genomic_position_ens.adj+1, as.character(runValue(seqnames(gen_pos_result[[1]]))))
    REF.seq <- REF.fullseq
    ALT.seq <- substr(REF.fullseq, 1,1)
    
    
    CHROM = as.character(runValue(seqnames(gen_pos_result$gen_pos[[1]])))
    POS = genomic_position_ens.adj
    REF = REF.seq
    ALT = ALT.seq
    print(paste(CHROM, POS, REF, ALT, sep = " "))
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
    
    #print(var_toAdd)
    return(var_toAdd)
  }
  
  # build delins VCF entryß
  build_VCF_entry_delins <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside delins function for ", qgene2[[2]]))
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    
    #note that for deletions, if on antisense / neg strand, must calc Gen Position from "end"
    if(exons.ens.list[[1]]$strand[1] == 1)
    {
      var_loc <- var_parsed$start
    }
    else if (exons.ens.list[[1]]$strand[1] == -1)
    {
      var_loc <- var_parsed$end
    }
    
    
    
    if(any(grepl(var_tx_rs_accnum,transcripts.ens$transcript_ccds_match)))
    {
      #tx ref accnum matches ens version by ccds
      #tx ref accnum matches ens version by ccds
      # ensembl way
      canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == var_tx_rs_accnum), "ensembl_transcript_id"]
      
      #possible that ccds_match will return multiple ens transcripts. take the one that matches MANE or otherwise the first one
      if(length(canon_transcript) > 1 & 
         transcripts.ens[which(transcripts.ens$transcript_is_canonical ==1), "transcript_mane_select"] == var_tx_rs_accnum)
      {
        canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_mane_select == var_tx_rs_accnum), "ensembl_transcript_id"]
      }
      else
      {
        canon_transcript <- canon_transcript[1]
      }
      
      genomic_position_ens_result <- calcGenomicPositionEns(var_loc, canon_transcript)
      genomic_position_ens <- start(genomic_position_ens_result[[1]])
      
      
      
      #pull ens_ref_transr fasta and cds interval
      
      ens_vers_transcr.fasta <- getFasta(var_tx_rs_accnum)
      ens_vers_transcr.cds_int <- getCDS(var_tx_rs_accnum)
      ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
      ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
      
      
      #cDNA position on the ens tx
      cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
      
      print("supplied var ref transcr matches ref transcr in ensembl table; proceed with esnembl way of deriving position")
      
      
    }
    else if(var_tx_rs_accnum != "NA")
    {
      #tx ref accnum does not match ens canon_tx
      
      #manual way
      #pull NCBI accession number linked by CCDS to accession of supplied reference transcript in esnembl db (ie latest version)
      accession <- sub("\\..*$", "", var_tx_rs_accnum)
      pattern <- paste0("^", accession)
      ens_vers_transcr <- unique(transcripts.ens[grep(pattern, transcripts.ens$transcript_ccds_cap), "transcript_ccds_match"])
      use_canon_tx<-0
      
      #if there is no match, assign canoncial to ens_ver_transcr
      if(is_empty(ens_vers_transcr))
      {
        ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_mane_select"]
        use_canon_tx<-1
        print("no match in ensembl, using ensembl canonical transcript")
      }
      
      #pull var_ref_transcr fasta and cds_interval
      var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
      var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
      var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
      var_tx_rs_accnum.cds <- paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = "")
      
      # var_fun <- "NM_001256789.3"
      # 
      # var_fun.fasta <- getFasta(var_fun)
      # var_fun.cds_int <- getCDS(var_tx_rs_accnum)
      # var_fun.cds_len <- (var_fun.cds_int[,3] - var_fun.cds_int[,2])
      # 
      
      #pull ens_ref_transr fasta and cds interval
      ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
      ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
      ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
      ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
      
      # #pairwise align full tx
      # ens_vs_var_transcr <- pairwiseAlignment(paste0(ens_vers_transcr.fasta[[1]], collapse = ""),
      #                                         paste0(var_tx_rs_accnum.fasta[[1]], collapse = ""))
      
      #pattern is ens
      #subject is var
      
      ens_vs_var_transcr.cds <- pairwiseAlignment(
        ens_vers_transcr.cds,
        var_tx_rs_accnum.cds
      )
      
      # check if CDS length equal
      
      if (var_tx_rs_accnum.cds_len == ens_vers_transcr.cds_len)
      {
        
        #cDNA position on the ens tx
        cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
        
        #check if percent identity is 100% i.e. CDS are identical between versions
        if (pid(ens_vs_var_transcr.cds) == 100 )
        {
          #standard ensembl way to calc genomic position
          canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
          if(use_canon_tx == 1)
          {
            canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
          }
          genomic_position_ens_result <- calcGenomicPositionEns(var_loc, canon_transcript)
          genomic_position_ens <- start(genomic_position_ens_result[[1]])
        }
        else
        {
          #check if var_parsed$start is not in mismatch table
          mmt <- mismatchTable(ens_vs_var_transcr.cds)
          if (checkVarVsMismatch(var_parsed, mmt))
          {
            #standard ensembl way to calc genomic position
            canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
            genomic_position_ens_result <- calcGenomicPositionEns(var_loc, canon_transcript)
            genomic_position_ens <- start(genomic_position_ens_result[[1]])
          }
          else
          {
            #error out
            print("ref seqeunce changed, variant is likely incorrect")
            
            #dump qgene to some sort of log
            return(NULL)
            
          }
        }
      }
      else 
      {
        #CDS length is not equal
        
        canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
        var_loc_adjusted <- adjustVarLoc(var_loc, ens_vs_var_transcr.cds)
        
        genomic_position_ens_result <- calcGenomicPositionEns(var_loc_adjusted, canon_transcript)
        genomic_position_ens <- start(genomic_position_ens_result[[1]])
        
        #new cDNA position on the ens tx
        cdna_position <- var_loc_adjusted + ens_vers_transcr.cds_int[1,"from"] - 1
        
        
        print("ref CDS sequence length changed, position adjusted for indels between transcripts")  
        
      }
      
    }
    else if(var_tx_rs_accnum == "NA")
    {
      #no tx ref accnum supplied, use ens canon_tx
      #no tx ref accnum supplied, use ens canon_tx
      canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
      ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_ccds_match"]
      
      
      genomic_position_ens_result <- calcGenomicPositionEns(var_loc, canon_transcript)
      genomic_position_ens <- start(genomic_position_ens_result[[1]])
      
      #get tx seq
      ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
      ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
      ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
      ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
      
      
      #cDNA position on the ens tx
      cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
      
      
      #manually extract genomic position
      genomic_position_man <- calcGenomicPosition(var_loc =  var_loc, 
                                                  exons.ens.list =  exons.ens.list, 
                                                  transcripts.ens =  transcripts.ens, 
                                                  canon_transcript = canon_transcript, 
                                                  ensembl = ensembl)
      
      if(genomic_position_man == genomic_position_ens) #matches; likely correct transcript chosen
      {
        print("No var trascript ref supplied, used canonical transcript - both ens and manual ways align")
      }
      else
      {
        print("no match between ens and man")
        return(NULL)
      }
    }
    
    
    #from delins old code
    
    #for deletions and duplications that have both start and end
    if(!is.null(var_parsed$end))
    {
      REF.seq <- toupper(paste0(ens_vers_transcr.fasta[[1]][cdna_position:(cdna_position+(var_parsed$end-var_parsed$start))],
                                collapse = ""))
      ALT.seq <- var_parsed$variant
    }
    
    
    #extract ref sequence at extracted genomic location - to check if matches variant
    #checkGenRef <- getSeq(BSgenome.Hsapiens.NCBI.GRCh38, genomic_position_ens_result)
    
    CHROM = unique(transcripts.ens[transcripts.ens$ensembl_transcript_id == canon_transcript,]$chromosome_name)
    POS = genomic_position_ens
    REF = checkSense(REF.seq, strand = exons.ens.list[[1]]$strand[1])
    ALT = checkSense(ALT.seq, strand = exons.ens.list[[1]]$strand[1])
    
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
    #print(paste(CHROM, POS, REF, ALT, sep = " "))
    #print(var_toAdd)
    return(var_toAdd)
  }
  
  # build duplication VCF entry
  build_VCF_entry_duplication <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside duplication function for ", qgene2[[2]]))
    
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    
    var_start <- var_parsed$start
    var_end <- var_parsed$end
    
    #note that for deletions, if on antisense / neg strand, must calc Gen Position from "end"
    if(exons.ens.list[[1]]$strand[1] == 1)
    {
      
      if(!is_empty(var_parsed$startIntron))
      {
        var_startIntron_shift <- var_parsed$startIntron
      }
      if(!is_empty(var_parsed$endIntron))
      {
        var_endIntron_shift <- var_parsed$endIntron
      }
      
    }
    else if (exons.ens.list[[1]]$strand[1] == -1)
    {
      
      if(!is_empty(var_parsed$startIntron))
      {
        var_startIntron_shift <- -(var_parsed$startIntron)
      }
      if(!is_empty(var_parsed$endIntron))
      {
        var_endIntron_shift <- -(var_parsed$endIntron)
      }
    }
    
    genomic_position_ens_result <- getGenomicPosition.del(var_tx_rs_accnum, var_loc_start = var_start,
                                                          var_loc_end = var_end, transcripts.ens, 
                                                          exons.ens.list, ensembl)
    
    start_coord <- min(start(genomic_position_ens_result$gen_pos[[1]])) #start of the sequence is there the dup is from (i.e. end of dup)
    end_coord <- max(end(genomic_position_ens_result$gen_pos[[1]])) #end of the sequence is where it is inserted i.e. pos
    
    #end coord+1 is where the duplication / insertion starts
    canon_transcript <- genomic_position_ens_result$canon_tx
    
    var_end_vcf <- end_coord
    var_start_vcf <- start_coord
    
    #####
    #construct output
    if(!is_empty(var_parsed$startIntron))
    {
      var_start_vcf <- start_coord + var_startIntron_shift
    }
    if(!is_empty(var_parsed$endIntron))
    {
      var_end_vcf <- end_coord + var_endIntron_shift
    }
    
    #for special cases where start_coord == end_coord
    if(start_coord == end_coord)
    {
      var_start_vcf <- min(start_coord + var_startIntron_shift, start_coord + var_endIntron_shift)
      var_end_vcf <- max(end_coord + var_startIntron_shift, end_coord + var_endIntron_shift)
    }
    
    REF.seq <- getRefSeq(server, var_start_vcf, var_end_vcf, as.character(runValue(seqnames(genomic_position_ens_result[[1]]))))
    
    #build VCF file entries
    
    #### NOTE: snpeff throws an error and calls <DUP> variant as intron_variant
    
    ## work around
    
    if(var_end_vcf - var_start_vcf < 20)
    {
      CHROM = as.character(runValue(seqnames(genomic_position_ens_result$gen_pos[[1]])))
      POS = var_end_vcf
      REF = last(strsplit(REF.seq, NULL)[[1]])
      #REF = strsplit(REF.seq, NULL)[[1]][1]
      ALT = paste0(REF,REF.seq)
      INFO.VCF <- paste0(
        "SVTYPE=INS;END=",var_end_vcf + 1
      )
      
      var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT, INFO = INFO.VCF)
      #print(paste(CHROM, POS, REF, ALT, sep = " "))
      #print(var_toAdd)
      return(var_toAdd)
    }
    
    CHROM = as.character(runValue(seqnames(genomic_position_ens_result$gen_pos[[1]])))
    POS = var_start_vcf
    #REF = last(strsplit(REF.seq, NULL)[[1]])
    REF = strsplit(REF.seq, NULL)[[1]][1]
    ALT = "<DUP>"
    INFO.VCF <- paste0(
      "SVTYPE=DUP;SVLEN=",
      length(strsplit(REF.seq, NULL)[[1]]),";END=",var_start_vcf + length(strsplit(REF.seq, NULL)[[1]]) -1
    )
    
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT, INFO = INFO.VCF)
    #print(paste(CHROM, POS, REF, ALT, sep = " "))
    #print(var_toAdd)
    return(var_toAdd)
  }
  
  # build deletion VCF entry
  build_VCF_entry_deletion <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside deletion function for ", qgene2[[2]]))
    
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    
    var_start <- var_parsed$start
    var_end <- var_parsed$end
    
    genomic_position_ens_result <- getGenomicPosition.del(var_tx_rs_accnum, var_loc_start = var_start,
                                                          var_loc_end = var_end, transcripts.ens, 
                                                          exons.ens.list, ensembl)
    
    start_coord <- min(start(genomic_position_ens_result$gen_pos[[1]]))
    end_coord <- max(end(genomic_position_ens_result$gen_pos[[1]]))
    canon_transcript <- genomic_position_ens_result$canon_tx
    
    var_start_vcf <- start_coord-1
    var_end_vcf <- end_coord-1
    
    #note that for deletions, if on antisense / neg strand, must calc Gen Position from "end"
    if(exons.ens.list[[1]]$strand[1] == 1)
    {
      
      if(!is_empty(var_parsed$startIntron))
      {
        var_startIntron_shift <- var_parsed$startIntron
      }
      if(!is_empty(var_parsed$endIntron))
      {
        var_endIntron_shift <- var_parsed$endIntron
      }
      
      
    }
    else if (exons.ens.list[[1]]$strand[1] == -1)
    {
      
      if(!is_empty(var_parsed$startIntron))
      {
        var_startIntron_shift <- -(var_parsed$startIntron)
      }
      if(!is_empty(var_parsed$endIntron))
      {
        var_endIntron_shift <- -(var_parsed$endIntron)
      }
    }
    
    
    #####
    
    
    #####
    #construct output
    if(!is_empty(var_parsed$startIntron) & start_coord != end_coord)
    {
      var_start_vcf <- var_start_vcf + var_startIntron_shift
    }
    if(!is_empty(var_parsed$endIntron) & start_coord != end_coord)
    {
      var_end_vcf <- var_end_vcf + var_endIntron_shift
    }
    
    #for special cases where start_coord == end_coord
    if(start_coord == end_coord & !is_empty(var_parsed$endIntron) & !is_empty(var_parsed$startIntron))
    {
      var_start_vcf <- min(var_start_vcf + var_startIntron_shift, var_end_vcf + var_endIntron_shift)
      var_end_vcf <- max(var_start_vcf + var_startIntron_shift, var_end_vcf + var_endIntron_shift)
    }
    
    #construct query to pull ref seq
    
    REF.seq <- getRefSeq(server, var_start_vcf, var_end_vcf, as.character(runValue(seqnames(genomic_position_ens_result[[1]]))))
    
    #build VCF file entries
    CHROM = as.character(runValue(seqnames(genomic_position_ens_result$gen_pos[[1]])))
    POS = var_start_vcf
    REF = strsplit(REF.seq, NULL)[[1]][1]
    ALT = "<DEL>"
    INFO.VCF <- paste0(
      "SVTYPE=DEL;END=",
      var_end_vcf+1
    )
    
    
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT, INFO = INFO.VCF)
    #print(paste(CHROM, POS, REF, ALT, sep = " "))
    #print(var_toAdd)
    return(var_toAdd)
  }
  
  # build duplication VCF entry
  build_VCF_entry_insertion <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside insertion function for ", qgene2[[2]]))
    
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    var_start <- var_parsed$start
    
    #note that for deletions, if on antisense / neg strand, must calc Gen Position from "end"
    if(exons.ens.list[[1]]$strand[1] == 1)
    {
      var_len <- length(strsplit(var_parsed$variant, NULL)[[1]])
      var_end <- 1
    }
    else if (exons.ens.list[[1]]$strand[1] == -1)
    {
      var_len <- -(length(strsplit(var_parsed$variant, NULL)[[1]]))
      var_end <- -1
    }
    
    gen_pos_result <- getGenomicPosition.simple(var_tx_rs_accnum, var_start, transcripts.ens, exons.ens.list ,ensembl)
    genomic_position_ens <- start(gen_pos_result$gen_pos[[1]])
    canon_transcript <- gen_pos_result$canon_tx
    
    #####
    # if(any(grepl(var_tx_rs_accnum,transcripts.ens$transcript_ccds_match)))
    # {
    #   #tx ref accnum matches ens version by ccds
    #   # ensembl way
    #   canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == var_tx_rs_accnum), "ensembl_transcript_id"]
    #   
    #   #possible that ccds_match will return multiple ens transcripts. take the one that matches MANE or otherwise the first one
    #   if(length(canon_transcript) > 1 & 
    #      transcripts.ens[which(transcripts.ens$transcript_is_canonical ==1), "transcript_mane_select"] == var_tx_rs_accnum)
    #   {
    #     canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_mane_select == var_tx_rs_accnum), "ensembl_transcript_id"]
    #   }
    #   else
    #   {
    #     canon_transcript <- canon_transcript[1]
    #   }
    #   
    #   genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #   print("supplied var ref transcr matches ref transcr in ensembl table; proceed with esnembl way of deriving position")
    #   
    # }
    # else if(var_tx_rs_accnum != "NA")
    # {
    #   #tx ref accnum does not match ens canon_tx
    #   #manual way
    #   #pull NCBI accession number linked by CCDS to accession of supplied reference transcript in esnembl db (ie latest version)
    #   accession <- sub("\\..*$", "", var_tx_rs_accnum)
    #   pattern <- paste0("^", accession)
    #   ens_vers_transcr <- unique(transcripts.ens[grep(pattern, transcripts.ens$transcript_ccds_cap), "transcript_ccds_match"])
    #   use_canon_tx<-0
    #   
    #   #if there is no match, assign canoncial to ens_ver_transcr
    #   if(is_empty(ens_vers_transcr))
    #   {
    #     ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_mane_select"]
    #     use_canon_tx<-1
    #     print("no match in ensembl, using ensembl canonical transcript")
    #   }
    #   
    #   #pull var_ref_transcr fasta and cds_interval
    #   var_tx_rs_accnum.fasta <- getFasta(var_tx_rs_accnum)
    #   var_tx_rs_accnum.cds_int <- getCDS(var_tx_rs_accnum)
    #   var_tx_rs_accnum.cds_len <- (var_tx_rs_accnum.cds_int[,3] - var_tx_rs_accnum.cds_int[,2])
    #   var_tx_rs_accnum.cds <- paste0(var_tx_rs_accnum.fasta[[1]][var_tx_rs_accnum.cds_int[,2]:var_tx_rs_accnum.cds_int[,3]], collapse = "")
    #   
    #   #pull ens_ref_transr fasta and cds interval
    #   ens_vers_transcr.fasta <- getFasta(ens_vers_transcr)
    #   ens_vers_transcr.cds_int <- getCDS(ens_vers_transcr)
    #   ens_vers_transcr.cds_len <- (ens_vers_transcr.cds_int[,3] - ens_vers_transcr.cds_int[,2])
    #   ens_vers_transcr.cds <- paste0(ens_vers_transcr.fasta[[1]][ens_vers_transcr.cds_int[,2]:ens_vers_transcr.cds_int[,3]], collapse = "")
    #   
    #   # #pairwise align full tx
    #   # ens_vs_var_transcr <- pairwiseAlignment(paste0(ens_vers_transcr.fasta[[1]], collapse = ""),
    #   #                                         paste0(var_tx_rs_accnum.fasta[[1]], collapse = ""))
    #   
    #   #pattern is ens
    #   #subject is var
    #   
    #   ens_vs_var_transcr.cds <- pairwiseAlignment(
    #     ens_vers_transcr.cds,
    #     var_tx_rs_accnum.cds
    #   )
    #   
    #   # check if CDS length equal
    #   
    #   if (var_tx_rs_accnum.cds_len == ens_vers_transcr.cds_len)
    #   {
    #     
    #     #cDNA position on the ens tx
    #     #cdna_position <- var_parsed$start + ens_vers_transcr.cds_int[1,"from"] - 1
    #     
    #     #check if percent identity is 100% i.e. CDS are identical between versions
    #     if (pid(ens_vs_var_transcr.cds) == 100 )
    #     {
    #       #standard ensembl way to calc genomic position
    #       canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #       if(use_canon_tx == 1)
    #       {
    #         canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
    #       }
    #       genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #       
    #     }
    #     else
    #     {
    #       #check if var_parsed$start is not in mismatch table
    #       mmt <- mismatchTable(ens_vs_var_transcr.cds)
    #       if (checkVarVsMismatch(var_parsed, mmt))
    #       {
    #         #standard ensembl way to calc genomic position
    #         canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #         genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #         
    #       }
    #       else
    #       {
    #         #error out
    #         print("ref seqeunce changed, variant is likely incorrect")
    #         
    #         #dump qgene to some sort of log
    #         return(NULL)
    #         
    #       }
    #     }
    #   }
    #   else 
    #   {
    #     #CDS length is not equal
    #     
    #     canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_ccds_match == ens_vers_transcr), "ensembl_transcript_id"][1]
    #     var_start_adjusted <- adjustVarLoc(var_parsed$start, ens_vs_var_transcr.cds)
    #     
    #     genomic_position_ens_result <- calcGenomicPositionEns(var_start_adjusted, canon_transcript)
    #     
    #     #new cDNA position on the ens tx
    #     # cdna_position <- var_loc_adjusted + ens_vers_transcr.cds_int[1,"from"] - 1
    #     
    #     
    #     print("ref CDS sequence length changed, position adjusted for indels between transcripts")  
    #     
    #   }
    #   
    # }
    # else if(var_tx_rs_accnum == "NA")
    # {
    #   #no tx ref accnum supplied, use ens canon_tx
    #   canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "ensembl_transcript_id"]
    #   #ens_vers_transcr <- transcripts.ens[which(transcripts.ens$transcript_is_canonical == 1), "transcript_ccds_match"]
    #   
    #   
    #   genomic_position_ens_result <- calcGenomicPositionEns(var_parsed$start, canon_transcript)
    #   
    #   genomic_position_man <- calcGenomicPosition(var_loc =  var_parsed$start, 
    #                                                     exons.ens.list =  exons.ens.list, 
    #                                                     transcripts.ens =  transcripts.ens, 
    #                                                     canon_transcript = canon_transcript, 
    #                                                     ensembl = ensembl)
    #   
    #   if(genomic_position_man == start(genomic_position_ens_result[[1]]) ) #matches; likely correct transcript chosen
    #   {
    #     print("No var trascript ref supplied, used canonical transcript - both ens and manual ways align")
    #   }
    #   else
    #   {
    #     print("no match between ens and man")
    #     return(NULL)
    #   }
    # }
    
    #####
    #construct query to pull ref seq
    
    #determine start and end coordinates based on strand orientation
    if(var_end == -1 )
    {
      var_start_vcf <- end(gen_pos_result$gen_pos[[1]])+var_end
      var_end_vcf <- start(gen_pos_result$gen_pos[[1]])
      
      REF.seq <- getRefSeq(server, var_start_vcf, var_end_vcf, as.character(runValue(seqnames(gen_pos_result[[1]]))))
      
    } else if (var_end == 1)
    {
      var_start_vcf <- start(gen_pos_result$gen_pos[[1]])-var_end
      var_end_vcf <- end(gen_pos_result$gen_pos[[1]])
      
      REF.seq <- getRefSeq(server, var_start_vcf, var_end_vcf, as.character(runValue(seqnames(gen_pos_result[[1]]))))
    }
    
    #build VCF file entries
    CHROM = as.character(runValue(seqnames(gen_pos_result[[1]])))
    POS = var_start_vcf
    REF = strsplit(REF.seq, NULL)[[1]][1]
    ALT = paste0(REF,checkSense(var_parsed$variant, exons.ens.list[[1]]$strand[1]))
    INFO.VCF <- paste0(
      "SVTYPE=INS;END=",
      var_end_vcf
    )
    
    
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT, INFO = INFO.VCF)
    
    return(var_toAdd)
  }
  
  # build Structural Variant VCF entry
  build_VCF_entry_SV <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print(paste0("inside Structural Variant function for ", qgene2, collapse = " "))
    #if reference transcript is known
    var_tx_rs_accnum <- qgene2[[4]]
    
    if(any(grepl(var_tx_rs_accnum,transcripts.ens$transcript_ccds_match)))
    {
      #tx ref accnum matches ens version by ccds
      
    }
    else if(var_tx_rs_accnum != "NA")
    {
      #tx ref accnum does not match ens canon_tx
      
    }
    else if(var_tx_rs_accnum == "NA")
    {
      #no tx ref accnum supplied, use ens canon_tx
    }
    
    return(NULL)
  }
  
  
  # build SNV VCF entry
  build_VCF_entry_SNV <- function(qgene2, var_parsed, transcripts.ens, exons.ens.list ,ensembl)
  {
    print("inside build entry SNV")
    # if the reference transcript is known
    #var_tx_rs_accnum <- "NM_005183.2" #will be extracted from table, if available - part of qgene2
    var_tx_rs_accnum <- qgene2[[4]]
    
    genomicPositionEns <- getGenomicPosition.simple(var_tx_rs_accnum, var_loc = var_parsed$start, 
                                                    transcripts.ens = transcripts.ens, exons.ens.list = exons.ens.list,ensembl = ensembl)
    
    genomic_position_ens <- start(genomicPositionEns$gen_pos[[1]])
    canon_transcript <- genomicPositionEns$canon_tx
    
    
    REF.seq <- getRefSeq(server, genomic_position_ens, 
                         genomic_position_ens, as.character(runValue(seqnames(genomicPositionEns$gen_pos[[1]]))))
    
    
    
    
    #####
    CHROM = as.character(runValue(seqnames(genomicPositionEns$gen_pos[[1]])))
    POS = genomic_position_ens
    REF = checkSense(var_parsed$ancestral, exons.ens.list[[1]]$strand[1])
    ALT = checkSense(var_parsed$variant, exons.ens.list[[1]]$strand[1])
    var_toAdd <- constructVCF_line(CHROM = CHROM, POS = POS, REF = REF, ALT = ALT)
    #print(paste(CHROM, POS, REF, ALT, sep = " "))
    #print(var_toAdd)
    return(var_toAdd)
    # fill if above returns false to find the right transcript via ccds
    
  }
  
  
  
  
  getCanonEnsTx <- function(vcf_table, ensembl, output_vcf_path)
  {
    
    genevars <- vcf_table[is.na(canon_tx),Gene_Variant]
    total <- length(genevars)
    i <- 1
    
    for (qgene in genevars)
    {
      qgene2 <- splitGeneVariant(qgene)
      print(paste0("Completed ", qgene2, ", ", i, " out of ", total))
      i<-i+1
      geneInfo <- retrieveGeneInfo(qgene2,ensembl)
      exons.ens <- geneInfo$exons
      transcripts.ens <- geneInfo$transcripts
      
      canon_transcript <- transcripts.ens[which(transcripts.ens$transcript_is_canonical ==1), "ensembl_transcript_id"]
      #print(qgene2)
      #print(canon_transcript)
      
      vcf_table[Gene_Variant == qgene, canon_tx := canon_transcript]
      
      
      write.table(vcf_table, file = output_vcf_path, row.names = FALSE, sep = "\t")
      
    }
    
    return(vcf_table)
  }
  
  
  
  #add PtID to vcf table
  
  #add in PtID to the start of INFO
  
  
  # Write the VCF file once verified working
  writeVCF <- function(vcf_table, file_path = "vcf_output.vcf")
  {
    vcf_var_names <- init_VCF_paste()
    
    
    vcf_table_to_paste <- as.data.table(vcf_table)[(CHROM != "NA"),.SD, .SDcols = names(vcf_var_names)]
    
    #create contigs chacarter
    contigs <- as.character()
    seq_rep <- fread("sequence_report.tsv", header = TRUE)
    for (i in 1:24)
    {
      contigs[i] <- paste0("##contig=<ID=",seq_rep[i,"Chromosome name"], ",length=", seq_rep[i,"Seq length"], ">")  
    }
    
    
    vcf_header <- c(
      "##fileformat=VCFv4.2",
      "##source=nERG dataset",
      "##reference=GRCh38.p14",
      '##INFO=<ID=SVTYPE,Number=.,Type=String,Description="Type of Structural Variant">',
      '##INFO=<ID=SVLEN,Number=.,Type=String,Description="Length of Structural Variant">',
      '##INFO=<ID=END,Number=.,Type=String,Description="End Position of Structural Variant">',
      '##INFO=<ID=PtID,Number=.,Type=String,Description="PatientID">',
      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
      '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
      '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
      contigs,
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA0001"
    )
    
    # File path for the VCF file
    
    # add contigs
    
    
    
    #vcfhead <- c(vcf_header,contigs)
    #writeLines(vcfhead, "vcf_head.vcf")
    
    
    
    # Write the header
    writeLines(vcf_header, file_path)
    
    #add pt ID to info part
    
    
    write.table(vcf_table_to_paste, file_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
    
  }
  
  
  #initialize dummy variable to populate with clinvar data
  init_VCF_paste <- function()
  {
    var_vcf <- data_frame(CHROM = character(),
                          POS = numeric(),
                          ID = character(),
                          REF = character(),
                          ALT = character(),
                          QUAL = character(),
                          FILTER = character(),
                          INFO = character(),
                          FORMAT = character(),
                          NA00001 = character())
    
  }
  
  # Function to continuously try to create VCF until successful
  run_create_vcf <- function(genevars, ensembl, melted_df, mdf2_dump, output_file) {
    success <- FALSE
    attempts <- 0
    max_attempts <- 10  # Set a reasonable limit to avoid infinite loops
    
    while (!success && attempts < max_attempts) {
      attempts <- attempts + 1
      tryCatch({
        # Attempt to load the latest saved VCF file
        saved_vcf_file <- fread("vcf_tbl_apr_16_2024.txt", header = TRUE)
        # Try running the function
        result <- CreateVCF(genevars, ensembl, melted_df, mdf2_dump = mdf2_dump, output_file, old_table = saved_vcf_file)
        
        # If CreateVCF finishes without error, set success to TRUE to exit the loop
        success <- TRUE
        return(result)
      }, error = function(e) {
        cat("Error encountered: ", e$message, "\n")
        cat("Attempt ", attempts, " failed. Retrying...\n")
        Sys.sleep(5)  # Wait for 5 seconds before retrying to avoid hammering the server
      })
    }
    
    if (!success) {
      stop("Failed to complete after ", max_attempts, " attempts.")
    }
  }
  
  
  
