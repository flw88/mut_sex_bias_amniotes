#!/usr/bin/env Rscript
rm(list=ls())

suppressMessages(library(stringr))
suppressMessages(library(ape))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

##### FUNCTIONS #####
IsEmptyChar <- function(x){
  return(length(x) == 0)
}

SplitMut <- function(mut.strs){
  # Split mutation type
  tmp <- tstrsplit(mut.strs, ">")
  return(cbind(tmp[[1]], tmp[[2]]))
}

### Complement sequence (e.g., TTC becomes AAG)
### NB This is NOT  reverse!
ComplementSeq <- function(in.seq){
  # Expects in.seq to be uppercase
  in.seq <- toupper(in.seq)
  in.seq <- gsub(in.seq, pattern="A", replacement="t")
  in.seq <- gsub(in.seq, pattern="C", replacement="g")
  in.seq <- gsub(in.seq, pattern="G", replacement="c")
  in.seq <- gsub(in.seq, pattern="T", replacement="a")
  return(toupper(in.seq))
}

### Reverse complement sequence (e.g., TTC becomes GAA)
RevComplementSeq <- function(in.seq){
  str.arr <- str_split(ComplementSeq(in.seq), "")[[1]]
  len <- length(str.arr)
  mid <- which(str.arr %in% ">")
  mid <- ifelse(length(mid) == 0, -1, mid)
  
  # Is in.seq is a valid mutation type
  is.mut.type <- (mid != 1) && (mid == (1+len)/2) && ((mid-1) == (len-mid))
  
  if(is.mut.type){
    out.str <- c(rev(str.arr[1:(mid-1)]), ">", rev(str.arr[(mid+1):len]))
    out.str <- str_c(out.str, collapse="")
  } else {
    out.str <- str_c(rev(str.arr), collapse="")
  }
  return( out.str )
}

GetRoot <- function(phy, force=FALSE){
  # Return NA if not rooted
  if(!is.rooted(phy) && !force){
    return( NA )
  } else {
    return( nodepath(phy)[[1]][1] )
  }
}

GetParent <- function(phy, node){
  r.node <- GetRoot(phy, force=TRUE)
  if(node == r.node){
    warning(sprintf("Node %d is top-most node (usually root) so no parent.\n", node))
    return(NA)
  } else {
    p.node <- nodepath(phy, from=node, to=r.node)[2]
    return(p.node)
  }
}

DistTips2Parent <- function(phy){
  dist.mat <- dist.nodes(phy) # Node distance matrix
  path.list <- nodepath(phy) # List of root to tip node paths
  
  # Get distance from tip to parent
  out.vals <- lapply(path.list, function(x){y <- tail(x, n=2); dist.mat[y[1], y[2]]})
  
  names(out.vals) <- phy$tip.label
  return(out.vals)
}

InitMultiPhylo <- function(base.phy, n=NULL, phy.names=NULL){
  if(is.null(n) && !is.null(phy.names)){
    n <- length(phy.names)
  } else if(!is.null(n) && !is.null(phy.names)){
    stopifnot(n == length(phy.names))
  } else if(is.null(n) && is.null(phy.names)){
    stop("One of n or phy.names arguments must be provided")
  }
  
  if(n < 1){
    stop("n must be integer greater than 0")
  }
  
  mult.phy <- base.phy; i <- 2
  while(i <= n){
    mult.phy <- c(mult.phy, base.phy)
    i <- i + 1
  }
  
  mult.phy <- .compressTipLabel(mult.phy)
  
  if(!is.null(phy.names)){
    names(mult.phy) <- phy.names
  }
  
  return(mult.phy)
}

SliceMultiPhylo <- function(in.phy, patt=c(), complement=TRUE){
  # Patterns should be something like "AT>CT"
  # Wildcards can be specified like "A.>AT" which should match
  #   any of "AA>AT", "AT>AT", "AC>AT", "AG>AT"
  # Anti-matching can be specified like "C[^G]>T[^G]" which should match
  #   any of "CA>TA", "CT>TT", "CC>TC"
  # complement: return the complement nucleotides for the given
  #   patterns as well
  nms <- names(in.phy)
  mask <- rep(FALSE, length(in.phy))
  for(p in patt){
    mask <- mask | grepl(nms, pattern=p)
    if(complement){
      mask <- mask | grepl(nms, pattern=RevComplementSeq(p))
    }
  }
  
  return(in.phy[nms[mask]])
}

SumMultiPhylo <- function(in.phy){
  # in.phy should be an object of class multiPhylo.
  # Sum the corresponding branch lengths of each tree
  # in the multiPhylo to produce a single tree
  out.phy <- NULL
  for(i in 1:length(in.phy)){
    cur.phy <- in.phy[[i]]
    if(is.null(out.phy)){
      out.phy <- cur.phy
    } else {
      out.phy$edge.length <- out.phy$edge.length + cur.phy$edge.length
    }
  }
  return(out.phy)
}

PhyloArith <- function(a.phy, b.phy, operand="+"){
  out.phy <- a.phy
  if(operand == "+"){
    out.phy$edge.length <- out.phy$edge.length + b.phy$edge.length
  } else if(operand == "-"){
    out.phy$edge.length <- out.phy$edge.length - b.phy$edge.length
  } else if(operand == "*"){
    out.phy$edge.length <- out.phy$edge.length * b.phy$edge.length
  } else if(operand == "/"){
    out.phy$edge.length <- out.phy$edge.length / b.phy$edge.length
  } else {
    stop(sprintf("operand '%s' not recognized\n", operand))
  }
  return(out.phy)
}

ReadExptotsub <- function(fn, base.phy){
  # cat(str_interp("Reading ${fn}\n"))
  ### INTERNAL FUNCTIONS
  ParseBranchDesc <- function(str.val, base.phy){
    # Branch description is usually of the form of either:
    #   "Branch above node 5 (leaf labeled 'hg38'):"
    #      or
    #   "Branch above node 4 (leaf labeled 'hg38-ponAbe2'):"
    # We want the corresponding edge in base.phy
    # Returns the index of the edge in base.phy
    
    split.str <- str_split(str.val, " ")[[1]]
    node.str <- split.str[which(split.str == "labeled") + 1]
    for(patt in c("\\)", "'", ":")){
      node.str <- gsub(node.str, pattern=patt, replacement="")
    }
    
    is.tip <- node.str %in% as.character(base.phy$tip.label) #!grepl(node.str, pattern="-")
    if(is.tip){ # Is tip
      node.n <- which(base.phy$tip.label == node.str)
    } else { # Is an inner node
      sp.pair <- str_split(node.str, "-")[[1]]
      if(all(sp.pair %in% as.character(base.phy$tip.label))){
        node.n <- getMRCA(base.phy, tip=sp.pair)
      } else { # Else it's a named node
        node.n <- which(base.phy$node.label == node.str) + Ntip(base.phy)
        if(length(node.n) == 0){
          stop(sprintf("Unable to determine node from string %s\n", node.str))
        }
      }
    }
    
    par.n <- GetParent(base.phy, node.n) # Parental node
    
    edge.index <- which( ((base.phy$edge[,1] == node.n) & (base.phy$edge[,2] == par.n)) |
                           ((base.phy$edge[,2] == node.n) & (base.phy$edge[,1] == par.n))   )
    return(edge.index)
  }
  
  IsBlockStart <- function(cur.line){
    block.start.patt <- "^Branch above"
    out.bool <- grepl(cur.line, pattern=block.start.patt)
    return(out.bool)
  }
  
  ReadBlock <- function(con){
    # Read the current table of mutations
    
    is.hdr <- TRUE
    cur.line <- ""
    all.lines <- c()
    while(!IsEmptyChar(cur.line) && !IsBlockStart(cur.line)){
      if(cur.line != ""){
        if(is.hdr){
          cur.line <- paste("rn", cur.line, sep=" ")
          is.hdr <- FALSE
        }
        # cat(cur.line, "\n", sep="", file=tmp.fn, append=TRUE)
        all.lines <- c(all.lines, cur.line)
      }
      
      cur.line <- readLines(con, n=1)
    }
    
    tmp <- fread(text=all.lines, sep="auto", header="auto")
    rn <- tmp[, rn]
    tmp[, rn := NULL]
    dat <- as.data.frame(tmp)
    rownames(dat) <- rn
    
    return(list("dat"=dat, "cur.line"=cur.line))
  }
  
  ### READ THROUGH FILE BLOCK BY BLOCK
  con <- file(fn, "r") # Read from stdin
  
  # Find first block start
  cur.line <- ""
  while(!IsEmptyChar(cur.line) && !IsBlockStart(cur.line)){
    cur.line <- readLines(con, n=1)
  }
  
  if(IsEmptyChar(cur.line)){
    close(con)
    stop(sprintf("File %s doesn't look like a phyloFit .exptotsub file.\n", fn))
  }
  
  # Compile substitution count tables for each edge
  sub.tabs <- vector(length=nrow(base.phy$edge), mode="list") # (List of data.tables)
  
  while(!IsEmptyChar(cur.line)){
    i <- ParseBranchDesc(cur.line, base.phy) # Which edge are we considering?
    x <- ReadBlock(con)
    sub.tabs[[i]] <- x$dat
    cur.line <- x$cur.line
  }
  
  close(con)
  
  ### Convert edge substitutions into multiphylo object
  phy.names <- apply(expand.grid(rownames(sub.tabs[[1]]), colnames(sub.tabs[[1]])), 1,
                     function(x) paste(x, collapse=">"))
  count.phy <- InitMultiPhylo(base.phy, phy.names=phy.names)
  
  for(m in names(count.phy)){ # For mutation type 'm'
    cur.edge.len <- rep(0.0, length(sub.tabs))
    for(i in 1:length(sub.tabs)){ # For edge 'i'
      m.pair <- SplitMut(m)
      x <- m.pair[1]; y <- m.pair[2]
      
      cur.edge.len[i] <- sub.tabs[[i]][x, y]
    }
    
    count.phy[[m]][["edge.length"]] <- cur.edge.len
  }
  
  return(count.phy)
}

CountSubstitutions <- function(base.phy, u2s.fn=NULL, unr.fn=NULL){
  CpgStatus <- function(s){
    if(str_detect(s, pattern="^noncpg_")){
      x <- 0
    } else if(str_detect(s, pattern="^cpg_")){
      x <- 1
    } else { # CpG status doesn't apply
      x <- -1
    }
    return(x)
  }
  
  
  if(!is.null(u2s.fn)){
    two.phy <- ReadExptotsub(u2s.fn, base.phy=base.phy) # Dinucleotide
  }
  if(! is.null(unr.fn)){
    one.phy <- ReadExptotsub(unr.fn, base.phy=base.phy) # Single nucleotide
  }
  
  if(is.null(u2s.fn) && is.null(unr.fn)){
    stop("At least one of 'u2s.fn' or 'unr.fn' must be provided.")
  }
  
  if(is.null(unr.fn)){ # U2S only
    # Denominators
    den.phy <- InitMultiPhylo(base.phy, phy.names = c("C", "cpg_C", "noncpg_C", "T"))
    patt <- c("C.>..", ".C>..", "G.>..", ".G>..")
    tmp <- SumMultiPhylo(SliceMultiPhylo(two.phy, patt[1], complement=FALSE))
    for(i in 2:length(patt)){
      tmp <- PhyloArith(tmp, SumMultiPhylo(SliceMultiPhylo(two.phy, patt[i], complement=FALSE)), "+")
    }
    den.phy[["C"]] <- tmp
    
    patt <- c("A.>..", ".A>..", "T.>..", ".T>..")
    tmp <- SumMultiPhylo(SliceMultiPhylo(two.phy, patt[1], complement=FALSE))
    for(i in 2:length(patt)){
      tmp <- PhyloArith(tmp, SumMultiPhylo(SliceMultiPhylo(two.phy, patt[i], complement=FALSE)), "+")
    }
    den.phy[["T"]] <- tmp
    
    # den.phy[["C"]] <- SumMultiPhylo(SliceMultiPhylo(one.phy, "C>."))
    # den.phy[["T"]] <- SumMultiPhylo(SliceMultiPhylo(one.phy, "T>."))
    cg.phy <- SliceMultiPhylo(two.phy, "CG>..")
    den.phy[["cpg_C"]] <- SumMultiPhylo(c(cg.phy, cg.phy))
    den.phy[["noncpg_C"]] <- PhyloArith(den.phy[["C"]], den.phy[["cpg_C"]], "-")
    
    den.phy[["S"]] <- den.phy[["C"]]
    den.phy[["W"]] <- den.phy[["T"]]
    den.phy[["cpg_S"]]    <- den.phy[["cpg_C"]]
    den.phy[["noncpg_S"]] <- den.phy[["noncpg_C"]]
    
    den.phy[["all"]] <- PhyloArith(den.phy[["C"]], den.phy[["T"]], "+")
    
    # Numerators
    sub.types <- c("C>A", "C>G", "C>T",
                   "T>A", "T>C", "T>G",
                   "cpg_C>A", #"noncpg_C>A",
                   "cpg_C>G", #"noncpg_C>G",
                   "cpg_C>T")#, "noncpg_C>T")
    
    num.phy <- InitMultiPhylo(base.phy, phy.names=sub.types)
    for(s in sub.types){
      a.base <- str_sub(str_extract(s, ".>."), 1, 1) # Ancestral base
      m.base <- str_sub(str_extract(s, ".>."), 3, 3) # Mutant base
      if(CpgStatus(s) == 1){ # CpG
        patt <- sprintf("CG>%sG", m.base)
        num.phy[[s]] <- SumMultiPhylo(SliceMultiPhylo(two.phy, patt, complement=TRUE))
      } else if(CpgStatus(s) == -1){
        patt <- c(sprintf("%s.>%s.", a.base, m.base), 
                  sprintf(".%s>.%s", a.base, m.base),
                  sprintf("%s.>%s.", ComplementSeq(a.base), ComplementSeq(m.base)), 
                  sprintf(".%s>.%s", ComplementSeq(a.base), ComplementSeq(m.base)))
        
        tmp <- SumMultiPhylo(SliceMultiPhylo(two.phy, patt[1], complement=FALSE))
        for(i in 2:length(patt)){
          tmp <- PhyloArith(tmp, SumMultiPhylo(SliceMultiPhylo(two.phy, patt[i], complement=FALSE)), "+")
        }
        num.phy[[s]] <- tmp
      }
    }
    for(s in c("C>A", "C>G", "C>T")){
      m.base <- str_sub(s, 3, 3)
      
      patt <- sprintf("CG>%s.", m.base)
      sub.phy <- SumMultiPhylo(SliceMultiPhylo(two.phy, patt, complement=TRUE))
      
      # sub.phy <- num.phy[[str_c("cpg_", s)]]
      
      num.phy[[str_c("noncpg_", s)]] <- PhyloArith(num.phy[[s]], sub.phy, "-")
    }
    
    num.phy[["cpg_S>S"]]    <- num.phy[["cpg_C>G"]]
    num.phy[["noncpg_S>S"]] <- num.phy[["noncpg_C>G"]]
    num.phy[["cpg_S>W"]]    <- PhyloArith(num.phy[["cpg_C>A"]], num.phy[["cpg_C>T"]], "+")
    num.phy[["noncpg_S>W"]] <- PhyloArith(num.phy[["noncpg_C>A"]], num.phy[["noncpg_C>T"]], "+")
    num.phy[["W>S"]] <-        PhyloArith(num.phy[["T>C"]], num.phy[["T>G"]], "+")
    num.phy[["W>W"]] <-        num.phy[["T>A"]]
    
    num.phy[["T_ts"]]        <- num.phy[["T>C"]]
    num.phy[["T_tv"]]        <- PhyloArith(num.phy[["T>A"]], num.phy[["T>G"]], "+")
    num.phy[["C_ts"]]        <- num.phy[["C>T"]]
    num.phy[["C_tv"]]        <- PhyloArith(num.phy[["C>A"]], num.phy[["C>G"]], "+")
    num.phy[["cpg_C_ts"]]    <- num.phy[["cpg_C>T"]]
    num.phy[["cpg_C_tv"]]    <- PhyloArith(num.phy[["cpg_C>A"]], num.phy[["cpg_C>G"]], "+")
    num.phy[["noncpg_C_ts"]] <- num.phy[["noncpg_C>T"]]
    num.phy[["noncpg_C_tv"]] <- PhyloArith(num.phy[["noncpg_C>A"]], num.phy[["noncpg_C>G"]], "+")
    
    tmp1                     <- PhyloArith(num.phy[["T_ts"]], num.phy[["T_tv"]], "+")
    tmp2                     <- PhyloArith(num.phy[["C_ts"]], num.phy[["C_tv"]], "+")
    num.phy[["all"]]         <- PhyloArith(tmp1, tmp2, "+")
    
  } else if(is.null(u2s.fn)){ # UNR only
    # Denominators
    den.phy <- InitMultiPhylo(base.phy, phy.names = c("C", "T"))
    den.phy[["C"]] <- SumMultiPhylo(SliceMultiPhylo(one.phy, "C>.", complement=TRUE))
    den.phy[["T"]] <- SumMultiPhylo(SliceMultiPhylo(one.phy, "T>.", complement=TRUE))
    
    # den.phy[["S"]] <- den.phy[["C"]]
    # den.phy[["W"]] <- den.phy[["T"]]
    
    # den.phy[["all"]] <- PhyloArith(den.phy[["C"]], den.phy[["T"]], "+")
    
    # Numerators
    sub.types <- c("C>A", "C>G", "C>T",
                   "T>A", "T>C", "T>G")
    
    num.phy <- InitMultiPhylo(base.phy, phy.names=sub.types)
    for(s in sub.types){ 
      num.phy[[s]] <- SumMultiPhylo(SliceMultiPhylo(one.phy, s, complement=TRUE))
    }
    
    ### FW: These types aren't really useful for the current use case
    # num.phy[["W>S"]]  <- PhyloArith(num.phy[["T>C"]], num.phy[["T>G"]], "+")
    # num.phy[["W>W"]]  <- num.phy[["T>A"]]
    # 
    # num.phy[["T_ts"]] <- num.phy[["T>C"]]
    # num.phy[["T_tv"]] <- PhyloArith(num.phy[["T>A"]], num.phy[["T>G"]], "+")
    # num.phy[["C_ts"]] <- num.phy[["C>T"]]
    # num.phy[["C_tv"]] <- PhyloArith(num.phy[["C>A"]], num.phy[["C>G"]], "+")
    # 
    # num.phy[["all"]] <- num.phy[[sub.types[1]]]
    # for(s in sub.types[-1]){
    #   num.phy[["all"]]  <- PhyloArith(num.phy[["all"]], num.phy[[s]], "+")  
    # }
    
  } else { # Both U2S and UNR files provided
    # Denominators
    den.phy <- InitMultiPhylo(base.phy, phy.names = c("C", "cpg_C", "noncpg_C", "T"))
    den.phy[["C"]] <- SumMultiPhylo(SliceMultiPhylo(one.phy, "C>."))
    den.phy[["T"]] <- SumMultiPhylo(SliceMultiPhylo(one.phy, "T>."))
    cg.phy <- SliceMultiPhylo(two.phy, "CG>..")
    den.phy[["cpg_C"]] <- SumMultiPhylo(c(cg.phy, cg.phy))
    den.phy[["noncpg_C"]] <- PhyloArith(den.phy[["C"]], den.phy[["cpg_C"]], "-")
    
    # den.phy[["S"]] <- den.phy[["C"]]
    # den.phy[["W"]] <- den.phy[["T"]]
    # den.phy[["cpg_S"]]    <- den.phy[["cpg_C"]]
    # den.phy[["noncpg_S"]] <- den.phy[["noncpg_C"]]
    # 
    # den.phy[["all"]] <- PhyloArith(den.phy[["C"]], den.phy[["T"]], "+")
    
    # Numerators
    sub.types <- c("C>A", "C>G", "C>T",
                   "T>A", "T>C", "T>G",
                   "cpg_C>A", "noncpg_C>A",
                   "cpg_C>G", "noncpg_C>G",
                   "cpg_C>T", "noncpg_C>T")
    
    num.phy <- InitMultiPhylo(base.phy, phy.names=sub.types)
    for(s in sub.types){
      x.base <- gsub(s, pattern="^.+>", replacement="") # Mutant base
      if(CpgStatus(s) == 1){
        x.patt <- sprintf("CG>%sG", x.base)
        num.phy[[s]] <- SumMultiPhylo(SliceMultiPhylo(two.phy, x.patt))
      } else if(CpgStatus(s) == 0){
        # Take total number of C>X and subtract
        #   (CG>X. + CG>.Y)
        #   where Y is complement base of X
        y.base <- ComplementSeq(x.base)
        x.patt <- sprintf("CG>%s.", x.base)
        y.patt <- sprintf("CG>.%s", y.base)
        
        tot.phy <- SumMultiPhylo(SliceMultiPhylo(one.phy, sprintf("C>%s", x.base)))
        sub.phy <- PhyloArith(SumMultiPhylo(SliceMultiPhylo(two.phy, x.patt, complement=FALSE)),
                              SumMultiPhylo(SliceMultiPhylo(two.phy, y.patt, complement=FALSE)), "+")
        num.phy[[s]] <- PhyloArith(tot.phy, sub.phy, "-")
      } else if(CpgStatus(s) == -1){
        num.phy[[s]] <- SumMultiPhylo(SliceMultiPhylo(one.phy, s))
      }
    }
  }
  
  return(list(num=num.phy, den=den.phy))
}

# Find the maximum replicate 
# Note that fmt.str should include a '${i}' somewhere to indicate where the 
MaxReplicate <- function(fmt.str, n=6){
  fns <- sapply(1:6, function(i) str_interp(fmt.str))
  if(any(!file.exists(fns))){
    warning(str_c(fmt.str, ":\nMissing mod files, returning NA."))
    return(NA)
  }
  max.lhood <- -Inf
  for(i in 1:n){
    cur.mod <- str_interp(fmt.str)
    lhood <- fread(cmd=str_interp("grep TRAINING_LNL ${cur.mod}"))$V2
    if(lhood > max.lhood){
      max.lhood <- lhood
      max.i <- i
    }
  }
  
  return(max.i)
}

##### MAIN #####
# main.dir <- "/Users/felix/mt_mp_lab/projects/male_mutation_bias_XA/"
main.dir <- "/moto/palab/projects/male_mutation_bias_XA/"

region.dat <- fread(str_c(main.dir, "mut_sex_bias_amniotes/scripts/region_files/Homo_sapiens.1Mb.bed"),
                    col.names=c("chrom", "start", "end"))


# Keep only autosomes
region.dat <- region.dat[(chrom != "chrX") & (chrom != "chrY")]

base.phy <- read.tree(str_interp("${main.dir}/mut_sex_bias_amniotes/trees/Mammals.nwk"))

counts.fn <- str_interp("${main.dir}/mut_sex_bias_amniotes/scripts/merged_features/microcebus_subrate.RData")

if(file.exists(counts.fn)){
  cat(str_interp("Loading counts from ${counts.fn}\n\n"))
  load(counts.fn)
} else {
  counts <- NULL
  fn <- c("unr"="", "u2s"="")

  for(j in 1:region.dat[,.N]){ # Iterate over chunks
    cur.chr <- region.dat[j, chrom]
    chunk <- str_c(region.dat[j, str_c(start, ".", end)])
    
    skip.chunk <- FALSE
    for(nm in names(fn)){
      if(nm == "unr"){
        midfix <- "."
      } else {
        midfix <- ".U2S."
      }
      fmt.str <- str_c(main.dir, "MAFs_2021/Homo_sapiens/", cur.chr, "/", cur.chr, ".", chunk,".Mammals.${i}.filtered.nonCGI", midfix, "mod")
      
      # Get the maximum replicate
      i <- MaxReplicate(fmt.str)
      if(is.na(i)){ skip.chunk <- TRUE; break } # Skip this chunk if not all mod files are present
          
      # Get exptotsub filename
      fn[nm] <- str_replace(str_interp(fmt.str), "\\.mod$", ".exptotsub")
    }
    if(skip.chunk){
      # Toss chunk if files don't exist
      next
    }

    cat(str_interp("Reading counts from ${fn[1]} and ${fn[2]}\n"))
    cur.counts <- CountSubstitutions(base.phy, unr.fn=fn["unr"], u2s.fn=fn["u2s"])

    # Sum chunks
    if(is.null(counts)){
      counts <- cur.counts
    } else {
      for(k in names(counts)){
        for(mt in names(counts[[k]])){
          counts[[k]][[mt]] <- PhyloArith(counts[[k]][[mt]], cur.counts[[k]][[mt]], "+")
        }
      }
    }
  }

  cat(str_interp("Saving counts to ${counts.fn}\n\n"))
  save(counts, file=counts.fn)
}




##### COMPARE SUM AND ALL #####
# plt.dat <- list()
# for(mt in names(dat$sum_windows$counts$num)){
#   plt.dat[[mt]] <- data.table("SUM" = dat$sum_windows$counts$num[[mt]]$edge.length,
#                               "CAT" = dat$all_windows$counts$num[[mt]]$edge.length,
#                               "TYPE"= mt)
# }
# plt.dat <- rbindlist(plt.dat)

# theme_set(theme_bw() +
#             theme(axis.text    = element_text(size=12), 
#                   panel.border = element_rect(size = 1.5)))

# p <- ggplot(aes(x=CAT, y=SUM), data=plt.dat) + geom_point() + 
#   geom_smooth(method="lm") + geom_abline(intercept=0, slope=1, linetype=2, color="gray") +
#   facet_wrap(~ TYPE, nrow=3)
# ggsave(p, filename=str_interp("${main.dir}/mut_sex_bias_amniotes/scripts/pdfs/microcebus.cat_vs_sum.pdf"), width=12, height=9)

##### GET SPLIT TIME DATA #####
split.sp <- c("Microcebus_murinus", "Callithrix_jacchus")
split.time <- unlist(DistTips2Parent(keep.tip(read.tree(str_interp("${main.dir}/mut_sex_bias_amniotes/trees/mammals241.TimeTree.nwk")), split.sp)))
split.time <- split.time * 1e6

##### CALCULATE YEARLY MUTATION RATE #####

# sub.all <- unlist(DistTips2Parent(keep.tip(dat$all_windows$unr$mod.phy, split.sp)))
# cat("Yearly mutation rate (all sub types):\n")
# print(sub.all / split.time[names(sub.all)])
# cat("\n")1

##### CALCULATE YEARLY C>T MUTATION RATE #####
sub.cpg <- unlist(DistTips2Parent(keep.tip(counts$num$`cpg_C>T`, split.sp))) /
              unlist(DistTips2Parent(keep.tip(counts$den$`cpg_C`, split.sp)))

sub.noncpg <- unlist(DistTips2Parent(keep.tip(counts$num$`noncpg_C>T`, split.sp))) /
                 unlist(DistTips2Parent(keep.tip(counts$den$`noncpg_C`, split.sp)))

cat("Yearly mutation rate (CpG C>T):\n")
print(sub.cpg    / split.time[names(sub.cpg)])
cat("\n")

cat("Yearly mutation rate (non-CpG C>T):\n")
print(sub.noncpg / split.time[names(sub.noncpg)])
cat("\n")

cat("Rate ratio of CpG C>T vs. non-CpG C>T):\n")
print(sub.cpg    / sub.noncpg)
cat("\n")



