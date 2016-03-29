# The following scripts are required for implementation of PyClone analysis
# scipt N.L.McGranahan 

# updates needed:
# include sex chromosomes

# script log:
# 27/04/15 scripting begins




identify.subclonal.mut.copy.number.ascat <- function(x,sub.mat.mut, sub.mat.copy,region,sample,sex='male')
{
  
  
  mut                  <- sub.mat.mut[x,,drop=FALSE]
  ww                   <- which(as.numeric(sub.mat.copy$chr)==as.numeric(mut$chr)
                                &as.numeric(sub.mat.copy$startpos)<=as.numeric(mut$start)
                                &as.numeric(sub.mat.copy$endpos)>=as.numeric(mut$stop))
  copy                 <- sub.mat.copy[ww,,drop=FALSE]
  
  
  mutation_id         <- paste(sample,mut$chr,mut$start,mut$ref,sep=":")
  ref_counts          <- mut[,gsub("-","\\.",paste(region, ".ref_count",sep=""))]
  var_counts          <- mut[,gsub("-","\\.",paste(region, ".var_count",sep=""))]    

  
  normal_cn           <- 2
  region              <- region
  Reference_Base      <- mut$ref
  Alternate_Base      <- mut$var
  
  if(nrow(copy)!=1)
  {
    minor_cn          <- NA
    major_cn          <- NA
    major_raw         <- NA
    minor_raw         <- NA
    fracA             <- NA
    nMaj_A            <- NA
    nMin_A            <- NA
    fracB             <- NA
    nMaj_B            <- NA
    nMin_B            <- NA
    fracC             <- NA
    nMaj_C            <- NA
    nMin_C            <- NA
    fracD             <- NA
    nMaj_D            <- NA
    nMin_D            <- NA  
    
    output            <-  data.frame(mutation_id
                                     , ref_counts
                                     , var_counts
                                     , normal_cn
                                     , minor_cn
                                     , major_cn
                                     , major_raw
                                     , minor_raw
                                     , region
                                     , Reference_Base
                                     , Alternate_Base
                                     , fracA            
                                     , nMaj_A           
                                     , nMin_A           
                                     , fracB            
                                     , nMaj_B           
                                     , nMin_B           
                                     , fracC            
                                     , nMaj_C           
                                     , nMin_C           
                                     , fracD            
                                     , nMaj_D           
                                     , nMin_D           
                                     , stringsAsFactors=FALSE)
    return(output)
  }
  
  minor_cn            <- min(c(copy$nMajor,copy$nMinor))
  major_cn            <- max(c(copy$nMajor,copy$nMinor))
  major_raw           <- c(copy$nAraw)
  minor_raw           <- c(copy$nBraw)
  fracA               <- as.numeric(copy$fracA)
  nMaj_A              <- as.numeric(copy$nMaj_A)
  nMin_A              <- as.numeric(copy$nMin_A)
  fracB               <- as.numeric(copy$fracB)
  nMaj_B              <- as.numeric(copy$nMaj_B)
  nMin_B              <- as.numeric(copy$nMin_B)
  fracC               <- as.numeric(copy$fracC)
  nMaj_C              <- as.numeric(copy$nMaj_C)
  nMin_C              <- as.numeric(copy$nMin_C)
  fracD               <- as.numeric(copy$fracD)
  nMaj_D              <- as.numeric(copy$nMaj_D)
  nMin_D              <- as.numeric(copy$nMin_D)
  
  output            <-  data.frame(mutation_id
                                   , ref_counts
                                   , var_counts
                                   , normal_cn
                                   , minor_cn
                                   , major_cn
                                   , major_raw
                                   , minor_raw
                                   , region
                                   , Reference_Base
                                   , Alternate_Base
                                   , fracA            
                                   , nMaj_A           
                                   , nMin_A           
                                   , fracB            
                                   , nMaj_B           
                                   , nMin_B           
                                   , fracC            
                                   , nMaj_C           
                                   , nMin_C           
                                   , fracD            
                                   , nMaj_D           
                                   , nMin_D           
                                   , stringsAsFactors=FALSE)
  
  return(output)
  
  
  
  
}



subclonalDissection <- function(region
                                ,complete.mutation.table
                                ,purity
                                ,region.seg.copy
                                ,order.by.pos=TRUE
                                ,min.subclonal=0.1)
{
  
  complete.mutation.table <- data.frame(complete.mutation.table,stringsAsFactors=FALSE)
  mut.order               <- complete.mutation.table$mutation_id
  if(order.by.pos)
  {
    chr <- do.call(rbind,strsplit(as.character(complete.mutation.table$mutation_id),split=":"))
    complete.mutation.table  <- complete.mutation.table[order(as.numeric(as.character(chr[,2])),as.numeric(as.character(chr[,3]))),]
    mut.order                <- complete.mutation.table$mutation_id
    
  }
  
  
  suppressPackageStartupMessages(library(sequenza))
  suppressPackageStartupMessages(library(bootstrap))
  suppressPackageStartupMessages(library(boot))
  
  # Get the table ready, with only information for specific patient
  pyClone.tsv  <- complete.mutation.table[complete.mutation.table$region==region,,drop=FALSE]
  row.names    <- pyClone.tsv$mutation_id
  cellularity  <- as.numeric(purity)
  
  
  f.function <- function (c,purity,local.copy.number)
  {
    
    return(min(c((purity*c) / (2*(1-purity) + purity*local.copy.number)),1))
    
  }
  
  get.mut.mult <- function(CNt,Vaf,cellularity,CNn)
  {
    
    return((Vaf *1/cellularity)*((cellularity*CNt)+CNn*(1-cellularity)))
    
  }
  
  pyClone.tsv$expected.VAF       <- f.function(rep(1,nrow(pyClone.tsv)),cellularity,unlist(pyClone.tsv$major_raw)+unlist(pyClone.tsv$minor_raw))
  
  
  # let's put it back together
  rownames(pyClone.tsv) <- pyClone.tsv$mutation_id
  pyClone.tsv           <- pyClone.tsv[unlist(mut.order),]
  pyClone.tsv$obs.VAF   <- as.numeric(pyClone.tsv$var_counts)/(as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts))
  
  get.mutCopyNum<- function(i)
  {
    
    # let's simplify things, and assume there are max two subclonal populations
    mutCopyNumber <- get.mut.mult(CNt=(as.numeric(pyClone.tsv$major_raw[i])+as.numeric(pyClone.tsv$minor_raw[i])),Vaf=pyClone.tsv$obs.VAF[i],cellularity=cellularity,CNn=unlist(pyClone.tsv$normal_cn[i]))
    return(mutCopyNumber)
    
    
  }
  get.absolute.ccf<- function(i)
  {
    #print(i)
    absolute.cancer.cell.fraction <- function(n.alt, depth, purity, local.copy.number)
    {
      f.function <- function (c,purity,local.copy.number)
      {
        
        return(min(c((purity*c) / (2*(1-purity) + purity*local.copy.number),1)))
        
      }
      x              <- dbinom(n.alt,depth, prob=sapply(seq(0.01,1,length.out=100),f.function,purity,local.copy.number))
      if(min(x)==0)
      {
        x[length(x)] <- 1
      }
      
      names(x)       <- seq(0.01,1,length.out=100)
      sub.cint <- function(x, prob = 0.95,n.alt,depth) {
        xnorm   <- x/sum(x)
        xsort   <- sort(xnorm, decreasing = TRUE)
        xcumLik <- cumsum(xsort)
        n = sum(xcumLik < prob) + 1
        LikThresh <- xsort[n]
        cint  <- x[xnorm >= LikThresh]
        all   <- as.numeric(names(x))
        cellu <- as.numeric(names(cint))
        l.t   <- cellu[1]
        r.t   <- cellu[length(cellu)]
        m     <- cellu[which.max(cint)]
        
        prob.subclonal <- sum(xnorm[1:90])# 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='less')$p.val
        prob.clonal    <- sum(xnorm[91:100]) # 1-prop.test(n.alt,depth,p=f.function(1,purity,local.copy.number),alternative='greater')$p.val
        
        data.frame(left = l.t, est = m, right = r.t,prob.subclonal=prob.subclonal,prob.clonal=prob.clonal)
      }
      
      
      return(sub.cint(x,n.alt=n.alt,depth=depth))
      
      
    }
    
    absolute.calc        <- absolute.cancer.cell.fraction(n.alt=unlist(pyClone.tsv$var_counts)[i],depth=(as.numeric(pyClone.tsv$ref_counts[i])+as.numeric(pyClone.tsv$var_counts[i])),purity=cellularity,local.copy.number=(as.numeric(pyClone.tsv$major_raw[i])+as.numeric(pyClone.tsv$minor_raw[i])))
    absolute.ccf.0.05    <- absolute.calc[1]
    absolute.ccf.0.95    <- absolute.calc[3]
    absolute.ccf         <- absolute.calc[2]
    prob.subclonal       <- absolute.calc[4]
    prob.clonal          <- absolute.calc[5]
    return(cbind(absolute.ccf,absolute.ccf.0.05,absolute.ccf.0.95,absolute.ccf,prob.subclonal,prob.clonal))
    
    
    
    
  }
  
  pyClone.tsv$mutCopyNum        <- sapply(1:nrow(pyClone.tsv),get.mutCopyNum)
  absolute.ccfs                 <- t(sapply(1:nrow(pyClone.tsv),get.absolute.ccf))
  absolute.ccfs                 <- data.frame(absolute.ccfs,stringsAsFactors=FALSE)
  
  pyClone.tsv$absolute.ccf        <- as.numeric(absolute.ccfs$est)
  pyClone.tsv$absolute.ccf.0.05   <- as.numeric(absolute.ccfs$left)
  pyClone.tsv$absolute.ccf.0.95   <- as.numeric(absolute.ccfs$right)
  pyClone.tsv$phyloCCF            <- as.numeric(pyClone.tsv$mutCopyNum)
  pyClone.tsv$phyloCCF.0.05       <- sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,CNtumor,cellularity,CNn,i){
    VAF <- prop.test(var.count[i]
                     ,depth[i]
                     ,e[i]
    )$conf.int[1]
    return(get.mut.mult(CNt = CNtumor[i],Vaf = VAF,cellularity = cellularity,CNn = CNn))
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e = as.numeric(pyClone.tsv$expected.VAF)
  ,CNtumor=as.numeric(pyClone.tsv$major_raw)+as.numeric(pyClone.tsv$minor_raw)
  ,cellularity=cellularity
  ,CNn=2)
  pyClone.tsv$phyloCCF.0.95      <- sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,CNtumor,cellularity,CNn,i){
    VAF <- prop.test(var.count[i]
                     ,depth[i]
                     ,e[i]
    )$conf.int[2]
    return(get.mut.mult(CNt = CNtumor[i],Vaf = VAF,cellularity = cellularity,CNn = CNn))
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e = as.numeric(pyClone.tsv$expected.VAF)
  ,CNtumor=as.numeric(pyClone.tsv$major_raw)+as.numeric(pyClone.tsv$minor_raw)
  ,cellularity=cellularity
  ,CNn=2)
  
  
  pyClone.tsv$no.chrs.bearing.mut <- 1
  
  # let's check which subclonal mutations can be explained by copy number
  #pyClone.tsv.subclonal         <- pyClone.tsv[which(pyClone.tsv$absolute.ccf.0.95<1&pyClone.tsv$mutCopyNum>0.01),]
  subclonal.mutations            <- which(pyClone.tsv$absolute.ccf.0.95<1&pyClone.tsv$mutCopyNum>0.01)
  subclonal.mutations            <- which(absolute.ccfs$prob.subclonal>0.5&pyClone.tsv$mutCopyNum>0.01)
  pyClone.tsv$whichFrac          <- NA
  
  
  
  if(length(subclonal.mutations)>0)
  {
    
    # assume subclonal muts are on one chromosome copy
    # therefore mutation copy number must be subclonal fraction 
    # of the higher CN subclone (i.e. lost in lower CN subclone) 
    # or 1 (i.e. present in both subclones)
    for (a in subclonal.mutations)
    {
      #print(a)
      mut.info <- pyClone.tsv[a,,drop=FALSE]
      # check whether the fracA is hundred%
      if (mut.info$fracA==1)
      {
        pyClone.tsv$whichFrac[a] <- c('A,B')
        next;
      }
      
      # let's see what happens when we only have A and B
      if (is.na(mut.info$fracC))
      {
        # determine which fraction has a copy number loss
        possible.subclonal.fractions <- c(1)
        if(unlist(mut.info$nMaj_A)>unlist(mut.info$nMaj_B))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracA))
        }
        if(unlist(mut.info$nMaj_B)>unlist(mut.info$nMaj_A))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracB))
        }
        if(unlist(mut.info$nMin_B)>unlist(mut.info$nMin_A))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracB))
        }
        if(unlist(mut.info$nMin_A)>unlist(mut.info$nMin_B))
        {
          possible.subclonal.fractions <- c(possible.subclonal.fractions,unlist(mut.info$fracA))
        }
        
        
        best.CN                      <- possible.subclonal.fractions[which.min(abs(mut.info$mutCopyNum/possible.subclonal.fractions - 1))]
        if(best.CN==1)
        {
          pyClone.tsv$whichFrac[a] <- c('A,B')
          next;
        }
        var.count       = as.numeric(mut.info$var_counts)
        depth.count     =  as.numeric(mut.info$var_counts)+as.numeric(mut.info$ref_counts)
        expected.prop   = pyClone.tsv$expected.VAF[a] * best.CN
        
        # check whether subclonal CN results in clonal mutation
        # otherwise subclonal CN doesn't explain subclonal MCN
        if(best.CN != 1 & prop.test(var.count,depth.count*purity,expected.prop/purity,alternative = 'less')$p.value > 0.01)
        {
          pyClone.tsv$phyloCCF[a]  = mut.info$mutCopyNum/best.CN
          pyClone.tsv$phyloCCF.0.05[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = 2)/best.CN 
          pyClone.tsv$phyloCCF.0.95[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[2],cellularity = cellularity,CNn = 2)/best.CN 
          pyClone.tsv$no.chrs.bearing.mut[a] = best.CN
          #pyClone.tsv$whichFrac[a] = c('A')
          pyClone.tsv$expected.VAF[a] = expected.prop
          
        }
      }
      if (!is.na(mut.info$fracC))
      {
        possible.subclonal.fractions <- c(1)
        # there are going to be quite a few options here, let's list them
        if((unlist(mut.info$nMaj_A)+unlist(mut.info$nMaj_B))>(unlist(mut.info$nMaj_C)+unlist(mut.info$nMaj_D)))
        {
          fracAB                       <- unlist(mut.info$fracA)+unlist(mut.info$fracB)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracAB)
        }
        if((unlist(mut.info$nMaj_A)+unlist(mut.info$nMaj_B))<(unlist(mut.info$nMaj_C)+unlist(mut.info$nMaj_D)))
        {
          fracCD                       <- unlist(mut.info$fracC)+unlist(mut.info$fracD)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracCD)
        }
        if((unlist(mut.info$nMin_A)+unlist(mut.info$nMin_D))>(unlist(mut.info$nMin_C)+unlist(mut.info$nMin_B)))
        {
          fracAD                       <- unlist(mut.info$fracA)+unlist(mut.info$fracD)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracAD)
        }
        if((unlist(mut.info$nMin_A)+unlist(mut.info$nMin_D))<(unlist(mut.info$nMin_C)+unlist(mut.info$nMin_B)))
        {
          fracBC                       <- unlist(mut.info$fracB)+unlist(mut.info$fracC)
          possible.subclonal.fractions <- c(possible.subclonal.fractions,fracBC)
        }
        
        best.CN                      <-  possible.subclonal.fractions[which.min(abs(mut.info$mutCopyNum/possible.subclonal.fractions - 1))]
        if(best.CN==1)
        {
          pyClone.tsv$whichFrac[a] <- c('A,B,C,D')
          next;
        }
        var.count       = as.numeric(mut.info$var_counts)
        depth.count     =  as.numeric(mut.info$var_counts)+as.numeric(mut.info$ref_counts)
        expected.prop   = pyClone.tsv$expected.VAF[a] * best.CN
        
        # check whether subclonal CN results in clonal mutation
        # otherwise subclonal CN doesn't explain subclonal MCN
        if(best.CN != 1 & prop.test(var.count,depth.count*purity,expected.prop/purity,alternative = 'less')$p.value > 0.01)
        {
          pyClone.tsv$phyloCCF[a]       = mut.info$mutCopyNum/best.CN
          pyClone.tsv$phyloCCF.0.05[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[1],cellularity = cellularity,CNn = 2)/best.CN 
          pyClone.tsv$phyloCCF.0.95[a]  = get.mut.mult(CNt = as.numeric(unlist(mut.info$major_raw))+as.numeric(unlist(mut.info$minor_raw)),Vaf = prop.test(var.count,depth.count,expected.prop)$conf.int[2],cellularity = cellularity,CNn = 2)/best.CN 
          pyClone.tsv$no.chrs.bearing.mut[a] = best.CN
          #pyClone.tsv$min.or.maj[a] = min.or.maj
          pyClone.tsv$expected.VAF[a] = expected.prop
          
        }
      }
      
      
      
    }
  }
  # Next, let's deal with potentially amplified mutations
  # convert MCN to subclonal fraction - tricky for amplified mutations
  # test for mutations in more than 1 copy
  
  
  
  p.vals = sapply(1:nrow(pyClone.tsv),function(var.count,depth,e,i){
    prop.test(var.count[i]
              ,depth[i]
              ,e[i]
              ,alternative="greater")$p.value 
  },var.count=as.numeric(pyClone.tsv$var_counts)
  ,depth=as.numeric(pyClone.tsv$var_counts)+as.numeric(pyClone.tsv$ref_counts)
  ,e=as.numeric(pyClone.tsv$expected.VAF))
  
  amplified.muts <- which(p.vals<=0.05&pyClone.tsv$mutCopyNum>1)
  
  if(length(amplified.muts)>0)
  {  	
    for(a in amplified.muts)
    {
      mut.info    = pyClone.tsv[a,,drop=FALSE]
      if(is.na(mut.info$fracC))
      {
        max.CN2     = 0
        max.CN1     = unlist(mut.info$nMaj_A)
        frac1.mut   = unlist(mut.info$fracA)
        frac2.mut   = 0
        
        
        #swap subclones, so that the one with the higher CN is first
        if(unlist(mut.info$nMaj_B)>max.CN1){
          max.CN2 = max.CN1
          max.CN1 = unlist(mut.info$nMaj_B)
          frac2.mut = frac1.mut
          frac1.mut = unlist(mut.info$fracB)   	
        }else{
          max.CN2   = unlist(mut.info$nMaj_B)
          frac2.mut = unlist(mut.info$fracB)
        }
        
        
        best.err = mut.info$mutCopyNum - 1
        best.CN=1
        for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN = j * frac1.mut + k * frac2.mut
            err = abs(mut.info$mutCopyNu/potential.CN-1)
            if(err<best.err){
              
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err=err
              best.CN = potential.CN
            }
          }
        }
        pyClone.tsv$phyloCCF[a] = pyClone.tsv$mutCopyNum[a] / best.CN
        pyClone.tsv$expected.VAF[a]   = pyClone.tsv$expected.VAF[a]*best.CN
        pyClone.tsv$phyloCCF.0.05[a]  = pyClone.tsv$phyloCCF.0.05[a]/best.CN
        pyClone.tsv$phyloCCF.0.95[a]  = pyClone.tsv$phyloCCF.0.95[a]/best.CN
        
      }
      if(!is.na(mut.info$fracC))
      {
        
        max.CN2         = 0
        max.CN1         = unlist(mut.info$nMaj_A)
        frac1.mut       = unlist(mut.info$fracA)+unlist(mut.info$fracB)
        frac2.mut       = 0 
        
        
        #swap subclones, so that the one with the higher CN is first
        if(mut.info$nMaj_C>max.CN1){
          max.CN2       = max.CN1
          max.CN1       = unlist(mut.info$nMaj_C)
          frac2.mut     = frac1.mut
          frac1.mut     = unlist(mut.info$fracC)+unlist(mut.info$fracD)
          
          
        }else{
          max.CN2       = unlist(mut.info$nMaj_C)
          frac2.mut     = unlist(mut.info$fracC)+unlist(mut.info$fracD)
          
        }
        
        
        best.err = mut.info$mutCopyNum - 1
        best.CN=1
        for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN  = j * frac1.mut + k * frac2.mut
            err           = abs(mut.info$mutCopyNu/potential.CN-1)
            
            
            if(err<best.err){
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err= err
              best.CN = potential.CN
            }
          }
        }
         
        # next, let's see whether we should also look at the minor allele
        max.CN2         = 0
        max.CN1         = unlist(mut.info$nMin_A)
        frac1.mut       = unlist(mut.info$fracA)+unlist(mut.info$fracD)
        frac2.mut       = 0 
        
        
        #swap subclones, so that the one with the higher CN is first
        if(unlist(mut.info$nMin_B)>max.CN1){
          max.CN2       = max.CN1
          max.CN1       = unlist(mut.info$nMin_B)
          frac2.mut     = frac1.mut
          frac1.mut     = unlist(mut.info$fracB)+unlist(mut.info$fracC)
          
          
        }else{
          max.CN2       = unlist(mut.info$nMin_B)
          frac2.mut     = unlist(mut.info$fracB)+unlist(mut.info$fracC)
          
        }
        
        
         for(j in 1:max.CN1){
          for(k in (j-1):min(j,max.CN2)){
            potential.CN  = j * frac1.mut + k * frac2.mut
            err           = abs(mut.info$mutCopyNu/potential.CN-1)
            
            
            if(err<best.err){
              pyClone.tsv$no.chrs.bearing.mut[a] = potential.CN
              best.err= err
              best.CN = potential.CN
            }
          }
        }
        pyClone.tsv$phyloCCF[a]       = pyClone.tsv$mutCopyNum[a] / best.CN
        pyClone.tsv$expected.VAF[a]   = pyClone.tsv$expected.VAF[a]*best.CN
        pyClone.tsv$phyloCCF.0.05[a]  = pyClone.tsv$phyloCCF.0.05[a]/best.CN
        pyClone.tsv$phyloCCF.0.95[a]  = pyClone.tsv$phyloCCF.0.95[a]/best.CN 
      }
    }
  }
  
  # finally, let's sort out 'missing' ones
  pyClone.tsv[pyClone.tsv$var_counts==0,'no.chrs.bearing.mut'] <- 0
  pyClone.tsv[pyClone.tsv$var_counts==0,'expected.VAF']        <- 0
  pyClone.tsv[pyClone.tsv$var_counts==0,'absolute.ccf']        <- 0
  return(pyClone.tsv)
  
}



create.subclonal.copy.number <- function(seg.mat.copy
                                         ,min.subclonal=0.1)
{
  seg.out           <- seg.mat.copy
  seg.out$nMaj1     <- floor(as.numeric(seg.out$nAraw))
  seg.out$nMaj2     <- ceiling(as.numeric(seg.out$nAraw))
  seg.out$nMin1     <- floor(as.numeric(seg.out$nBraw))
  seg.out$nMin2     <- ceiling(as.numeric(seg.out$nBraw))
  
  
  seg.out$fracMaj1  <- as.numeric(seg.out$nMaj2)-as.numeric(seg.out$nAraw)
  seg.out$fracMaj2  <- 1-as.numeric(seg.out$fracMaj1)
  seg.out$fracMin1  <- as.numeric(seg.out$nMin2)-as.numeric(seg.out$nBraw)
  seg.out$fracMin2  <- 1-as.numeric(seg.out$fracMin1)
  
  
  # next, let's deal with the minimum subclonal 
  seg.out$fracMaj1 <- ifelse(seg.out$fracMaj1<min.subclonal,0
                             ,ifelse(seg.out$fracMaj1>(1-min.subclonal),1,seg.out$fracMaj1))
  
  seg.out$fracMaj2 <- ifelse(seg.out$fracMaj2<min.subclonal,0
                             ,ifelse(seg.out$fracMaj2>(1-min.subclonal),1,seg.out$fracMaj2))
  
  seg.out$fracMin1 <- ifelse(seg.out$fracMin1<min.subclonal,0
                             ,ifelse(seg.out$fracMin1>(1-min.subclonal),1,seg.out$fracMin1))
  
  seg.out$fracMin2 <- ifelse(seg.out$fracMin2<min.subclonal,0
                             ,ifelse(seg.out$fracMin2>(1-min.subclonal),1,seg.out$fracMin2))
  
  # how much of the genome for each tumour region is subject to subclonal copy number
  for (region in unique(seg.mat.copy$SampleID))
  {
    seg.region <- seg.out[seg.out$SampleID==region,,drop=FALSE]
    prop.aber  <- 1-length(which((seg.region$fracMaj1==1|seg.region$fracMaj1==0)&(seg.region$fracMin1==1|seg.region$fracMin1==0)))/nrow(seg.region)
    print(prop.aber)
  }
  
  # which ones work?
  
  seg.sorted      <- seg.out[(seg.out$fracMaj1==seg.out$fracMin1|seg.out$fracMaj1==seg.out$fracMin2)|(seg.out$fracMaj1==1)|(seg.out$fracMaj1==0)|(seg.out$fracMin1==1)|(seg.out$fracMin1==0),,drop=FALSE]
  seg.problem     <- seg.out[which(!rownames(seg.out)%in%rownames(seg.sorted)),]
  
  
  
  # let's divide the problem further
  seg.out.final  <- c()
  
  for (a in 1:nrow(seg.problem))
  {
    seg.info         <- seg.problem[a,,drop=FALSE]
    seg.info$fracA   <- seg.info$fracMaj1*seg.info$fracMin1  
    seg.info$fracB   <- seg.info$fracMaj1-seg.info$fracA
    seg.info$fracC   <- seg.info$fracMaj2*seg.info$fracMin2
    seg.info$fracD   <- seg.info$fracMaj2-seg.info$fracC
    
    # let's see how this works. 
    seg.info$nMaj_A <- seg.info$nMaj1
    #seg.info$nMaj_A <- seg.info$nMaj1
    #seg.info$nMin_A <- seg.info$nMin1
    seg.info$nMin_A <- seg.info$nMin1
    
    seg.info$nMaj_B <- seg.info$nMaj1
    #seg.info$nMaj2_B <- seg.info$nMaj1
    #seg.info$nMin1_B <- seg.info$nMin2
    seg.info$nMin_B <- seg.info$nMin2
    
    seg.info$nMaj_C <- seg.info$nMaj2
    #seg.info$nMaj2_C <- seg.info$nMaj2
    seg.info$nMin_C <- seg.info$nMin2
    #seg.info$nMin2_C <- seg.info$nMin2
    
    seg.info$nMaj_D <- seg.info$nMaj2
    #seg.info$nMaj2_D <- seg.info$nMaj2
    seg.info$nMin_D <- seg.info$nMin1
    #seg.info$nMin2_D <- seg.info$nMin1
    
    
    
    seg.out.final <- rbind(seg.out.final,seg.info)
    
    
  }
  
  # let's make the sorted easier to read
  seg.sorted$fracA            <- NA
  seg.sorted$fracB            <- NA
  seg.sorted$fracC            <- NA
  seg.sorted$fracD            <- NA
  seg.sorted$nMaj_A           <- NA
  seg.sorted$nMaj_B           <- NA
  seg.sorted$nMaj_C           <- NA
  seg.sorted$nMaj_D           <- NA
  seg.sorted$nMin_A           <- NA
  seg.sorted$nMin_B           <- NA
  seg.sorted$nMin_C           <- NA
  seg.sorted$nMin_D           <- NA
  
  # put the columns we do need in the right order
  cols.maj <- apply(cbind(seg.sorted$fracMaj1,seg.sorted$fracMaj2),1,which.max)
  cols.min <- apply(cbind(seg.sorted$fracMin1,seg.sorted$fracMin2),1,which.max)
  for (i in 1:nrow(seg.sorted))
  {
    poss.frac.major         <- c(seg.sorted$fracMaj1[i],seg.sorted$fracMaj2[i])
    poss.frac.minor         <- c(seg.sorted$fracMin1[i],seg.sorted$fracMin2[i])
    poss.cpn.major          <- c(seg.sorted$nMaj1[i],seg.sorted$nMaj2[i])
    poss.cpn.minor          <- c(seg.sorted$nMin1[i],seg.sorted$nMin2[i])
    cols.major         <- cols.maj[i]
    cols.minor         <- cols.min[i]
    
    fracA.major         <- poss.frac.major[cols.major]
    fracA.minor         <- poss.frac.minor[cols.minor]
    fracB.major         <- poss.frac.major[which(!1:2%in%cols.major)]
    fracB.minor         <- poss.frac.minor[which(!1:2%in%cols.minor)]
    
    seg.sorted$nMaj_A[i]              <- poss.cpn.major[cols.major]
    seg.sorted$nMin_A[i]              <- poss.cpn.minor[cols.minor]
    seg.sorted$nMaj_B[i]              <- poss.cpn.major[which(!1:2%in%cols.major)]
    seg.sorted$nMin_B[i]              <- poss.cpn.minor[which(!1:2%in%cols.minor)]
    
    seg.sorted$fracA[i]               <- fracA.major
    seg.sorted$fracB[i]               <- fracB.major
    if(fracA.major==1)
    {
      seg.sorted$fracA[i]  <- fracA.minor
      seg.sorted$fracB[i]  <- fracB.minor
      seg.sorted$nMaj_B[i] <- poss.cpn.major[cols.major]
    }
    
    if(fracA.minor==1)
    {
      seg.sorted$fracA[i]  <- fracA.major
      seg.sorted$fracB[i]  <- fracB.major
      seg.sorted$nMin_B[i] <- poss.cpn.minor[cols.minor]
    }
    
    
    
  }
  seg.final <- rbind(seg.out.final,seg.sorted) 
  
  #let's order this correctly
  seg.final <- seg.final[order(seg.final$SampleID,seg.final$chr,seg.final$startpos),,drop=FALSE]
  # finally, let's choose the columns we want and the order we want
  colnames(seg.final)
  col.names <- c("SampleID","chr","startpos","endpos","n.het","cnTotal","nMajor","nMinor","Ploidy"
                 ,"ACF","nAraw","nBraw","fracA","nMaj_A","nMin_A","fracB","nMaj_B","nMin_B"
                 ,"fracC", "nMaj_C",  "nMin_C","fracD","nMaj_D","nMin_D")
  seg.final  <- seg.final[,col.names]
  return(seg.final)
  
}





PasteVector <- function(v,sep=""){
  
  vt <- v[1];
  if(length(v) > 1){
    for(g in 2:length(v)){
      vt <- paste(vt,v[g],sep=sep)
      
    }
  }
  vt <- paste(vt," EnD",sep="");
  out.v <- sub(" EnD","",vt);
  out.v <- sub("NA , ","",out.v);
  out.v <- sub(" , NA","",out.v);
  out.v <- sub(" , NA , "," , ",out.v);
  return(out.v);
}

#### Plotting ######

plot.EarlyOrLate.raw    <- function( seg.mat.patient
                                     , TCGA.earlyLate
                                     , TCGA.purity
                                     , TCGA.barcode=""
                                     #, prob.early=0.75
                                     #, prob.late=0.75
                                     , max.cpn = 5
                                     , min.probes = 10
                                     , sub.clonal = 1
                                     , min.vaf.present = 0.01
)
{
  
  
  seg.mat.patient$Chromosome    <- as.numeric(as.character(seg.mat.patient$Chromosome))
  seg.mat.patient$StartPosition <- as.numeric(as.character(seg.mat.patient$StartPosition))
  seg.mat.patient$EndPosition   <- as.numeric(as.character(seg.mat.patient$EndPosition))
  seg.mat.patient$nr.probes     <- as.numeric(as.character(seg.mat.patient$nr.probes))
  seg.mat.patient$major         <- seg.mat.patient$nAraw
  seg.mat.patient$min.allele    <- seg.mat.patient$nBraw
  
  
  chrom.length.copy   <- fun.chrom.length(seg.mat.patient[,2:4])
  chrom.segs          <- fun.add.chrom(seg.mat.patient[,2:4],chrom.length.copy)
  seg.mat.plot        <- seg.mat.patient
  seg.mat.plot[,2:4]  <- chrom.segs
  major               <- seg.mat.plot$nAraw
  seg.mat.plot        <- cbind(seg.mat.plot,major)
  #seg.mat.plot        <- seg.mat.plot[seg.mat.plot$nr.probes>=min.probes,]
  
  
  # Order correctly
  TCGA.plot          <- TCGA.earlyLate
  Chromosome         <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,2])
  Start_pos          <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,3])
  TCGA.plot          <- cbind(TCGA.plot,Chromosome,Start_pos)
  TCGA.plot          <- data.frame(apply(TCGA.plot,2,unlist),stringsAsFactors=FALSE)
  
  min.x               <- 0
  min.y               <- -0.25
  max.y               <- max.cpn +1
  max.x               <- as.numeric(max(seg.mat.plot$EndPosition))
  
  # Make sure any copy numbers that are too high are still included
  seg.mat.plot$major       <- ifelse(seg.mat.plot$major>max.cpn,max.cpn,seg.mat.plot$major)
  seg.mat.plot$min.allele  <- ifelse(seg.mat.plot$min.allele>max.cpn,max.cpn,seg.mat.plot$min.allele)
  TCGA.plot$mutCopyNum      <- ifelse(TCGA.plot$mutCopyNum>max.cpn,max.cpn,TCGA.plot$mutCopyNum)
  
  #layout(rbind(1,2))
  #par(mar=c(5,5,5,5))
  plot(1
       ,xlim=c(min.x,max.x)
       ,ylim=c(min.y,max.y)
       ,type='n'
       ,xaxt='n'
       ,yaxs='i'
       ,yaxt='n'
       ,xaxs='i'
       ,ylab=""
       ,xlab=""
       ,lwd=2
       # ,main=region
       ,bty="n")
  
  box(lwd=0.5) 
  
  axis(side=2
       ,at=seq(0,max.cpn,by=1)
       ,labels = c(seq(0,max.cpn-1,by=1),paste('>',max.cpn,sep=""))
       ,las=2
       ,cex.axis=0.7
       ,lwd=0.5)
  
  mtext(text=TCGA.barcode
        ,side=3
        ,line=2
  )
  
  
  
  # Now plot each segment
  fun.plot.segment <- function(x,seg.mat.plot)
  {
    #print(x)
    #start with major allele
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- seg.mat.plot[x,'major']+0.05
    y1 <- y0
    #     if(y1>=1)
    #     {
    #       for (i in 1:(y1-1))
    #         segments(x0,i+0.05,x1,i+0.05
    #                  ,col="#00000090"
    #                  ,lwd=1.4
    #                  ,lty='dotted')
    #     }
    #     
    
    segments(x0,y0,x1,y1,col='black',lwd=3)
    
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- as.numeric(seg.mat.plot[x,'min.allele'])-0.05
    y1 <- y0
    
    #     if(y1>=1)
    #     {
    #       for (i in 1:(y1-1))
    #         segments(x0,i-0.05,x1,i-0.05
    #                  ,col="#009E7395"
    #                  ,lwd=1.4
    #                  ,lty='dotted')
    #     }  
    
    #  
    segments(x0,y0,x1,y1,col='#009E73',lwd=3)
    
  }
  
  sapply(1:nrow(seg.mat.plot),fun.plot.segment,seg.mat.plot)
  # Draw lines to separate the different chromosomes
  for (i in sort(unique(as.numeric(seg.mat.plot$Chromosome))))
  {
    abline(v=max(as.numeric(seg.mat.plot[seg.mat.plot$Chromosome==i,'EndPosition']))
           ,lwd=0.5)
    #,lty='dashed')
  }
  
  # We need to update the 
  TCGA.plot <- TCGA.plot[order(as.numeric(TCGA.plot$Chromosome),as.numeric(TCGA.plot$Start_pos)),]
  
  TCGA.plot$Start_pos <- as.numeric(fun.add.chrom(cbind(TCGA.plot$Chromosome,TCGA.plot$Start_pos,TCGA.plot$Start_pos),chrom.length.copy)[,2])
  
  # # Let's add the mutations
  absent.muts    <- TCGA.plot[as.numeric(TCGA.plot$var_counts)/(as.numeric(TCGA.plot$ref_counts)+as.numeric(TCGA.plot$var_counts))<min.vaf.present,,drop=FALSE]
  
  
  # Start with early mutations
  early.muts      <- TCGA.plot[TCGA.plot$phyloCCF.0.95>=sub.clonal&!TCGA.plot$mutation_id%in%absent.muts$mutation_id,]
  #early.muts      <- early.muts[as.numeric(early.muts$Exp.Cpn.Likelihood)>prob.early,]
  early.clonal    <- early.muts[early.muts$absolute.ccf.0.95>=sub.clonal,]
  early.subclonal <- early.muts[early.muts$absolute.ccf.0.95<sub.clonal,]
  
  late.muts  <- TCGA.plot[TCGA.plot$phyloCCF.0.95<sub.clonal&!TCGA.plot$mutation_id%in%absent.muts$mutation_id,]
  #late.muts  <- late.muts[as.numeric(late.muts$Exp.Cpn.Likelihood)>prob.late,]
  late.clonal <- late.muts[late.muts$absolute.ccf.0.95>=sub.clonal,]
  late.subclonal <- late.muts[late.muts$absolute.ccf.0.95<sub.clonal,]
  
  nontimed.muts       <- TCGA.plot[!TCGA.plot$mutation_id%in%c(early.muts$mutation_id,late.muts$mutation_id,absent.muts$mutation_id),]
  nontimed.clonal     <- nontimed.muts[nontimed.muts$ccf.btstr.0.95>=sub.clonal,]
  nontimed.subclonal  <- nontimed.muts[nontimed.muts$ccf.btstr.0.95<sub.clonal,]
  
  
  # let's determine early or late or not possible
  plot.earlyORlateCol <- function(x,timed.muts,nontimed=FALSE)
  {
    #print(x)
    mut <- timed.muts[x,,drop=FALSE]
    
    #determine cell multiplicity
    mut.multiplicity <- mut$mutCopyNum
    #mut.multiplicity <- mut.multiplicity*(as.numeric(mut$minor_cn)+as.numeric(mut$major_cn))
    
    if(nontimed)
    {
      if(as.numeric(mut$absolute.ccf.0.95)<sub.clonal)
      {
        points(mut$Start_pos
               ,mut.multiplicity
               ,cex=0.7
               ,pch=17
               ,col='#99999965') # cannot be determined
        
      } 
      
      if(as.numeric(mut$absolute.ccf.0.95)>=sub.clonal)
      {
        points(mut$Start_pos
               ,mut.multiplicity
               ,cex=0.8
               ,pch=16
               ,col='#99999995') # cannot be determined
        
      } 
      
    }
    
    if(!nontimed)
    {
      # is it an early event?
      if(mut$phyloCCF>=0.75)
      {
        # is it subclonal
        if(as.numeric(mut$absolute.ccf.0.95)<sub.clonal)
        {
          points(mut$Start_pos
                 ,mut.multiplicity
                 ,cex=0.7
                 ,pch=17
                 ,col="#0072B245") # early event
        }
        
        # is it clonal
        if(as.numeric(mut$absolute.ccf.0.95)>=sub.clonal)
        {
          points(mut$Start_pos
                 ,mut.multiplicity
                 ,cex=0.9
                 ,pch=16
                 ,col="#0072B295") # early event
        }
        
        
      }
      
      if(mut$phyloCCF<0.75)
      {
        if(as.numeric(mut$absolute.ccf.0.95)<sub.clonal)
        {
          points(mut$Start_pos
                 ,mut.multiplicity
                 ,cex=0.7
                 ,pch=17
                 ,col="#D55E0045") # late event
        }
        
        if(as.numeric(mut$absolute.ccf.0.95)>=sub.clonal)
        {
          points(mut$Start_pos
                 ,mut.multiplicity
                 ,cex=0.9
                 ,pch=16
                 ,col="#D55E0095") # late event
        }
        
        
        
      }
      
    }
    
  }
  
  
  if(nrow(early.muts)>=1)
  {
    sapply(1:nrow(early.muts),plot.earlyORlateCol,early.muts)  
  }
  if(nrow(late.muts)>=1)
  {
    sapply(1:nrow(late.muts),plot.earlyORlateCol,late.muts)  
  }
  if(nrow(nontimed.muts)>=1)
  {
    sapply(1:nrow(nontimed.muts),plot.earlyORlateCol,nontimed.muts,nontimed=TRUE)
  }
  
  if(nrow(absent.muts)>=1)
  {
    points(absent.muts$Start_pos
           ,absent.muts$mutCopyNum
           ,cex=0.7
           ,pch=25
           ,col='#99999980'
           ,bg='#99999980') # cannot be determined
    
  }
  
  
  
  mtext(side=3
        ,at=fun.chrom.mean(seg.mat.plot[,2:4])
        ,text=sort(unique(seg.mat.plot[,2]))
        ,cex=seq(0.6,0.4,length.out=length(unique(seg.mat.plot[,2])))
        ,line=-1
        #,lwd=0.5
  )
  sapply(1:nrow(seg.mat.plot),fun.plot.segment,seg.mat.plot)
  
}


##
plot.region.mutCopyNum  <- function(phylo.region.list
                                    , seg.mat.copy
                                    , mostLikelyClusters = most.likely.cluster
                                    , plot.separate.clusters=FALSE
)
{
  if(!plot.separate.clusters)
  {
    regions.to.use <- names(phylo.region.list)
    lyout <- c()
    for (i in 1:length(regions.to.use))
    {
      lyout <- rbind(lyout,rbind(c(rep(i,9),length(regions.to.use)+1)))
    }
    layout(lyout)
    for (region in regions.to.use)
    {
      
      
      # add extra column to 
      sub.mat.copy               <- seg.mat.copy[seg.mat.copy$SampleID%in%region,,drop=FALSE]
      sub.mat.copy$Copy.Number   <- sub.mat.copy$cnTotal
      sub.mat.copy$min.allele    <- apply(cbind(sub.mat.copy$nMinor,sub.mat.copy$nMajor),1,min)
      colnames(sub.mat.copy)[2]  <- 'Chromosome'
      colnames(sub.mat.copy)[3]  <- 'StartPosition'
      colnames(sub.mat.copy)[4]  <- 'EndPosition'
      colnames(sub.mat.copy)[5]  <- 'nr.probes'
      region.phyloCCF            <- phylo.region.list[[region]]
      
      #pdf(early.late.pdf)
      par(mar=c(5,4,5,0.2))
      if(nrow(lyout)==1)
      {
        par(mar=c(5,4,5,0.2))
      }
      
      par(lend=1)
      #layout(rbind(c(rep(1,9),2)))
      plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                              ,TCGA.earlyLate=region.phyloCCF
                              #,TCGA.purity=sample.purity
                              #,TCGA.barcode=region
                              #,prob.early=prob.early
                              #,prob.late=prob.late
                              ,sub.clonal=1
                              ,most.likely.cluster=mostLikelyClusters
      )
      mtext(unlist(strsplit(region,split="\\."))[2], side=2,cex=1,line=2,las=2)
      #its hard to distinguish more than 8 different colours
    }
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    require(RColorBrewer)
    cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    par(mar=c(5,0.5,5,0))
    if(nrow(lyout)==1)
    {
      par(mar=c(5,0.5,5,0))
    }
    plot(1, type="n"
         , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
    for (k in 1:no.optima)
    {
      points(0
             ,k
             ,col=cols.opac[k]
             ,bg = cols[k]
             ,pch=pch.vals[k]
             ,cex=1.5)
      text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
           ,pos=4)
    }
    
    text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
  }
  if(plot.separate.clusters)
  {
    regions.to.use <- names(phylo.region.list)
    lyout <- c()
    for (i in 1:length(regions.to.use))
    {
      lyout <- rbind(lyout,rbind(c(rep(i,9),length(regions.to.use)+1)))
    }
    layout(lyout)
    for (region in regions.to.use)
    {
      
      
      # add extra column to 
      sub.mat.copy               <- seg.mat.copy[seg.mat.copy$SampleID%in%region,,drop=FALSE]
      sub.mat.copy$Copy.Number   <- sub.mat.copy$cnTotal
      sub.mat.copy$min.allele    <- apply(cbind(sub.mat.copy$nMinor,sub.mat.copy$nMajor),1,min)
      colnames(sub.mat.copy)[2]  <- 'Chromosome'
      colnames(sub.mat.copy)[3]  <- 'StartPosition'
      colnames(sub.mat.copy)[4]  <- 'EndPosition'
      colnames(sub.mat.copy)[5]  <- 'nr.probes'
      region.phyloCCF            <- phylo.region.list[[region]]
      
      #pdf(early.late.pdf)
      par(mar=c(0.5,4,0.5,0.2))
      #par(mar=c(0.5,4,0.5,0.2))
      if(nrow(lyout)==1)
      {
        par(mar=c(5,4,5,0.2))
      }
      par(lend=1)
      #layout(rbind(c(rep(1,9),2)))
      plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                              ,TCGA.earlyLate=region.phyloCCF
                              #,TCGA.purity=sample.purity
                              #,TCGA.barcode=region
                              #,prob.early=prob.early
                              #,prob.late=prob.late
                              ,sub.clonal=1
                              ,most.likely.cluster=mostLikelyClusters
      )
      mtext(unlist(strsplit(region,split="\\."))[2], side=2,cex=1,line=2,las=2)
      #its hard to distinguish more than 8 different colours
    }
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    require(RColorBrewer)
    cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    par(mar=c(2,0.5,2,0))
    if(nrow(lyout)==1)
    {
      par(mar=c(5,0.5,5,0))
    }
    
    plot(1, type="n"
         , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
    for (k in 1:no.optima)
    {
      points(0
             ,k
             ,col=cols.opac[k]
             ,bg = cols[k]
             ,pch=pch.vals[k]
             ,cex=1.5)
      text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
           ,pos=4)
    }
    
    text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
    
    lyout <- c()
    for (i in seq(1,length(regions.to.use)*2,by=2))
    {
      lyout <- rbind(lyout,rbind(c(rep(i,9),i+1)))
    }
    
    layout(lyout)
    
    # next, let's plot each cluster separately
    for (cluster in 1:max(mostLikelyClusters))
    {
      for (region in regions.to.use)
      {
        
        
        # add extra column to 
        sub.mat.copy               <- seg.mat.copy[seg.mat.copy$SampleID%in%region,,drop=FALSE]
        sub.mat.copy$Copy.Number   <- sub.mat.copy$cnTotal
        sub.mat.copy$min.allele    <- apply(cbind(sub.mat.copy$nMinor,sub.mat.copy$nMajor),1,min)
        colnames(sub.mat.copy)[2]  <- 'Chromosome'
        colnames(sub.mat.copy)[3]  <- 'StartPosition'
        colnames(sub.mat.copy)[4]  <- 'EndPosition'
        colnames(sub.mat.copy)[5]  <- 'nr.probes'
        region.phyloCCF            <- phylo.region.list[[region]]
        region.phyloCCF.cluster    <- region.phyloCCF[names(which(mostLikelyClusters==cluster)),,drop=FALSE]
        #pdf(early.late.pdf)
        par(mar=c(0.5,4,0.5,0.2))
        if(nrow(lyout)==1)
        {
          par(mar=c(5,4,5,0.2))
        }
        par(lend=1)
        #layout(rbind(c(rep(1,9),2)))
        plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                                ,TCGA.earlyLate=region.phyloCCF
                                #,TCGA.purity=sample.purity
                                #,TCGA.barcode=region
                                #,prob.early=prob.early
                                #,prob.late=prob.late
                                ,sub.clonal=1
                                ,most.likely.cluster=mostLikelyClusters
                                ,cluster=cluster
        )
        mtext(unlist(strsplit(region,split="\\."))[2], side=2,cex=1,line=2,las=2)
        
        if(nrow(region.phyloCCF.cluster)==1)
        {
          plot(NULL, type = "n"
               , xlim = c(0, max(A$density))
               , ylim=c(-0.25,6)
               ,bty='n'
               ,xaxs='i'
               ,xaxt='n'
               ,yaxt='n'
               ,yaxs='i'
               ,xlab=""
               ,main=""
               ,ylab=""
          )
          next;
        }
        ds <- density(ifelse(as.numeric(region.phyloCCF.cluster$mutCopyNum)>5,5,as.numeric(region.phyloCCF.cluster$mutCopyNum)))
        ds1 <- ds
        ds1$x <- ds$y
        ds1$y <- ds$x
        par(mar=c(0.5,0,0.5,4))
        if(nrow(lyout)==1)
        {
          par(mar=c(5,0,5,4))
        }
        A <- hist(ifelse(as.numeric(region.phyloCCF.cluster$mutCopyNum)>5,5,as.numeric(region.phyloCCF.cluster$mutCopyNum)),breaks=seq(-0.25,6,by=0.1),plot=FALSE)
        plot(NULL, type = "n"
             , xlim = c(0, max(A$density))
             , ylim=c(-0.25,6)
             ,bty='n'
             ,xaxs='i'
             ,xaxt='n'
             ,yaxt='n'
             ,yaxs='i'
             ,xlab=""
             ,main=""
             ,ylab=""
        )
        rect(0, A$breaks[1:(length(A$breaks) - 1)], A$density, A$breaks[2:length(A$breaks)]
             ,border=TRUE,col="#CC6666")
        
        lines(ds1)
        
        
      }
      
      
    }
    
    
  }
  
}


plot.region  <- function(phylo.region.list
                         , seg.mat.copy
                         , mostLikelyClusters = most.likely.cluster
)
{
  regions.to.use <- names(phylo.region.list)
  lyout <- c()
  for (i in 1:length(regions.to.use))
  {
    lyout <- rbind(lyout,rbind(c(rep(i,9),length(regions.to.use)+1)))
  }
  layout(lyout)
  for (region in regions.to.use)
  {
    
    
    # add extra column to 
    sub.mat.copy               <- seg.mat.copy[seg.mat.copy$SampleID%in%region,,drop=FALSE]
    sub.mat.copy$Copy.Number   <- sub.mat.copy$cnTotal
    sub.mat.copy$min.allele    <- apply(cbind(sub.mat.copy$nMinor,sub.mat.copy$nMajor),1,min)
    colnames(sub.mat.copy)[2]  <- 'Chromosome'
    colnames(sub.mat.copy)[3]  <- 'StartPosition'
    colnames(sub.mat.copy)[4]  <- 'EndPosition'
    colnames(sub.mat.copy)[5]  <- 'nr.probes'
    region.phyloCCF            <- phylo.region.list[[region]]
    
    #pdf(early.late.pdf)
    par(mar=c(0.5,4,0.5,0.2))
    par(lend=1)
    #layout(rbind(c(rep(1,9),2)))
    plot.simpleClusters.raw(seg.mat.patient=sub.mat.copy
                            ,TCGA.earlyLate=region.phyloCCF
                            #,TCGA.purity=sample.purity
                            #,TCGA.barcode=region
                            #,prob.early=prob.early
                            #,prob.late=prob.late
                            ,sub.clonal=1
                            ,most.likely.cluster=mostLikelyClusters
    )
    mtext(unlist(strsplit(region,split="\\."))[2], side=2,cex=1,line=2,las=2)
    #its hard to distinguish more than 8 different colours
  }
  no.optima = length(unique(mostLikelyClusters))
  max.cols = 12
  require(RColorBrewer)
  cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
  cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
  cols.opac      = paste(cols,'99',sep="")
  
  
  
  pch.vals       = c()
  k <- 21
  m <- 1  
  while(m<=no.optima)
  {
    pch.vals    = c(pch.vals,k)
    if (m==max.cols)
    {
      k <- k+1
    }
    m <- m+1
  }
  
  pch.vals     = pch.vals
  par(mar=c(2,0.5,2,0))
  plot(1, type="n"
       , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
  for (k in 1:no.optima)
  {
    points(0
           ,k
           ,col=cols.opac[k]
           ,bg = cols[k]
           ,pch=pch.vals[k]
           ,cex=1.5)
    text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
         ,pos=4)
  }
  
  text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
  
  
  
}



plot.simpleClusters    <- function( seg.mat.patient
                                    , TCGA.earlyLate
                                    #, TCGA.purity
                                    , TCGA.barcode=""
                                    #, prob.early=0.75
                                    #, prob.late=0.75
                                    , max.cpn = 5
                                    , min.probes = 10
                                    , sub.clonal = 1
                                    , min.vaf.present = 0.01
                                    , most.likely.cluster
                                    ,cluster=NA
                                    
)
{
  
  
  seg.mat.patient$Chromosome    <- as.numeric(as.character(seg.mat.patient$Chromosome))
  seg.mat.patient$StartPosition <- as.numeric(as.character(seg.mat.patient$StartPosition))
  seg.mat.patient$EndPosition   <- as.numeric(as.character(seg.mat.patient$EndPosition))
  seg.mat.patient$nr.probes     <- as.numeric(as.character(seg.mat.patient$nr.probes))
  
  chrom.length.copy   <- fun.chrom.length(seg.mat.patient[,2:4])
  chrom.segs          <- fun.add.chrom(seg.mat.patient[,2:4],chrom.length.copy)
  seg.mat.plot        <- seg.mat.patient
  seg.mat.plot[,2:4]  <- chrom.segs
  major               <- seg.mat.plot$Copy.Number - seg.mat.plot$min.allele
  seg.mat.plot        <- cbind(seg.mat.plot,major)
  #seg.mat.plot        <- seg.mat.plot[seg.mat.plot$nr.probes>=min.probes,]
  
  
  # Order correctly
  TCGA.plot          <- TCGA.earlyLate
  Chromosome         <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,2])
  Start_pos          <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,3])
  TCGA.plot          <- cbind(TCGA.plot,Chromosome,Start_pos)
  TCGA.plot          <- data.frame(apply(TCGA.plot,2,unlist),stringsAsFactors=FALSE)
  
  min.x               <- 0
  min.y               <- -0.25
  max.y               <- max.cpn +1
  max.x               <- as.numeric(max(seg.mat.plot$EndPosition))
  
  # Make sure any copy numbers that are too high are still included
  seg.mat.plot$major       <- ifelse(seg.mat.plot$major>max.cpn,max.cpn,seg.mat.plot$major)
  seg.mat.plot$min.allele  <- ifelse(seg.mat.plot$min.allele>max.cpn,max.cpn,seg.mat.plot$min.allele)
  TCGA.plot$mut.multi      <- ifelse(TCGA.plot$mut.multi>max.cpn,max.cpn,TCGA.plot$mut.multi)
  
  #layout(rbind(1,2))
  #par(mar=c(5,5,5,5))
  plot(1
       ,xlim=c(min.x,max.x)
       ,ylim=c(min.y,max.y)
       ,type='n'
       ,xaxt='n'
       ,yaxs='i'
       ,yaxt='n'
       ,xaxs='i'
       ,ylab=""
       ,xlab=""
       ,lwd=2
       # ,main=region
       ,bty="n")
  
  box(lwd=0.5) 
  
  axis(side=2
       ,at=seq(0,max.cpn,by=1)
       ,labels = c(seq(0,max.cpn-1,by=1),paste('>',max.cpn,sep=""))
       ,las=2
       ,cex.axis=0.7
       ,lwd=0.5)
  
  mtext(text=TCGA.barcode
        ,side=3
        ,line=2
  )
  
  
  
  # Now plot each segment
  fun.plot.segment <- function(x,seg.mat.plot)
  {
    #print(x)
    #start with major allele
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- seg.mat.plot[x,'major']+0.05
    y1 <- y0
    #     if(y1>=1)
    #     {
    #       for (i in 1:(y1-1))
    #         segments(x0,i+0.05,x1,i+0.05
    #                  ,col="#00000090"
    #                  ,lwd=1.4
    #                  ,lty='dotted')
    #     }
    #     
    
    segments(x0,y0,x1,y1,col='black',lwd=2.4)
    
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- as.numeric(seg.mat.plot[x,'min.allele'])-0.05
    y1 <- y0
    
    #     if(y1>=1)
    #     {
    #       for (i in 1:(y1-1))
    #         segments(x0,i-0.05,x1,i-0.05
    #                  ,col="#009E7395"
    #                  ,lwd=1.4
    #                  ,lty='dotted')
    #     }  
    
    #  
    segments(x0,y0,x1,y1,col='#009E73',lwd=2.4)
    
  }
  
  sapply(1:nrow(seg.mat.plot),fun.plot.segment,seg.mat.plot)
  # Draw lines to separate the different chromosomes
  for (i in sort(unique(as.numeric(seg.mat.plot$Chromosome))))
  {
    abline(v=max(as.numeric(seg.mat.plot[seg.mat.plot$Chromosome==i,'EndPosition']))
           ,lwd=0.5)
    #,lty='dashed')
  }
  
  # We need to update the 
  TCGA.plot          <- TCGA.plot[order(as.numeric(TCGA.plot$Chromosome),as.numeric(TCGA.plot$Start_pos)),]
  TCGA.plot$Start_pos <- as.numeric(fun.add.chrom(cbind(TCGA.plot$Chromosome,TCGA.plot$Start_pos,TCGA.plot$Start_pos),chrom.length.copy)[,2])
  
  
  max.cols = 12
  require(RColorBrewer)
  no.optima      = length(unique(most.likely.cluster))
  cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
  cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
  cols.opac      = paste(cols,'99',sep="")
  
  pch.vals       = c()
  k <- 21
  m <- 1  
  while(m<=no.optima)
  {
    pch.vals    = c(pch.vals,k)
    if (m==max.cols)
    {
      k <- k+1
    }
    m <- m+1
  }
  
  pch.vals     = rev(pch.vals)
  
  
  if (is.na(cluster[1]))
  {
    for (i in 1:length(unique(most.likely.cluster)))
    {
      
      cluster.muts <- TCGA.plot[TCGA.plot$mutation_id%in%names(which(most.likely.cluster==i)),,drop=FALSE]
      points(cluster.muts$Start_pos
             ,cluster.muts$mut.multi
             ,cex=0.7
             ,pch=pch.vals[i]
             ,col=cols[i]
             ,bg =cols.opac[i] ) # ca
      
      
      
    }
    
  }
  if (!is.na(cluster[1]))
  {
    for (i in cluster)
    {
      
      cluster.muts <- TCGA.plot[TCGA.plot$mutation_id%in%names(which(most.likely.cluster==i)),,drop=FALSE]
      points(cluster.muts$Start_pos
             ,cluster.muts$mut.multi
             ,cex=0.7
             ,pch=pch.vals[as.numeric(i)]
             ,col=cols[i]
             ,bg =cols.opac[as.numeric(i)] ) # ca
      
      
      
    }
    
  }
  
  
  
  
  
  
  
  mtext(side=3
        ,at=fun.chrom.mean(seg.mat.plot[,2:4])
        ,text=sort(unique(seg.mat.plot[,2]))
        ,cex=seq(0.6,0.4,length.out=length(unique(seg.mat.plot[,2])))
        ,line=-1
        #,lwd=0.5
  )
  
  
}


plot.simpleClusters.raw    <- function( seg.mat.patient
                                        , TCGA.earlyLate
                                        #, TCGA.purity
                                        , TCGA.barcode=""
                                        #, prob.early=0.75
                                        #, prob.late=0.75
                                        , max.cpn = 5
                                        , min.probes = 10
                                        , sub.clonal = 1
                                        , min.vaf.present = 0.01
                                        , most.likely.cluster
                                        ,cluster=NA
                                        
)
{
  
  
  seg.mat.patient$Chromosome    <- as.numeric(as.character(seg.mat.patient$Chromosome))
  seg.mat.patient$StartPosition <- as.numeric(as.character(seg.mat.patient$StartPosition))
  seg.mat.patient$EndPosition   <- as.numeric(as.character(seg.mat.patient$EndPosition))
  seg.mat.patient$nr.probes     <- as.numeric(as.character(seg.mat.patient$nr.probes))
  seg.mat.patient$major         <- seg.mat.patient$nAraw
  seg.mat.patient$min.allele    <- seg.mat.patient$nBraw
  
  
  chrom.length.copy   <- fun.chrom.length(seg.mat.patient[,2:4])
  chrom.segs          <- fun.add.chrom(seg.mat.patient[,2:4],chrom.length.copy)
  seg.mat.plot        <- seg.mat.patient
  seg.mat.plot[,2:4]  <- chrom.segs
  major               <- seg.mat.plot$nAraw
  seg.mat.plot        <- cbind(seg.mat.plot,major)
  #seg.mat.plot        <- seg.mat.plot[seg.mat.plot$nr.probes>=min.probes,]
  
  
  # Order correctly
  TCGA.plot          <- TCGA.earlyLate
  Chromosome         <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,2])
  Start_pos          <- as.numeric(do.call(rbind,strsplit(unlist(TCGA.earlyLate$mutation_id),split=":"))[,3])
  TCGA.plot          <- cbind(TCGA.plot,Chromosome,Start_pos)
  TCGA.plot          <- data.frame(apply(TCGA.plot,2,unlist),stringsAsFactors=FALSE)
  
  min.x               <- 0
  min.y               <- -0.25
  max.y               <- max.cpn +1
  max.x               <- as.numeric(max(seg.mat.plot$EndPosition))
  
  # Make sure any copy numbers that are too high are still included
  seg.mat.plot$major       <- ifelse(seg.mat.plot$major>max.cpn,max.cpn,seg.mat.plot$major)
  seg.mat.plot$min.allele  <- ifelse(seg.mat.plot$min.allele>max.cpn,max.cpn,seg.mat.plot$min.allele)
  TCGA.plot$mutCopyNum      <- ifelse(TCGA.plot$mutCopyNum>max.cpn,max.cpn,TCGA.plot$mutCopyNum)
  
  #layout(rbind(1,2))
  #par(mar=c(5,5,5,5))
  plot(1
       ,xlim=c(min.x,max.x)
       ,ylim=c(min.y,max.y)
       ,type='n'
       ,xaxt='n'
       ,yaxs='i'
       ,yaxt='n'
       ,xaxs='i'
       ,ylab=""
       ,xlab=""
       ,lwd=2
       # ,main=region
       ,bty="n")
  
  box(lwd=0.5) 
  
  axis(side=2
       ,at=seq(0,max.cpn,by=1)
       ,labels = c(seq(0,max.cpn-1,by=1),paste('>',max.cpn,sep=""))
       ,las=2
       ,cex.axis=0.7
       ,lwd=0.5)
  
  mtext(text=TCGA.barcode
        ,side=3
        ,line=2
  )
  
  
  # Now plot each segment
  fun.plot.segment <- function(x,seg.mat.plot)
  {
    #print(x)
    #start with major allele
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- seg.mat.plot[x,'major']+0.05
    y1 <- y0
    #     if(y1>=1)
    #     {
    #       for (i in 1:(y1-1))
    #         segments(x0,i+0.05,x1,i+0.05
    #                  ,col="#00000090"
    #                  ,lwd=1.4
    #                  ,lty='dotted')
    #     }
    #     
    
    segments(x0,y0,x1,y1,col='black',lwd=2.4)
    
    x0 <- as.numeric(seg.mat.plot[x,'StartPosition'])
    x1 <- as.numeric(seg.mat.plot[x,'EndPosition'])
    y0 <- as.numeric(seg.mat.plot[x,'min.allele'])-0.05
    y1 <- y0
    
    #     if(y1>=1)
    #     {
    #       for (i in 1:(y1-1))
    #         segments(x0,i-0.05,x1,i-0.05
    #                  ,col="#009E7395"
    #                  ,lwd=1.4
    #                  ,lty='dotted')
    #     }  
    
    #  
    segments(x0,y0,x1,y1,col='#009E73',lwd=2.4)
    
  }
  
  sapply(1:nrow(seg.mat.plot),fun.plot.segment,seg.mat.plot)
  # Draw lines to separate the different chromosomes
  for (i in sort(unique(as.numeric(seg.mat.plot$Chromosome))))
  {
    abline(v=max(as.numeric(seg.mat.plot[seg.mat.plot$Chromosome==i,'EndPosition']))
           ,lwd=0.5)
    #,lty='dashed')
  }
  
  # We need to update the 
  TCGA.plot          <- TCGA.plot[order(as.numeric(TCGA.plot$Chromosome),as.numeric(TCGA.plot$Start_pos)),]
  TCGA.plot$Start_pos <- as.numeric(fun.add.chrom(cbind(TCGA.plot$Chromosome,TCGA.plot$Start_pos,TCGA.plot$Start_pos),chrom.length.copy)[,2])
  
  
  max.cols = 12
  require(RColorBrewer)
  no.optima      = length(unique(most.likely.cluster))
  cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
  cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
  cols.opac      = paste(cols,'99',sep="")
  
  pch.vals       = c()
  k <- 21
  m <- 1  
  while(m<=no.optima)
  {
    pch.vals    = c(pch.vals,k)
    if (m==max.cols)
    {
      k <- k+1
    }
    m <- m+1
  }
  
  pch.vals     = rev(pch.vals)
  
  
  if (is.na(cluster[1]))
  {
    for (i in 1:length(unique(most.likely.cluster)))
    {
      
      cluster.muts <- TCGA.plot[TCGA.plot$mutation_id%in%names(which(most.likely.cluster==i)),,drop=FALSE]
      points(cluster.muts$Start_pos
             ,cluster.muts$mutCopyNum
             ,cex=0.7
             ,pch=pch.vals[i]
             ,col=cols[i]
             ,bg =cols.opac[i] ) # ca
      
      
      
    }
    
  }
  if (!is.na(cluster[1]))
  {
    for (i in cluster)
    {
      
      cluster.muts <- TCGA.plot[TCGA.plot$mutation_id%in%names(which(most.likely.cluster==i)),,drop=FALSE]
      points(cluster.muts$Start_pos
             ,cluster.muts$mutCopyNum
             ,cex=0.7
             ,pch=pch.vals[as.numeric(i)]
             ,col=cols[i]
             ,bg =cols.opac[as.numeric(i)] ) # ca
      
      
      
    }
    
  }
  
  
  
  
  
  
  
  mtext(side=3
        ,at=fun.chrom.mean(seg.mat.plot[,2:4])
        ,text=sort(unique(seg.mat.plot[,2]))
        ,cex=seq(0.6,0.4,length.out=length(unique(seg.mat.plot[,2])))
        ,line=-1
        #,lwd=0.5
  )
  
  
}

determinePhylogeny <- function(regionList
                               ,mostLikelyClusters
                               ,driverCat=c(1:3)
                               ,ccf.type ='phylo.ccf'
                               ,mutTable
                               #,pyclone.results
                               #,min.thresh=0.05
)
{
  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(wordcloud))
  no.optima = length(unique(mostLikelyClusters))
  max.cols = 12
  require(RColorBrewer)
  cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
  cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
  cols.opac      = paste(cols,'99',sep="")
  
  pch.vals       = c()
  k <- 21
  m <- 1  
  while(m<=no.optima)
  {
    pch.vals    = c(pch.vals,k)
    if (m==max.cols)
    {
      k <- k+1
    }
    m <- m+1
  }
  
  
  
  line.type       = c()
  k <- 1
  m <- 1  
  while(m<=no.optima)
  {
    line.type    = c(line.type,k)
    if (m==max.cols)
    {
      k <- k+1
    }
    m <- m+1
  }
  
  pch.vals       = pch.vals
  

  if(ccf.type=='phylo.ccf.PyClone')
  {
    for (region in names(regionList))
    {
      region.mut          <- regionList[[region]]
      layout(rbind(1,2))
      par(mar=c(0.2,5,5,5))
      
      for (cl in unique(mostLikelyClusters))
      {
        
        cl.nr             <-  length(which(unlist(region.mut[,'phyloCCF_PyClone'])!=0))
        clust.ccf         <-  as.numeric(unlist(region.mut[unlist(region.mut$mutation_id)%in%names(which(mostLikelyClusters==cl)),'phyloCCF_PyClone']))
        cl.size           <-  log((length(clust.ccf)/cl.nr)*1000)
        
        if(length(clust.ccf)==1)
        {
          next;
        }
        
        if(max(clust.ccf)==0)
        {
          next;
        }
        ds.cl               <- density(clust.ccf,from=0,to=1.15,kernel='gaussian',bw=0.05)
        ds.cl$y             <- (ds.cl$y/max(ds.cl$y))*cl.size
        ds.cl$x             <- rev(ds.cl$x)
        
        if(cl==1)
        {
          plot(ds.cl
               ,xlim=c(0.1,1)
               ,ylim=c(0,log(1000))
               ,xaxt='n'
               ,yaxt='n'
               ,xaxs='i'
               ,yaxs='i'
               ,xlab=""
               ,ylab='Density (a.u.)'
               ,main=region
               ,lwd=1
               ,col=cols[cl]
               ,cex=0.85
               ,cex.lab=0.85
               ,cex.main=0.9
          )
        }
        
        if(cl!=1)
        {
          lines(ds.cl
                ,xlim=c(0.1,1)
                #,ylim=c(0,ylim)
                ,xaxt='n'
                ,yaxt='n'
                #,xaxs='i'
                #,yaxs='i'
                ,xlab=""
                ,ylab='Density (a.u.)'
                ,main=region
                ,lwd=1
                ,col=cols[cl]
                ,cex.lab=0.85
                ,lty=line.type[cl]
                ,cex=0.85
          )
        }
        
        
      }
      
      #abline(v=0.15)
      
      par(mar=c(5,5,0.2,5))
      rownames(region.mut) <- region.mut$mutation_id
      
      plot(1-(unlist(region.mut[names(mostLikelyClusters),'phyloCCF_PyClone']))
           ,log(as.numeric(unlist(region.mut[names(mostLikelyClusters),'var_counts']))
                +as.numeric(unlist(region.mut[names(mostLikelyClusters),'ref_counts'])))
           ,ylim=c(log(5),log(1000))
           ,xlim=c(-0.05,0.85)
           ,xaxt='n'
           ,yaxt='n'
           ,xaxs='i'
           ,yaxs='i'
           ,xlab="Cancer Cell Fraction"
           ,ylab="Tumour Coverage"
           ,cex.axis=0.85
           ,cex.lab=0.75
           ,main=""
           ,lwd=3
           ,cex=0.85
           ,pch=pch.vals[mostLikelyClusters]
           ,col=cols.opac[mostLikelyClusters]
           ,bg=cols.opac[mostLikelyClusters]
      )
      
      #abline(v=-0)
      
      driver.table <- mutTable[mutTable$driverCategory%in%driverCat,,drop=FALSE]
      driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
      driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
      driver.table$driverCategoryNumeric<-as.numeric(driver.table$driverCategory)
      if(nrow(driver.table)>1)
      {
        textplot(1-(unlist(region.mut[rownames(driver.table),'phyloCCF_PyClone']))
                 ,log(as.numeric(unlist(region.mut[rownames(driver.table),'var_counts']))
                      +as.numeric(unlist(region.mut[rownames(driver.table),'ref_counts'])))
                 ,words=c(driver.table$Gene.refGene)
                 ,new=FALSE
                 ,col=gray(driver.table$driverCategoryNumeric/7)
                 ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                 ,cex=0.7
                 ,show.lines=TRUE
        )
      }
      
      if(nrow(driver.table)==1)
      {
        textplot(c(1-(unlist(region.mut[rownames(driver.table),'phyloCCF_PyClone'])),0)
                 ,c(log(as.numeric(unlist(region.mut[rownames(driver.table),'var_counts']))
                        +as.numeric(unlist(region.mut[rownames(driver.table),'ref_counts']))),0)
                 ,words=c(c(driver.table$Gene.refGene),'')
                 ,new=FALSE
                 ,col=gray(driver.table$driverCategoryNumeric/7)
                 ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                 ,cex=0.7
                 ,show.lines=TRUE
        )
      }
      
      
      
      
      
      axis(side=1,at=seq(0,1,by=0.2),labels=rev(seq(0,1,by=0.2)),cex.axis=0.75)
      cov.seq <- c(5,10,20,50,100,200,500,1000,cex.axis=0.75)
      axis(side=2,at=log(cov.seq),labels=cov.seq,las=2,cex.axis=0.75)
      
    }
  }
  
  
  
}

plot.pyclone.clusters <- function(regionList
                                  ,mostLikelyClusters
                                  ,driverCat
                                  ,ccf='sanger'
                                  ,colour.choice=NA
                                  ,plot.median=TRUE)
  
{
  suppressPackageStartupMessages(library(plyr))
  if (ccf=='absolute')
  {
    
    #its hard to distinguish more than 8 different colours
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    require(RColorBrewer)
    cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    
    
    for(i in 1:(length(regionList)-1)){
      for(j in (i+1):length(regionList)){
        regionX = regionList[[i]]
        regionY = regionList[[j]]
        #regionX$region <- regionList[i]
        #regionY$region <- regionList[j]
        rownames(regionX) <- regionX$mutation_id
        rownames(regionY) <- regionY$mutation_id
        # let's make sure we're looking at the same thing
        tmp     <- intersect(rownames(regionX),rownames(regionY))
        tmp     <- intersect(names(mostLikelyClusters),tmp)
        regionX <- regionX[tmp,,drop=FALSE]
        regionY <- regionY[tmp,,drop=FALSE]
        mostLikelyClusters <- mostLikelyClusters[tmp]
        
        layout(cbind(rbind(c(1,1,1,1,1)),2))
        par(mar=c(5,5,7,1))
        plot(regionX$absolute.ccf,regionY$absolute.ccf
             ,type = "n"
             , main = patient
             , xlab = paste(names(regionList)[i]," CCF",sep="")
             , ylab = paste(names(regionList)[j]," CCF",sep="")
             , xlim = c(-0.05,1.25)
             , ylim = c(-0.05,1.25)
             ,xaxt='n'
             ,yaxt='n'
             ,cex.axis=1.3)
        
        abline(h=1,lty='dashed')
        abline(v=1,lty='dashed')
        axis(side = 1,at = seq(0,1,by=0.2))
        axis(side = 2,at = seq(0,1,by=0.2),las=2) 
        # c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.optima){
          
          points(unlist(regionX$absolute.ccf)[mostLikelyClusters==n],unlist(regionY$absolute.ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        for(n in 1:no.optima){
          
          points(unlist(regionX$absolute.ccf)[mostLikelyClusters==n],unlist(regionY$absolute.ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        
        
        x.vals <- c()
        y.vals <- c()
        colors <- c()
        labs   <- c()
        for (n in 1:no.optima)
        {
          
          x.vals <- c(x.vals,median(unlist(regionX$absolute.ccf)[mostLikelyClusters==n]))
          y.vals <- c(y.vals,median(unlist(regionY$absolute.ccf)[mostLikelyClusters==n]))
          labs   <- c(labs,n)
          points(median(unlist(regionX$absolute.ccf)[mostLikelyClusters==n])
                 ,median(unlist(regionY$absolute.ccf)[mostLikelyClusters==n])
                 ,pch=21
                 ,col=cols[n]
                 ,bg=cols[n]
                 ,cex=1.2*2.5)
          
          
          #           text(median(unlist(regionX$absolute.ccf)[mostLikelyClusters==n])
          #                ,median(unlist(regionY$absolute.ccf)[mostLikelyClusters==n])
          #                ,labels=n
          #                ,col='white',offset=0,cex=0.7)
          
        }
        
        require(wordcloud)
        
        if(length(x.vals)>1)
        {
          textplot(x.vals,y.vals,words = labs,new=FALSE
                   ,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white')
                   ,show.lines=TRUE) 
        }
        
        if(length(x.vals)==1)
        {
          text(x.vals,y.vals,labels = labs)
        }
        
        
        
        
        
        # let's add the drivers
        driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
        driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
        driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
        driver.table$driverCategoryNumeric<-as.numeric(revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2,"2"=3, "3"=4))) 
        if(nrow(driver.table)>1)
        {
          #         text(unlist(regionX[rownames(driver.table),'absolute.ccf'])
          #              ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
          #              ,labels=c(driver.table$Gene.refGene)
          #              ,cex=0.7
          #              ,pos=4
          #              ,offset=0.5)
          
          textplot(unlist(regionX[rownames(driver.table),'absolute.ccf'])
                   ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
                   ,words=c(driver.table$Gene.refGene)
                   ,new=FALSE
                   ,col=gray(driver.table$driverCategoryNumeric/7)
                   ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                   ,cex=0.7
                   ,show.lines=FALSE
                   ,pos=4)
        }
        if(nrow(driver.table)==1)
        {
          #         text(unlist(regionX[rownames(driver.table),'absolute.ccf'])
          #              ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
          #              ,labels=c(driver.table$Gene.refGene)
          #              ,cex=0.7
          #              ,pos=4
          #              ,offset=0.5)
          
          text(x = unlist(regionX[rownames(driver.table),'absolute.ccf']),
               ,y=unlist(regionY[rownames(driver.table),'absolute.ccf'])
               ,labels=c(driver.table$Gene.refGene)
               #,new=FALSE
               ,col=gray(driver.table$driverCategoryNumeric/7)
               ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
               ,cex=0.7
               #,show.lines=FALSE
               ,pos=4)
        }
        
        
        
        
        #plot.new()
        par(mar=c(3,0.5,7,0))
        plot(1, type="n"
             , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
        for (k in 1:no.optima)
        {
          points(0
                 ,k
                 ,col=cols.opac[k]
                 ,bg = cols[k]
                 ,pch=pch.vals[k]
                 ,cex=1.5)
          text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
               ,pos=4)
        }
        
        text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
        
      }
    }
  }
  
  if (ccf=='sanger')
  {
    #its hard to distinguish more than 8 different colours
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    require(RColorBrewer)
    cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    
    
    for(i in 1:(length(regionList)-1)){
      for(j in (i+1):length(regionList)){
        regionX = regionList[[i]]
        regionY = regionList[[j]]
        regionX$region <- regions.to.use[i]
        regionY$region <- regions.to.use[j]
        rownames(regionX) <- regionX$mutation_id
        rownames(regionY) <- regionY$mutation_id
        # let's make sure we're looking at the same thing
        tmp     <- intersect(rownames(regionX),rownames(regionY))
        tmp     <- intersect(names(mostLikelyClusters),tmp)
        regionX <- regionX[tmp,,drop=FALSE]
        regionY <- regionY[tmp,,drop=FALSE]
        mostLikelyClusters <- mostLikelyClusters[tmp]
        
        layout(cbind(rbind(c(1,1,1,1,1)),2))
        par(mar=c(5,5,7,1))
        plot(regionX$ccf,regionY$ccf
             ,type = "n"
             , main = patient
             , xlab = paste(names(regionList)[i]," CCF",sep="")
             , ylab = paste(names(regionList)[j]," CCF",sep="")
             , xlim = c(-0.05,1.5)
             , ylim = c(-0.05,1.5)
             ,xaxt='n'
             ,yaxt='n'
             ,cex.axis=1.3)
        
        abline(h=1,lty='dashed')
        abline(v=1,lty='dashed')
        axis(side = 1,at = seq(0,1,by=0.2))
        axis(side = 2,at = seq(0,1,by=0.2),las=2) 
        # c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.optima){
          
          points(unlist(regionX$ccf)[mostLikelyClusters==n],unlist(regionY$ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        for(n in 1:no.optima){
          
          points(unlist(regionX$ccf)[mostLikelyClusters==n],unlist(regionY$ccf)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        
        
        x.vals <- c()
        y.vals <- c()
        colors <- c()
        labs   <- c()
        for (n in 1:no.optima)
        {
          
          x.vals <- c(x.vals,median(unlist(regionX$ccf)[mostLikelyClusters==n]))
          y.vals <- c(y.vals,median(unlist(regionY$ccf)[mostLikelyClusters==n]))
          labs   <- c(labs,n)
          points(median(unlist(regionX$ccf)[mostLikelyClusters==n])
                 ,median(unlist(regionY$ccf)[mostLikelyClusters==n])
                 ,pch=21
                 ,col=cols[n]
                 ,bg=cols[n]
                 ,cex=1.2*2.5)
          
          
          #           text(median(unlist(regionX$absolute.ccf)[mostLikelyClusters==n])
          #                ,median(unlist(regionY$absolute.ccf)[mostLikelyClusters==n])
          #                ,labels=n
          #                ,col='white',offset=0,cex=0.7)
          
        }
        
        require(wordcloud)
        textplot(x.vals,y.vals,words = labs,new=FALSE
                 ,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white')
                 ,show.lines=TRUE)
        
        
        
        
        
        # let's add the drivers
        driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
        driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
        driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
        driver.table$driverCategoryNumeric<-as.numeric(revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2,"2"=3, "3"=4)))
        
        
        if(nrow(driver.table)>=1)
        {
          #         text(unlist(regionX[rownames(driver.table),'absolute.ccf'])
          #              ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
          #              ,labels=c(driver.table$Gene.refGene)
          #              ,cex=0.7
          #              ,pos=4
          #              ,offset=0.5)
          
          textplot(unlist(regionX[rownames(driver.table),'absolute.ccf'])
                   ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
                   ,words=c(driver.table$Gene.refGene)
                   ,new=FALSE
                   ,col=gray(driver.table$driverCategoryNumeric/7)
                   ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                   ,cex=0.7
                   ,show.lines=FALSE
                   ,pos=4)
        }
        
        
        
        
        #plot.new()
        par(mar=c(3,0.5,7,0))
        plot(1, type="n"
             , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
        for (k in 1:no.optima)
        {
          points(0
                 ,k
                 ,col=cols.opac[k]
                 ,bg = cols[k]
                 ,pch=pch.vals[k]
                 ,cex=1.5)
          text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
               ,pos=4)
        }
        
        text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
        
      }
    }
  }
  
  if (ccf=='phylo')
  {
    #its hard to distinguish more than 8 different colours
    no.optima = length(unique(mostLikelyClusters))
    max.cols = 12
    require(RColorBrewer)
    cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
    cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
    cols.opac      = paste(cols,'99',sep="")
    
    if(!is.na(colour.choice))
    {
      cols        = rep(colour.choice,length(cols)) 
      cols.opac   = rep(colour.choice,length(cols)) 
    }
    
    pch.vals       = c()
    k <- 21
    m <- 1  
    while(m<=no.optima)
    {
      pch.vals    = c(pch.vals,k)
      if (m==max.cols)
      {
        k <- k+1
      }
      m <- m+1
    }
    
    pch.vals     = pch.vals
    
    
    for(i in 1:(length(regionList)-1)){
      for(j in (i+1):length(regionList)){
        regionX = regionList[[i]]
        regionY = regionList[[j]]
        regionX$region <- regions.to.use[i]
        regionY$region <- regions.to.use[j]
        rownames(regionX) <- regionX$mutation_id
        rownames(regionY) <- regionY$mutation_id
        # let's make sure we're looking at the same thing
        tmp     <- intersect(rownames(regionX),rownames(regionY))
        tmp     <- intersect(names(mostLikelyClusters),tmp)
        regionX <- regionX[tmp,,drop=FALSE]
        regionY <- regionY[tmp,,drop=FALSE]
        mostLikelyClusters <- mostLikelyClusters[tmp]
        
        layout(cbind(rbind(c(1,1,1,1,1)),2))
        par(mar=c(5,5,7,1))
        plot(regionX$phyloCCF,regionY$phyloCCF
             ,type = "n"
             , main = patient
             , xlab = paste(names(regionList)[i]," CCF",sep="")
             , ylab = paste(names(regionList)[j]," CCF",sep="")
             , xlim = c(-0.05,1.5)
             , ylim = c(-0.05,1.5)
             ,xaxt='n'
             ,yaxt='n'
             ,cex.axis=1.3)
        
        abline(h=1,lty='dashed')
        abline(v=1,lty='dashed')
        axis(side = 1,at = seq(0,1,by=0.2))
        axis(side = 2,at = seq(0,1,by=0.2),las=2) 
        # c(0,max(plot.data[,i])*1.25))
        for(n in 1:no.optima){
          
          points(unlist(regionX$phyloCCF)[mostLikelyClusters==n],unlist(regionY$phyloCCF)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        for(n in 1:no.optima){
          
          points(unlist(regionX$phyloCCF)[mostLikelyClusters==n],unlist(regionY$phyloCCF)[mostLikelyClusters==n]
                 ,bg=cols.opac[n]
                 ,col = cols[n]
                 ,pch=pch.vals[n]
                 ,cex=1.2)
          
        }
        
        if(plot.median)
        {
          
          
          x.vals <- c()
          y.vals <- c()
          colors <- c()
          labs   <- c()
          for (n in 1:no.optima)
          {
            
            x.vals <- c(x.vals,median(unlist(regionX$phyloCCF)[mostLikelyClusters==n]))
            y.vals <- c(y.vals,median(unlist(regionY$phyloCCF)[mostLikelyClusters==n]))
            labs   <- c(labs,n)
            points(median(unlist(regionX$phyloCCF)[mostLikelyClusters==n])
                   ,median(unlist(regionY$phyloCCF)[mostLikelyClusters==n])
                   ,pch=21
                   ,col=cols[n]
                   ,bg=cols[n]
                   ,cex=1.2*2.5)
            
            
            #           text(median(unlist(regionX$absolute.phyloCCF)[mostLikelyClusters==n])
            #                ,median(unlist(regionY$absolute.phyloCCF)[mostLikelyClusters==n])
            #                ,labels=n
            #                ,col='white',offset=0,cex=0.7)
            
          }
          
          require(wordcloud)
          if(length(x.vals)>1)
          {        textplot(x.vals,y.vals,words = labs,new=FALSE
                            ,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white')
                            ,show.lines=TRUE)
                   
          }
          if(length(x.vals)>1)
          {        
            text(x.vals,y.vals,labels = labs,col=ifelse(x.vals==0.01&y.vals==0.01,'black','white'))
            
          }
          
          
        }
        
        
        
        # let's add the drivers
        driver.table <- mut.table[mut.table$driverCategory%in%driverCat,,drop=FALSE]
        driver.table <- driver.table[driver.table$mutation_id%in%names(mostLikelyClusters),,drop=FALSE]
        driver.table$Gene.refGene <- do.call(rbind,strsplit(driver.table$Gene.refGene,split="\\("))[,1]
        driver.table$driverCategoryNumeric<-as.numeric(revalue(as.character(driver.table$driverCategory), c("1A"=1,"1"=2,"2"=3, "3"=4)))
        
        
        if(nrow(driver.table)>1)
        {
          #         text(unlist(regionX[rownames(driver.table),'absolute.phyloCCF'])
          #              ,unlist(regionY[rownames(driver.table),'absolute.phyloCCF'])
          #              ,labels=c(driver.table$Gene.refGene)
          #              ,cex=0.7
          #              ,pos=4
          #              ,offset=0.5)
          
          textplot(unlist(regionX[rownames(driver.table),'phyloCCF'])
                   ,unlist(regionY[rownames(driver.table),'phyloCCF'])
                   ,words=c(driver.table$Gene.refGene)
                   ,new=FALSE
                   ,col=gray(driver.table$driverCategoryNumeric/7)
                   ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
                   ,cex=0.7
                   ,show.lines=FALSE
                   ,pos=4)
        }
        
        if(nrow(driver.table)==1)
        {
          #         text(unlist(regionX[rownames(driver.table),'absolute.ccf'])
          #              ,unlist(regionY[rownames(driver.table),'absolute.ccf'])
          #              ,labels=c(driver.table$Gene.refGene)
          #              ,cex=0.7
          #              ,pos=4
          #              ,offset=0.5)
          
          text(x = unlist(regionX[rownames(driver.table),'phyloCCF']),
               ,y=unlist(regionY[rownames(driver.table),'phyloCCF'])
               ,labels=c(driver.table$Gene.refGene)
               #,new=FALSE
               ,col=gray(driver.table$driverCategoryNumeric/7)
               ,font=ifelse(driver.table$driverCategoryNumeric==1|driver.table$driverCategoryNumeric==2,2,1)
               ,cex=0.7
               #,show.lines=FALSE
               ,pos=4)
        }
        
        
        
        #plot.new()
        par(mar=c(3,0.5,7,0))
        plot(1, type="n"
             , axes=F, xlab="", ylab="",ylim=c(-0.5,no.optima+1),xlim=c(-0.25,5))
        for (k in 1:no.optima)
        {
          points(0
                 ,k
                 ,col=cols.opac[k]
                 ,bg = cols[k]
                 ,pch=pch.vals[k]
                 ,cex=1.5)
          text(0,k,cex=0.65, labels=paste(k," (",length(which(mostLikelyClusters==k))," SNVs)",sep="")
               ,pos=4)
        }
        
        text(-0.25,no.optima+0.75,pos=4,labels='Clusters:',cex=0.8)
        
      }
    }
  }
  
}


##
mtexti <- function(text, side, off = 0.25,
                   srt = if(side == 2) 90  else
                     if(side == 4) 270 else 0, ...) {
  # dimensions of plotting region in user units
  usr <- par('usr')
  # dimensions of plotting region in inches
  pin <- par('pin')
  # user units per inch
  upi <- c(usr[2]-usr[1],
           usr[4]-usr[3]) / pin
  # default x and y positions
  xpos <- (usr[1] + usr[2])/2
  ypos <- (usr[3] + usr[4])/2
  if(1 == side)
    ypos <- usr[3] - upi[2] * off
  if(2 == side)
    xpos <- usr[1] - upi[1] * off
  if(3 == side)
    ypos <- usr[4] + upi[2] * off
  if(4 == side)
    xpos <- usr[2] + upi[1] * off
  text(x=xpos, y=ypos, text, xpd=NA, srt=srt, ...)
}

##
fun.chrom.length <- function(seg){
  # This function returns the length of chromosome 1,...,22
  # seg: first three columns of the minimum conistent region matrix
  
  chrom.length <- c()
  for(i in sort(as.numeric(unique(seg[,1])))){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.length <- c(chrom.length, sub.chrom[nrow(sub.chrom), 3])
    
  }
  return(chrom.length)
  
}

##   

fun.add.chrom <- function(seg, chrom.length){
  # This function adds the length of chromosome 1,...,(i - 1) to chromsomes 2,...,22
  # seg: first three columns of the minimum conistent region matrix
  # chrom.length: output of fun.chrom.length(...)
  if(1 %in% seg[, 1])
    seg.tmp <-  subset(seg, seg[, 1] == 1)
  else
    seg.tmp <- matrix(NA, nr = 0, nc = ncol(seg))  
  for(i in 2:max(max(as.numeric(seg[,1])),3)){
    if(!(i %in% seg[, 1]))
      next()
    sub.chrom <- subset(seg, seg[, 1] == i)
    sub.chrom[, 2] <- as.numeric(sub.chrom[, 2]) + sum(as.numeric(chrom.length[1:(i - 1)]))
    sub.chrom[, 3] <- as.numeric(sub.chrom[, 3]) + sum(as.numeric(chrom.length[1:(i - 1)]))
    seg.tmp <- rbind(seg.tmp, sub.chrom)
  }
  return(seg.tmp)
  
} 

##

fun.smooth.score <- function(score, seg){
  # This function applies smoothing (rollmean) to chromosome 1,...,22 separetely
  # score: score to be smoothed (e.g. correlation)
  # seg: Corresponding positions
  # window: window size for rollmean
  # The R package "zoo" is required
  score.tmp <- c()
  for(i in 1:22)
    
    score.tmp <- c(score.tmp, rollmean(score[seg[, 1] == i], k = round(length(score[seg[, 1] == i])/10), na.pad = T))
  
  
  return(score.tmp)
  
  
}  

##

fun.chrom.mean <- function(seg){
  # This function finds the middle of the positions for each chromosome.
  # seg: Matrix with three columns, containing the positions.
  chrom.mean <- c()
  for(i in sort(as.numeric(unique(seg[,1])))){
    sub.chrom <- subset(seg, seg[, 1] == i)
    chrom.mean <- c(chrom.mean, (min(sub.chrom[, 2]) + max(sub.chrom[, 3]))/2)
    
  }
  return(chrom.mean)
  
}



