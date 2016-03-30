# This script is designed to estimate the cancer cell fraction using PyClone

####################################################################################################### Load libraries
######################################################################################################################
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(gdata))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(wordcloud))
suppressPackageStartupMessages(library(sequenza))
suppressPackageStartupMessages(library(bootstrap))
suppressPackageStartupMessages(library(boot))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(wordcloud))

options(expressions = 500000)

############################################################################################### Command line arguments
######################################################################################################################
cmdArgs                  <- commandArgs(trailingOnly = TRUE);
patient                  <- cmdArgs[1];
saveDir                  <- cmdArgs[2];
ascatDir                 <- cmdArgs[3];
snvDir                   <- cmdArgs[4];
PyClone                  <- cmdArgs[5];
template.config.yaml     <- cmdArgs[6];
PyCloneFunctions         <- cmdArgs[7]

interactive          <- FALSE
run.pyclone          <- TRUE

# IMPORTANT #
# The following script is fully dependent on PyClone
# For information on PyClone, see https://bitbucket.org/aroth85/pyclone/wiki/Home


if(interactive)
{
  patient              <- "LMS025"
  saveDir              <- "/farm/home/lr-tct-lif/mcgran01/clonalDissectionPipeLine/" 
  ascatDir             <- "/farm/home/lr-tct-lif/mcgran01/clonalDissectionPipeLine/ExampleFiles/ASCAT/"
  snvDir               <- "/farm/home/lr-tct-lif/mcgran01/clonalDissectionPipeLine/ExampleFiles/SNV/"  
  PyClone              <- "/farm/babs/redhat6/software/python-2.7.3/bin/PyClone"
  template.config.yaml <- "/farm/home/lr-tct-lif/mcgran01/clonalDissectionPipeLine/ExampleFiles/template.config.yaml"
  PyCloneFunctions     <- "/farm/home/lr-tct-lif/mcgran01/clonalDissectionPipeLine/clonalDissectionFunctions.R"
}

#######  Source helper functions ######

source(PyCloneFunctions)


########  Define your own parameters  ##########
# it is assumed that you have applied filters to your mutation calls
# e.g. depth and p.val thresholds
# The following script essentially assumes all mutations are correct

sample             <- patient
new.dir            <- paste(saveDir,"/",patient,"_PyClone_phylo/",sep="")
driver.cat         <- 1:3

options(expressions = 10000)
burn.in              <- 1000


###### Run script #####################

cat('\nModified PyClone analysis of the following patient:\n')
print(patient)
cat('\n')


if( !file.exists(new.dir))
{
  if( !dir.create(new.dir, recursive = TRUE) )
  {
    stop("Unable to create root directory.\n")
  }
}

# first of all, let's load the data we need ####
# mutation.table
mut.table.loc     <- list.files(snvDir,full.names=TRUE)
if(length(mut.table.loc)!=1)
{
  stop("Unable to find one mut.table.loc.\n")  
}

# Load the mutation table
mut.table         <- read.table(mut.table.loc
                                ,sep="\t"
                                ,stringsAsFactors=FALSE
                                ,header=TRUE)

# load the segmented copy number table
ascatsegRData     <- grep('ascat.seg.RData$',list.files(ascatDir,full.names=TRUE),value=TRUE)
segRData          <- load(ascatsegRData)
seg.mat.copy      <- get(segRData)


# at present we can only reliably use autosomal mutations, restrict to chromosomes 1 to 22 ####
mut.table      <- mut.table[mut.table$chr%in%1:22,,drop=FALSE]
seg.mat.copy   <- seg.mat.copy[seg.mat.copy$chr%in%1:22,,drop=FALSE]
seg.mat.copy$sample <- gsub("-","\\.",seg.mat.copy$sample)
colnames(seg.mat.copy)[1]  <- c('SampleID')
seg.mat.phylo      <- create.subclonal.copy.number(seg.mat.copy = seg.mat.copy,min.subclonal = 0.01)
mut.table$mutation_id  <- paste(sample,mut.table$chr,mut.table$start,mut.table$ref,sep=":")
rownames(mut.table)    <- mut.table$mutation_id
regions.to.use         <- seg.mat.phylo$SampleID[1]
region                 <- regions.to.use

if(nrow(mut.table)==0)
{
  stop('No mutations passed filtering, stopping PyClone phylo clustering')
}

patient.list          <- list()
phylo.region.list     <- list()
cellularity           <- rep(NA,length(regions.to.use))


# Let's plot the copy number. 
patient.list      <- list()
phylo.region.list <- list()
cellularity       <- rep(NA,length(regions.to.use))

pdf(paste(new.dir, sample, ".subclonal.mut.cpn.pdf",sep=""),width=6,height=3)
{
  lyout <- rbind(c(rep(1,9),2))
  layout(lyout)
  
  region.mut.table <- mut.table
  region.seg.copy  <- seg.mat.copy[seg.mat.copy$SampleID%in%region,,drop=FALSE]
  region.seg.phylo <- seg.mat.phylo[seg.mat.phylo$SampleID%in%region,,drop=FALSE]
  pyclone.table    <- data.frame(t(sapply(1:nrow(region.mut.table),identify.subclonal.mut.copy.number.ascat,region.mut.table,region.seg.phylo,region,sample))
                                 ,stringsAsFactors=FALSE)
  pyclone.table    <- pyclone.table[!is.na(pyclone.table$minor_cn),]
  pyclone.table    <- pyclone.table[!is.na(pyclone.table$ref_counts),]
  pyclone.table    <- pyclone.table[!duplicated(pyclone.table$mutation_id),]
  
  
  # let's load the purity estimate from VAF purity
   sample.purity   <- region.seg.copy$ACF[1]

  pyclone.table   <- pyclone.table[as.numeric(pyclone.table$ref_counts)+as.numeric(pyclone.table$var_counts)>=1,,drop=FALSE]
  region.phyloCCF                                                          <- subclonalDissection(region=region,complete.mutation.table=pyclone.table,purity=sample.purity,order.by.pos = TRUE)
  phylo.region.list[[region]]                                              <- region.phyloCCF
  
  cellularity[region]                                                      <- sample.purity
  
  
  # add extra column to 
  sub.mat.copy               <- region.seg.copy
  sub.mat.copy$Copy.Number   <- sub.mat.copy$cnTotal
  sub.mat.copy$min.allele    <- apply(cbind(sub.mat.copy$nMinor,sub.mat.copy$nMajor),1,min)
 
  colnames(sub.mat.copy)[2]  <- 'Chromosome'
  colnames(sub.mat.copy)[3]  <- 'StartPosition'
  colnames(sub.mat.copy)[4]  <- 'EndPosition'
  colnames(sub.mat.copy)[5]  <- 'nr.probes'
  
  
  
  par(mar=c(5,4,5,0.2))

  
  par(lend=1)
  plot.EarlyOrLate.raw(seg.mat.patient=sub.mat.copy
                       ,TCGA.earlyLate=region.phyloCCF
                       ,TCGA.purity=sample.purity
                       ,sub.clonal=1
  )
  mtext(unlist(strsplit(region,split="\\."))[2], side=2,cex=1,line=2,las=2)
  
  ds <- density(ifelse(as.numeric(region.phyloCCF$mutCopyNum)>5,5,as.numeric(region.phyloCCF$mutCopyNum)))
  ds1 <- ds
  ds1$x <- ds$y
  ds1$y <- ds$x
  par(mar=c(5,0,5,4))

  A <- hist(ifelse(as.numeric(region.phyloCCF$mutCopyNum)>5,5,as.numeric(region.phyloCCF$mutCopyNum)),breaks=seq(-0.25,6,by=0.1),plot=FALSE)
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
dev.off()


# let's run PyClone, but correct for copy number before running

region.mut.table <- mut.table
region.seg.copy  <- seg.mat.phylo[seg.mat.phylo$SampleID%in%region,,drop=FALSE]
pyclone.table    <- data.frame(t(sapply(1:nrow(region.mut.table),identify.subclonal.mut.copy.number.ascat,region.mut.table,region.seg.copy,region,sample))
                                 ,stringsAsFactors=FALSE)
  
na.mutations           <- pyclone.table[is.na(pyclone.table$minor_cn),,drop=FALSE]
pyclone.table          <- pyclone.table[!is.na(pyclone.table$minor_cn),,drop=FALSE]
loss.mutations         <- pyclone.table[as.numeric(pyclone.table$major_cn)==0|(as.numeric(pyclone.table$var_counts)+as.numeric(pyclone.table$ref_counts)==0),]
error.muts             <- rbind(na.mutations,loss.mutations)
error.muts             <- unlist(na.mutations$mutation_id,loss.mutations$mutation_id)
error.muts.table       <- paste(new.dir,"/",region,".error.muts.tsv",sep="")
  
  
  
if (length(error.muts)>=1)
{
    write.table  (error.muts
                  ,sep="\t"
                  ,quote=FALSE
                  ,col.names=TRUE
                  ,row.names=FALSE
                  ,file=error.muts.table
    )
    
  }
  
  
  # Make sure everything in PyClone table is there. 
  pyclone.table    <- pyclone.table[!is.na(pyclone.table$minor_cn),]
  pyclone.table    <- pyclone.table[!is.na(pyclone.table$ref_counts),]
  pyclone.table    <- pyclone.table[!duplicated(pyclone.table$mutation_id),]
  pyclone.table    <- pyclone.table[as.numeric(pyclone.table$major_cn)>=1,]
  pyclone.table    <- pyclone.table[!is.na(pyclone.table$minor_cn),]
  
  
  # Now, let's check what the cancer cell fraction estimates are for this region
  region.ccf               <- phylo.region.list[[region]]
  region.ccf               <- data.frame(region.ccf,stringsAsFactors=FALSE)
  rownames(region.ccf)     <- region.ccf$mutation_id
  tmp                      <- intersect(unlist(pyclone.table$mutation_id),unlist(region.ccf$mutation_id))
  rownames(pyclone.table)  <- pyclone.table$mutation_id
  pyclone.table            <- pyclone.table[tmp,,drop=FALSE]
  region.ccf               <- region.ccf[tmp,,drop=FALSE]
  
  tmp                      <- round(((unlist(pyclone.table$var_counts))/(unlist(region.ccf$phyloCCF)/2))-unlist(pyclone.table$var_counts))
  tmp[is.na(tmp)]          <- unlist(pyclone.table$ref_counts[(is.na(tmp))])
  pyclone.table$ref_counts <- tmp
  pyclone.table$minor_cn   <- 0
  pyclone.table$major_cn   <- 2
  pyclone.table$ref_counts <- apply(cbind(pyclone.table$ref_counts,2),1,max)
  
  sample.purity    <- 0.5
  
  
  pyclone.tsv   <- paste(new.dir,"/",region,".tsv",sep="")
  pyclone.out   <- apply(pyclone.table,2,as.character)
  write.table  (pyclone.out
                ,sep="\t"
                ,quote=FALSE
                ,col.names=TRUE
                ,row.names=FALSE
                ,file=pyclone.tsv
  )
  
  #   Run PyClone build_mutations_file TSV_FILE where TSV_FILE is the input file you have created.
  pyclone.yaml   <- paste(new.dir,"/",region,".yaml",sep="")
  
  cmd <- paste(PyClone
               ," build_mutations_file "
               ,pyclone.tsv
               ," "
               ,pyclone.yaml
               ," --var_prior BB"
               ," --ref_prior normal"
               ,sep="")
  cat('\n')
  
  if(run.pyclone)
  {
    cat(cmd)
    system(cmd)
    
  }
  
  cat('\n')
  


# Next let's create a configuration file
pyclone.config.yaml <- paste(new.dir,"/",sample,".config.yaml",sep="")
pyclone.config      <- readLines(template.config.yaml)
start.samples       <- (grep("samples",pyclone.config)+1)
end.samples         <- length(pyclone.config)

sample.lines        <- pyclone.config[start.samples:end.samples]
pyclone.config      <- pyclone.config[-c(start.samples:end.samples)]
pyclone.config      <- gsub("working.directory.location", new.dir, pyclone.config)

write.table(pyclone.config
            , file=pyclone.config.yaml
            ,col.names=FALSE
            ,row.names=FALSE
            ,quote=FALSE)


for (region in names(phylo.region.list))
{
  
  sample.config       <- gsub("TCGA.barcode", region, sample.lines)
  pyclone.yaml        <- paste(region,".yaml",sep="")
  sample.config       <- gsub("mutations.yaml", pyclone.yaml, sample.config)
  region.purity       <- 0.5
  
  sample.config       <- gsub("value: 1.0", paste("value: ", signif(region.purity,3), sep=""), sample.config)
  sample.config       <- sample.config[1:8]
  
  write.table(sample.config,file=pyclone.config.yaml
              ,append=TRUE
              ,quote=FALSE
              ,col.names=FALSE
              ,row.names=FALSE)
  
}

# Run the PyClone analysis using the PyClone analyse 
cmd <- paste(PyClone
             ," analyse "
             ,pyclone.config.yaml
             ,sep=""
)
cat('\n')

if(run.pyclone)
{
  cat(cmd)
  system(cmd)
}

cat('\n')

# Load the results
sample.results <- paste(new.dir,"/",sample,'.results.tsv',sep="")
cmd <- paste(PyClone
             ," build_table "
             , pyclone.config.yaml
             , " "
             , sample.results
             , " --burnin 1000"
             , sep="")
cat('\n')

if(run.pyclone)
{
  cat(cmd)
  system(cmd)
  
}

cat('\n')


pyclone.results <- read.table(sample.results
                              ,sep="\t"
                              ,header=TRUE
                              ,stringsAsFactors=FALSE
)

# let's make sure the same mutations are being used. 
rownames(pyclone.results)  <- pyclone.results$X
pyclone.results            <- pyclone.results[rownames(pyclone.results)%in%mut.table$mutation_id,,drop=FALSE]

most.likely.cluster        <- pyclone.results$cluster_id
names(most.likely.cluster) <- pyclone.results$X

no.optima = length(unique(most.likely.cluster))
max.cols = 12

cols           = paste(brewer.pal(min(max.cols,no.optima),name="Paired"),sep="")
cols           = rep(cols,ceiling(no.optima/max.cols))[1:no.optima]
cols.opac      = paste(cols,'99',sep="")


# let's plot the copy number clusters
pdf(paste(new.dir,"/",patient,"_pyclone_cluster_assignment_copynumber",".pdf",sep=""),width=6,height=3)
plot.region.mutCopyNum(phylo.region.list = phylo.region.list,seg.mat.copy = seg.mat.copy,mostLikelyClusters = most.likely.cluster,plot.separate.clusters = TRUE)
dev.off()

# first, let's import the trace files
region.trace              <- list()
region.preClustPosterior  <- list()
region.postClustPosterior <- list()

pdf(paste(new.dir,sample,'.GlobalClusteringResults.pdf',sep=""))
for (region in regions.to.use)
{
  trace.region <- read.table(paste(new.dir,"trace/",region, ".cellular_frequencies.tsv.bz2",sep="")
                             ,stringsAsFactors=FALSE
                             ,sep="\t"
                             ,header=TRUE)
  trace.region             <- trace.region[(burn.in+1):nrow(trace.region),]
  colnames(trace.region)   <- gsub("\\.",":",(colnames(trace.region)))
  region.trace[[region]]   <- trace.region
  phylo.region.list[[region]]$phyloCCF_PyClone          <- NA
  phylo.region.list[[region]]$phyloCCF_PyClone.0.05     <- NA
  phylo.region.list[[region]]$phyloCCF_PyClone.0.95     <- NA
  phylo.region.list[[region]]$phyloCCF_PyClone.cluster  <- NA
  
  mutation_ids    <- unlist(phylo.region.list[[region]]$mutation_id)
  tmp             <- intersect(mutation_ids,colnames(trace.region))
  pyclone.est     <- pyclone.results[tmp,region]
  quants          <- apply(trace.region,2,quantile,prob=c(0.05,0.95))
  
  phylo.region.list[[region]]$phyloCCF_PyClone.0.05[mutation_ids%in%tmp]      <- quants[1,tmp]
  phylo.region.list[[region]]$phyloCCF_PyClone.0.95[mutation_ids%in%tmp]      <- quants[2,tmp]
  phylo.region.list[[region]]$phyloCCF_PyClone[mutation_ids%in%tmp]           <- pyclone.est
  phylo.region.list[[region]]$phyloCCF_PyClone.cluster[mutation_ids%in%tmp]   <- most.likely.cluster[tmp]
  
  # let's plot every mutation, all together and separately
  # load the pyclone.tsv file
  pyclone.tsv <- read.table(paste(new.dir,region, '.tsv',sep="")
                            ,stringsAsFactors=FALSE
                            ,header=TRUE)
  
  #pyclone.tsv$var_counts/(pyclone.tsv$ref_counts+pyclone.tsv$var_counts)
  n.alt <- pyclone.tsv$var_counts
  depth <- (pyclone.tsv$ref_counts+pyclone.tsv$var_counts)
  get.posterior <- function(n.alterntive,depth.all)
  {
    x <- dbinom(n.alterntive,depth.all, prob=seq(0.0,0.5,length.out=100))
    return(x/sum(x))
  }
  posterior.preCluster <- c()
  for (i in 1:length(n.alt))
  {
    posterior.preCluster <- cbind(posterior.preCluster,get.posterior(n.alt[i],depth[i]))
  }
  
  colnames(posterior.preCluster) <- pyclone.tsv$mutation_id
  region.preClustPosterior[[region]] <- posterior.preCluster
  
  par(mfrow=c(2,1))
  par(mar=c(0.5,5,5,3))
  plot(0,xlim=c(-0.01,1.01),ylim=c(0,1),xaxt='n',yaxt='n',type='n',xlab='PhyloCCF',ylab='density')
  text(x = 0.4,y=0.8,labels = paste('PreClustering',region))
  #axis(side = 1,at = seq(0,1,by=0.1),labels = seq(0,1,by=0.1))
  axis(side = 2,at = seq(0,1,by=0.1),labels = seq(0,1,by=0.1),las=2)
  for (i in 1:ncol(posterior.preCluster))
  {
    lines(y=c(0,posterior.preCluster[,i],0)
          ,x=c(-0.00005,seq(0,1,length.out=100),1.00005)
          ,col='#2b8cbe99')
    polygon(y=c(0,posterior.preCluster[,i],0),border = FALSE
            ,x=c(-0.00005,seq(0,1,length.out=100),1.00005)
            ,col='#2b8cbe85'
    )
  }
  
  
  posterior.postCluster  <- c()
  par(mar=c(5,5,0.5,3))
  
  plot(0,xlim=c(-0.01,1.01),ylim=c(0,1),xaxt='n',yaxt='n',type='n',xlab='PhyloCCF',ylab='density')
  text(x = 0.4,y=0.8,labels = paste('PostClustering',region))
  axis(side = 1,at = seq(0,1,by=0.1),labels = seq(0,1,by=0.1))
  axis(side = 2,at = seq(0,1,by=0.1),labels = seq(0,1,by=0.1),las=2)
  for (i in 1:ncol(trace.region))
  {
    trace.i <- trace.region[,i]
    y       <- density(trace.i,from=0,to=1,n=100)
    y$y     <- y$y/sum(y$y)
    posterior.postCluster <- cbind(posterior.postCluster,y$y)
    lines(y=c(0,y$y,0)
          ,x=c(-0.00005,seq(0,1,length.out=length(y$y)),1.00005)
          ,col='#de2d2699')
    polygon(y=c(0,y$y,0),border = FALSE
            ,x=c(-0.00005,seq(0,1,length.out=length(y$y)),1.00005)
            ,col='#de2d2685'
    )
    
  }
  colnames(posterior.postCluster)    <- colnames(trace.region)
  region.postClustPosterior[[region]] <- posterior.postCluster
  
  
}
dev.off()

pdf(paste(new.dir,"/",patient,"_pyclone_clustersPatient",".pdf",sep=""),
    ,useDingbats=FALSE
    ,width=4.5,height=3.85
)
par(mgp=c(2,0.5,0))
determinePhylogeny(regionList = phylo.region.list
                   ,mostLikelyClusters = most.likely.cluster
                   ,mutTable = mut.table
                   ,driverCat = c('1','2','3')
                   ,ccf.type = 'phylo.ccf.PyClone'
)

dev.off()

# Finally, let's add the results back to the mutation table
mut.table$PreClusterCCF   <- NA
mut.table$PostClusterCCF  <- NA
mut.table$PyCloneCluster  <- NA
mut.table$MutCopyNum      <- NA

PreClusterCCF <- rep(NA,nrow(mut.table))
names(PreClusterCCF) <- mut.table$mutation_id
PostClusterCCF       <- PreClusterCCF
PyCloneCluster       <- PreClusterCCF
MutCopyNum           <- PreClusterCCF

PrephyloCCF         <- phylo.region.list[[1]]$phyloCCF
names(PrephyloCCF)  <- phylo.region.list[[1]]$mutation_id
PostphyloCCF        <- phylo.region.list[[1]]$phyloCCF_PyClone
names(PostphyloCCF) <- phylo.region.list[[1]]$mutation_id
mutcopyNum          <- phylo.region.list[[1]]$mutCopyNum
names(mutcopyNum)   <- phylo.region.list[[1]]$mutation_id

PreClusterCCF[names(PrephyloCCF)] <- PrephyloCCF
PostClusterCCF[names(PostphyloCCF)] <- PostphyloCCF
PyCloneCluster[names(most.likely.cluster)] <- most.likely.cluster
MutCopyNum[names(mutcopyNum)] <- mutcopyNum


mut.table$PreClusterCCF  <- PreClusterCCF
mut.table$PostClusterCCF <- PostClusterCCF
mut.table$PyCloneCluster <- PyCloneCluster
mut.table$MutCopyNum     <- MutCopyNum

# let's write the new mutation table to the new directory

write.table(mut.table, file=paste(new.dir,patient,".clonalDissection.txt",sep=""),quote=FALSE,sep="\t",col.names=NA)
cat(paste('\n', 'output File written to:',paste(new.dir,patient,".clonalDissection.txt\n",sep=""),sep=""))

# Finally, let's put this into a megatable, and write this to an appropriate place
# save(phylo.region.list,file=paste(new.dir,sample,'.PhyloRegionList.RData',sep=""))
sessInfo <- sessionInfo()
save.image(file=paste(new.dir, sample, ".PyClone.RData",sep=""))


