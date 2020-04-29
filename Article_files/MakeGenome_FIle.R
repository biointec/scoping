#import library
library(GSSimTPUpdate) 
library(hypred)

#import data
data("CAP.haploids")
data("CAP.markers")

if (!file.exists("Genome")){
  dir.create(Genome)
}

#choose parameters for QTL and heritability
n.QTL <- 100 
h2 <- 0.5
n.iterations<-250 #number of iterations 

n.chr.snps <- tapply(X = CAP.markers$rs, INDEX = CAP.markers$chrom, length)
chr.len <- as.numeric(tapply(X = CAP.markers$pos, INDEX = CAP.markers$chrom, FUN = max))
genetic.map.list <- tapply(X = CAP.markers$pos, INDEX = CAP.markers$chrom, FUN = function(chr) list(chr))

#construct basic genome
hv.genome <- make.genome( n.chr = length(chr.len), 
                          chr.len = chr.len, 
                          n.chr.snps = n.chr.snps,
                          genetic.map = genetic.map.list)

#pre alloceer 
Genome<-list(data=c(),genome=c())
#save genomic simulation parameter
Genome$data<-list(n.iterations=n.iterations,
                  heritability=h2, 
                  number_QTL=n.QTL, 
                  n.chr.snps=n.chr.snps,
                  chr.len=chr.len, 
                  genetic.map.list=genetic.map.list)

for(i in seq(n.iterations)){
  #create genome architecture
  hv.genome <- trait.architecture(genome = hv.genome,
                                  n.QTL = n.QTL, 
                                  qtl.index = NULL, 
                                  qtl.dom.index = NULL, 
                                  qtl.perf.index = NULL, 
                                  qtl.add.eff = "geometric", 
                                  qtl.dom.eff = NULL)
  
  #save genome in list
  Genome$genome<-c(Genome$genome, list(hv.genome))
}

filename<-paste0("Genome/CAP_Genome_QTL", n.QTL,"_h",h2,"_",n.iterations,"_Iterations.RData")
save("Genome", file = filename)
