#' @genome the genome of the breeding population
#' @candidate.af.i the frequency of each marker 
#' @candidate.haploid.i the haploid of the individuals in the breeding pool

# calculate the maximum reachable genetic value, maximum genetic value, fixed value, max fixed value of the breeding pool
# used in the baseline, back-crossing, scoping and back-crossing with scoping methods
QTL.max<-function(genome=hv.genome, candidate.af.i=candidate.af.i, haploid = candidate.haploid.i){

  #pre allocation
  geno_fixed<-c()
  geno_total<-c()
  geno_max<-c()
  max_notFixed<-c()
  n_fixed<-c()
  max_fixed<-c()
  names<-c()
  mean_notFixed_value<-c()
  haploid.mat<-do.call("rbind",haploid)

  n.chr<-length(genome)
  loci.per.chr <- sapply(X = genome, FUN = function(chr) chr@pos.snp %>% length())
  # Split the haploid.mat into chromosomes
  haploid.genos.split <- lapply(X = split(x = seq(ncol(haploid.mat)), rep(seq(n.chr), loci.per.chr)), 
                                FUN = function(chr) haploid.mat[,chr])
  
 
  #Calculate the genotype of the alleles
  for (i in  seq(n.chr)){
    
    #QTL ID on chromosome i
    qtl.index<-genome[[i]]@pos.add.qtl$ID
    
    #additive effect of each QTL
    add.eff<-genome[[i]]@add.and.dom.eff$add
    
    #haploids of chromosome 1
    chr.haploid.genos<-haploid.genos.split[[i]]
    
    #haploids of QTLs  
    qtl.alleles <-   chr.haploid.genos[,qtl.index]
    f=colMeans(qtl.alleles)  
    
    names<-cbind(names, t(colnames(qtl.alleles)))
    # Split the gametes into pairs
    haploid.line.split <- split(seq(nrow(qtl.alleles)), rep(seq(nrow(qtl.alleles)/2), each = 2))
                                 
    
    # lapply to calculate genotypic value
    geno <- sapply(X = haploid.line.split, FUN = function(pair) {
    # Pull out the allele states for an individual
    line.gametes <- qtl.alleles[pair,]
    # Sum and subtract 1 to get the multiplier for a
    line.genos <- t(as.matrix(colSums(line.gametes) - 1))
                                   
    # Mulitple by the qtl allele effect to get the value at each QTL
    # Then add the domiance effect for hets
    return(line.genos)
  
    })
  
  fixed<-(f==1 | f==0)
 
  qtl.geno<-t(as.matrix(geno))
  
  geno_fixed<-cbind(geno_fixed, sum(qtl.geno[1,fixed]*add.eff[fixed]))
  max_fixed<-cbind(max_fixed, sum(abs(add.eff[fixed])))
  geno_total<-cbind(geno_total, qtl.geno%*%add.eff)
  geno_max<-cbind(geno_max, sum(abs(add.eff)))
  max_notFixed<-cbind(max_notFixed, sum(abs(add.eff[!fixed])))
  n_fixed<-cbind(n_fixed, sum(fixed))
  #mean_notFixed_value<-cbind( mean_notFixed_value, mean(qtl.geno[!fixed]%*%add.eff[!fixed]))
  

}


output<-c(maximum = as.numeric(sum(geno_max)), n.fixed = as.numeric(sum(n_fixed)), fixed.value = as.numeric(sum(geno_fixed)), maximum.notFixed = as.numeric(sum(max_notFixed)), maximum.fixed=as.numeric(sum(max_fixed)), max_reachable=sum(max_notFixed)+sum(geno_fixed), mean_NotFixed_value=mean_notFixed_value)
return(output)
}
 
