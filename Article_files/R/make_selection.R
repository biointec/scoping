#' @haploid the haploid of the breeding pool
#' @genome the genome of the breeding pool
#' @N number of individuals which need to be selected 
# Selects parents to maximize the frequency of the favorable QTL alleles

make.selection<-function(haploids=parent.haploids,genome=hv.genome, N){
  
#genotype parents
  Genos.i <- genotype.loci(haploid.genos = haploids, 
                                     genome = genome, 
                                     include.QTL = T)
  
  
  n.chr<-length(genome)
  loci.per.chr <- sapply(X = genome, FUN = function(chr) chr@pos.snp %>% length())
  # Split the haploid.mat into chromosomes
  length_genome<-0
  qtl.index<-c()
  add.effect<-c()
  #Calculate the genotype of the alleles
  for (i in  seq(n.chr)){
    
    #QTL ID on chromosome i
    qtl.index<-c(qtl.index, genome[[i]]@pos.add.qtl$ID+length_genome)
    length_genome<-length_genome+genome[[i]]@num.snp.chr+genome[[i]]@num.add.qtl.chr
    #additive effect of each QTL
    add.effect<-c(add.effect, genome[[i]]@add.and.dom.eff$add>0)
  }
  

  qtl.candidates<-Genos.i[,qtl.index]
  #find the best candidate to start
  ideal.qtl<-matrix(-1,1,100)
  ideal.qtl[add.effect]=1
  score.matrix=matrix(0,dim(Genos.i)[1])
  
  for(j in seq(dim(Genos.i)[1])){
  score.matrix[j]<-sum(qtl.candidates[j,]==ideal.qtl)
  }
  
  selectie<-which.max(score.matrix)
  score<-max(score.matrix)
 
  sel.qtl<-matrix(F,1,100)
  sel.qtl[qtl.candidates[selectie,]==ideal.qtl]=T
  
  
  for (sel in seq(N-1)){
    vrije<-seq(dim(Genos.i)[1])[-selectie]
    
   temp.choose<-0
    for(candidate in vrije){
      temp<-sel.qtl
      temp[qtl.candidates[candidate,]==ideal.qtl]=T
      if (sum(temp)>score){
        score<-sum(temp)
        temp.choose<-candidate
      }
    }
   if(temp.choose!=0){
     selectie<-c(selectie, temp.choose)
     sel.qtl[qtl.candidates[temp.choose,]==ideal.qtl]=T
     
   }else{
     temp.choose<-which.max(score.matrix[vrije])
     temp.choose<-vrije[temp.choose]
     selectie<-c(selectie, temp.choose)
   }
  }
  return(rownames(Genos.i)[selectie])
  
}
