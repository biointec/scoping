#' @Description
#' @Parents contains a list of names wich we can select the next ofspring with.
#' @Haploids contains the haploid data of all those Parents
#' @N is the number of crosses needed

#' @output: Selection is a list with N arrays containing 2 Parents wich will be crossed with each other

#' @Goal: making crossing blocks that maximize the genetic variability

make.crossing.block.scoping <- function(Parents= parent.selections.i,
                                        Haploids = parent.haploids.i, 
                                        N = n.crosses, 
                                        genome = hv.genome) {
  
  library(resample)
  #-------------------------------------------------------------------------------------------------------------------------------------------
  #SCOPING METHODOLOGY------------------------------------------------------------------------------------------------------------------------
  #-------------------------------------------------------------------------------------------------------------------------------------------
  Parent.genos <- genotype.loci(haploid.genos =  Haploids, 
                                     genome = genome, 
                                     include.QTL = T)
  #If T: QTL alleles will be included in the genotype, simulating the case when markers also covers the QTL 
  #if F: QTL alleles will not be included in the genotype, simulating the case when the the markers are only in LD with QTL
  
  order_id<-order(Parents$value.sel,decreasing = T)
  
  marker_allele=matrix(0,2,dim(Haploids)[2])
  marker_seq=matrix(TRUE,1,dim(Haploids)[2])
  
  Genotype<-Parent.genos[order_id,]
  Parents_available<-Parents$lines.sel[order_id]
  free.position<-seq(length(Parents_available))
  all_parents<-c()
  
  #select first parent------------------------------------------------------------
  name1<-Parents_available[free.position[1]] 
  gen1<-Genotype[free.position[1],]
  chosen<-free.position[1]
  all_parents<-rbind(all_parents, Parent.genos[name1,])
  p1<- name1
  
  score<-0
  sel<-0
    for (i in free.position[-chosen]){
    temp.name<- Parents_available[i]
    temp.geno<- as.matrix(Genotype[temp.name,])
    temp.score<-sum(colVars(rbind(all_parents, t(temp.geno)))[marker_seq])
    if (temp.score == 0 ){ #Indien de huidige combo geen nieuwe info inbrengen of alle info is al opgenomen, neem dan de beste kandidaten
      temp.score<-sum(colVars(rbind(all_parents, t(temp.geno))))
    }
    
    if (temp.score>score){
      score<-temp.score
      sel<-i
    }
  }  
  chosen<-c(chosen, sel)
  p2<-c(Parents_available[sel])
  gen2<-Genotype[sel,]
  
  all_parents<-rbind(all_parents, Parent.genos[p2,])
  
  marker_allele[1,(gen1+gen2)==-2]<-marker_allele[1,(gen1+gen2)==-2]+1
  marker_allele[2,(gen1+gen2)== 2]<-marker_allele[2,(gen1+gen2)==2]+1
  #nul wordt niet meegeteld want er bestaat slechts een kleine kans dat beide allelen worden ingebouwd
  marker_allele[,abs(gen1+gen2)== 1]<-marker_allele[,abs(gen1+gen2)== 1]+1
  marker_seq[marker_allele[1,]>0 & marker_allele[2,]>0]<- FALSE
  
  
  for (sel.round in seq(2,N)){

    name1<-Parents_available[free.position[-chosen][1]]
    gen1<-Genotype[free.position[-chosen][1],]
    chosen<-c(chosen, free.position[-chosen][1])
    all_parents<-rbind(all_parents, Parent.genos[name1,])
    
    
    p1<-c(p1, name1)
    geno1<- as.matrix(Genotype[name1,])
    score<-0
    sel<-0
    for (i in free.position[-chosen]){
      temp.name<-Parents_available[i]
      temp.geno<- as.matrix(Genotype[temp.name,])
      
      temp.score<-sum(colVars(rbind(all_parents, t(temp.geno)))[marker_seq])
      if (temp.score == 0 ){ #if no new marker alleles can be found, select the indidivuals greedily
        temp.score<-sum(colVars(rbind(all_parents, t(temp.geno))))
      }
      
      if (temp.score>score){
            score<-temp.score
            sel<-i
       }
     }
            
      chosen<-c(chosen, sel)
      p2<-c(p2, Parents_available[sel])
      gen2<-Genotype[sel,]
      
      marker_allele[1,(gen1+gen2)==-2]<-marker_allele[1,(gen1+gen2)==-2]+1
      marker_allele[2,(gen1+gen2)== 2]<-marker_allele[2,(gen1+gen2)==2]+1
      
      marker_allele[,abs(gen1+gen2)== 1]<-marker_allele[,abs(gen1+gen2)== 1]+1
            
      marker_seq[marker_allele[1,]>0 & marker_allele[2,]>0]<- FALSE
            
      all_parents<-rbind(all_parents, Parent.genos[Parents_available[sel],])
      score<-temp.score
      sel<-i
      
  }

output<-data.frame(Parent1=p1,Parent2=p2)
colnames(output)<-c(c("Parent1", "Parent2"))
return(output)

}
