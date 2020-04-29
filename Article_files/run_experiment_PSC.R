## Genomic selection simulations
#-------------------------------------------------------------------------
#Copyright (c) 2016 Jeff Neyhart
#-------------------------------------------------------------------------

load("Genome/CAP_Genome_QTL100_h0.5_250_Iterations.RData")#------------------------------Set Path To Genome
#Penalty weight of PSC method (often negative)
D_w<--0.15
#source hypred map, imported library as tar.gz file
library(hypred)
library(rrBLUP)
library(GSSimTPUpdate)
library(rrBLUP)
library(parallel)
library(stringr)

save.dir<-"own_results"
# First and only argument is the pop makeup (MN, ND, or MNxND)
# Second argument is how the TP should be combined after each cycle (cumulative or window)
# Third argument is the number of QTL and the heritability
pop.makeup <- "MNxND" #wijst naar de datasets die gebruikt worden, MN en ND
tp.formation <- "window"
h2 <- Genome$data$heritability #heritability
tp.change <- "tails"


# Verify that the TP changes in the arguments are acceptable
if (!all(tp.change %in% c("best", "worst", "random", "nochange", "PEVmean", "CDmean", "tails")) )
  stop("The TP change arguments are not acceptable.")

# Load the datasets
data("CAP.haploids")
data("CAP.markers")

# Set the number of cores by detection
max.cores <- detectCores()
n.cores <- 1#----------------------------------------------------------------------------NumberOfCores

#if number of selected cores exceed the available number of cores
if (n.cores>max.cores){
  n.cores<-max.cores
}

# Other simulation parameters
n.QTL <- Genome$data$number_QTL
# Make sure heritability is numeric
h2 <- as.numeric(h2)

# How many cycles?
n.cycles =10 #----------------------------------------------------------------------------------------Cycles

# Number of phenotyping environments and reps
n.env = 3
n.rep = 1

# Minor allele frequency cut-off for markers
min.maf = 0.03

# Barley population genetics data
mutation.rate.snp = 0
mutation.rate.qtl = 0

# Selection intensity and the number of crosses
parents.sel.intensity = 100
n.crosses = 50

# The number of lines to add the TP after each cycle
tp.update.increment = 150
# Size of the TP to maintain - this is the same as the starting TP
tp.size <- nrow(CAP.haploids) / 2

# Parent selection and crossing parameters
ind.per.cross = 20
cycle.candidate.size = n.crosses * ind.per.cross

# Standardized selection intensity
std.sel.intensity = parents.sel.intensity / cycle.candidate.size

# Computation parameters
n.iterations = Genome$data$n.iterations

date <- format(Sys.time(), "%d%m%y-%H%M%S")

# Save the metadata to a list
metadata <- list(h2 = h2,
                 n.cycles = n.cycles,
                 n.QTL = n.QTL,
                 min.marker.maf = min.maf,
                 parents.sel.intensity = parents.sel.intensity,
                 n.env = n.env, 
                 n.rep = n.rep,
                 mutation.rate.qtl = mutation.rate.qtl,
                 mutation.rate.snp = mutation.rate.snp,
                 tp.update.increment = tp.update.increment,
                 tp.size = tp.size,
                 n.crosses = n.crosses,
                 ind.per.cross = ind.per.cross,
                 cycle.candidate.size = cycle.candidate.size,
                 std.sel.intensity = std.sel.intensity,
                 n.iterations = n.iterations,
                 pop.makeup = pop.makeup,
                 date = date)


#### Define genome characteristics ####
# Find the snps per chromsome   rs contains the name of markers, index the number of chromosome (factor), function  length
# tAPPLY will split cap.markers$rs into the different factors and give the length of each part
n.chr.snps <- Genome$data$n.chr.snps

# Find the chromsome lengths
#the cromosome length is calculated as the max distance we find a marker
chr.len <- Genome$data$chr.len

# Create a list of loci positions
genetic.map.list <- Genome$data$genetic.map.list


for (change in tp.change) {
  
  # Split iterations into cores
  if (n.cores > 1) {
    iters.per.core <- split(x = 1:n.iterations, factor(cut(x = 1:n.iterations, breaks = n.cores)))
    names(iters.per.core) <- paste("set", seq_along(iters.per.core), sep = "")
  } else {
    iters.per.core <- 1:n.iterations
  }    
  
  # Apply the iterations over cores
  experiment.sub.results <- mclapply(X = iters.per.core, FUN = function(iter.set) {
    
    # Create a vector of rep names
    reps <- str_c("rep", iter.set)
    
    # Iterate over reps
    sapply(X = reps, FUN = function(r) {
      
      rep.iter<-as.numeric(substring(r,4))
 
      # All code below this line is different in each iteration of the simulation

      #### Define trait parameters ####
      hv.genome <- Genome$genome[[rep.iter]]
      
      TP.haploids.i <- CAP.haploids
      # Convert the gametes to genotypes
      TP.genos <- genotype.loci(haploid.genos = TP.haploids.i, genome = hv.genome)
      
      # Find the allele freq of all marker snps
      marker.af <- measure.af(genome = hv.genome, haploid.genos = TP.haploids.i)$snp
      
      loci.af<-measure.af(genome = hv.genome, haploid.genos = TP.haploids.i)$loci
      Fixed<-sum( loci.af==0 | loci.af == 1)
      # Calculate maf minor allele freq
      marker.maf <- sapply(marker.af, FUN = function(freq) min(freq, 1 - freq))
      
      # Determine which are below 
      markers.below.maf <- which(marker.maf < min.maf)
      
      
      # Calculate the frequency of the 1 allele in the base training population
      p_i <- apply(X = TP.genos + 1, MARGIN = 2, FUN = mean) / 2
      # Calculate the P matrix
      P = matrix(2 * (p_i - 0.5))
      # Calculate the normalization constant
      c = 2 * sum(p_i * (1 - p_i))
      
      ## Select the TP lines for use as parents
      # First separate MN and ND lines
      line.names <- row.names(TP.genos)
      ND.lines <- str_subset(string = line.names, pattern = "^ND")
      MN.lines <- setdiff(x = line.names, y = ND.lines)
      
      ## Set the inital variances for the heritability
      # True genetic variance
      TP.V_g <- genotypic.value(genome = hv.genome, haploid.genos = TP.haploids.i) %>%
        var()
      G_G0 <- genotypic.value(genome= hv.genome, haploid.genos = TP.haploids.i) %>% 
        mean()
      G_G0_top <- sort(genotypic.value(genome= hv.genome, haploid.genos = TP.haploids.i),decreasing=TRUE)[1:10]%>%
        mean()
      
      # Environmental variance (scale * 8 as in Bernardo 2015)
      V_E = TP.V_g * 8
      # Residual variance scaled to achieve desired h2
      V_e = n.rep * n.env * ((TP.V_g / h2) - TP.V_g)
      
      
      # Phenotype the training population
      
      TP.values <- phenotype.population(genome = hv.genome,
                                        haploid.genos = TP.haploids.i,
                                        V_E = V_E,
                                          V_e = V_e,
                                        n.env = n.env,
                                        n.rep = n.rep,
                                        run.anova = T)
      
      
      
      TP.phenos <- TP.values$mean.pheno.values
      
      # Next select the top MN and top ND
      top.MN.lines <- names(sort(TP.phenos[MN.lines,], decreasing = T)[1:parents.sel.intensity])
      top.ND.lines <- names(sort(TP.phenos[ND.lines,], decreasing = T)[1:parents.sel.intensity])
      
      # Create a list to store the line names
      pop.makeup.list <- list(
        MN = list(p1 = top.MN.lines, p2 = top.MN.lines),
        ND = list(p1 = top.ND.lines, p2 = top.ND.lines),
        MNxND = list(p1 = top.MN.lines[seq(parents.sel.intensity / 2)], p2 = top.ND.lines[seq(parents.sel.intensity / 2)])
      )
      
      # Set the parent gamete data input
      parent.lines.list <- pop.makeup.list[[pop.makeup]]
      parent.haploids <- select.haploids(haploid.genos = TP.haploids.i, line.names = unique(unlist(parent.lines.list)))
      
      
      # Set dummy variables for the phenos and genos
      TP.phenos.i <- TP.phenos
      TP.genos.i <- TP.genos
      
      # Create an initial data list
      simulation.results <- list()
      
      # Loop over the number of cycles
      for (breeding.cycle in seq(n.cycles)) {

        ##### Start the Cycle Executions #####

        ##### Step 1 - Crossing and inbreeding
        # Make a crossing block
        crossing.block.i <- make.crossing.block(parent1.lines = parent.lines.list$p1, 
                                                parent2.lines = parent.lines.list$p2, 
                                                n.crosses = n.crosses, 
                                                use.parents.once = T)

  
        # According to the crossing block, parents are crossed to form F1s, then
        ## the progeny are inbred to the F3 generation. Since each F1 plant is inbred
        ## individually, the resulting families consist of F1:3 lines
        library(dplyr)
        library(stringr)
        library(Rcpp)

        candidate.haploid.i <- make.population(genome = hv.genome, 
                                               parental.haploids = parent.haploids,
                                               crossing.block = crossing.block.i,
                                               N = ind.per.cross,
                                               cycle.number = breeding.cycle,
                                               generations = 2,
                                               pop.type = "inbred",
                                               mutation.rate.snp = mutation.rate.snp,
                                               mutation.rate.qtl = mutation.rate.qtl)
        

        
        #--------------------------------------------------------------------------------------------------------------
        ##### Step 2 - Genotype
        # Find the genotypes of the markers and QTL
      
        candidate.genos.i <- genotype.loci(haploid.genos = candidate.haploid.i, 
                                           genome = hv.genome, 
                                           include.QTL = T)
        # Just the marker genotypes
        candidate.marker.genos.i <- genotype.loci(haploid.genos = candidate.haploid.i, 
                                                  genome = hv.genome, 
                                                  include.QTL = F)
        
        
        ##### Step 3 - Genotypic Summary Statistics
        # Measure the frequency of the 1 allele in each of the haploids set
        candidate.af.i <- measure.af(genome = hv.genome, 
                                     haploid.genos = candidate.haploid.i)
        
        TP.af.i <- measure.af(genome = hv.genome, 
                              haploid.genos = TP.haploids.i)
        
        source("R/QTL_maxv2.R")
        qtl.out<-QTL.max(genome=hv.genome, candidate.af.i=candidate.af.i, haploid = candidate.haploid.i)
        if (breeding.cycle==1){
        source("R/QTL_maxv2_base.R")
        qtl.base<-QTL.max.base(genome=hv.genome, candidate.af.i = TP.af.i,haploid= TP.haploids.i)
        }
        ### LD measures
        # Candidates
        
        candidate.LD.window <- measure.LD(genome = hv.genome, 
                                          genos = candidate.genos.i, 
                                          Morgan.window = 0.25)
        
        # No window
        candidate.LD.genome <- measure.LD(genome = hv.genome,
                                          genos = candidate.genos.i)
        
        # For the whole genome, find the mean LD value across those
        ## max LD values per QTL
        
        #if the candidate.LD.genome is NA, next calculations cannot be calculated.
        if (!any(is.na(candidate.LD.genome))){
          
        candidate.mean.max.LD.genome <- apply(X = candidate.LD.genome, MARGIN = 1, FUN = function(qtl)
          max(qtl^2) ) %>%
          mean(na.rm = T)

        
        # Measure genomic LD on the TP
        TP.LD.genome <- measure.LD(genome = hv.genome, 
                                   genos = genotype.loci(haploid.genos = TP.haploids.i,
                                                         genome = hv.genome,
                                                         include.QTL = T) )
        
        # Measure the mean LD value across those
        ## max LD values per QTL in the TP
        TP.mean.max.LD.genome <- apply(X = TP.LD.genome, MARGIN = 1, FUN = function(qtl)
          max(qtl^2) ) %>%
          mean(na.rm = T)
        
        ## Persistance of LD phase
        # First find the common polymorphic QTL
        common.poly.QTL <- intersect( row.names(TP.LD.genome), row.names(candidate.LD.genome) )
        common.poly.markers <- intersect( colnames(TP.LD.genome), colnames(candidate.LD.genome) )
        
        # Subset the TP and candidates for those markers and QTL, then create a
        # data.frame
        TP.candidate.LD <- data.frame(TP = TP.LD.genome[common.poly.QTL, common.poly.markers] %>%
                                        as.vector(),
                                      candidates = candidate.LD.genome[common.poly.QTL, common.poly.markers] %>%
                                        as.vector())
        
        
        # Correlate
        TP.candidate.persistance.of.phase <- cor(TP.candidate.LD) %>% 
          .[upper.tri(.)]
        
        # Determine the number of loci used to calculate LD
        n.marker.LD <- min(ncol(TP.LD.genome), ncol(candidate.LD.genome))
        n.qtl.LD <- min(nrow(TP.LD.genome), nrow(candidate.LD.genome))
        
        
        # Create a list to save
        qtl.marker.LD.i <- list(sc.mean.max.genome = candidate.mean.max.LD.genome,
                                tp.mean.max.genome = TP.mean.max.LD.genome,
                                persistance.of.phase = TP.candidate.persistance.of.phase,
                                n.qtl.LD = n.qtl.LD,
                                n.marker.LD = n.marker.LD)
        } else {
          qtl.marker.LD.i<-NA #indicates that the qtl.marker.LD.i cannot be calculated.
        }
        ### Measure the average relationship between the TP and the candidates
        # Assign M
        M <- rbind(TP.genos.i, candidate.marker.genos.i)
        # Subtract P to make Z (need to convert P into a repeated matrix)
        W = M - matrix(P, nrow(M), length(P), byrow = T)
        # Calculate the relationship matrix
        A = tcrossprod(W) / c
        
        # Subset the relationship matrix for the selection candidates
        A.sc <- A[row.names(candidate.marker.genos.i), row.names(candidate.marker.genos.i)]
        
        # Calculate the mean relationship between the TP and the candidates
        mu.relationship <- A[row.names(TP.genos.i), row.names(candidate.marker.genos.i)] %>%
          mean()
        
        ## Find the average inbreeding coefficient among the selection candidates
        candidate.inbreeding <- A.sc %>%
          diag() %>%
          - 1 %>%
          mean()
        
        
        ##### Step 4 - Prediction
        # Remove the markers with maf below the threshold (set at the start of the sim)
        ## also remove monomorphic markers

        
        markers.to.remove <- c( which(candidate.af.i$snp == 0), 
                                which(TP.af.i$snp == 0),
                                markers.below.maf ) %>%
          unique() %>%
          sort()

        
        # Filter the TP and candidate marker matrices for those markers
        TP.genos.use <- TP.genos.i[,-markers.to.remove]
        candidate.genos.use <- candidate.marker.genos.i[,-markers.to.remove]
        
        # Estimate marker effects
        predictions.out <- try(make.predictions(pheno.train = TP.phenos.i,
                                            geno.train = TP.genos.use,
                                            geno.pred = candidate.genos.use))
                                            
        # If it errors out, just return NA for this iteration
        if (class(predictions.out) == "try-error")
          return(NA)
        
        ##### Step 5 - Phenotype the population
        # We use the haploid genotypes from the F1:3 generation to measure the
        ## the phenotypes
        
        # Measure the phenotype and true genotypic values of all selection candidates
        candidate.values.i <- phenotype.population( genome = hv.genome,
                                                    haploid.genos = candidate.haploid.i,
                                                    V_E = V_E,
                                                    V_e = V_e,
                                                    n.env = n.env,
                                                    n.rep = n.rep )
        
        # Validate the predictions
        # Find the correlation between the GEBVs and the true genotypic value
        pred.accuracy.i <- cor(candidate.values.i$geno.values,
                               predictions.out$GEBV) %>%
          as.numeric()


        
        ##### Step 6 - Select the parents of the next generation --------------------------------------- 6 Select Parents
       
        #Use the GEBVs to guide the parental selection  
        value.mat <- predictions.out$GEBV
        l.names <- rownames(value.mat) #Line names
        
        sort.id <- order(value.mat, decreasing=T) #Sort from low to high GEBVs
        #First parental selection (top-100 GEBVs)
        select.id <- sort.id[seq(parents.sel.intensity)] 
        
        #Keep track of the selected and available individuals
        id <- seq(length(value.mat))
        available.mask <- matrix(T,length(value.mat),1)
        available.mask[select.id] <- F
        
        #Calculate mean genetic value of the parental population
        g_w <- mean(value.mat[select.id])
        
        #Genotype current population
        M_pop <- candidate.marker.genos.i

        #Find the inbreeding coefficients of the population
        F_pop<-A.sc %>%
          diag() %>%
          -1
        
        ## Calculate the average inbreeding coefficient of the parental population
        F_inbr <- F_pop[select.id] %>%
          mean()
        
        #Calculate the PSC
        B_w <- g_w + D_w*F_inbr
        
        #Initiate B_w_old
        B_w_old<-0
        
        #Maximize B_w
        while (B_w!=B_w_old){
          B_w_old=B_w
          #Iterate over the selected parents
          for (i in seq(parents.sel.intensity)){
            available.id<-id[available.mask]
            #iterate over the available individuals
            for(j in seq(length(available.id))){
              #Replace the i-th parent with the j-th available individual
              select.temp <- select.id[-i]
              select.temp <- c(select.temp, available.id[j])
            
              #Calculate the new mean genetic value of the parental population
              g_w_temp <- mean(value.mat[select.temp])
              
              
              #Calculate the average inbreeding coefficient of the parental population
              F_temp <- F_pop[select.temp]%>%
                mean()
              #Calculate the new PSC 
              B_w_temp <- g_w_temp + D_w*F_temp
              
              
              #If the replacement of the i-th parent with the j-th individual increases the PSC, the parents are replaced 
              if (B_w_temp>B_w){
                #Update the PSC
                B_w <- B_w_temp
                #The j-th individual becomes unavailable, the i-th parent becomes available again
                available.mask[available.id[j]]<-F
                available.mask[select.id[i]]<-T
                select.id[i]<-available.id[j]
              }
            }
          }
        }
        
        #Process the parental selection into a list
        sel.val<-as.matrix(value.mat[select.id])
        rownames(sel.val)<-names(value.mat[select.id,])
        parent.selections.i<-list(value.sel=sel.val ,lines.sel=rownames(sel.val))
                                   
        parent.lines.list <- list(p1 = parent.selections.i$lines.sel, 
                                  p2 = parent.selections.i$lines.sel)
        
        # The parents are selected and crossed at the F3 stage, so subset the haploid genotpyes from the F1:3
        parent.haploids <- select.haploids(haploid.genos = candidate.haploid.i,
                                          line.names = parent.selections.i$lines.sel)
        
        parent.values <- select.values(pheno.values.list = candidate.values.i, 
                                       line.names = parent.selections.i$lines.sel)
        
        ##### Step 7 - Update the TP --------------------------------------------------------------------7 Update TP
        
        # Skip this step if not called
        if (change != "nochange") {
          
          # If the TP change is best, worst, or random, simply subset the population.
          if (change %in% c("best", "worst", "random", "tails")) {
            
            TP.addition.list <- 
              list(TP.addition.lines = select.population(value.mat = predictions.out$GEBV,
                                                         sel.intensity = tp.update.increment,
                                                         selection = change)$lines.sel )
          }

          if (change %in% c("PEVmean", "CDmean")) {
            
            # Analyze using PEVmean or CDmean
            # We want to see what optimized TP is best for the parents, so we will optimize the training set based
            ## on the lines from the whole candidate set, including the parents
            phenotyped.index <- which( row.names(candidate.marker.genos.i) %in%
                                         row.names(A.sc) )
            unphenotyped.index <- which( parent.selections.i$lines.sel %in%
                                           row.names(A.sc) )

            # V_e is estimated from maximum likelihood
            V_e.i <- predictions.out$solve.out$Ve
            # V_a is estimated as the variance among marker effects * the number of markers
            V_a.i <- predictions.out$solve.out$Vu * ncol(TP.genos.use)
            
            # Run the optimization algorithm
            optimized.TP <- optimize.subset(method = change, A = A.sc,
                                            n.training = tp.update.increment,
                                            phenotyped.index = phenotyped.index,
                                            unphenotyped.index = unphenotyped.index,
                                            V_a = V_a.i, V_e = V_e.i, max.iter = 500)
            
            # Using the optimized index of "phenotyped" entries, determine
            # which candidates should be added to the TP
            TP.addition.list <- list(TP.addition.lines = row.names(A.sc)[optimized.TP$OTP],
                                     optimization.info = optimized.TP)
            
          } # Close the tp optimization algorithm if statement
          
          ### Collect information on the new TP lines
          
          # TP additions
          TP.addition.lines <- TP.addition.list$TP.addition.lines
          
          TP.addition.haploids <- select.haploids(haploid.genos = candidate.haploid.i,
                                                  line.names = TP.addition.lines)
          
          # Subset the geno matrix for these lines
          TP.addition.genos <- candidate.marker.genos.i[TP.addition.lines,]

          # Gather genotypic and phenotypic values of the TP additions
          TP.addition.values <- select.values(pheno.values.list = candidate.values.i, 
                                              line.names = TP.addition.lines)
          # Separate the phenotypes
          TP.addition.phenos <- TP.addition.values$mean.pheno.values
          
          # Measure the average relationship among these lines
          TP.addition.mu.relationship <- A[TP.addition.lines, TP.addition.lines] %>%
            .[upper.tri(.)] %>% 
            mean()
          
          # Measure inbreeding among the additions
          TP.addition.inbreeding <- A[TP.addition.lines, TP.addition.lines] %>%
            diag() %>% 
            - 1 %>%
            mean()
          
          # Combine the new data to the TP
          TP.phenos.i <- rbind(TP.phenos.i, TP.addition.phenos)
          TP.genos.i <- rbind(TP.genos.i, TP.addition.genos)
          TP.haploids.i <- rbind(TP.haploids.i, TP.addition.haploids)
          
          ## Measure the expected heterozygosity of the additions, using all
          ## loci including QTL
          TP.addition.list[["Exp.het"]] <- 
            measure.expected.het(genome = hv.genome, 
                                 haploid.genos = TP.addition.haploids)
          
          
          # If the TP formation calls for a sliding window, use only the ~750 most recent training individuals
          if (tp.formation == "window") {
            
            # If the breeding cycle * tp addition size is less than the starting tp size, 
            ## include the most recent 150 additions, then randomly sample from the
            ## remaining TP
            if ((breeding.cycle * tp.update.increment) < tp.size) {
              
              # Find the index of the most recent additions
              tp.recent.index <- tail(1:nrow(TP.genos.i), (breeding.cycle * tp.update.increment))
              # Randomly select among the remaining index to maintain the tp.size
              tp.random.index <- sort( sample(setdiff(1:nrow(TP.genos.i), tp.recent.index), size = (tp.size - length(tp.recent.index)) ) )
              # Combine
              tp.keep.index <- sort(c(tp.recent.index, tp.random.index))
              
            } else { # Otherwise just take the last tp.size individuals added to the TP
              tp.keep.index <- tail(1:nrow(TP.genos.i), tp.size)
              
            }
              # Set the TP.pheno and TP.genos
              TP.phenos.i <- as.matrix(TP.phenos.i[tp.keep.index,])
              TP.genos.i <- as.matrix(TP.genos.i[tp.keep.index,])
              TP.haploids.i <- select.haploids(haploid.genos = TP.haploids.i,
                                               line.names = row.names(TP.genos.i))
          }
          
        } else {
          
          TP.addition.list <- list(TP.addition.lines = NA,
                                   Exp.het = NA)
          
          TP.addition.inbreeding <- NA
          TP.addition.mu.relationship <- NA
          
          
        } # Close the tp.change if statement
                              
        
        print( paste("Cycle", breeding.cycle, "complete.") )
        
        cycle.name <- paste("cycle", breeding.cycle, sep = "")
        
        # Gather data for analysis
        simulation.results[[cycle.name]] <- 
          list(geno.summary.stats = list(candidate.af = candidate.af.i,
                                         TP.af = TP.af.i,
                                         qtl.marker.LD = qtl.marker.LD.i,
                                         mu.TP.candidate.rel = mu.relationship),
               MM.solve = predictions.out$solve.out,
               candidate.values = candidate.values.i,
               QTL_info = qtl.out,
               QTL_Gen0=qtl.base,
               Genotype_Gen0=G_G0,
               Genotype_Gen0_top=G_G0_top,
               candidate.values=candidate.values.i,
               candidate.GEBV=predictions.out$GEBV,
               prediction.accuracy = pred.accuracy.i,
               parents <- parent.lines.list,
               tp.update = TP.addition.list,
               inbreeding = list(candidates = candidate.inbreeding,
                                 TP.additions = TP.addition.inbreeding),
               relationship = list(TP.candidates = mu.relationship,
                                   TP.additions = TP.addition.mu.relationship))
      }          
      
      # Return a list
      list(sim.results = simulation.results, genome = hv.genome)
      
      # Close the sapply over replications
    }, simplify = F)
    
    # End mclapply
  }, mc.cores = n.cores)
  
  # Save the tp.change data
  filename <- file.path(save.dir, paste("simulation_results_Population_Inbreeding", "_D_",D_w ,"_rep_", n.iterations,"_", date, ".RData", sep = "") )
  
  
  save(list = c("experiment.sub.results", "change", "metadata"), file = filename)
  
  
  
} # Close the tp.change for loop


