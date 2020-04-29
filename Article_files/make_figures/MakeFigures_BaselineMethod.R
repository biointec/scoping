library('dplyr')
library(ggplot2)
library(egg)


map<-"data/Baseline"

if (file.exists(paste0(map,"/Results"))){
  unlink(paste0(map,"/Results"), recursive=TRUE)
}

files<-dir(map, full.names = T)
names<-dir(map, full.names = F)

dir.create(paste0(map,"/Results"))

for(expir in seq(length(files))){
  load(files[expir])
  n.set<-length(experiment.sub.results)
  results<-c()
  acc_fixed<-c()

  for(s in seq(n.set)){
    set <- experiment.sub.results[[s]]
    if(!is.null(set)  & (class(set) != "try-error")){
    nrep<-length(set)
    for(r in seq(nrep)){ 
      rep<-set[[r]]
      
      if (!any(is.na(rep))){
      genome<-rep$genome
      sim<-rep$sim.results
      out<-sim$cycle1$QTL_Gen0    
      geneticValue<-sim$cycle1$Genotype_Gen0/sim$cycle1$QTL_Gen0[['maximum']]
      geneticValue.top<-sim$cycle1$Genotype_Gen0_top/sim$cycle1$QTL_Gen0[['maximum']]
      
      
      for(br.cycle in seq(metadata$n.cycles)){
          gen<-sim[[br.cycle]]
          candidate.af.i<-gen$geno.summary.stats$candidate.af
          geneticValue<-cbind(geneticValue, gen$candidate.values$mu.g/gen$QTL_info[['maximum']])
          geneticValue.top<-cbind(geneticValue.top, mean(sort(gen$candidate.values$geno.values,decreasing = TRUE)[1:10])/gen$QTL_info[['maximum']])
          out<-cbind(out, gen$QTL_info)
      }
      

      Data<-data.frame(t(out))
      Data$cycle<-seq(0,metadata$n.cycles)
     
      
      acc_fixed<-rbind(acc_fixed, data.frame(cbind(Data$n.fixed, Data$cycle))) #----------------------------------------------------------------aanpasssing
     
      # nodige data in 1 column zetten met factor for ggplot legend
      cycle<-Data$cycle
      m1<-cbind(cycle, Data$max_reachable/Data$maximum , 1)
      m2<-cbind(cycle, t(geneticValue), 2)
      m3<-cbind(cycle, t(geneticValue.top), 3)
      m<-rbind(m1,m2,m3)
      m<-data.frame(m)
      colnames(m)<-c("cycle","Genetic_value","Legend")
      m$Legend<-factor(m$Legend, levels = c("1","2","3"), labels=c("Maximum Reachable \nGenetic Value","Genetic Value", "Genetic Value Top"))
      results<-rbind(results, m)
      }
    }
  }
  }
  colnames(acc_fixed)<-c("n.fixed","cycle")

  #plot figure
  library('ggplot2')
  library(gridExtra)

  #bereken de gemiddelde per generatie en de std om deze dan te plotten
  sum.result<-data.frame(cycle = seq(metadata$n.cycles))
  
  result.top<-c()
  result.mean<-c()
  result.reach<-c()
  result.fixed<-c()
  result.cycle<-c()
 
 
  
  for (b.c in seq(0,metadata$n.cycles)){
    
    results %>%
      filter(cycle == b.c) %>%
      filter(Legend=='Maximum Reachable \nGenetic Value') %>%
      summarise(mean_Genetic_Value = mean(Genetic_value), std_Genetic_Value = sd(Genetic_value))->temp.cycle
    temp.cycle$cycle<-b.c
    temp.cycle$Parameter<-'Maximum Reachable \nGenetic Value'
    temp.cycle$Dataset<-'Maximum Reachable GV'
    result.reach<-temp.cycle
    
    results %>%
      filter(cycle == b.c) %>%
      filter(Legend=='Genetic Value') %>%
      summarise(mean_Genetic_Value = mean(Genetic_value), std_Genetic_Value = sd(Genetic_value))->temp.cycle
    temp.cycle$cycle<-b.c
    temp.cycle$Parameter<-'Mean Genetic Value'
    temp.cycle$Dataset<-'Mean GV of Population'
    result.mean<-temp.cycle
    
    results %>%
      filter(cycle == b.c) %>%
      filter(Legend=='Genetic Value Top') %>%
      summarise(mean_Genetic_Value = mean(Genetic_value), std_Genetic_Value = sd(Genetic_value))->temp.cycle
    temp.cycle$cycle<-b.c
    temp.cycle$Parameter<-'Mean Genetic Value'
    temp.cycle$Dataset<-'Mean GV of Top 10 Individuals'
    result.top<-temp.cycle
    
    result.cycle<-rbind(result.cycle, temp.cycle,result.reach, result.mean, result.top)
    
    
    acc_fixed %>%
      filter(cycle == b.c) %>%
      summarise(mean_fixed = mean(n.fixed), std_fixed = sd(n.fixed))->temp.cycle
      temp.cycle$cycle<-b.c
    
    result.fixed<-rbind(result.fixed, temp.cycle)
  }

  
  title<-names[expir]
  title<-substr(title, 26, nchar(title)-20)
  title<-paste(title, "basic breeding strategy")
#  f1<-ggplot(data=result.cycle, aes(x=cycle, y=mean_Genetic_Value,colour=Legend)) +
    #geom_errorbar(aes(ymin=mean_Genetic_Value-std_Genetic_Value, ymax=mean_Genetic_Value+std_Genetic_Value)) +
 #   geom_line()  +
#    ggtitle(title) +
 #   ylab('Genetic Value')
  
  
  f1<-ggplot(data=result.cycle, aes(x=cycle, y=mean_Genetic_Value)) +
    geom_line(aes(linetype=Dataset),size=1.2)  +
    scale_linetype_manual(values=c("twodash", "solid","dotted"))+
    ylab('Genetic Value (GV)')+
    xlab('Breeding Cycle')+
    theme_bw()+
    scale_y_continuous(breaks = seq(0.0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
    scale_x_continuous(breaks = seq(0, 50, 10),expand=c(0,0),limits=c(1,50),sec.axis = dup_axis(name = NULL, labels = NULL))+
    theme(legend.title =  element_blank(),text = element_text(size=30) , legend.position= c(0.55,0.16),legend.direction = "vertical", panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme( axis.text = element_text(colour="black"))+
    theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
    theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"),legend.text = element_text(size=25))+
    theme(axis.text.x =element_text(margin=unit(c(t = 0.5, r = 0, b = 0, l = 0), "cm")))+
    theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.5, b = 0, l = 0), "cm")))+
    theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
  
  
  f2<-ggplot(data=result.fixed, aes(x=cycle, y=mean_fixed)) +

    geom_line(size=1.2)+
    theme_bw()+
    ylab("Percentage of \n Fixed QTL")+
    xlab("Breeding Cycle")+
    theme(text = element_text(size=30),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    scale_y_continuous(breaks = seq(0, 100, 20), limits=c(0,100),expand=c(0,0),sec.axis = dup_axis(name = NULL, labels = NULL))+
    scale_x_continuous(breaks = seq(0, 50, 10), limits=c(1,50),expand=c(0,0),sec.axis = dup_axis(name = NULL, labels = NULL))+
    theme(axis.text = element_text(colour="black"))+
    theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
    theme(axis.text.x =element_text(margin=unit(c(t = 0.5, r = 0, b = 0, l = 0), "cm")))+
    theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.5, b = 0, l = 0), "cm")))+
    theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
  

  
  figure<-ggarrange(f1,f2,nrow=2)
  titel<-names[expir]
  titel<-substring(titel,25,nchar(titel)-6) 
  filename<-paste0(map,"/Results/Baseline",titel,"_1.eps")
  ggsave(filename, height=8 ,width= 10, f1,device="eps")
  filename<-paste0(map,"/Results/Baseline",titel,"_2.eps")
  ggsave(filename, height=8 ,width= 10, f2,device="eps")

  filename<-paste0(map,"/Results/",titel,"_1.jpeg")
  ggsave(filename, height=8 ,width= 10, f1,device="jpeg")
  filename<-paste0(map,"/Results/",titel,"_2.jpeg")
  ggsave(filename, height=8 ,width= 10,dpi=1000, f2,device="jpeg")
  print(expir)
}


