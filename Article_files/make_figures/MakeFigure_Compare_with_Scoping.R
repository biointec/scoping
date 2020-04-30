
library('dplyr')
library(ggplot2)
library(egg)
map<-"data/Compare_merrit"

if (file.exists(paste0(map,"/Results"))){
  unlink(paste0(map,"/Results"), recursive=TRUE)
}
#Voeg alle scopings toe aan dezelfde figuur
files<-dir(map, full.names = T)
names<-dir(map, full.names = F)
map<-paste0(map,"/Results")
dir.create(map)

results<-c()
results_v<-c()
level<-c()
fac<-c()

for(expir in seq(length(files))){
  a<-names[expir]
  load(files[expir])
  n.set<-length(experiment.sub.results)
  for(s in seq(n.set)){
    set <- experiment.sub.results[[s]]
    if(!is.null(set)  & (class(set) != "try-error")){
      nrep<-length(set)
      for(r in seq(nrep)){ #veiligheid indien er een fout gebeurd tijdens een repeat wordt deze overgeslagen
        rep<-set[[r]]
        
        if (!any(is.na(rep))){
          genome<-rep$genome
          sim<-rep$sim.results
          
         

          reachable<-c()
          fixed<-c()
          geneticValue_top<-c()
          geneticValue<-c()
          geneticVar<-c()
          acc<-c()
          
          for(br.cycle in seq(metadata$n.cycles)){
            gen<-sim[[br.cycle]]
            candidate.af.i<-gen$geno.summary.stats$candidate.af
            Candidate.value<-gen$candidate.values$geno.values
            Candidate.value.sort<-sort(Candidate.value,decreasing = TRUE)
            geneticValue_top<-cbind(geneticValue_top, mean(Candidate.value.sort[1:10])/gen$QTL_info[[1]]) #top 10
            geneticValue<-cbind(geneticValue, mean(Candidate.value)/gen$QTL_info[[1]])           #mean
            reachable<-cbind(reachable,   gen$QTL_info[['max_reachable']]/gen$QTL_info[[1]])
            geneticVar<-cbind(geneticVar, gen$candidate.values$true.var.components$V_g )
            fixed<-cbind(fixed, gen$QTL_info[[2]])
            acc<-cbind(acc, gen$prediction.accuracy)
           # fixed<-cbind(fixed, (gen$QTL_info[[3]]-(gen$QTL_info[[1]]-gen$QTL_info[[6]]))/gen$QTL_info[[1]])
          }
          # nodige data in 1 column zetten met factor for ggplot legend
          cycle<-seq(metadata$n.cycles)
          m1<-cbind(cycle, t(geneticValue),t(geneticValue_top), t(reachable),t(fixed), t(geneticVar), t(acc), expir)
          m<-rbind(m1)
          m<-data.frame(m)
          colnames(m)<-c("cycle","Genetic_value","Genetic_value_top","Reachable","Fixed","Genetic_variance", "Genetic_acc","Legend")
          results<-rbind(results, m)

        
        }
      }
    }
  }
  level<-c(level, toString(expir))
  
  fac<-c(fac, substring(a,1,nchar(names[expir])-6) )#voor de scoping_rep150
  
}
  results$Legend<-factor(results$Legend, labels = fac)
 
  #plot figure
  library('ggplot2')
  library(gridExtra)
  result.cycle<-c()
  result.cycle2<-c()
  result.cycle.top<-c()
  result.fixation<-c()
  result.variance<-c()
  result.acc<-c()
  
for (b.c in seq(metadata$n.cycles)){
   results %>%
   filter(cycle == b.c) %>%
   group_by(Legend) %>%
   summarise(mean_Genetic_Value = mean(Genetic_value), mean_reachable_Value = mean(Reachable), mean_variance = mean(Genetic_variance))->temp.cycle
   temp.cycle$cycle<-b.c
   result.cycle<-rbind(result.cycle, temp.cycle)
  
   results %>% 
   filter(cycle == b.c) %>%
   group_by(Legend) %>%
   summarise(mean_Genetic_Value = mean(Genetic_value))->temp.cycle
   temp.cycle$type<-paste("Mean GV Population" ,temp.cycle$Legend)
   temp.cycle$cycle<-b.c
   result.cycle2<-rbind(result.cycle2, temp.cycle)
   
   
   results %>% 
    filter(cycle == b.c) %>%
    group_by(Legend) %>%
    summarise(mean_Genetic_Value = mean(Genetic_value_top))->temp.cycle
    temp.cycle$type<-paste("Mean GV Top 10" ,temp.cycle$Legend)
    temp.cycle$cycle<-b.c
  result.cycle.top<-rbind(result.cycle.top, temp.cycle)
   
   results %>% 
    filter(cycle == b.c) %>%
   group_by(Legend) %>%
    summarise(mean_Genetic_Value = mean(Reachable))->temp.cycle
    temp.cycle$type<-paste("Maximum Reachable GV" ,temp.cycle$Legend)
    temp.cycle$cycle<-b.c
   #result.cycle2<-rbind(result.cycle2, temp.cycle)
   result.cycle.top<-rbind(result.cycle.top, temp.cycle)
   
   results %>% 
     filter(cycle == b.c) %>%
     group_by(Legend) %>%
     summarise(mean_Genetic_Value = mean(Fixed))->temp.cycle
   temp.cycle$type<-paste("QTL fixation" ,temp.cycle$Legend)
   temp.cycle$cycle<-b.c
   result.fixation<-rbind(result.fixation , temp.cycle)
   
   results %>% 
     filter(cycle == b.c) %>%
     group_by(Legend) %>%
     summarise(mean_Genetic_Value = mean(Genetic_variance))->temp.cycle
   temp.cycle$type<-paste("Genetic Variance" ,temp.cycle$Legend)
   temp.cycle$cycle<-b.c
   result.variance<-rbind(result.variance , temp.cycle)
   
   results %>% 
     filter(cycle == b.c) %>%
     group_by(Legend) %>%
     summarise(mean_Genetic_Value = mean(Genetic_acc,na.rm=T))->temp.cycle
   temp.cycle$type<-paste("Genetic Variance" ,temp.cycle$Legend)
   temp.cycle$cycle<-b.c
   result.acc<-rbind(result.acc , temp.cycle)
   
}
result.cycle2$type<-factor(result.cycle2$type)
result.cycle.top$type<-factor(result.cycle.top$type)

title<-'Genetic Value for fixed scoping factors'


  

fig1<-ggplot(data=result.cycle.top, aes(x=cycle, y=mean_Genetic_Value,colour=type,linetype=type)) +
  geom_line(size=1)  +
  scale_linetype_manual(values=c(1,2,3,4,5,6,7,8))+
  scale_color_manual(values=c("#000000","#0009bd","#00aeff","#ffd20a","#000000","#0009bd","#00aeff","#ffd20a"))+ 
  #scale_shape_manual(values=c(5,0,1,2,19,5,0,1,2,19))+
  ylab('Genetic Value \n(GV)') +
  xlab('Breeding Cycle')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(breaks = seq(0, 50, 10),expand=c(0,0),limits=c(1,50),sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme(legend.title =  element_blank(),text = element_text(size=9) , legend.position =  c(0.55,0.27),legend.direction = "vertical", panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.text = element_text(size=9), axis.text = element_text(colour="black"))+
  theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.15,"cm"))+
  theme(legend.key.width = unit(0.5,"cm"),legend.key.height = unit(0.1,"cm"), legend.key = element_rect(fill = NA))+
  theme(axis.text.x =element_text(margin=unit(c(t = 0.2, r = 0, b = 0, l = 0), "cm")))+
  theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.2, b = 0, l = 0), "cm")))+
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.),"cm"))+
  guides(linetype=guide_legend(ncol=1))+
  theme(legend.background = element_rect(fill=NA))

fig2<-ggplot(data=result.variance, aes(x=cycle, y=mean_Genetic_Value,colour=type,linetype=type)) +
  geom_line(size=1)  +
  scale_linetype_manual(values=c(6,7,8,9))+
  scale_color_manual(values=c("#000000","#0009bd","#00aeff","#ffd20a"))+ 
  scale_shape_manual(values=c(5,0,1,2,19))+
  ylab('Genetic Variance') +
  xlab('Breeding Cycle')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 10, 5),expand=c(0,0),limits=c(0,10),sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(breaks = seq(0, 50, 10),expand=c(0,0),limits=c(1,50),sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme(legend.title =  element_blank(),text = element_text(size=9) , legend.position =  c(0.50,0.84),legend.direction = "vertical", panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.text = element_text(size=9), axis.text = element_text(colour="black"))+
  theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.15,"cm"))+
  theme(legend.key.width = unit(0.5,"cm"),legend.key.height = unit(0.1,"cm"), legend.key = element_rect(fill = NA))+
  theme(axis.text.x =element_text(margin=unit(c(t = 0.2, r = 0, b = 0, l = 0), "cm")))+
  theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.2, b = 0, l = 0), "cm")))+
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.),"cm"))+
  guides(linetype=guide_legend(ncol=1))+
  theme(legend.background = element_rect(fill=NA))

fig3<-ggplot(data=result.acc, aes(x=cycle, y=mean_Genetic_Value,colour=type,linetype=type)) +
  geom_line(size=1)  +
  scale_linetype_manual(values=c(6,7,8,9,10))+
  scale_color_manual(values=c("#000000","#0009bd","#00aeff","#ffd20a"))+ 
  ylab('Accuracy') +
  xlab('Breeding Cycle')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(breaks = seq(0, 50, 10),expand=c(0,0),limits=c(1,50),sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme(legend.title =  element_blank(),text = element_text(size=9) , legend.position =  c(0.50,0.84),legend.direction = "vertical", panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.text = element_text(size=9), axis.text = element_text(colour="black"))+
  theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.15,"cm"))+
  theme(legend.key.width = unit(0.5,"cm"),legend.key.height = unit(0.1,"cm"), legend.key = element_rect(fill = NA))+
  theme(axis.text.x =element_text(margin=unit(c(t = 0.2, r = 0, b = 0, l = 0), "cm")))+
  theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.2, b = 0, l = 0), "cm")))+
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.),"cm"))+
  guides(linetype=guide_legend(ncol=1))+
  theme(legend.background = element_rect(fill=NA))

filename<-paste0(map,"/Compared_GV.eps")
ggsave(filename,heigh=14.5,width=20, fig1,device="eps")

filename<-paste0(map,"/Compared_Variance.eps")
ggsave(filename,heigh=14.5,width=20, fig2,device="eps")

filename<-paste0(map,"/Compared_Acc.eps")
ggsave(filename,heigh=14.5,width=20, fig3,device="eps")

#save the results to make latex figure
save("result.cycle.top", file=paste0(map,"/MakeFigure_CompareGV_withScoping.RData"))
save("result.variance", file=paste0(map,"/MakeFigure_CompareVar_withScoping.RData"))
save("result.acc", file=paste0(map,"/MakeFigure_CompareAcc_withScoping.RData"))
print("Program succesfully finished")


