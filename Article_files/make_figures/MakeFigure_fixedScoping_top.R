
library('dplyr')
library(ggplot2)
library(egg)
map<-"data/Scoping"

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
 

  a<-names[[expir]]
  
  if(substr(a,15,17)=="0.0" | substr(a,15,17)=="0.1" | substr(a,15,17)=="0.3" | substr(a,15,17)=="0.6" ){
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
          
          if( substr(a,15,17)=="0.0"){ #Data is lichtjes anders opgeslagen in de baseline metho
            reachable<-sim$cycle1$QTL_Gen0[[6]]/sim$cycle1$QTL_Gen0[[1]]
            fixed<-sim$cycle1$QTL_Gen0[[2]]
            geneticValue_top<-sim$cycle1$Genotype_Gen0_top/sim$cycle1$QTL_Gen0[[1]]
            geneticValue<-sim$cycle1$Genotype_Gen0/sim$cycle1$QTL_Gen0[[1]]
            geneticVar<-sim$cycle1$Genotype_Gen0/sim$cycle1$QTL_Gen0[[1]]
            
            
          }else{

          reachable<-sim$cycle1$QTL_base[[6]]/sim$cycle1$QTL_base[[1]]
          fixed<-sim$cycle1$QTL_base[[2]]
          geneticValue_top<-sim$cycle1$gen0_top/sim$cycle1$QTL_info[[1]]
          geneticValue<-sim$cycle1$gen0/sim$cycle1$QTL_info[[1]]
          geneticVar<-sim$cycle1$Generatie_nul$true.var.components$V_g
          }
          
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
           # fixed<-cbind(fixed, (gen$QTL_info[[3]]-(gen$QTL_info[[1]]-gen$QTL_info[[6]]))/gen$QTL_info[[1]])
          }
          # nodige data in 1 column zetten met factor for ggplot legend
          cycle<-seq(0,metadata$n.cycles)
          m1<-cbind(cycle, t(geneticValue),t(geneticValue_top), t(reachable),t(fixed), t(geneticVar), expir)
          m<-rbind(m1)
          m<-data.frame(m)
          colnames(m)<-c("cycle","Genetic_value","Genetic_value_top","Reachable","Fixed","Genetic_variance","Legend")
          results<-rbind(results, m)

        
        }
      }
    }
  }
  level<-c(level, toString(expir))
  if ( substr(a,15,17)=="0.0"){
    fac<-c(fac, paste("Baseline Method") )#voor de scoping_rep150
  }else{
  fac<-c(fac, paste0("Scoping Method (SR ",substring(names[expir],15,nchar(names[expir])-19)  ,")"))#voor de scoping_rep150
  }
#fac<-c(fac, substring(names[expir],0,nchar(names[expir])-22) ) #voor fixed_scoping_150rep
}#haakje van de if scoping == 0.1 0.3 0.5
}
  results$Legend<-factor(results$Legend, labels = fac)
 
  #plot figure
  library('ggplot2')
  library(gridExtra)
  result.cycle<-c()
  result.cycle2<-c()
  result.cycle.top<-c()
  result.fixation<-c()
for (b.c in seq(0,metadata$n.cycles)){
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
   
}
result.cycle2$type<-factor(result.cycle2$type)
result.cycle.top$type<-factor(result.cycle.top$type)

title<-'Genetic Value for fixed scoping factors'

fig1<-ggplot(data=result.cycle2, aes(x=cycle, y=mean_Genetic_Value,colour=type,linetype=type)) +
    geom_line(size=1.2)  +
   scale_linetype_manual(values=c("solid","solid","solid", "solid"))+
    scale_color_manual(values=c("gray80","gray60","gray40","black"))+ 
    ylab('Genetic Value \n(GV)') +
    xlab('Breeding Cycle')+
    theme(legend.position="none")+
    theme_bw()+
    scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
    scale_x_continuous(breaks = seq(0, 50, 10),expand=c(0,0),limits=c(1,50),sec.axis = dup_axis(name = NULL, labels = NULL))+
    theme(legend.title =  element_blank(),text = element_text(size=30) , legend.position = c(0.58,0.18),legend.direction = "vertical", panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(legend.text = element_text(size=25), axis.text = element_text(colour="black"))+
    theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
    theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"))+
    theme(axis.text.x =element_text(margin=unit(c(t = 0.7, r = 0, b = 0, l = 0), "cm")))+
    theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.7, b = 0, l = 0), "cm")))+
    theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))+
    guides(linetype=guide_legend(ncol=3)) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),legend.background = element_rect(fill=NA))

  

fig2<-ggplot(data=result.cycle.top, aes(x=cycle, y=mean_Genetic_Value,colour=type,linetype=type)) +
  geom_line(size=1.2)  +
  scale_linetype_manual(values=c("twodash", "twodash","twodash","twodash", "solid","solid","solid", "solid"))+
  scale_color_manual(values=c("gray80","gray60","gray40","black","gray80","gray60","gray40","black"))+ 
  ylab('Genetic Value \n(GV)') +
  xlab('Breeding Cycle')+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(breaks = seq(0, 50, 10),expand=c(0,0),limits=c(1,50),sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme(legend.title =  element_blank(),text = element_text(size=30) , legend.position = c(0.55,0.23),legend.direction = "vertical", panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.text = element_text(size=25), axis.text = element_text(colour="black"))+
  theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
  theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"))+
  theme(axis.text.x =element_text(margin=unit(c(t = 0.7, r = 0, b = 0, l = 0), "cm")))+
  theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.7, b = 0, l = 0), "cm")))+
  theme(plot.margin=unit(c(0.7,0.7,0.7,0.7),"cm"))+
  guides(linetype=guide_legend(ncol=2))+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),legend.background = element_rect(fill=NA))


fig3<-ggplot(data=result.fixation, aes(x=cycle, y=mean_Genetic_Value,colour=type,linetype=type)) +
  geom_line(size=1.2)  +
  scale_linetype_manual(values=c("twodash", "twodash","twodash","twodash", "solid","solid","solid", "solid"))+
  scale_color_manual(values=c("gray80","gray60","gray40","black","gray80","gray60","gray40","black"))+ 
  ylab('Percentage \n of Fixed QTLs') +
  xlab('Breeding Cycle')+
  theme(legend.position="none")+
  theme_bw()+
  scale_y_continuous(breaks = seq(0, 100, 20),expand=c(0,0),limits=c(0,100),sec.axis = dup_axis(name = NULL, labels = NULL))+
  scale_x_continuous(breaks = seq(0, 50, 10),expand=c(0,0),limits=c(1,50),sec.axis = dup_axis(name = NULL, labels = NULL))+
  theme(legend.title =  element_blank(),text = element_text(size=30) , legend.position = c(0.65,0.17),legend.direction = "vertical", panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
  theme(legend.text = element_text(size=25), axis.text = element_text(colour="black"))+
  theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
  theme(legend.key.width = unit(1,"cm"),legend.key.height = unit(0.8,"cm"))+
  theme(axis.text.x =element_text(margin=unit(c(t = 0.7, r = 0, b = 0, l = 0), "cm")))+
  theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.7, b = 0, l = 0), "cm")))+
  theme(plot.margin=unit(c(0.7,0.7,0.1,0.5),"cm"), legend.background = element_rect(fill=NA))+
  guides(linetype=guide_legend(ncol=2))

#fig3<-ggplot(data=result.cycle, aes(x=cycle, y=mean_reachable_Value,colour=Legend)) +
#  geom_line(size=1.2)  +
#  theme(plot.title = element_text(hjust = 0.5))+
#  ylab('Maximum reachable \n genetic value') +
#  theme(legend.position="bottom")+
#  guides(col=guide_legend(nrow=2))+
#  scale_y_continuous(breaks = seq(0, 1, 0.2),expand=c(0,0),limits=c(0,1),sec.axis = dup_axis(name = NULL, labels = NULL))+
#  scale_x_continuous(breaks = seq(0, 50, 10),expand=c(0,0),limits=c(0,50),sec.axis = dup_axis(name = NULL, labels = NULL))+
# theme(legend.title =  element_blank(),text = element_text(size=25) , legend.position= c(0.565,0.2),legend.direction = "vertical", legend.title.align  = 0.5,panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
#  theme(legend.text = element_text(size=19), axis.text = element_text(colour="black"))+
#  theme(axis.line = element_line(size=0.5),axis.ticks.length = unit(-0.25,"cm"))+
#  theme(legend.key.width = unit(1,"cm"))+    theme(axis.text.x =element_text(margin=unit(c(t = 0.5, r = 0, b = 0, l = 0), "cm")))+
#  theme(axis.text.y =element_text(margin=unit(c(t = 0, r = 0.5, b = 0, l = 0), "cm")))+
#  theme(plot.margin=unit(c(1,1,0.5,0.5),"cm"))
 

#fig3<-ggplot(data=result.cycle, aes(x=cycle, y=mean_variance, colour=Legend)) +
#  geom_line()+
#  theme(text = element_text(size=20))+
#  ylab('Genetic variance')
layout<-rbind(c(1),c(1),c(2),c(3))
figure<-grid.arrange(fig2,fig1,fig3, nrow=3,layout_matrix=layout)

filename<-paste0(map,"/Fixed_Scoping.eps")
ggsave(filename,heigh=14.5,width=20, figure,device="eps")

filename<-paste0(map,"/Fixed_Scoping.jpeg")
ggsave(filename,heigh=14,width=20, dpi=1000, fig1,device="jpeg")
#ggsave(paste0(map,"/Fixed_Scoping_variance.pdf"),heigh=4,width=8,fig3,device="pdf")
print("Program succesfully finished")


