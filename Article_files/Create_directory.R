#Create directories to store simulation results



if (!dir.exists("own_results")){
  succes<-dir.create("own_results")
  if (!succes){
    print("error in create directory")
  }
}

if (!dir.exists("data/Backcrossing")){
  succes<-dir.create("data/Backcrossing")
  if (!succes){
    print("error in create directory")
  }
}

if (!dir.exists("data/BackCrossing_withscoping")){
  succes<-dir.create("data/BackCrossing_withscoping")
  if (!succes){
    print("error in create directory")
  }
}


if (!dir.exists("data/Baseline")){
  succes<-dir.create("data/Baseline")
  if (!succes){
    print("error in create directory")
  }
}

if (!dir.exists("data/Oracle")){
  succes<-dir.create("data/Oracle")
  if (!succes){
    print("error in create directory")
  }
}

if (!dir.exists("data/Scoping")){
  succes<-dir.create("data/Scoping")
  if (!succes){
    print("error in create directory")
  }
}


if (!dir.exists("data/Scoping_TP")){
  succes<-dir.create("data/Scoping_TP")
  if (!succes){
    print("error in create directory")
  }
}