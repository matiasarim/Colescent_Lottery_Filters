
#############################
## THE FUNCTION ITERATION.MODEL SIMULATE A COALESCENT FOLLOWED BY A LOTTERY METACOMMUNITY DYNAMIC IN A GIVEN LANDSCAPE
#
# replicas: is the number of times diversity is estimated in this pondscape
# Meta.pool: is the vector of species abundance in the pool a log normal distribution is assumed
# m.pool: is the migration from the species pool to a local community 
# Js: is the vector of total community size (numner of individuals) in each local community
# id.module: is the vector of module membership of each local community
# filter.env: is a matrix with ncol= number of local communities and nrow= number of species. 
#      filter.env elements are the relative probability [0,1] of species success in an environment
# M.dist: is a community x community matrix with distance between communities as elements


iteration.model<-function(replicas, Meta.pool, m.pool, Js, id.module, filter.env,
                          M.dist, D50, m.max,
                          id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                          Lottery, it, prop.dead.by.it, id.obs){                   
  library(doParallel)
  registerDoParallel(cores=detectCores()-1) #poner 2 para que divida los 2 nucleos, entonces divide esos dos nucleos, pero poner it bajos para que no se sature la maquina
  sub.it<-detectCores()-1
  tandas<-ceiling(replicas/sub.it)
  out.all<-NULL
  for(ttt in 1:tandas){
    #print(paste("iteration  ", ttt, " of ", tandas), sep="   ")    
    
    out<-iteration.model.par(sub.it=sub.it, tandas=tandas, Meta.pool=Meta.pool, m.pool=m.pool, Js=Js, id.module=id.module, filter.env=filter.env,
                             M.dist=M.dist, D50=D50, m.max=m.max,
                             id.fixed =id.fixed, D50.fixed=D50.fixed, 
                             m.max.fixed=m.max.fixed, comm.fixed=comm.fixed,
                             Lottery=Lottery, it=it, prop.dead.by.it=prop.dead.by.it, 
                             id.obs=id.obs)
    out.all<-rbind(out.all,out)
  }
  out.all
}

iteration.model.par<-function(sub.it, tandas,replicas, Meta.pool, m.pool, Js, id.module,id.env, filter.env,
                              M.dist, D50, m.max,
                              id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                              
                              Lottery, it, prop.dead.by.it, id.obs){                   
  #  library(doParallel)
  #  registerDoParallel(cores=detectCores()-2) ASÍ DECÍA MATIAS QUE LO DEJARAMOS
  #  sub.it<-cores=detectCores()
  #  tandas<-ceiling(replicas/sub.it)
  
  #    print(c("iteration  ", i, " of ", tandas))    
  out<-foreach(sss = 1:sub.it, .combine=rbind)%dopar% {
    H2020_Coalescent.and.lottery.exp.Kernel.J(Meta.pool=Meta.pool, m.pool=m.pool, Js=Js, id.module=id.module, filter.env=filter.env,
                                              M.dist=M.dist, D50=D50, m.max=m.max,
                                              id.fixed =id.fixed, D50.fixed=D50.fixed, 
                                              m.max.fixed=m.max.fixed, comm.fixed=comm.fixed,
                                              id.fixed2=id.fixed2, D50.fixed2=D50.fixed2, m.max.fixed2=m.max.fixed2, comm.fixed2=comm.fixed2,
                                              Lottery=Lottery, it=it, prop.dead.by.it=prop.dead.by.it, 
                                              id.obs=id.obs)
  }
  out
}

resume.out<-function(out){
  out2<-list()
  out2[[1]]<-apply(out,2,quantile, 0.5, na.rm=T)    # NA originates if only a single module is involved
  out2[[2]]<-apply(out,2,sd, na.rm=T)
  out2[[3]]<-apply(out,2,quantile, 0.975, na.rm=T)
  out2[[4]]<-apply(out,2,quantile, 0.025, na.rm=T)
  names(out2)<-c("Median", "Standard Deviation", "out.IC.up","out.IC.inf")
  out2
  
}

H2020_Coalescent.and.lottery.exp.Kernel.J<-function(Meta.pool, m.pool, Js, id.module, filter.env,
                                                    M.dist, D50, m.max,
                                                    id.fixed, D50.fixed, m.max.fixed, comm.fixed,
                                                    id.fixed2, D50.fixed2, m.max.fixed2, comm.fixed2,
                                                    Lottery, it, prop.dead.by.it, id.obs){
  
  if(length(Meta.pool)==1)Meta.pool<-round(rlnorm(n = Meta.pool, meanlog = log(100), sdlog = log(10)))
  library(vegan)
  Meta.pool<-Meta.pool/sum(Meta.pool)
  comm.fixed<-comm.fixed/sum(comm.fixed)
  comm.fixed2<-comm.fixed2/sum(comm.fixed2)
  Meta<-NULL
  
  # Function to update: M.migra from Graph to a function based in distance matrix
  M.migra<-H2020_migration.matrix.kernel.all(M.dist=M.dist, m.pool=m.pool, D50=D50, m.max=m.max,             # Funciton defined above. It estimates Migration matrix
                                             id.fixed=id.fixed, D50.fixed=D50.fixed, m.max.fixed=m.max.fixed,
                                             id.fixed2=id.fixed2, D50.fixed2=D50.fixed2,m.max.fixed2=m.max.fixed2)
  
  for(i in 1:ncol(M.migra)){
    Meta<-cbind(Meta, rmultinom(1,1,Meta.pool*filter.env[,i])) #Acá me da error
  }
  for (ii in 2:max(Js)){
    id.j<-which(Js>=ii)
    cat("coalescent construction in J: ", ii," de" ,max(Js),"\n")  
    if(length(id.fixed)>0)Meta[,id.fixed]<-comm.fixed*(ii-1)# scale vector of abundances in the fixed community to the abundance of all other communities
    if(length(id.fixed2)>0)Meta[,id.fixed2]<-comm.fixed2*(ii-1)
    Pool.neighbor<-(Meta%*%M.migra)         # estimates potential reclutants including immigrants for all communities weighted by local abundances
    Pool.neighbor<-Pool.neighbor*filter.env # IMPORTANT: element by element adjustment of species abundances to local filters
    if(length(id.j)>1){new<-apply(Pool.neighbor[,id.j],2,born,dead.by.it = 1, M.pool = Meta.pool, m.pool = m.pool)   # random selection of new individuals from reclutants pool 
    Meta[,id.j]<-Meta[,id.j]+new} else {
      Meta[,id.j]<-Meta[,id.j]+born(n = Pool.neighbor[,id.j],dead.by.it = 1, M.pool = Meta.pool, m.pool = m.pool) 
    }                          # upadate communities 
  }
  
  
  if(Lottery==T){                                        # START LOTTERY ################################################
    dead.by.it<-round(prop.dead.by.it*Js,0)    # estiamte individual to remove in each iteration and local community
    dead.by.it<-ifelse(dead.by.it<1,1,dead.by.it)
    if(length(id.fixed)>0) dead.by.it[id.fixed]<-0       # fixed communities are not updated
    if(length(id.fixed2)>0) dead.by.it[id.fixed2]<-0
    max.dead.by.it<-max(dead.by.it)
    if(length(id.fixed)>0)Meta[,id.fixed]<-round(comm.fixed*max(Js),0)# update abundances of the fixed community to community size
    if(length(id.fixed2)>0)Meta[,id.fixed2]<-round(comm.fixed2*max(Js),0)
    for(iteration in 1:it){                              # start lottery iterations  
      print(c(iteration, " of ", it))
      for(dead in 1:max.dead.by.it) {
        id.no.dead<-which(dead.by.it>=dead)               # identify communities to remove individuals
        if(length(id.no.dead)>1)Meta[,id.no.dead]<-Meta[,id.no.dead]-apply(Meta[,id.no.dead]*(1*filter.env[,id.no.dead]),2,FUN = change, change=1) # remove individuals along all communities IMPORTANT:dead is inverselly proportional to filter matrix (because matrix elements are performance)
        if(length(id.no.dead)==1)Meta[,id.no.dead]<-Meta[,id.no.dead]-change(Meta[,id.no.dead]*(1*filter.env[,id.no.dead]), change=1) # remove individuals along all communities
      }
      Pool.neighbor<-(Meta%*%M.migra)                    # estimates potential reclutant from all communities
      Pool.neighbor<-Pool.neighbor*filter.env # IMPORTANT: element by element adjustment of species abundances to local filters
      id.comm<-0
      for(reclutants in dead.by.it) {                    # dead.by.it is the vector of number of indiviuals to update in each iteration (a fixed fraction of J)
        id.comm<-id.comm+1
        Meta[,id.comm]<-Meta[,id.comm]+rmultinom(1,reclutants,
                                                 prob = (1-m.pool)*(Pool.neighbor[,id.comm]/sum(Pool.neighbor[,id.comm]))+m.pool*Meta.pool) # random selection of reclutants from neighbours, local or external pool
      }
    }
  }
  BB<-as.matrix(vegdist(t(Meta[,id.obs]), method = "jaccard"))
  Bett.all<-apply(BB,2,mean)
  Bett.intra<-rep(NA,length(id.module))
  Bett.inter<-rep(NA,length(id.module))
  for (bb in unique(id.module)){
    id.i<-which(id.module==bb)
    b.intra<-apply(BB[id.i,id.i],2,mean)
    Bett.intra[id.i]<-b.intra
    b.inter<-apply(BB[-id.i,id.i],2,mean)
    Bett.inter[id.i]<-b.inter
  }
  #  Meta.out<-c("m.pool"=m.pool, "Js.max"=max(Js),"Js.min"=min(Js), "D50"=D50, "m.max"=m.max,
  #             "D50.fixed"=D50.fixed, "m.max.fixed"=m.max.fixed,
  #            "Lottery"=ifelse(Lottery==TRUE,1,0), "it"=it, "Gamma"=Gamma,
  #            "S.loc"=apply(ifelse(Meta[,id.obs]>0,1,0),2,sum),
  #            "B.loc.all"=Bett.all,
  #            "B.loc.intra.module"=Bett.intra,
  #            "B.loc.inter.module"=Bett.inter)
  
  #Meta.out
  #Meta2<- rbind("B.loc.all"=Bett.all,
  #              "B.loc.intra.module"=Bett.intra,
  #             "B.loc.inter.module"=Bett.inter)
  Meta<-rbind(#Meta,
              "m.pool"=m.pool, "Js.max"=max(Js),"Js.min"=min(Js), "D50"=D50, "m.max"=m.max,
                      "D50.fixed"=D50.fixed, "m.max.fixed"=m.max.fixed, "D50.fixed2"=D50.fixed2, "m.max.fixed2"=m.max.fixed2,
                     "Lottery"=ifelse(Lottery==TRUE,1,0), "it"=it, 
                    "S.loc"=apply(ifelse(Meta[,id.obs]>0,1,0),2,sum),
              "B.loc.all"=Bett.all,
              "B.loc.intra.module"=Bett.intra,
              "B.loc.inter.module"=Bett.inter
  )
  Meta
  #  Meta2
}

H2020_migration.matrix.kernel.all<-function(M.dist, m.pool, D50,m.max, 
                                            id.fixed, D50.fixed, m.max.fixed,
                                            id.fixed2, D50.fixed2, m.max.fixed2){
  diag(M.dist)=NA
  M.dist = M.dist-min(M.dist, na.rm=T) # min distance is the distance between neighbour cells. 
  # connected cells has distance zero and migration m.max
  b = -log(0.5)/D50         # b is estimated | m(D50)=m.max*0.5
  M.migra = m.max*exp(-b*M.dist) 
  if(length(id.fixed)>0){
    b.fixed= -log(0.5)/D50.fixed 
    M.migra[id.fixed,]<- m.max.fixed*exp(-b.fixed*M.dist[id.fixed,]) # migration from outlet
  }                                                          # M.migra is the potential migration between communities, 
  diag(M.migra)<-1                                             # selfrecruitment is considered as 1=m.intra community
  M.migra<-apply(M.migra,2,m_to_1,m.pool)                      # standirize migrations to 1: (m.intra+m.pool+m.neigh=1)
  M.migra
  
  if(length(id.fixed2)>0){
    b.fixed2= -log(0.5)/D50.fixed2 
    M.migra[id.fixed2,]<- m.max.fixed2*exp(-b.fixed2*M.dist[id.fixed2,]) # migration from outlet
  }                                                          # M.migra is the potential migration between communities, 
  diag(M.migra)<-1                                             # selfrecruitment is considered as 1=m.intra community
  M.migra<-apply(M.migra,2,m_to_1,m.pool)                      # standirize migrations to 1: (m.intra+m.pool+m.neigh=1)
  M.migra
}

############
m_to_1<-function(m, m.pool) (1-m.pool)*m/sum(m) # standarize a vertor of migration to add 1 
# also considering migration from an external pool

born<-function(n, dead.by.it, M.pool, m.pool)rmultinom(1,dead.by.it,(1-m.pool)*(n/sum(n))+m.pool*M.pool)
change<-function(n,change)rmultinom(1,change,n)

