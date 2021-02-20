
plot_trace_MPG_circular<-function(WorkSpIndex=1,dataIndex=1 ){
folderNames=c("01", "05", "10", "15")
folder_loc=paste0("/Users/subhadippal/Desktop/BCUP/Sim Circular/kappa",folderNames,"/")
fileNamePattern=c("", "kappa_5_","kappa_10_", "kappa_15_")

WorkSpaceName=paste0("MC_OBject_data_", fileNamePattern[WorkSpIndex],"id_",1:20,".RData")
WorkSpaceName_folder=paste0(folder_loc[WorkSpIndex], WorkSpaceName)

print(WorkSpaceName_folder[dataIndex])
load(WorkSpaceName_folder[dataIndex])
kappa_sample=MC_OBject$McSample_kappa


file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/Circular_MPG_trace_kappa_"

file=paste0(file,fileNamePattern[WorkSpIndex], "_data_",dataIndex, ".pdf"  )
pdf(file = file,width = 6.5, height =3.5)
plot(kappa_sample, type="l", xlab="Iteration Number", ylab=expression(paste("Posterior Samples of ",kappa,"") ))
dev.off()
rm(MC_OBject)
}

#setwd("/Users/subhadippal/Desktop/BCUP/Sim Circular/")

#dda=VioplotMCAcfAll(WorkSpaceName,MAXLAG = 25,Plot_type = "kappa", MC_CHOICE = "MPG")

circular_data_MPG_acf<-function(MAXLAG=25 ){


    dd=data.frame( cat=(as.vector(replicate(20, 1:(MAXLAG)))))
    # NUMofData=length(AllList_3D)
    NUMofData=20
     #browser()

    folderNames=c("01", "05", "10", "15")
    folder_loc=paste0("/Users/subhadippal/Desktop/BCUP/Sim Circular/kappa",folderNames,"/")
    fileNamePattern=c("", "kappa_5_","kappa_10_", "kappa_15_")
    for( WorkSpIndex in 1:4 ){
      #MC_OBject_data_kappa_5_id_2
      WorkSpaceName=paste0("MC_OBject_data_", fileNamePattern[WorkSpIndex],"id_",1:40,".RData")
      WorkSpaceName_folder=paste0(folder_loc[WorkSpIndex], WorkSpaceName)
      acf_value=NULL
      for( dataIndex in 1:NUMofData){
        #WorkSpaceName=paste0("kappa_",FolderName[WorkSpIndex],"/MC_kappa",FolderName[WorkSpIndex],"_data_",dataIndex ,".RData")
        dataIndex1=dataIndex
        if(WorkSpIndex==3){dataIndex1=dataIndex+20}

        print(WorkSpaceName_folder[dataIndex1])
        load(WorkSpaceName_folder[dataIndex1])
        eval(parse(text=paste0("kappa_sample=MC_OBject$McSample_kappa[1:5000]")))
        print(summary(MC_OBject$McSample_kappa))
       # if(is.null(kappa_sample)){ eval(parse(text=paste0("kappa_sample=MC_VONV$MC$kappa")))}
         acf_val=(acf(kappa_sample, lag.max = MAXLAG,plot=FALSE))$acf
             #eval(parse(text=paste0("kappa_sample=AllList_3D[[dataIndex]]$lst_",MC_CHOICE,"_3D$McSample_kappa")))
        acf_value=cbind(acf_value,(acf(kappa_sample,lag.max = MAXLAG,plot=FALSE))$acf)
      }

      eval(parse(text=paste0("dd$kappa_",folderNames[WorkSpIndex],"=as.vector(acf_value[-1,])")))
      rm(acf_value)
    }

  return(dd)
 }



CalculateTable_MPG_circular<-function(df,lag_ind=5){
  out_list=list()
  Out_table=NULL
  TableFormat1=NULL
  for(iii in 1:lag_ind){

    current_out_table=cbind(Mean=apply(df[df$cat==iii, ], 2, mean), SD=apply(df[df$cat==iii, ], 2, sd))[-1,]
    Out_table=rbind(Out_table,current_out_table)
    TableFormat1=cbind(TableFormat1, paste0( round(current_out_table[,1],3),"(", round(current_out_table[,2],3), ")") )


  }
  #browser()
  out_list$Out_table=Out_table
  num_of_cases=dim(df)[2]-1
  TableFormat=cbind(c(replicate(lag_ind, 1:num_of_cases)), paste0( round(Out_table[,1],3),"(",round(1000*Out_table[,2],3), ")"))
  out_list$TableFormat= TableFormat
  out_list$TableFormat1=TableFormat1

  return(out_list)

}



plot_density_MCMC_kappa_all<-function(WorkSpIndex=1,dataIndex=1, sel_section=c(2000:5000), choose=1 ){
  folderNames=c("01", "05", "10", "15")
  folder_loc=paste0("/Users/subhadippal/Desktop/BCUP/Sim Circular/kappa",folderNames,"/")
  fileNamePattern=c("", "kappa_5_","kappa_10_", "kappa_15_")

  WorkSpaceName=paste0("MC_OBject_data_", fileNamePattern[WorkSpIndex],"id_",1:20,".RData")
  WorkSpaceName_folder=paste0(folder_loc[WorkSpIndex], WorkSpaceName)
  kappaAll=NULL
  kappa_sample_all=NULL;Grp_index=NULL
  for(dataIndex in 1: 20){
  print(WorkSpaceName_folder[dataIndex])
  load(WorkSpaceName_folder[dataIndex])
  #eval(parse(text=paste0("kappa_sample",i ," =MC_OBject$McSample_kappa[sel_section]")))
  kappa_sample_all=c(kappa_sample_all,MC_OBject$McSample_kappa[sel_section])
  Grp_index=c(Grp_index, replicate(length(sel_section), dataIndex))
}

  df=data.frame(kappa=kappa_sample_all, Grp=Grp_index)
  file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/Circular_MPG_trace_kappa_"

  file=paste0(file,fileNamePattern[WorkSpIndex], "_data_",dataIndex, ".pdf" )
  #pdf(file = file,width = 6.5, height =3.5)
  #plot(kappa_sample, type="l", xlab="Iteration Number", ylab=expression(paste("Posterior Samples of ",kappa,"") ))
   #browser()
  #plot(density(kappa_sample), xlab="Support of the Distribution")
 # p1 <- ggplot(data=diamonds, aes(x=price, group=cut, fill=cut)) +
  #p1 <- ggplot(data=data.frame(kappa_sample), aes(x=kappa_sample)) +
   # geom_density(adjust=1.5)
    #theme_ipsum()


  p1 <- ggplot(data=df, aes(x=kappa, group=Grp, fill=Grp)) +
    geom_density(adjust=1.5, alpha=.4)


  #dev.off()
  rm(MC_OBject)
  return(p1)
}



plot_density_MCMC_kappa_all_circular_mpg<-function(WorkSpIndex=1,dataIndex=1, sel_section=c(2000:5000) ){
  folderNames=c("01", "05", "10", "15")
  folder_loc=paste0("/Users/subhadippal/Desktop/BCUP/Sim Circular/kappa",folderNames,"/")
  fileNamePattern=c("", "kappa_5_","kappa_10_", "kappa_15_")

  WorkSpaceName=paste0("MC_OBject_data_", fileNamePattern[WorkSpIndex],"id_",1:20,".RData")
  WorkSpaceName_folder=paste0(folder_loc[WorkSpIndex], WorkSpaceName)
  kappaAll=NULL
  kappa_sample_all=NULL;Grp_index=NULL
  for(dataIndex in 1: 20){
    print(WorkSpaceName_folder[dataIndex])
    load(WorkSpaceName_folder[dataIndex])
    #eval(parse(text=paste0("kappa_sample",i ," =MC_OBject$McSample_kappa[sel_section]")))
    kappa_sample_all=c(kappa_sample_all,MC_OBject$McSample_kappa[sel_section])
    Grp_index=c(Grp_index, replicate(length(sel_section), dataIndex))
  }

  df=data.frame(kappa=kappa_sample_all, Grp=Grp_index)
  file="/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WriteUp/V6.0/figures/Circular_MPG_trace_kappa_"

  file=paste0(file,fileNamePattern[WorkSpIndex], "_data_",dataIndex, ".pdf" )
  #pdf(file = file,width = 6.5, height =3.5)
  #plot(kappa_sample, type="l", xlab="Iteration Number", ylab=expression(paste("Posterior Samples of ",kappa,"") ))
  #browser()
  #plot(density(kappa_sample), xlab="Support of the Distribution")
  # p1 <- ggplot(data=diamonds, aes(x=price, group=cut, fill=cut)) +
  #p1 <- ggplot(data=data.frame(kappa_sample), aes(x=kappa_sample)) +
  # geom_density(adjust=1.5)
  #theme_ipsum()


  p1 <- ggplot(data=df, aes(x=kappa, group=Grp, fill=Grp)) +
    geom_density(size=.05,  alpha=.4)


  #dev.off()
  rm(MC_OBject)
  return(p1)
}





