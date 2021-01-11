


#' plot_Acf_violiin_Slice
#' @param data :  A dataframe of size n*lag X number of different kappa seq.
#' @examples
#' data(ACF_Slice)
#' plot_Acf_violiin_Slice(ACF_Slice)
#' @export
plot_Acf_violiin_Slice<-function(data, MAXLAG=15,scale_type="width"){
dd=data
  #
# dd=data.frame( cat=(as.vector(replicate(20, 1:(MAXLAG)))))
# #NUMofData=length(AllList_3D)
# NUMofData=20
# #browser()
# FolderName=c("1","5","10","15")
# WorkSpIndex=1;dataIndex=1
#
# for( WorkSpIndex in 1:4 ){
#   acf_value=NULL
#   for( dataIndex in 1:NUMofData){
#     WorkSpaceName=paste0("kappa_",FolderName[WorkSpIndex],"/MC_kappa",FolderName[WorkSpIndex],"_data_",dataIndex ,".RData")
#     print(WorkSpaceName)
#     load(WorkSpaceName)
#     eval(parse(text=paste0("kappa_sample=MC_VONV$MC$MC$kappa")))
#     if(is.null(kappa_sample)){ eval(parse(text=paste0("kappa_sample=MC_VONV$MC$kappa")))}
#     #acf_val=(acf(kappa_sample, lag.max = MAXLAG,plot=FALSE))$acf
#         #eval(parse(text=paste0("kappa_sample=AllList_3D[[dataIndex]]$lst_",MC_CHOICE,"_3D$McSample_kappa")))
#     acf_value=cbind(acf_value,(acf(kappa_sample,lag.max = MAXLAG,plot=FALSE))$acf)
#   }
#
#   eval(parse(text=paste0("dd$val",WorkSpIndex,"=as.vector(acf_value[-1,])")))
#   rm(acf_value)
# }
#
 #browser()
 #eval(parse( text=paste0("data$",var_names[1])))
  var_names=names(data)
p <- ggplot(data, aes(factor(lag), kappa_1))
p <- p + geom_violin(aes(fill = "kappa=01"), alpha = 0.7, scale= scale_type)
#q <- p + geom_violin(aes(y = val2, fill = "kappa=01"), alpha = 0.7,scale = scale_type)
q<-p+  geom_violin(aes(y = kappa_5, fill = "kappa=05"), alpha = 0.7,scale =  scale_type)
q<-q+  geom_violin(aes(y = kappa_10, fill = "kappa=10"), alpha = 0.7,scale =  scale_type)
q<-q+  geom_violin(aes(y = kappa_15, fill = "kappa=15"), alpha = 0.7,scale =  scale_type)
#q + scale_fill_brewer(palette="Dark2") + theme_minimal()
#q + scale_fill_brewer(palette="Blues") + theme_classic()
q=q + theme( panel.grid.minor = element_blank())+ theme(legend.title=element_blank())+theme(legend.position="bottom")
q=q+xlab("Lag Numbers")+ylab("Autocorrelation")+ ylim(.95, 1)



return( plt=q)
}







CalculateTable_Slice<-function(df,lag_ind=5){
  out_list=list()
  Out_table=NULL
  TableFormat1=NULL
  for(iii in 1:lag_ind){

    current_out_table=cbind(Mean=apply(df[df$cat==iii, ], 2, mean), SD=apply(df[df$cat==iii, ], 2, sd))[-1,]
    Out_table=rbind(Out_table,current_out_table)
    TableFormat1=cbind(TableFormat1, paste( round(current_out_table[,1],3),"(", round(current_out_table[,2],2), ")") )


  }
  #browser()
  out_list$Out_table=Out_table
  num_of_cases=dim(df)[2]-1
  TableFormat=cbind(c(replicate(lag_ind, 1:num_of_cases)), paste( round(Out_table[,1],3),"(", round(1000*Out_table[,2],2), ")"))
  out_list$TableFormat= TableFormat
  out_list$TableFormat1=TableFormat1

  return(out_list)

}





#' plot_Acf_violiin_Slice
#' @param data :  A dataframe of size n*lag X number of different kappa seq.
#' @examples
#' data(ACF_Slice)
#' plot_Acf_violiin(ACF_Slice, Ylimit=c(.95, 1))
#' data(ACF_Sim1_3D_MPG_DA_vioplot_data)
#' plot_Acf_violiin(Vio_ACF_MPG, size=.2, ifGray=TRUE)
#' data(ACF_Sim1_3D_GEO_DA_vioplot_data)
#' plot_Acf_violiin(Vio_ACF_MPG, size=.2, ifGray=TRUE, GridLines="Vertical")
#' data(ACF_sim2_pD_MPG_DA_vioplot_data)
#' plot_Acf_violiin(Vio_ACF_MPG_pd,size=.2, ifGray=TRUE, GridLines="none")
#' @export
plot_Acf_violiin<-function(data, MAXLAG=15,scale_type="width",size=.5, Ylimit=NULL, ifGray=TRUE, GridLines="Vertical"){
  dd=data
  #
  # dd=data.frame( cat=(as.vector(replicate(20, 1:(MAXLAG)))))
  # #NUMofData=length(AllList_3D)
  # NUMofData=20
  # #browser()
  # FolderName=c("1","5","10","15")
  # WorkSpIndex=1;dataIndex=1
  #
  # for( WorkSpIndex in 1:4 ){
  #   acf_value=NULL
  #   for( dataIndex in 1:NUMofData){
  #     WorkSpaceName=paste0("kappa_",FolderName[WorkSpIndex],"/MC_kappa",FolderName[WorkSpIndex],"_data_",dataIndex ,".RData")
  #     print(WorkSpaceName)
  #     load(WorkSpaceName)
  #     eval(parse(text=paste0("kappa_sample=MC_VONV$MC$MC$kappa")))
  #     if(is.null(kappa_sample)){ eval(parse(text=paste0("kappa_sample=MC_VONV$MC$kappa")))}
  #     #acf_val=(acf(kappa_sample, lag.max = MAXLAG,plot=FALSE))$acf
  #         #eval(parse(text=paste0("kappa_sample=AllList_3D[[dataIndex]]$lst_",MC_CHOICE,"_3D$McSample_kappa")))
  #     acf_value=cbind(acf_value,(acf(kappa_sample,lag.max = MAXLAG,plot=FALSE))$acf)
  #   }
  #
  #   eval(parse(text=paste0("dd$val",WorkSpIndex,"=as.vector(acf_value[-1,])")))
  #   rm(acf_value)
  # }
  #
  #browser()
  #eval(parse( text=paste0("data$",var_names[1])))
  var_names=names(data)
  n_var=length(var_names)
  p <- ggplot(data, aes(factor(get(var_names[1])), get(var_names[2])))
  p <- p + geom_violin(aes(fill = var_names[2]), size=size,alpha = 0.9, scale= scale_type)


   #browser()
  if(n_var>=3){
   code_str=paste0("p<-p+geom_violin(aes(y =", (var_names[3:n_var]),", fill = '",var_names[3:n_var],"'),size=size, alpha = 0.7,scale =  scale_type)")
   eval(parse(text=code_str))
   }

  #q<-q+  geom_violin(aes(y = get(var_names[4]), fill = "kappa=10"), size=size,alpha = 0.7,scale =  scale_type)
  #q<-q+  geom_violin(aes(y = get(var_names[5]), fill = "kappa=15"), size=size,alpha = 0.7,scale =  scale_type)
  #q + scale_fill_brewer(palette="Dark2") + theme_minimal()
  #q + scale_fill_brewer(palette="Blues") + theme_classic()
  p=p + theme( panel.grid.minor = element_blank())+ theme(legend.title=element_blank())+theme(legend.position="bottom")
  p=p+xlab("Lag Numbers")+ylab("Autocorrelation")
  if( GridLines=="none"){p=p+theme(panel.grid.major = element_blank())}
  if( GridLines=="Vertical"){p=p+theme(panel.grid.major.y = element_blank())}


  if(!is.null(Ylimit)){
  p=p+ylim(Ylimit[1], Ylimit[2])
  }
  if(ifGray){p=p+scale_fill_grey(start = 1, end = 0)}
  return( plt=p)
}







CalculateTable_Slice<-function(df,lag_ind=5){
  out_list=list()
  Out_table=NULL
  TableFormat1=NULL
  for(iii in 1:lag_ind){

    current_out_table=cbind(Mean=apply(df[df$cat==iii, ], 2, mean), SD=apply(df[df$cat==iii, ], 2, sd))[-1,]
    Out_table=rbind(Out_table,current_out_table)
    TableFormat1=cbind(TableFormat1, paste( round(current_out_table[,1],3),"(", round(current_out_table[,2],2), ")") )


  }
  #browser()
  out_list$Out_table=Out_table
  num_of_cases=dim(df)[2]-1
  TableFormat=cbind(c(replicate(lag_ind, 1:num_of_cases)), paste( round(Out_table[,1],3),"(", round(1000*Out_table[,2],2), ")"))
  out_list$TableFormat= TableFormat
  out_list$TableFormat1=TableFormat1

  return(out_list)

}




######################
#' @export
VioplotMCAcfAll<-function(WorkSpaceName,MC_CHOICE="MPG",MAXLAG=15,IF_PLOT=TRUE,Plot_type="p",  scale_type="width"){

  dd=data.frame( cat=(as.vector(replicate(20, 1:(MAXLAG)))))
  #NUMofData=length(AllList_3D)
  NUMofData=20



  for( WorkSpIndex in 1:length(WorkSpaceName) ){
    load(WorkSpaceName[WorkSpIndex])
    eval(parse(text=paste0("kappa_sample=AllList_3D[[1]]$lst_",MC_CHOICE,"_3D$McSample_kappa")))
    acf_value=(acf(kappa_sample, lag.max = MAXLAG,plot=FALSE))$acf
    for( dataIndex in 2:NUMofData){
      eval(parse(text=paste0("kappa_sample=AllList_3D[[dataIndex]]$lst_",MC_CHOICE,"_3D$McSample_kappa")))
      acf_value=cbind(acf_value,(acf(kappa_sample,lag.max = MAXLAG,plot=FALSE))$acf)
    }



    eval(parse(text=paste0("dd$val",WorkSpIndex,"=as.vector(acf_value[-1,])")))
    rm(acf_value)
  }

  if(IF_PLOT){
    if(Plot_type=="p"){


      library(ggplot2)
      p <- ggplot(dd, aes(factor(cat), val1))
      p <- p + geom_violin(aes(fill = "p=010"), alpha = 0.7,scale = scale_type)
      q <- p + geom_violin(aes(y = val2, fill = "p=020"), alpha = 0.7,scale =  scale_type)
      q<-q+  geom_violin(aes(y = val3, fill = "p=050"), alpha = 0.7,scale =  scale_type)
      q<-q+  geom_violin(aes(y = val4, fill = "p=100"), alpha = 0.7,scale =  scale_type)
      #q + scale_fill_brewer(palette="Dark2") + theme_minimal()
      #q + scale_fill_brewer(palette="Blues") + theme_classic()
      q=q + theme( panel.grid.minor = element_blank())+ theme(legend.title=element_blank())+theme(legend.position="bottom")
      q=q+xlab("Lag Numbers")+ylab("Autocorrelation")
      (q)

    }
    if(Plot_type=="kappa"){


      library(ggplot2)
      p <- ggplot(dd, aes(factor(cat), val1))
      p <- p + geom_violin(aes(fill = "kappa=.01"), alpha = 0.7, scale= scale_type)
      q <- p + geom_violin(aes(y = val2, fill = "kappa=01"), alpha = 0.7,scale = scale_type)
      q<-q+  geom_violin(aes(y = val3, fill = "kappa=05"), alpha = 0.7,scale =  scale_type)
      q<-q+  geom_violin(aes(y = val4, fill = "kappa=10"), alpha = 0.7,scale =  scale_type)
      q<-q+  geom_violin(aes(y = val5, fill = "kappa=15"), alpha = 0.7,scale =  scale_type)
      #q + scale_fill_brewer(palette="Dark2") + theme_minimal()
      #q + scale_fill_brewer(palette="Blues") + theme_classic()
      q=q + theme( panel.grid.minor = element_blank())+ theme(legend.title=element_blank())+theme(legend.position="bottom")
      q=q+xlab("Lag Numbers")+ylab("Autocorrelation")
      (q)

    }


  }


  return(list(dd=dd, plt=q,WorkSpaceName=WorkSpaceName))
}



#' @export
Sim_workSpaceSummary<-function(){
library(ggplot2)
setwd("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim1/")
WorkSpaceName=paste0("AllList_3D_Kappa_eq_",c("0.01","1" ,"5","10","15"  ),".RData")
dda=VioplotMCAcfAll(WorkSpaceName,MAXLAG = 25,Plot_type = "kappa", MC_CHOICE = "MPG")
dda$plt+scale_fill_brewer(palette="Spectral")

Vio_ACF_MPG=dda$dd
names(Vio_ACF_MPG)=c( "lag", paste0("kappa",c("0.01","01" ,"05","10","15"  )))
save(Vio_ACF_MPG, file="ACF_3D_MPG_DA_vioplot_data.RData")
save(dda, file="ACF_3D_MPG_DA_vioplot_data_and_plot.RData" )



setwd("/Users/subhadippal/Dropbox/projects/Data Augmentation for Vonmises Distribution/WorkSpaces/Sim2/")
WorkSpaceName=paste0("AllList_",c(10,20,50,100),"D_Kappa_eq_15.RData")
dda=VioplotMCAcfAll(WorkSpaceName,MAXLAG = 25, Plot_type = "p")
dda$plt+scale_fill_brewer(palette="Spectral")

Vio_ACF_MPG_pd=dda$dd
names(Vio_ACF_MPG_pd)=c( "lag", paste0("dim_",c(10,20,50,100)))
save(Vio_ACF_MPG_pd, file="ACF_sim2_pD_MPG_DA_vioplot_data.RData")
save(dda, file="ACF_sim2_pD_MPG_DA_vioplot_data_and_plot.RData" )
}

