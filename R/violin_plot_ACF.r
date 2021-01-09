


#' plot_Acf_violiin_Slice
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


