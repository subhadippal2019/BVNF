

fPlot_MPG_1<-function(m, type="Re"){
  str_code=paste0(type,"(ch_MPG(x*1i))")
  val=apply(matrix(m, ncol = 1), MARGIN = 1, FUN = function(x){ eval(parse(text= str_code )) });
  return(val)
}
#fdata=fPlot_MPG_1(t, type = "Im")






ch_MPG<-function(m){
  ch_MPG_one<-function(m1){
    #characteristic and mgf for MPG
    val1= alpha^2-m1
    val=BesselI(z = alpha, nu =nu )*(alpha^2-m1)^(nu/2)/(alpha^nu*BesselI(z = sqrt(val1), nu = nu))
    return(val)
  }
  val=apply(matrix(m, ncol = 1), MARGIN = 1, FUN = function(x){ (ch_MPG_one(x*1i)) });
  return(val)
  return(val)
}


ch_GIG<-function(m){
  ch_GIG_one<-function(m1){
    #characteristic and mgf for GIG
    val1= sqrt(alpha^2-m1)
    val=((val1/alpha)^(nu+1))*BesselK(z = val1,nu = -nu-1)/BesselK(z = alpha, nu = -nu-1)
    return(val)
  }
  val=apply(matrix(m, ncol = 1), MARGIN = 1, FUN = function(x){ (ch_GIG_one(x*1i)) });
  return(val)
}



#' comp_characteristic_function_alt
#'
#' @param t=c( -20000:20000);
#' @param alpha=10;
#' @param  nu=100
#' @examples
#' t=c( -20000:20000);  alpha=10; nu=100
#' comp_characteristic_function_alt(t,fun = c("ch_GIG", "ch_MPG"), c("GIG ","MPG"), showBcGround=FALSE, colr = c("red","blue") )
#' comp_characteristic_function_alt(t,fun = c("ch_GIG", "ch_MPG"), c("GIG ","MPG"), showBcGround=TRUE )
#' @export
comp_characteristic_function_alt<-function(t, fun=NULL, Name=NULL, showBcGround=TRUE,colr=c("white","black")){
  library(plotly)

  if(length(fun)>2){return(NULL)}

  f1= eval(parse(text=paste0( fun[1] ,"(t)")))
  x=Re(f1);y=Im(f1)
  data=data.frame(x=x, y=y, t=t/10)
  if(length(fun)>1){
    f2= eval(parse(text=paste0( fun[2] ,"(t)")))
    x1=Re(f2);y1=Im(f2)
    data$x1=x1
    data$y1=y1
  }





  ###### using plotly #####
  fig=plot_ly(data, x = ~x, y = ~y, z = ~t, type = 'scatter3d', mode = 'lines',
              line = list(color=colr[1], width = 1.2), name=Name[1])
  #rgba(11, 56, 128, 1)
  if(length(fun)==2){
    #fig= fig %>% add_trace(x = ~x1, y = ~y1, z = ~t,
    # line = list(color = 'rgba(207, 54, 37, .8)', width = 15))
    fig= fig %>% add_trace(x = ~x1, y = ~y1, z = ~t,
                           line = list(color = colr[2], width = 1.2), name=Name[2])
  }
  #rgba(207, 54, 37, 1)


  axx <- list(
    backgroundcolor="rgb(169,169,169)",
    gridcolor="rgb(255,255,255)",
    showbackground=showBcGround,
    zerolinecolor="rgb(255,255,255)",
    nticks=1,
    showticklabels=FALSE,
    title = "Real Axis"
  )

  axy <- list(
    backgroundcolor="rgb(169,169,169)",
    gridcolor="rgb(255,255,255)",
    showbackground=showBcGround,
    zerolinecolor="rgb(255,255,255)",
    nticks=1,
    showticklabels=FALSE,
    title = "Imaginary Axis"
  )

  axz <- list(
    backgroundcolor="rgb(169,169,169)",
    gridcolor="rgb(255,255,255)",
    showbackground=showBcGround,
    zerolinecolor="rgb(255,255,255)",
    nticks=1,
    showticklabels=FALSE,
    title = "t"
  )


  fig <- fig %>%layout(
    scene = list(
      xaxis = axx,
      yaxis = axy,
      zaxis = axz
    ))


  return(fig)


}






#########
#' comp_characteristic_functions_2d
#'
#' @examples
#' t=c( -20000:20000);  alpha=10; nu=100
#' aa=comp_characteristic_functions_2d(t,fun = c("ch_GIG", "ch_MPG"), c("GIG ","MPG") )
#' plot(aa)

#' t=c( -20000:20000)/80;  alpha=5; nu=0
#' aa=comp_characteristic_functions_2d(t,fun = c("ch_GIG", "ch_MPG"), c("GIG ","MPG") )
#' plot(aa)
#' @export
comp_characteristic_functions_2d<-function(t, fun=NULL, Name=NULL){

  if(length(fun)>2){return(NULL)}
  #browser()
  f1= eval(parse(text=paste0( fun[1] ,"(t)")))
  grp=rep(0, length(t) )
  if(length(fun)==1){
    data=data.frame(x=t, real=Re(f1), imng=Im(f1), grp=grp)
  }

  if(length(fun)>1){
    f2= eval(parse(text=paste0( fun[2] ,"(t)")))
    data=data.frame(
      t=c(t,t),
      grp=c(grp,rep(1, length(t))),
      real=c(Re(f1),  Re(f2)),
      imng=c(Im(f1), Im(f2))
    )
  }


  ###### using plotly #####
  fig=ggplot(data=data, aes(x=t, y=real, group=grp)) +
    geom_line(aes(color=grp), size=.7)+scale_color_gradient(low="black", high="white")+
    ylim(-1,1)+
    xlab("Real Axis")+
    ylab("Characteristic Functions ")+
    labs(color=("Group"))

  fig<-fig+theme(
    panel.background = element_rect(fill = "gray",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.002, linetype = 'solid',
                                    colour = "white"),
    panel.grid.minor = element_line(size = 0.001, linetype = 'solid',
                                    colour = "lightgray")
  )

  # fig <- fig %>%layout(
  #scene = list(
  #  xaxis = list(title = "t"),
  #  yaxis = list(title = "Real Axis"))
  # )


  #########
  fig1=ggplot(data=data, aes(x=t, y=imng, group=grp)) +
    #geom_line(aes(color=alpha), size=.5)+scale_color_gradient(low="skyblue", high="#9999CC")+
    #geom_line(aes(color=alpha), size=.5)+scale_color_gradient(low="white", high="black")+
    geom_line(aes(color=grp), size=.7)+scale_color_gradient(low="black", high="white")+
    ylim(-1,1)+
    xlab("Imaginary Axis")+
    ylab(" Characteristic Functions ") +
    labs(color=("Group"))

  fig1<-fig1+theme(
    panel.background = element_rect(fill = "gray",
                                    colour = "gray",
                                    size = 0.5, linetype = "solid"),
    panel.grid.major = element_line(size = 0.002, linetype = 'solid',
                                    colour = "white"),
    panel.grid.minor = element_line(size = 0.001, linetype = 'solid',
                                    colour = "lightgray")
  )



  fig_all=grid.arrange(fig, fig1, layout_matrix =rbind(c(1),c(2)) )
  list(fig,fig1)
  return(fig_all)


}





#' plot_char_ggplot_ray

#' @example
#' alpha=10; nu=100
#'  fun=c("ch_GIG", "ch_MPG")
#'  t<-seq(0, 20000, length.out=10000)
#'  plot_char_ggplot_ray(t = t, fun ="ch_GIG" )
#' @export
plot_char_ggplot_ray<-function(t, fun){

  #alpha=10; nu=100

  f1= eval(parse(text=paste0( fun[1] ,"(t)")))
  x=Re(f1);y=Im(f1)
  data=data.frame(x=x, y=y, t=t)

  mtcars_gg = ggplot(data) +
    geom_path(aes(x=x,color=t,y=y),size=1) +
    scale_color_continuous(limits=c(0,max(t))) +
    ggtitle(" ") +
    theme(title = element_text(size=8),
          text = element_text(size=12))


  plt=plot_gg(mtcars_gg, height=3, width=3.5, multicore=FALSE, pointcontract = 0.8, soliddepth=-20)
  return(plt)
}



