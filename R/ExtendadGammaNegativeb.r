  #samplerG for negative x values

rExtendedGamma_negative_b<-function(alpha,a,b){
  
  if(1-(a>0)*(b>0)*(alpha>0)){  print("The parameters a or b and alpha have to be positive."); return(NULL) } 
  cut_off_variable=1
  # cut_off_variable: decides which functional approximations to be used to crate the upper bound for the rejection sampler. 
  ACCEPT=FALSE
  
  
  
  if(alpha<=cut_off_variable){ # Applly Strategy1
    #print("Inside rExtendedGamma_negative_b")
  counter=0;
  while(!ACCEPT){
    counter=counter+1
    shape_par=(a+b)/(2*a+b)*alpha; rate_par=a+b; Exponent_power=(a+b)/(2*a+b);
    
    Exp_g<-(  rgamma(1,shape=shape_par,rate=rate_par) )
    g=Exp_g^(Exponent_power)	
    if(     log(runif(1)) <=  (  -a*g^2 -  b*g +  (a+b)*Exp_g     )  ){
      ACCEPT=T
     # print(paste("STEP1=",counter))
      return(g)
    }
    if(counter>100000){print("NO acceptence after 100000 iteration Procedure1");ACCEPT=T;}
  }
  }
  
  
  
  if(alpha>cut_off_variable){ # Applly Strategy1
  #x^(alpha-1)exp(-ax^2-bx). a>0,b>0, alpha>1
  x_mode=(sqrt(b^2+8*(alpha-1)*a)-b)/(4*a)
  scale=x_mode
  #Scaled_rv= original_rv/scale; original_rv=Scaled_rv*scale
  new_a=a* (scale)^2;new_b=b*(scale);
    ACCEPT=F
    counter=0;
    while(!ACCEPT){
      counter=counter+1
      shape_par=(new_a+new_b)/(2*new_a+new_b)*alpha; rate_par=new_a+new_b; Exponent_power=(new_a+new_b)/(2*new_a+new_b);
      
      Exp_g<-(  rgamma(1,shape=shape_par,rate=rate_par) )
      g=Exp_g^(Exponent_power)	
      if(     log(runif(1)) <=  (  -new_a*g^2 -  new_b*g +  (new_a+new_b)*Exp_g)  ){
        ACCEPT=T
       # print(paste("STEP2=",counter))
        x_sampled=g*scale
        return(x_sampled)
      }
      if(counter>100000){print("NO acceptence after 100000 iteration Procedure1");ACCEPT=T;}
    }
  }
return(NULL)
}








