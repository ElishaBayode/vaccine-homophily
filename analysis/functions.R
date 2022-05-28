library(deSolve)
library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)                             
library(plotly) 
library(ggpubr)




sir_homophily <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    N_0 <- S_0 + I_0 + R_0 #total vax status 1
    N_1 <- S_1 + I_1 + R_1 #total vax status 2
    N_2 <- S_2 + I_2 + R_2 #total vax status 3
    N_3 <- S_3 + I_3 + R_3 #total vax status 4
    
    
    N = N_0 + N_1 + N_2 + N_3 #total population 
    #lambda_ab means force of infection of group b by infected in group a 
    
    #force of infection = (total vax status specific contacts) (probability of transmission per contact)
    #(infectivity)  (susceptibility)  (proportion infected)  
    
    
    #((1-v*0)/1) =1 
    #v_i vaccine efficacy per number of doses 
    #l_i level of adherence to control measures 
    
    lambda_00 <-   prop_contact_00*total_0*beta*v_0*l_0*(I_0/N_0)
    lambda_01 <-   prop_contact_01*total_0*beta*v_1*l_0*(I_0/N_0)
    lambda_02 <-   prop_contact_02*total_0*beta*v_2*l_0*(I_0/N_0)
    lambda_03 <-   prop_contact_03*total_0*beta*v_3*l_0*(I_0/N_0)
    
    lambda_10 <-   prop_contact_10*total_1*beta*v_0*l_1*(I_1/N_1)
    lambda_11 <-   prop_contact_11*total_1*beta*v_1*l_1*(I_1/N_1)
    lambda_12 <-   prop_contact_12*total_1*beta*v_2*l_1*(I_1/N_1)
    lambda_13 <-   prop_contact_13*total_1*beta*v_3*l_1*(I_1/N_1)
    
    lambda_20 <-   prop_contact_20*total_2*beta*v_0*l_2*(I_2/N_2)
    lambda_21 <-   prop_contact_21*total_2*beta*v_1*l_2*(I_2/N_2)
    lambda_22 <-   prop_contact_22*total_2*beta*v_2*l_2*(I_2/N_2)
    lambda_23 <-   prop_contact_23*total_2*beta*v_3*l_2*(I_2/N_2)
    
    
    lambda_30 <-   prop_contact_30*total_3*beta*v_0*l_3*(I_3/N_3)
    lambda_31 <-   prop_contact_31*total_3*beta*v_1*l_3*(I_3/N_3)
    lambda_32 <-   prop_contact_32*total_3*beta*v_2*l_3*(I_3/N_3)
    lambda_33 <-   prop_contact_33*total_3*beta*v_3*l_3*(I_3/N_3)
    
    dS_0 <- -(lambda_00 + lambda_10 +lambda_20 + lambda_30)*S_0 + sigma_0*R_0 - f_0
    dI_0 <-  (lambda_00 + lambda_10 +lambda_20 + lambda_30)*(S_0 + R_0*(1-nu_0))   - gamma_0*I_0 + f_0
    dR_0 <-  gamma_0*I_0 - ((lambda_00 + lambda_10 +lambda_20 + lambda_30)*(1-nu_0) + sigma_0)*R_0
    
    dS_1 <- -(lambda_01+lambda_11+lambda_21 +lambda_31)*S_1 + sigma_1*R_1 - f_1
    dI_1 <-  (lambda_01+lambda_11+lambda_21 +lambda_31)*(S_1 + R_1*(1-nu_1))   - gamma_1*I_1 + f_1
    dR_1 <-  gamma_1*I_1 - ((lambda_01+lambda_11+lambda_21 +lambda_31)*(1-nu_1) + sigma_1)*R_1
    
    dS_2 <- -(lambda_02+lambda_12+lambda_22 +lambda_32)*S_2 + sigma_2*R_2 - f_2
    dI_2 <-  (lambda_02+lambda_12+lambda_22 +lambda_32)*(S_2 + R_2*(1-nu_2))   - gamma_2*I_2 + f_2
    dR_2 <-  gamma_2*I_2 - ((lambda_02+lambda_12+lambda_22 +lambda_32)*(1-nu_2) + sigma_2)*R_2
    
    dS_3 <- -(lambda_03+lambda_13+lambda_23 +lambda_33)*S_3 + sigma_3*R_3 - f_3
    dI_3 <-  (lambda_03+lambda_13+lambda_23 +lambda_33)*(S_3 + R_3*(1-nu_3))   - gamma_3*I_3 + f_3
    dR_3 <-  gamma_3*I_3 - ((lambda_03+lambda_13+lambda_23 +lambda_33)*(1-nu_3) + sigma_3)*R_3
    
    
    return(list(c(dS_0, dI_0, dR_0,dS_1, dI_1, dR_1,dS_2, dI_2, dR_2,dS_3, dI_3, dR_3)))
  })
}





make_initial_sate <- function(total_pop=total_pop, 
                              prop_unvax=prop_unvax, #87.6 at least 1 dose
                              prop_dose1= prop_dose1, #first dose 
                              prop_dose2=prop_dose2,
                              prop_dose3= prop_dose3, #from vax data
                              prop_wan_0=prop_wan_0,
                              prop_wan_1= prop_wan_1,
                              prop_wan_2=prop_wan_2,
                              prop_wan_3=prop_wan_3,
                              initial_inf=initial_inf,
                              prop_initial_inf_0=prop_initial_inf_0,
                              prop_initial_inf_1=prop_initial_inf_1,
                              prop_initial_inf_2=prop_initial_inf_2,
                              prop_initial_inf_3=prop_initial_inf_3){
  
  N_0 = total_pop - total_pop*(prop_dose1 + prop_dose2 + prop_dose3)
  N_1 = total_pop*prop_dose1
  N_2 = total_pop*prop_dose2
  N_3 = total_pop - (total_pop*(prop_dose1+prop_dose2) + N_0)
  R_0 = N_0*0.6 #60% recovered 
  I_0 = initial_inf*prop_initial_inf_0
  S_0 = N_0 - (R_0 + I_0)
  R_1 = N_1*0.6 #60% recovered  
  I_1 = initial_inf*prop_initial_inf_1
  S_1 = N_1 - (R_1 + I_1)
  R_2 = N_2*0.6 #60% recovered  
  I_2 = initial_inf*prop_initial_inf_2
  S_2 = N_2 - (R_2 + I_2)
  R_3 = N_3*0.6 #60% recovered  
  I_3 = initial_inf*prop_initial_inf_3
  S_3 = N_3 - (R_3 + I_3)
  
  
  init_state = c(S_0,I_0,R_0, S_1,I_1,R_1, S_2,I_2,R_2, S_3,I_3,R_3 )
  
  
  return(c(init_state))
}


#weighted average 
totat_0 <- last(output$S_0) + last(output$I_0) + last(output$R_0)
totat_1 <- last(output$S_1) + last(output$I_1) + last(output$R_1)
totat_2 <- last(output$S_2) + last(output$I_2) + last(output$R_2)
totat_3 <- last(output$S_3) + last(output$I_3) + last(output$R_3)

prop_totat_0 <- totat_0 / (totat_0 + totat_1 + totat_2 + totat_3)
prop_totat_1 <- totat_1 / (totat_0 + totat_1 + totat_2 + totat_3)
prop_totat_2 <- totat_2 / (totat_0 + totat_1 + totat_2 + totat_3)
prop_totat_3 <- totat_3 / (totat_0 + totat_1 + totat_2 + totat_3)

sum_weight <- prop_totat_0+prop_totat_1+prop_totat_2+prop_totat_3

average_weighted_conctact <- (23*prop_totat_0 + 25*prop_totat_1 + 
                                22*prop_totat_2 + 10*prop_totat_3)/sum_weight 

#Level of compliance by doese:  from Kiffer's data 
compliance_0 <- 1-0.134
compliance_1= 1-0.174
compliance_2= 1-0.349  
compliance_3= 1-0.817

average_weighted_complience <- (compliance_0*prop_totat_0 + compliance_1*prop_totat_1 + 
                                  compliance_2*prop_totat_2 + compliance_3*prop_totat_3)/sum_weight


parameters <- c(
  #force of infection parameters 
  prop_contact_00=0.39,
  prop_contact_01=0.04,
  prop_contact_02=0.1552,
  prop_contact_03=0.4148,
  prop_contact_10=0.16,
  prop_contact_11=0.18,
  prop_contact_12=0.18,
  prop_contact_13=0.48,
  prop_contact_20=0.1293,
  prop_contact_21=0.0107,
  prop_contact_22=0.62,
  prop_contact_23=0.24,
  prop_contact_30=0.0462,
  prop_contact_31=0.0038,
  prop_contact_32=0.27,
  prop_contact_33=0.68,
  #contact per week a1=23   a2=25  a3=22   a4=10
  total_0=23/7, #contact per day (from data)
  total_1=25/7,
  total_2=22/7,  
  total_3=10/7,
  beta = 0.21,# probability of transmission given contact
  k=1,
  #Vaccine efficacy
  #https://www.nejm.org/doi/full/10.1056/NEJMoa2119451
  v_0=1-0, #vax efficacy by number of doses. These are made up for now   
  v_1=1-0.1,
  v_2=1-0.148, #%65.5 after after 4 weeks, 10% after 25 weeks  
  v_3=1-0.74,
  # m=0.5, # effectiveness of control measure (I can mak ethis a function)
  l_0=compliance_0, # compliance level for vax 0 
  l_1=compliance_1, # compliance level for vax 1 
  l_2=compliance_2, # compliance level for vax 2 
  l_3=compliance_3, # compliance level for vax 3+ 
  
  #other parameters 
  
  sigma_0=1/(0.5*365),   # 1 year duration of immunity 
  sigma_1=1/(0.5*365),
  sigma_2=1/(0.5*365),
  sigma_3=1/(0.5*365),
  
  gamma_0=1/4,  #1/4 rocovery rate 
  gamma_1=1/4,
  gamma_2=1/4,
  gamma_3=1/4,
  
  # nu=0.95,  # protection for newly recovered 
  nu_0=0.35,
  nu_1=0.65,
  nu_2=0.68,
  nu_3=0.83,
  #importation parameter
  f_0 = 0,
  f_1 = 0,
  f_2 = 150,
  f_3 = 150
)





plot_output <- function(data=data){
  ggplot(data, aes(x=time, y=value, fill=key, col=key)) + geom_line() + theme_light() + 
    geom_text(data = . %>% filter(key == unique(key)),
              aes(label = round(as.numeric(key), 2)), size=2, col="black") +
    theme(axis.text=element_text(size=15),plot.title = element_text(size=15, face="bold"),
          legend.position = "none", legend.title = element_text(size=15),
          legend.text = element_text(size=15),
          axis.title=element_text(size=12,face="bold")) +
    labs(y="Ratio of Infections (with homophily to without homophily)", x ="Time (days)") + #ylim(c(0.01,10)) #+ 
    geom_hline(yintercept=1, linetype="solid", color = "black", size=0.5) 
  
  
  
  
}
