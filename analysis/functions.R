library(deSolve)
library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyr)



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
    
    lambda_00 <-   prop_contact_00*total_0*beta*((1-v*0)/1)*(1-m*l_0)*(I_0/N_0)
    lambda_01 <-   prop_contact_01*total_0*beta*((1-v*k)/2)*(1-m*l_0)*(I_0/N_0)
    lambda_02 <-   prop_contact_02*total_0*beta*((1-v*k)/3)*(1-m*l_0)*(I_0/N_0)
    lambda_03 <-   prop_contact_03*total_0*beta*((1-v*k)/4)*(1-m*l_0)*(I_0/N_0)
    
    lambda_10 <-   prop_contact_10*total_1*beta*((1-v*0)/1)*(1-m*l_1)*(I_1/N_1)
    lambda_11 <-   prop_contact_11*total_1*beta*((1-v*k)/2)*(1-m*l_1)*(I_1/N_1)
    lambda_12 <-   prop_contact_12*total_1*beta*((1-v*k)/3)*(1-m*l_1)*(I_1/N_1)
    lambda_13 <-   prop_contact_13*total_1*beta*((1-v*k)/4)*(1-m*l_1)*(I_1/N_1)
    
    lambda_20 <-   prop_contact_20*total_2*beta*((1-v*0)/1)*(1-m*l_2)*(I_2/N_2)
    lambda_21 <-   prop_contact_21*total_2*beta*((1-v*k)/2)*(1-m*l_2)*(I_2/N_2)
    lambda_22 <-   prop_contact_22*total_2*beta*((1-v*k)/3)*(1-m*l_2)*(I_2/N_2)
    lambda_23 <-   prop_contact_23*total_2*beta*((1-v*k)/4)*(1-m*l_2)*(I_2/N_2)
    
    
    lambda_30 <-   prop_contact_30*total_3*beta*((1-v*0)/1)*(1-m*l_3)*(I_3/N_3)
    lambda_31 <-   prop_contact_31*total_3*beta*((1-v*k)/2)*(1-m*l_3)*(I_3/N_3)
    lambda_32 <-   prop_contact_32*total_3*beta*((1-v*k)/3)*(1-m*l_3)*(I_3/N_3)
    lambda_33 <-   prop_contact_33*total_3*beta*((1-v*k)/4)*(1-m*l_3)*(I_3/N_3)
    
    dS_0 <- -(lambda_00 + lambda_10 +lambda_20 + lambda_30)*S_0 + sigma_0*R_0 
    dI_0 <-  (lambda_00 + lambda_10 +lambda_20 + lambda_30)*(S_0 + R_0*(1-nu))   - gamma_0*I_0
    dR_0 <-  gamma_0*I_0 - ((lambda_00 + lambda_10 +lambda_20 + lambda_30)*(1-nu) + sigma_0)*R_0
    
    dS_1 <- -(lambda_01+lambda_11+lambda_21 +lambda_31)*S_1 + sigma_1*R_1 
    dI_1 <-  (lambda_01+lambda_11+lambda_21 +lambda_31)*(S_1 + R_1*(1-nu))   - gamma_1*I_1
    dR_1 <-  gamma_1*I_1 - ((lambda_01+lambda_11+lambda_21 +lambda_31)*(1-nu) + sigma_1)*R_1
    
    dS_2 <- -(lambda_02+lambda_12+lambda_22 +lambda_32)*S_2 + sigma_2*R_2 
    dI_2 <-  (lambda_02+lambda_12+lambda_22 +lambda_32)*(S_2 + R_2*(1-nu))   - gamma_2*I_2
    dR_2 <-  gamma_2*I_2 - ((lambda_02+lambda_12+lambda_22 +lambda_32)*(1-nu) + sigma_2)*R_2
    
    dS_3 <- -(lambda_03+lambda_13+lambda_23 +lambda_33)*S_3 + sigma_3*R_3
    dI_3 <-  (lambda_03+lambda_13+lambda_23 +lambda_33)*(S_3 + R_3*(1-nu))   - gamma_3*I_3
    dR_3 <-  gamma_3*I_3 - ((lambda_03+lambda_13+lambda_23 +lambda_33)*(1-nu) + sigma_3)*R_3
    
    
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
                              prop_initial_inf_0=prop_initial_inf_0,
                              prop_initial_inf_1=prop_initial_inf_1,
                              prop_initial_inf_2=prop_initial_inf_2,
                              prop_initial_inf_3=prop_initial_inf_3){
  
  N_0 = total_pop - total_pop*(prop_dose1 + prop_dose2 + prop_dose3)
  N_1 = total_pop*prop_dose1
  N_2 = total_pop*prop_dose2
  N_3 = total_pop - (total_pop*(prop_dose1+prop_dose2) + N_0)
  R_0 = N_0*0.6 #60% recovered 
  I_0 = N_0*prop_initial_inf_0
  S_0 = N_0 - (R_0 + I_0)
  R_1 = N_1*0.6 #60% recovered  
  I_1 = N_1*prop_initial_inf_1
  S_1 = N_1 - (R_1 + I_1)
  R_2 = N_2*0.6 #60% recovered  
  I_2 = N_2*prop_initial_inf_2
  S_2 = N_2 - (R_2 + I_2)
  R_3 = N_3*0.6 #60% recovered  
  I_3 = N_3*prop_initial_inf_3
  S_3 = N_3 - (R_3 + I_3)
  
  
  init_state = c(S_0,I_0,R_0, S_1,I_1,R_1, S_2,I_2,R_2, S_3,I_3,R_3 )
  
  
  return(c(init_state))
}

