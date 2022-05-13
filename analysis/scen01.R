setwd("~/Desktop/vaccine-homophily/vaccine-homophily")
source("analysis/functions.R")


init_state <- make_initial_sate(total_pop=5.07e6, 
                                prop_unvax=(1-0.876), #87.6 at least 1 dose
                                prop_dose1=0.0292856, #first dose 
                                prop_dose2=0.3134464,
                                prop_dose3= 1-(0.3134464+0.0292856+0.124), #from vax data
                                prop_wan_0=0.5,
                                prop_wan_1=0.3,
                                prop_wan_2=0.3,
                                prop_wan_3=0.2,
                                prop_initial_inf_0=0.01,
                                prop_initial_inf_1=7e-5,
                                prop_initial_inf_2=5e-5,
                                prop_initial_inf_3=1e-5)

init <- c(S_0=init_state[1], I_0=init_state[2], R_0=init_state[3],
          S_1=init_state[4], I_1=init_state[5], R_1=init_state[6], S_2=init_state[7],
          I_2=init_state[8], R_2=init_state[9], S_3=init_state[10], I_3=init_state[11], R_3=init_state[12]) #initial condition 

##sanity check for force of infection
#prop_contact_00 <- seq(0.05,0.8, length=20)
# total_0 <- 2.5
# beta <- 0.4
# v <- 0.8
# m <- 0.5
#  l_0 <- 0
#  fI <-NULL
# # 
#  FI_0 <- function(i){
# #   fI <- prop_contact_00*total_0*beta*((1-v)/2)*(1-m*l_0)  
#    fI <- (1-v*1)/i
#    return(fI)
#  }
# # 
# i<- seq(1,5, length=50)
# FI_0(i)
# # 
#  plot(i,1-FI_0(i))

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
  beta = 0.7,# probability of transmission given contact
  k=1,
  v_0=1-0, #vax efficacy by number of doses. These are made up for now   
  v_1=1-0.3,
  v_2=1-0.6,
  v_3=1-0.85,
 # m=0.5, # effectiveness of control measure (I can mak ethis a function)
  l_0=1-0.1, # compliance level for vax 0 
  l_1=1-0.5, # compliance level for vax 1 
  l_2=1-0.7, # compliance level for vax 2 
  l_3=1-0.7, # compliance level for vax 3+ 
  
  #other parameters 
                
  sigma_0=1/(1*365),   # 1 year duration of immunity 
  sigma_1=1/(1*365),
  sigma_2=1/(1*365),
  sigma_3=1/(1*365),
  
  gamma_0=1/4,  #1/4 rocovery rate 
  gamma_1=1/4,
  gamma_2=1/4,
  gamma_3=1/4,
  
  nu=0.95  # protection for newly recovered 
                
)




times= seq(1, 365, 1)

#With homophily 

output <- as.data.frame(ode(y = init, times = times, func = sir_homophily , parms = parameters))



#Without homophily 

parameters_no_phil <- parameters 


parameters_no_phil["prop_contact_00"] <- 0.25
parameters_no_phil["prop_contact_01"] <- 0.25
parameters_no_phil["prop_contact_02"] <- 0.25
parameters_no_phil["prop_contact_03"] <- 0.25
parameters_no_phil["prop_contact_10"] <- 0.25
parameters_no_phil["prop_contact_11"] <- 0.25
parameters_no_phil["prop_contact_12"] <- 0.25
parameters_no_phil["prop_contact_13"] <- 0.25
parameters_no_phil["prop_contact_20"] <- 0.25
parameters_no_phil["prop_contact_21"] <- 0.25
parameters_no_phil["prop_contact_22"] <- 0.25
parameters_no_phil["prop_contact_23"] <- 0.25
parameters_no_phil["prop_contact_30"] <- 0.25
parameters_no_phil["prop_contact_31"] <- 0.25
parameters_no_phil["prop_contact_32"] <- 0.25
parameters_no_phil["prop_contact_33"] <- 0.25


parameters_no_phil["total_0"] <- 20/7
parameters_no_phil["total_1"] <- 20/7
parameters_no_phil["total_2"] <- 20/7
parameters_no_phil["total_3"] <- 20/7

parameters_no_phil["l_0"] <-  1-0.5
parameters_no_phil["l_1"] <-  1-0.5
parameters_no_phil["l_2"]  <- 1-0.5
parameters_no_phil["l_3"] <-  1-0.5


output_no_phil <- as.data.frame(ode(y = init, times = times, func = sir_homophily , parms = parameters_no_phil))


#make plot 

cols <- c("with homophily" = "orange",
        "without homophily"="darkgreen")

 scen01_plot <- ggplot() + geom_line(data=output,aes(x=times,y= (I_0+I_1+I_2+I_3)
           , colour = "with homophily") ,size=1.2,alpha=0.4) +
          geom_line(data=output_no_phil,aes(x=times,y= (I_0+I_1+I_2+I_3)
       ,  colour = "without homophily") ,size=1.2,alpha=0.4) +
        labs(y="Infections",x="Time (days)") + ylim(c(0,80000))  +theme_light() +
         scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
         plot.title = element_text(size=15, face="bold"),
         legend.position = "bottom", legend.title = element_text(size=15),
        legend.text = element_text(size=15),
        axis.title=element_text(size=15,face="bold")) +labs(color = " ") 

#make bar plot  
 final_size <- sum(output$I_0+output$I_1+output$I_2+output$I_3)
 peak_size <- max(output$I_0+output$I_1+output$I_2+output$I_3)

 final_size_no <- sum(output_no_phil$I_0+output_no_phil$I_1+
                        output_no_phil$I_2+output_no_phil$I_3)
 peak_size_no <- max(output_no_phil$I_0+output_no_phil$I_1+
                       output_no_phil$I_2+output_no_phil$I_3)

 bar_plot <- data.frame(mixing_pattern = c("with homophily","without homophily"), 
                        Total_Infection = c(final_size,final_size_no))
 
 
 scen01_bar_plot <- ggplot( bar_plot, aes(x=mixing_pattern, y=Total_Infection, colour=mixing_pattern, 
                       fill=mixing_pattern),size=0.2) + 
   geom_bar(stat = "identity", size=0.02, width = 0.2) + scale_fill_hue(c = 40) + 
   ylab("Total infection in 1 year") 
 
 #save plot 
 ggsave(file="figures/scen01_plot.png",  scen01_plot, width = 10, height = 8)
 ggsave(file="figures/scen01_bar_plot.png", scen01_bar_plot, width = 10, height = 8)

 
 
 
 
 
 
 
 
 
# 
# total_contact_1 <- pref_contact_1*(b_11*b_11*pref_contact_1*N_1 + b_12*b_21*pref_contact_2*N_2 +
#                                      b_13*b_31*pref_contact_3*N_3) / (pref_contact_1*N_1 + pref_contact_2*N_2 + pref_contact_3*N_3)
# 
# total_contact_2 <- pref_contact_2*(b_21*b_12*pref_contact_1*N_1 + b_22*b_22*pref_contact_2*N_2 +
#                                      b_23*b_32*pref_contact_3*N_3) / (pref_contact_1*N_1 + pref_contact_2*N_2 + pref_contact_3*N_3)
# 
# total_contact_3 <- pref_contact_3*(b_31*b_13*pref_contact_1*N_1 + b_32*b_23*pref_contact_2*N_2 +
#                                      b_33*b_33*pref_contact_3*N_3) / (pref_contact_1*N_1 + pref_contact_2*N_2 + pref_contact_3*N_3)
# 
# 
# 
# lambda_11 <-   total_contact_1 * (susceplity_1*infectivity_11*(1-exp(trasmEvent*contactDuratn_11)))* (I_1/N_1)
# lambda_12 <-   total_contact_1 * (susceplity_2*infectivity_12*(1-exp(trasmEvent*contactDuratn_12)))* (I_1/N_1)
# lambda_13 <-   total_contact_1 * (susceplity_3*infectivity_13*(1-exp(trasmEvent*contactDuratn_13)))* (I_1/N_1)
# 
# lambda_21 <-   total_contact_2 * (susceplity_1*infectivity_21*(1-exp(trasmEvent*contactDuratn_21)))* (I_2/N_2)
# lambda_22 <-   total_contact_2 * (susceplity_2*infectivity_22*(1-exp(trasmEvent*contactDuratn_22)))* (I_2/N_2)
# lambda_23 <-   total_contact_2 * (susceplity_3*infectivity_23*(1-exp(trasmEvent*contactDuratn_23)))* (I_2/N_2)
# 
# lambda_31 <-   total_contact_3 * (susceplity_1*infectivity_31*(1-exp(trasmEvent*contactDuratn_31)))* (I_3/N_3)
# lambda_32 <-   total_contact_3 * (susceplity_2*infectivity_32*(1-exp(trasmEvent*contactDuratn_32)))* (I_3/N_3)
# lambda_33 <-   total_contact_3 * (susceplity_3*infectivity_33*(1-exp(trasmEvent*contactDuratn_33)))* (I_3/N_3)
# 
# 
# dS_1 <- -(lambda_11+lambda_21+lambda_31)*S_1 + sigma_1*R_1 
# dI_1 <-  (lambda_11+lambda_21+lambda_31)*(S_1 + R_1*(1-nu))   - gamma_1*I_1
# dR_1 <-  gamma_1*I_1 - ((lambda_11+lambda_21+lambda_31) + sigma_1)*(1-nu)*R_1
# 
# dS_2 <- -(lambda_12+lambda_22+lambda_32)*S_2 + sigma_2*R_2 
# dI_2 <-  (lambda_12+lambda_22+lambda_32)*(S_2 + R_2*(1-nu))   - gamma_2*I_2
# dR_2 <-  gamma_2*I_2 - ((lambda_12+lambda_22+lambda_32) + sigma_2)*R_2*(1-nu)
# 
# dS_3 <- -(lambda_13+lambda_23+lambda_33)*S_3 + sigma_3*R_3 
# dI_3 <-  (lambda_13+lambda_23+lambda_33)*(S_3 + R_3*(1-nu))   - gamma_3*I_3
# dR_3 <-  gamma_3*I_3 - ((lambda_13+lambda_23+lambda_33) + sigma_3)*R_3*(1-nu)
# 
# dS_3 <- -(lambda_14+lambda_24+lambda_44)*S_4 + sigma_4*R_4 
# dI_4 <-  (lambda_14+lambda_24+lambda_44)*(S_4 + R_4*(1-nu))   - gamma_4*I_4
# dR_4 <-  gamma_4*I_4 - ((lambda_14+lambda_24+lambda_44) + sigma_4)*R_4*(1-nu)



