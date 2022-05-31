setwd("~/Desktop/vaccine-homophily/vaccine-homophily")
source("analysis/functions.R")
#dat = readRDS("data/BC-dat.rds")
#start_date <- ymd("2022-02-16")
#stop_date <- ymd("2022-03-03")
#dat <- dat %>% filter(date >= start_date &  date <= stop_date)
#import BC COVID-19 reported cases for the study period (Feb 16th to March 3rd 2022)

init <- c(S_0=last(output[,2]), I_0=last(output[,3]), R_0=last(output[,4]),
          S_1=last(output[,5]), I_1=last(output[,6]), R_1=last(output[,7]), S_2=last(output[,8]),
          I_2=last(output[,9]), R_2=last(output[,10]), S_3=last(output[,11]), I_3=last(output[,12]), 
          R_3=last(output[,13])) #initial condition 



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




times= 1:90

#With homophily 

output_homophily <- as.data.frame(ode(y = init, times = times, func = sir_homophily , parms = parameters)) %>%
          mutate(incid = (I_0+I_1+I_2+I_3)*parameters[["gamma_3"]]*asc_frac) %>% 
  mutate(date=seq.Date(ymd(start_date), ymd(start_date)-1+length(times), 1))


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


parameters_no_phil["total_0"] <- average_weighted_conctact/7
parameters_no_phil["total_1"] <- average_weighted_conctact/7
parameters_no_phil["total_2"] <- average_weighted_conctact/7
parameters_no_phil["total_3"] <- average_weighted_conctact/7

parameters_no_phil["l_0"] <-  average_weighted_complience
parameters_no_phil["l_1"] <-  average_weighted_complience
parameters_no_phil["l_2"]  <- average_weighted_complience
parameters_no_phil["l_3"] <-  average_weighted_complience


parameters_no_phil["nu_0"] <-  average_weighted_protection
parameters_no_phil["nu_1"] <-  average_weighted_protection
parameters_no_phil["nu_2"]  <- average_weighted_protection
parameters_no_phil["nu_3"] <-  average_weighted_protection




output_no_phil <- as.data.frame(ode(y = init, times = times, func = sir_homophily ,
                                    parms = parameters_no_phil))

#make plot 

cols <- c("with homophily" = "orange",
        "without homophily"="darkgreen")


#to do
#make a grid plot for ratio infection with and without homophily
#importation rate and transmission rate 
# 
# scen01_plot <- ggplot() + geom_line(data= output_homophily,aes(x=times, y= (I_0+I_1+I_2+I_3)
#                ,colour = "with homophily") ,size=1.2,alpha=0.4) +
#                 geom_line(data=output_no_phil,aes(x=times,y= (I_0+I_1+I_2+I_3)
#                ,colour = "without homophily") ,size=1.2,alpha=0.4) +
#                labs(y="Infections",x="Time (days)") + ylim(c(0,3500))  +theme_light() +
#               scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
#                plot.title = element_text(size=15, face="bold"),
#                legend.position = "bottom", legend.title = element_text(size=15),
#                legend.text = element_text(size=15),
#                axis.title=element_text(size=15,face="bold")) +labs(color = " ") 
# 
# ratio_homophily <- (output_no_phil$I_0+output_no_phil$I_1+ output_no_phil$I_2+output_no_phil$I_3)/
#   (output_homophily$I_0+output_homophily$I_1+output_homophily$I_2+output_homophily$I_3)
# 
# plot(ratio_homophily)
# 



beta = data.frame(beta=seq(0.1,0.24, length.out=20))

vary_parameter <- function(x, parameters=parameters,init1=init){
  parameters["beta"] <- x
 
  out_var <- as.data.frame(deSolve::ode(y=init1,time=times,func=sir_homophily ,
                                        parms=c(parameters))) # parameters[["p"]]*
  prev_var = (out_var$I_0 + out_var$I_1 + out_var$I_2 + out_var$I_3)
}
prev_var <- data.frame(apply(beta, 1, vary_parameter, parameters=parameters))


vary_parameter_no_phil <- function(x, parameters=parameters_no_phil,init1=init){
  parameters_no_phil["beta"] <- x

  out_var <- as.data.frame(deSolve::ode(y=init1,time=times,func=sir_homophily ,
                                        parms=parameters_no_phil)) # parameters[["p"]]*
  prev_var_no_phil = (out_var$I_0 + out_var$I_1 + out_var$I_2 + out_var$I_3)
}


prev_var_no_phil <- data.frame(apply(beta, 1, vary_parameter_no_phil, 
                                     parameters=parameters_no_phil))

colnames(prev_var_no_phil) <- beta[,1]
colnames(prev_var) <- beta[,1]

long_form <- gather(prev_var_no_phil)
long_form_phil <- gather(prev_var)

ratio_phil_no_phil  <- long_form

ratio_phil_no_phil$value <- long_form_phil$value/long_form$value

long_form$time <- times
long_form_phil$time <- times
ratio_phil_no_phil$time <- times

plot_ratio <- plot_output(data= ratio_phil_no_phil )
plot_no_phil <- plot_output(data=long_form)
plot_phil <- plot_output(data=long_form_phil)

ggarrange(plot_ratio , plot_no_phil, plot_phil, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)



ggsave(file="figures/var_noinport_plot.png",  var_inport_plot, width = 10, height = 8)


plot(NA,NA, xlim=c(0,180), ylim=c(0,5))
lines(prev_var$X1/prev_var_no_phil$X1, col="yellow")
lines(prev_var$X2/prev_var_no_phil$X2, col="green")
lines(prev_var$X3/prev_var_no_phil$X3, col="red")
lines(prev_var$X4/prev_var_no_phil$X4, col="blue")

plot(NA,NA, xlim=c(0,1800), ylim=c(0,8000))
lines(prev_var$X1)
lines(prev_var_no_phil$X1, col="yellow")
lines(prev_var$X2)
lines(prev_var_no_phil$X2, col="yellow")
lines(prev_var$X3) 
lines(prev_var_no_phil$X3, col="yellow")
lines(prev_var$X4)
lines(prev_var_no_phil$X4, col="yellow")


#at low infection levels high homophily is beneficial 
# at high infection levels high homophily is worse (with importation + mainly in vaccinated) 







scen01_plot <- ggplot() + geom_line(data= output_homophily,aes(x=times, y= (I_0+I_1+I_2+I_3)
                ,colour = "with homophily") ,size=1.2,alpha=0.4) +
  geom_line(data=output_no_phil,aes(x=times,y= (I_0+I_1+I_2+I_3)
                                    ,colour = "without homophily") ,size=1.2,alpha=0.4)
scen01_plot

plot(output_homophily$I_0)
lines(output_homophily$I_1)
lines(output_homophily$I_2)
lines(output_homophily$I_3)

plot(output_no_phil$I_0)
lines(output_no_phil$I_1)
lines(output_no_phil$I_2)
lines(output_no_phil$I_3)

#make bar plot  
 final_size <- sum(output_homophily$I_0+output_homophily$I_1+output_homophily$I_2+
                     output_homophily$I_3)
 peak_size <- max(output_homophily$I_0+output_homophily$I_1+output_homophily$I_2+output_homophily$I_3)

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
 scen01_bar_plot
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



