source("analysis/functions.R")
source("analysis/setup.R")

times= 1:60

init <- c(S_0=last(output[,2]), I_0=last(output[,3]), R_0=last(output[,4]),
          S_1=last(output[,5]), I_1=last(output[,6]), R_1=last(output[,7]), S_2=last(output[,8]),
          I_2=last(output[,9]), R_2=last(output[,10]), S_3=last(output[,11]), I_3=last(output[,12]), 
          R_3=last(output[,13])) 
# 
# #With homophily 
# 
# output_homophily <- as.data.frame(ode(y = init, times = times, func = sir_homophily , parms = parameters)) %>%
#           mutate(incid = (I_0+I_1+I_2+I_3)*parameters[["gamma_3"]]*asc_frac) %>% 
#   mutate(date=seq.Date(ymd(start_date), ymd(start_date)-1+length(times), 1))
# 
# 
# #Without homophily 
# 
# 
# 
# output_no_phil <- as.data.frame(ode(y = init, times = times, func = sir_homophily ,
#      parms = parameters_no_phil)) %>%
#     mutate(incid = (I_0+I_1+I_2+I_3)*parameters[["gamma_3"]]*asc_frac) %>% 
#     mutate(date=seq.Date(ymd(start_date), ymd(start_date)-1+length(times), 1))


 

beta = data.frame(beta=c(seq(0.05,0.15, length.out=15)))

vary_parameter <- function(x, parameters=parameters,init1=init){
  parameters["beta"] <- x
  
  parameters["v_1"]  = 1 - (0.1)*6
  parameters["v_2"]  = 1 - (0.148)*6
  parameters["v_3"]  = 1 - (0.74)*1.25
  parameters["f_2"]  =  0
  parameters["f_3"]  =  0
  parameters["nu_0"]  <- 0.8
  parameters["nu_1"]  <- 0.85
  parameters["nu_2"]  <- 0.9
  parameters["nu_3"]  <- 0.97
  
  
  out_var <- as.data.frame(deSolve::ode(y=init1,time=times,func=sir_homophily ,
                                        parms=parameters))
  prev_var = (out_var$I_0 + out_var$I_1 + out_var$I_2 + out_var$I_3)
}

prev_var <- data.frame(apply(beta, 1, vary_parameter, parameters=parameters))


vary_parameter_no_phil <- function(x, parameters=parameters_no_phil,init1=init){
  parameters_no_phil["beta"] <- x
  
  parameters_no_phil["v_1"]  <- 1 - (0.1)*6
  parameters_no_phil["v_2"]  <- 1 - (0.148)*6
  parameters_no_phil["v_3"]  <- 1 - (0.74)*1.25
  parameters_no_phil["f_2"]  <- 0
  parameters_no_phil["f_3"]  <- 0
  parameters_no_phil["nu_0"]  <- 0.8
  parameters_no_phil["nu_1"]  <- 0.85
  parameters_no_phil["nu_2"]  <- 0.9
  parameters_no_phil["nu_3"]  <- 0.97

  out_var <- as.data.frame(deSolve::ode(y=init1,time=times,func=sir_homophily ,
                                        parms=parameters_no_phil))
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

plot_ratio <- plot_output(data= ratio_phil_no_phil ) + 
labs(y="Ratio of Infections (with homophily to without homophily)", x ="Time (days)")+
geom_hline(yintercept=1, linetype="solid", color = "black", size=0.5) 
plot_no_phil <- plot_output(data=long_form) + labs(y="Infections (without  homophily)", x ="Time (days)")
plot_phil <- plot_output(data=long_form_phil) + labs(y="Infections (with  homophily)", x ="Time (days)")

high_efficacy_low_trasm <- ggarrange(plot_ratio , plot_phil, plot_no_phil,  
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)


ggsave(file="figures/high_efficacy_low_trasm.png",  high_efficacy_low_trasm, width = 12, height = 6)



low_vary_parameter <- function(x, parameters=parameters,init1=init){
parameters["beta"] <- x
  
 parameters["v_1"]  <- 1 - (0.1)*0.1
 parameters["v_2"]  <- 1 - (0.148)*0.1
 parameters["v_3"]  <- 1 - (0.74)*0.1
 parameters["f_2"]  <- 0
 parameters["f_3"]  <- 0
 parameters_no_phil["nu_0"]  <- 0.2
 parameters_no_phil["nu_1"]  <- 0.4
 parameters_no_phil["nu_2"]  <- 0.65
 parameters_no_phil["nu_3"]  <- 0.8
  
  
  out_var <- as.data.frame(deSolve::ode(y=init1,time=times,func=sir_homophily ,
                                        parms=c(parameters))) # parameters[["p"]]*
  low_prev_var = (out_var$I_0 + out_var$I_1 + out_var$I_2 + out_var$I_3)
}

low_prev_var <- data.frame(apply(beta, 1, low_vary_parameter, parameters=parameters))


low_vary_parameter_no_phil <- function(x, parameters = parameters_no_phil,init1=init){
  parameters_no_phil["beta"] <- x
  
  
  parameters_no_phil["v_1"]  <- 1 - (0.1)*0.1
  parameters_no_phil["v_2"]  <- 1 - (0.148)*0.1
  parameters_no_phil["v_3"]  <- 1 - (0.74)*0.1
  parameters_no_phil["f_2"]  <- 0
  parameters_no_phil["f_3"]  <- 0
  parameters_no_phil["nu_0"]  <- 0.2
  parameters_no_phil["nu_1"]  <- 0.4
  parameters_no_phil["nu_2"]  <- 0.65
  parameters_no_phil["nu_3"]  <- 0.8
  
  
  out_var <- as.data.frame(deSolve::ode(y=init1,time=times,func=sir_homophily ,
                                        parms=parameters_no_phil)) 
  low_prev_var_no_phil = (out_var$I_0 + out_var$I_1 + out_var$I_2 + out_var$I_3)
}


low_prev_var_no_phil <- data.frame(apply(beta, 1, low_vary_parameter_no_phil, 
                                     parameters=parameters_no_phil))


colnames(low_prev_var_no_phil) <- beta[,1]
colnames(low_prev_var) <- beta[,1]

low_long_form <- gather(low_prev_var_no_phil)
low_long_form_phil <- gather(low_prev_var)

low_ratio_phil_no_phil  <- low_long_form

low_ratio_phil_no_phil$value <- low_long_form_phil$value/low_long_form$value

low_long_form$time <- times
low_long_form_phil$time <- times
low_ratio_phil_no_phil$time <- times

low_plot_ratio <- plot_output(data= low_ratio_phil_no_phil ) + 
labs(y="Ratio of Infections (with homophily to without homophily)", x ="Time (days)") +
geom_hline(yintercept=1, linetype="solid", color = "black", size=0.5) 

low_plot_no_phil <- plot_output(data=low_long_form) + labs(y="Infections (without homophily)", x ="Time (days)")
low_plot_phil <- plot_output(data=low_long_form_phil) + labs(y="Infections (with  homophily)", x ="Time (days)")


#Note--- High homophily in low vax settings acts differently to high homophily in high vax settings 

#high vaccine efficacy should make the effect of homophily stronger 



low_efficacy_low_trasm <- ggarrange(low_plot_ratio , low_plot_phil,  low_plot_no_phil, 
                           labels = c("A", "B", "C"),ncol = 3, nrow = 1)

ggsave(file="figures/low_efficacy_low_trasm.png",  low_efficacy_low_trasm, width = 12, height = 6)

