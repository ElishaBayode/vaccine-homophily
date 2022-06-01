source("analysis/Model_validation.R")
times= 1:60
init <- c(S_0=last(output[,2]), I_0=last(output[,3]), R_0=last(output[,4]),
          S_1=last(output[,5]), I_1=last(output[,6]), R_1=last(output[,7]), S_2=last(output[,8]),
          I_2=last(output[,9]), R_2=last(output[,10]), S_3=last(output[,11]), I_3=last(output[,12]), 
          R_3=last(output[,13])) 



output_no_phil <- as.data.frame(ode(y = init, times = times, func = sir_homophily ,
                                    parms = parameters_no_phil))


output_homophily <- as.data.frame(ode(y = init, times = times, func = sir_homophily , parms = parameters)) %>%
mutate(incid = (I_0+I_1+I_2+I_3)*parameters[["gamma_3"]]*asc_frac) %>% 
  mutate(date=seq.Date(ymd(start_date), ymd(start_date)-1+length(times), 1))


#make plot 

cols <- c("with homophily" = "orange",
          "without homophily"="darkgreen")



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
