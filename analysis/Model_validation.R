source("analysis/functions.R")
source("analysis/setup.R")

#import BC COVID-19 reported cases for the study period (Feb 16th to March 3rd 2022)
asc_frac <- 0.25


#vaccine uptake set to Feb 16th levels 



init <- c(S_0=init_state[1], I_0=init_state[2], R_0=init_state[3],
          S_1=init_state[4], I_1=init_state[5], R_1=init_state[6], S_2=init_state[7],
          I_2=init_state[8], R_2=init_state[9], S_3=init_state[10], I_3=init_state[11], R_3=init_state[12]) #initial condition 


times= 1:nrow(dat)

#last(output[,2])
#With homophily 

output <- as.data.frame(ode(y = init, times = times, func = sir_homophily , parms = parameters_validate)) %>%
  mutate(incid = (I_0+I_1+I_2+I_3)*parameters_validate[["gamma_3"]]*asc_frac) %>% 
  mutate(date=seq.Date(ymd(start_date), ymd(start_date)-1+length(times), 1))

cols <- c("model" = "darkgreen",
          "data"="black")

Model_validation <- ggplot(data=dat,aes(x=date,y= value,  colour = "data")) + geom_point(alpha=0.6)  + 
  geom_line(alpha=0.2) +
  geom_line(data=output,aes(x=date,y= incid, colour = "model"), size =2, alpha=0.8) + 
  scale_x_date(date_breaks = "1 days", date_labels = "%b-%d-%y") +
  labs(y="Reported cases",x="Date") +
  theme_light() + scale_color_manual(values = cols) +  theme(axis.text=element_text(size=15),
                                                             plot.title = element_text(size=15, face="bold"),
                                                             legend.position = "bottom", legend.title = element_text(size=15),
                                                             axis.text.x = element_text(angle = 45, hjust = 1),
                                                             legend.text = element_text(size=15),
                                                             axis.title=element_text(size=15,face="bold")) +labs(color = " ") 
Model_validation

#ggsave(file="figures/Model_validation.png",  Model_validation, width = 10, height = 8)
