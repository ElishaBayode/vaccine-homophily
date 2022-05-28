
source("analysis/functions.R")
dat = readRDS("data/BC-dat.rds")
start_date <- ymd("2022-02-16")
stop_date <- ymd("2022-03-03")
dat <- dat %>% filter(date >= start_date &  date <= stop_date)
#import BC COVID-19 reported cases for the study period (Feb 16th to March 3rd 2022)
asc_frac <- 0.25

parameters_validate <- parameters

parameters_validate["f_0"] <- 0
parameters_validate["f_1"] <- 0
parameters_validate["f_2"] <- 25
parameters_validate["f_3"] <- 25

initial_inf <- (first(dat$value)/asc_frac)/parameters_validate[["gamma_3"]]


#vaccine uptake set to Feb 16th levels 

init_state <- make_initial_sate(total_pop=5.07e6, 
                                prop_unvax=(1-0.85), #87.6 at least 1 dose
                                prop_dose1=0.05, #first dose 
                                prop_dose2=0.36,#second dose
                                prop_dose3= 1-(0.15+0.05+0.36), #from vax data
                                prop_wan_0=0.5,
                                prop_wan_1=0.3,
                                prop_wan_2=0.3,
                                prop_wan_3=0.2,
                                initial_inf =initial_inf,
                                prop_initial_inf_0=0.15,
                                prop_initial_inf_1=0.05,
                                prop_initial_inf_2=0.36,
                                prop_initial_inf_3=0.44)

init <- c(S_0=init_state[1], I_0=init_state[2], R_0=init_state[3],
          S_1=init_state[4], I_1=init_state[5], R_1=init_state[6], S_2=init_state[7],
          I_2=init_state[8], R_2=init_state[9], S_3=init_state[10], I_3=init_state[11], R_3=init_state[12]) #initial condition 


times= 1:nrow(dat)

last(output[,2])
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

ggsave(file="figures/Model_validation.png",  Model_validation, width = 10, height = 8)
