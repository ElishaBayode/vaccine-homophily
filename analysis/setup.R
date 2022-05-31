setwd("~/Desktop/vaccine-homophily/vaccine-homophily")
source("analysis/functions.R")
dat = readRDS("data/BC-dat.rds")
start_date <- ymd("2022-02-16")
stop_date <- ymd("2022-03-03")
dat <- dat %>% filter(date >= start_date &  date <= stop_date)

#dat = readRDS("data/BC-dat.rds")
#start_date <- ymd("2022-02-16")
#stop_date <- ymd("2022-03-03")
#dat <- dat %>% filter(date >= start_date &  date <= stop_date)
#import BC COVID-19 reported cases for the study period (Feb 16th to March 3rd 2022)




recov_rate <- 1/4
initial_inf <- (first(dat$value)/asc_frac)/recov_rate



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


#weighted average 
totat_0 <- init_state[1] + init_state[2] + init_state[3]
totat_1 <- init_state[4] + init_state[5] + init_state[6]
totat_2 <- init_state[7] + init_state[8] + init_state[9]
totat_3 <- init_state[10] + init_state[11] + init_state[12]

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

average_weighted_protection <- (0.35*prop_totat_0 + 0.65*prop_totat_1 + 
                                  0.68*prop_totat_2 + 0.83*prop_totat_3)/sum_weight



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

times= 1:60





# 
# #With homophily 
# 
# output_homophily <- as.data.frame(ode(y = init, times = times, func = sir_homophily , parms = parameters)) %>%
#           mutate(incid = (I_0+I_1+I_2+I_3)*parameters[["gamma_3"]]*asc_frac) %>% 
#   mutate(date=seq.Date(ymd(start_date), ymd(start_date)-1+length(times), 1))
# 
# 
# #Without homophily 

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




 
 #save plot 
 #ggsave(file="figures/scen01_plot.png",  scen01_plot, width = 10, height = 8)
 #ggsave(file="figures/scen01_bar_plot.png", scen01_bar_plot, width = 10, height = 8)

 
 
