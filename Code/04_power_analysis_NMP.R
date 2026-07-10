
# loading packages ----

library("dplyr")
library("ggplot2")
library("metafor")
library("ggpubr")


# load MP data ----

log_means = read.csv2("Data/NMP_log_means.csv", stringsAsFactors = T, fileEncoding = "latin1")


# calculating mean measurement values and corresponding variance ----

## within endpoints (ROS markers) ----

vars_control = data.frame(log_means %>% 
                            group_by(endpoint) %>% 
                            summarize(mean_value = mean(c_mean), mean_sd = sd(c_mean)))

vars_treatment = data.frame(log_means %>% 
                              group_by(endpoint) %>% 
                              summarize(mean_value = mean(t_mean), mean_sd = sd(t_mean)))

vars_control$endpoint = as.character(vars_control$endpoint)
vars_treatment$endpoint = as.character(vars_treatment$endpoint)


#100 t-tests per sample size

nrows = 98*nrow(vars_control)
alpha = 0.05

NMP_power = data.frame(
  sample_size = rep(NA, nrows),
  power_t_test = rep(NA, nrows),
  endpoint = rep(NA, nrows))

m = 1
set.seed(1234)
for(a in 1:nrow(vars_control)){
  control_mean_response = vars_control[a,2]
  control_sd = vars_control[a,3]
  treatment_mean_response = vars_treatment[a,2]
  treatment_sd = vars_treatment[a,3]
  name = vars_treatment[a,1]
  print(name)
  # loop over sample sizes
  for(k in seq(3,100)){
    sample_size = k
    for(l in 1){
      p.ttest = rep(NA, nrows)
      # repeat data simulation and analysis to get power estimates
      for(i in 1:100){
        simcontrol = rnorm(sample_size, mean = control_mean_response, sd = control_sd)
        simtreatment = rnorm(sample_size, mean = treatment_mean_response, sd = treatment_sd)
        # t test
        test = t.test(simcontrol, simtreatment)
        p.ttest[i] = test$p.value
      }
      NMP_power$sample_size[m] = k
      NMP_power$power_t_test[m] = mean(p.ttest <= alpha, na.rm = TRUE)
      NMP_power$endpoint[m] = name
      m = m + 1
    }
  }
}


over_80_percent_power = subset(NMP_power, power_t_test >= 0.8)
over_80_percent_power = data.frame(over_80_percent_power %>% group_by(endpoint) %>% slice_min(sample_size, n = 1))
power_at_n3 = data.frame(NMP_power %>% group_by(endpoint) %>% filter(sample_size %in% 3))
power_at_n8 = data.frame(NMP_power %>% group_by(endpoint) %>% filter(sample_size %in% 8))
power_at_n100 = data.frame(NMP_power %>% group_by(endpoint) %>% filter(sample_size %in% 100))
NMP_endpoints_power = rbind(power_at_n3, power_at_n8, power_at_n100, over_80_percent_power)

data.frame(power_at_n3 %>% summarize(mean(power_t_test), sd(power_t_test)))
data.frame(power_at_n8 %>% summarize(mean(power_t_test), sd(power_t_test)))

ros_power = subset(NMP_power, endpoint == "ROS")
cat_power = subset(NMP_power, endpoint == "CAT")
sod_power = subset(NMP_power, endpoint == "SOD")
gst_power = subset(NMP_power, endpoint == "GST")
pod_power = subset(NMP_power, endpoint == "POD")
sod_mRNA_power= subset(NMP_power, endpoint == "SOD mRNA expression")
cat_mRNA_power = subset(NMP_power, endpoint == "CAT mRNA expression")
gst_mRNA_power = subset(NMP_power, endpoint == "GST mRNA expression")
TAC_power = subset(NMP_power, endpoint == "TAC")
gsh_power = subset(NMP_power, endpoint == "GSH")
mda_power = subset(NMP_power, endpoint == "MDA")
OhdG_power = subset(NMP_power, endpoint == "8-OHdG")

ROS_plot = ggplot(data = ros_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "ROS concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
CAT_plot = ggplot(data = cat_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "CAT activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
SOD_plot = ggplot(data = sod_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "SOD activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
GST_plot = ggplot(data = gst_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "GST activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
POD_plot = ggplot(data = pod_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "POD activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
SOD_mRNA_plot = ggplot(data = sod_mRNA_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "SOD mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 83, linetype = "dashed", color = "red")
CAT_mRNA_plot = ggplot(data = cat_mRNA_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "CAT mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 10, linetype = "dashed", color = "red")
GST_mRNA_plot = ggplot(data = gst_mRNA_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "GST mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 7, linetype = "dashed", color = "red")
TAC_plot = ggplot(data = TAC_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "TAC") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 13, linetype = "dashed", color = "red")
GSH_plot = ggplot(data = gsh_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "GSH concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
MDA_plot = ggplot(data = mda_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "MDA concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
OhdG_plot = ggplot(data = OhdG_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "8-OhdG concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 

ggarrange(ROS_plot, CAT_plot, SOD_plot,
          GST_plot, POD_plot, SOD_mRNA_plot, 
          CAT_mRNA_plot, GST_mRNA_plot, TAC_plot,
          GSH_plot,  MDA_plot, OhdG_plot,
          nrow = 4,
          ncol = 3,
          legend = "none")


## within all combinations ----

vars_control = data.frame(log_means %>% 
                            group_by(endpoint,combined_concentration,polymer_type,shape,combined_size_mm,exposure_duration_effectsize_days) %>% 
                            summarize(mean_value = mean(c_mean), mean_sd = mean(c_SD))) # 867 combinations

vars_treatment = data.frame(log_means %>% 
                              group_by(endpoint,combined_concentration,polymer_type,shape,combined_size_mm, exposure_duration_effectsize_days) %>% 
                              summarize(mean_value = mean(t_mean), mean_sd = mean(t_SD)))

temporary_data = subset(vars_control, mean_sd == 0)
temporary_data2 = subset(vars_treatment, mean_sd == 0)

vars_control = vars_control[-as.numeric(row.names(temporary_data2)),]
vars_treatment = vars_treatment[-as.numeric(row.names(temporary_data)),]
vars_control = subset(vars_control, mean_sd != 0)
vars_treatment = subset(vars_treatment, mean_sd != 0) # 54 values with sd = 0 were removed
rm(temporary_data, temporary_data2)

vars_control$endpoint = as.character(vars_control$endpoint)
vars_treatment$endpoint = as.character(vars_treatment$endpoint)


nrows = 48*nrow(vars_control)
alpha = 0.05

NMP_power = data.frame(
  sample_size = rep(NA, nrows),
  power_t_test = rep(NA, nrows),
  endpoint = rep(NA, nrows))

m = 1
set.seed(1234)
for(a in 1:nrow(vars_control)){
  control_mean_response = vars_control[a,7]
  control_sd = vars_control[a,8]
  treatment_mean_response = vars_treatment[a,7]
  treatment_sd = vars_treatment[a,8]
  name = vars_treatment[a,1]
  # loop over sample sizes
  for(k in seq(3,50)){
    sample_size = k
    for(l in 1){
      p.ttest = rep(NA, nrows)
      # repeat data simulation and analysis to get power estimates
      for(i in 1:100){
        simcontrol = rnorm(sample_size, mean = control_mean_response, sd = control_sd)
        simtreatment = rnorm(sample_size, mean = treatment_mean_response, sd = treatment_sd)
        # t test
        test = t.test(simcontrol, simtreatment)
        p.ttest[i] = test$p.value
      }
      NMP_power$sample_size[m] = k
      NMP_power$power_t_test[m] = mean(p.ttest <= alpha, na.rm = TRUE)
      NMP_power$endpoint[m] = name
      m = m + 1
    }
  }
}


NMP_power = data.frame(NMP_power %>%
                         group_by(endpoint, sample_size) %>% 
                         summarize(mean_power = mean(power_t_test), sd_power = sd(power_t_test)))

over_80_percent_power = subset(NMP_power, mean_power >= 0.8)
over_80_percent_power = data.frame(over_80_percent_power %>% group_by(endpoint) %>% slice_min(sample_size, n = 1))
power_at_n3 = data.frame(NMP_power %>% group_by(endpoint) %>% filter(sample_size %in% 3))
power_at_n8 = data.frame(NMP_power %>% group_by(endpoint) %>% filter(sample_size %in% 8))
NMP_endpoints_power = rbind(power_at_n3, power_at_n8, over_80_percent_power)

ros_power = subset(NMP_power, endpoint == "ROS")
cat_power = subset(NMP_power, endpoint == "CAT")
sod_power = subset(NMP_power, endpoint == "SOD")
gst_power = subset(NMP_power, endpoint == "GST")
pod_power = subset(NMP_power, endpoint == "POD")
sod_mRNA_power= subset(NMP_power, endpoint == "SOD mRNA expression")
cat_mRNA_power = subset(NMP_power, endpoint == "CAT mRNA expression")
gst_mRNA_power = subset(NMP_power, endpoint == "GST mRNA expression")
TAC_power = subset(NMP_power, endpoint == "TAC")
gsh_power = subset(NMP_power, endpoint == "GSH")
mda_power = subset(NMP_power, endpoint == "MDA")
OhdG_power = subset(NMP_power, endpoint == "8-OHdG")


ROS_plot = ggplot(data = ros_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "ROS concentration") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 12, linetype = "dashed", color = "red") 
CAT_plot = ggplot(data = cat_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "CAT activity") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 16, linetype = "dashed", color = "red")
SOD_plot = ggplot(data = sod_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "SOD activity") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 8, linetype = "dashed", color = "red")
GST_plot = ggplot(data = gst_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "GST activity") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 9, linetype = "dashed", color = "red")
POD_plot = ggplot(data = pod_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "POD activity") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 12, linetype = "dashed", color = "red")
SOD_mRNA_plot = ggplot(data = sod_mRNA_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "SOD mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 5, linetype = "dashed", color = "red")
CAT_mRNA_plot = ggplot(data = cat_mRNA_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "CAT mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 3, linetype = "dashed", color = "red")
GST_mRNA_plot = ggplot(data = gst_mRNA_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "GST mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 6, linetype = "dashed", color = "red")
TAC_plot = ggplot(data = TAC_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "TAC") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 15, linetype = "dashed", color = "red")
GSH_plot = ggplot(data = gsh_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "GSH concentration") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 9, linetype = "dashed", color = "red")
MDA_plot = ggplot(data = mda_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "MDA concentration") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 15, linetype = "dashed", color = "red")
OhdG_plot = ggplot(data = OhdG_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + 
  labs(x = "sample size", y = "statistical power", title = "8-OhdG concentration") + theme_classic() + scale_x_continuous(breaks = c(3,5,10,15,20,25,30,35,40,45,50), labels = c(3,5,10,15,20,25,30,35,40,45,50), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 5, linetype = "dashed", color = "red")

ggarrange(ROS_plot, CAT_plot, SOD_plot,
          GST_plot, POD_plot, SOD_mRNA_plot, 
          CAT_mRNA_plot, GST_mRNA_plot, TAC_plot,
          GSH_plot,  MDA_plot, OhdG_plot,
          nrow = 4,
          ncol = 3,
          legend = "none")
