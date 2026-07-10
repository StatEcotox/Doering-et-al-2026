
# loading packages ----

library("dplyr")
library("ggplot2")
library("metafor")
library("ggpubr")


# load Neonic data ----

log_means = read.csv2("Data/Neonics_log_means.csv", stringsAsFactors = T, fileEncoding = "latin1")

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


# 100 t-tests per sample size

nrows = 98*nrow(vars_control)
alpha = 0.05

Neonic_power = data.frame(
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
        stacked = stack(as.data.frame(cbind(simcontrol, simtreatment)))
      }
      Neonic_power$sample_size[m] = k
      Neonic_power$power_t_test[m] = mean(p.ttest <= alpha, na.rm = TRUE)
      Neonic_power$endpoint[m] = name
      m = m + 1
    }
  }
}

over_80_percent_power = subset(Neonic_power, power_t_test >= 0.8)
over_80_percent_power = data.frame(over_80_percent_power %>% group_by(endpoint) %>% slice_min(sample_size, n = 1))
power_at_n3 = data.frame(Neonic_power %>% group_by(endpoint) %>% filter(sample_size %in% 3))
power_at_n8 = data.frame(Neonic_power %>% group_by(endpoint) %>% filter(sample_size %in% 8))
power_at_n100 = data.frame(Neonic_power %>% group_by(endpoint) %>% filter(sample_size %in% 100))
Neonic_endpoints_power = rbind(power_at_n3, power_at_n8, power_at_n100, over_80_percent_power)

data.frame(power_at_n3 %>% summarize(mean(power_t_test), sd(power_t_test)))
data.frame(power_at_n8 %>% summarize(mean(power_t_test), sd(power_t_test)))

ros_power = subset(Neonic_power, endpoint == "ROS")
hydroxyl_power = subset(Neonic_power, endpoint == "OH")
cat_power = subset(Neonic_power, endpoint == "CAT")
sod_power = subset(Neonic_power, endpoint == "SOD")
gst_power = subset(Neonic_power, endpoint == "GST")
pod_power = subset(Neonic_power, endpoint == "POD")
care_power = subset(Neonic_power, endpoint == "CarE")
gr_power = subset(Neonic_power, endpoint == "GR") 
sod_mRNA_power = subset(Neonic_power, endpoint == "SOD mRNA expression")
cat_mRNA_power = subset(Neonic_power, endpoint == "CAT mRNA expression")
gst_mRNA_power = subset(Neonic_power, endpoint == "GST mRNA expression")
gsh_power = subset(Neonic_power, endpoint == "GSH")
mda_power = subset(Neonic_power, endpoint == "MDA")
OhdG_power = subset(Neonic_power, endpoint == "8-OHdG")
pc_power = subset(Neonic_power, endpoint == "PC")

ROS_plot = ggplot(data = ros_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "ROS concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
Hydroxyl_plot = ggplot(data = hydroxyl_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "OH concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 11, linetype = "dashed", color = "red")
CAT_plot = ggplot(data = cat_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "CAT activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
SOD_plot = ggplot(data = sod_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "SOD activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
GST_plot = ggplot(data = gst_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "GST activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
POD_plot = ggplot(data = pod_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "POD activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
CarE_plot = ggplot(data = care_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "CarE activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
GR_plot = ggplot(data = gr_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "GR activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
SOD_mRNA_plot = ggplot(data = sod_mRNA_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "SOD mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 28, linetype = "dashed", color = "red")
CAT_mRNA_plot = ggplot(data = cat_mRNA_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "CAT mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 60, linetype = "dashed", color = "red")
GST_mRNA_plot = ggplot(data = gst_mRNA_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "GST mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
GSH_plot = ggplot(data = gsh_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "GSH concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 97, linetype = "dashed", color = "red")
MDA_plot = ggplot(data = mda_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "MDA concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
OhdG_plot = ggplot(data = OhdG_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "8-OhdG concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 11, linetype = "dashed", color = "red")
PC_plot = ggplot(data = pc_power, aes(x = sample_size, y = power_t_test)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = power_t_test-sd_power, ymax = power_t_test+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "PC concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(10, 100, 10)), labels = c(3,seq(10, 100, 10)), expand = c(0,0), limits = c(1,102)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 

ggarrange(ROS_plot, Hydroxyl_plot, CAT_plot,
          SOD_plot, GST_plot, POD_plot, 
          CarE_plot, GR_plot, SOD_mRNA_plot,
          CAT_mRNA_plot, GST_mRNA_plot, GSH_plot,
          MDA_plot, OhdG_plot, PC_plot,
          nrow = 5,
          ncol = 3,
          legend = "none")

## within all combinations ----

vars_control = data.frame(log_means %>% 
                            group_by(endpoint,neonic,nominal_dose_mass,exposure_duration_effectsize_days) %>% 
                            summarize(mean_value = mean(c_mean), mean_sd = mean(c_SD)))

vars_treatment = data.frame(log_means %>% 
                              group_by(endpoint,neonic,nominal_dose_mass,exposure_duration_effectsize_days) %>% 
                              summarize(mean_value = mean(t_mean), mean_sd = mean(t_SD)))

temporary_data = subset(vars_control, mean_sd == 0)
temporary_data2 = subset(vars_treatment, mean_sd == 0)

vars_control = vars_control[-as.numeric(row.names(temporary_data2)),]
vars_treatment = vars_treatment[-as.numeric(row.names(temporary_data)),]
vars_control = subset(vars_control, mean_sd != 0)
vars_treatment = subset(vars_treatment, mean_sd != 0) # 12 values with sd = 0 were removed
rm(temporary_data, temporary_data2)

vars_control$endpoint = as.character(vars_control$endpoint)
vars_treatment$endpoint = as.character(vars_treatment$endpoint)


# 100 t-tests per sample size

nrows = 48*nrow(vars_control)
alpha = 0.05

Neonic_power = data.frame(
  sample_size = rep(NA, nrows),
  power_t_test = rep(NA, nrows),
  endpoint = rep(NA, nrows))

m = 1
set.seed(1234)
for(a in 1:nrow(vars_control)){
  control_mean_response = vars_control[a,5]
  control_sd = vars_control[a,6]
  treatment_mean_response = vars_treatment[a,5]
  treatment_sd = vars_treatment[a,6]
  name = vars_treatment[a,1]
  #print(name)
  # loop over sample sizes
  for(k in seq(3,50)){
    sample_size = k
    for(l in 1){
      p.ttest = rep(NA, nrows)
      # repeat data simulation and analysis to get power estimates
      for(i in 1:100){
        simcontrol = rnorm(sample_size, mean = control_mean_response, sd = control_sd)
        simtreatment = rnorm(sample_size, mean = treatment_mean_response, sd = treatment_sd)
        test = t.test(simcontrol, simtreatment)
        p.ttest[i] = test$p.value
      }
      Neonic_power$sample_size[m] = k
      Neonic_power$power_t_test[m] = mean(p.ttest <= alpha, na.rm = TRUE)
      Neonic_power$endpoint[m] = name
      m = m + 1
    }
  }
}


Neonic_power = data.frame(Neonic_power %>%
                         group_by(endpoint, sample_size) %>% 
                         summarize(mean_power = mean(power_t_test), sd_power = sd(power_t_test)))

over_80_percent_power = subset(Neonic_power, mean_power >= 0.8)
over_80_percent_power = data.frame(over_80_percent_power %>% group_by(endpoint) %>% slice_min(sample_size, n = 1)) 
power_at_n3 = data.frame(Neonic_power %>% group_by(endpoint) %>% filter(sample_size %in% 3))
power_at_n8 = data.frame(Neonic_power %>% group_by(endpoint) %>% filter(sample_size %in% 8))
Neonic_endpoints_power = rbind(power_at_n3, power_at_n8, over_80_percent_power)

ros_power = subset(Neonic_power, endpoint == "ROS")
hydroxyl_power = subset(Neonic_power, endpoint == "OH")
cat_power = subset(Neonic_power, endpoint == "CAT")
sod_power = subset(Neonic_power, endpoint == "SOD")
gst_power = subset(Neonic_power, endpoint == "GST")
pod_power = subset(Neonic_power, endpoint == "POD")
care_power = subset(Neonic_power, endpoint == "CarE")
gr_power = subset(Neonic_power, endpoint == "GR") 
sod_mRNA_power = subset(Neonic_power, endpoint == "SOD mRNA expression")
cat_mRNA_power = subset(Neonic_power, endpoint == "CAT mRNA expression")
gst_mRNA_power = subset(Neonic_power, endpoint == "GST mRNA expression")
gsh_power = subset(Neonic_power, endpoint == "GSH")
mda_power = subset(Neonic_power, endpoint == "MDA")
OhdG_power = subset(Neonic_power, endpoint == "8-OHdG")
pc_power = subset(Neonic_power, endpoint == "PC")

ROS_plot = ggplot(data = ros_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "ROS concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 11, linetype = "dashed", color = "red")
Hydroxyl_plot = ggplot(data = hydroxyl_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "OH concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 22, linetype = "dashed", color = "red")
CAT_plot = ggplot(data = cat_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "CAT activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 14, linetype = "dashed", color = "red")
SOD_plot = ggplot(data = sod_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "SOD activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 12, linetype = "dashed", color = "red")
GST_plot = ggplot(data = gst_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "GST activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 24, linetype = "dashed", color = "red")
POD_plot = ggplot(data = pod_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "POD activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 23, linetype = "dashed", color = "red")
CarE_plot = ggplot(data = care_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "CarE activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
GR_plot = ggplot(data = gr_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "GR activity") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
SOD_mRNA_plot = ggplot(data = sod_mRNA_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "SOD mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 10, linetype = "dashed", color = "red")
CAT_mRNA_plot = ggplot(data = cat_mRNA_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "CAT mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 17, linetype = "dashed", color = "red")
GST_mRNA_plot = ggplot(data = gst_mRNA_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "GST mRNA expression") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 13, linetype = "dashed", color = "red")
GSH_plot = ggplot(data = gsh_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "GSH concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) 
MDA_plot = ggplot(data = mda_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "MDA concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 13, linetype = "dashed", color = "red")
OhdG_plot = ggplot(data = OhdG_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "8-OhdG concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 8, linetype = "dashed", color = "red")
PC_plot = ggplot(data = pc_power, aes(x = sample_size, y = mean_power)) + geom_point(shape = 1, size = 2) + #geom_errorbar(aes(ymin = mean_power-sd_power, ymax = mean_power+sd_power), width = 0.25) +
  labs(x = "sample size", y = "statistical power", title = "PC concentration") + theme_classic() + scale_x_continuous(breaks = c(3,seq(5, 50, 5)), labels = c(3,seq(5, 50, 5)), expand = c(0,0), limits = c(1,52)) +
  theme(axis.text = element_text(size = 10), axis.title = element_text(size = 11), title = element_text(size = 14)) + geom_vline(xintercept = 25, linetype = "dashed", color = "red")

ggarrange(ROS_plot, Hydroxyl_plot, CAT_plot,
          SOD_plot, GST_plot, POD_plot, 
          CarE_plot, GR_plot, SOD_mRNA_plot,
          CAT_mRNA_plot, GST_mRNA_plot, GSH_plot,
          MDA_plot, OhdG_plot, PC_plot,
          nrow = 5,
          ncol = 3,
          legend = "none")
