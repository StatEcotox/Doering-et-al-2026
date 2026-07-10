
# loading packages ----

library("ggplot2")
library("metafor")
library("orchaRd")
library("dplyr")


# load Neonic data ----

Neonics_contact_test = read.csv2("Data/Neonics_contact_test.csv", stringsAsFactors = T)

# transform to correct data type ----

Neonics_contact_test$c_unit = as.character(Neonics_contact_test$c_unit)
Neonics_contact_test$t_unit = as.character(Neonics_contact_test$t_unit)
Neonics_contact_test$c_n_replicates = as.integer(Neonics_contact_test$c_n_replicates)
Neonics_contact_test$t_n_replicates = as.integer(Neonics_contact_test$t_n_replicates)
Neonics_contact_test$c_mean = as.numeric(Neonics_contact_test$c_mean)
Neonics_contact_test$t_mean = as.numeric(Neonics_contact_test$t_mean)
Neonics_contact_test$c_SD = as.numeric(Neonics_contact_test$c_SD)
Neonics_contact_test$t_SD = as.numeric(Neonics_contact_test$t_SD)
Neonics_contact_test$c_SE = as.numeric(Neonics_contact_test$c_SE)
Neonics_contact_test$t_SE = as.numeric(Neonics_contact_test$t_SE)
(table(Neonics_contact_test$t_unit, Neonics_contact_test$endpoint))



# unify units ----

# Transform "µmol/g_protein" in "nmol/mg_protein"
Neonics_contact_test[Neonics_contact_test == "µmol/g_protein"] = "nmol/mg_protein"

# Transform "nmol/min/mg_protein" in "U/mg_protein"
temporary_data = subset(Neonics_contact_test, t_unit == "nmol/min/mg_protein")
temporary_data2 = subset(Neonics_contact_test, t_unit != "nmol/min/mg_protein")
temporary_data$c_mean = temporary_data$c_mean/1000
temporary_data$c_SD = temporary_data$c_SD/1000
temporary_data$c_SE = temporary_data$c_SE/1000
temporary_data$t_mean = temporary_data$t_mean/1000
temporary_data$t_SD = temporary_data$t_SD/1000
temporary_data$t_SE = temporary_data$t_SE/1000
temporary_data[temporary_data == "nmol/min/mg_protein"] = "U/mg_protein"
Neonics_contact_test = rbind(temporary_data, temporary_data2)

# Transform "U/g_protein" in "U/mg_protein"
temporary_data = subset(Neonics_contact_test, t_unit == "U/g_protein")
temporary_data2 = subset(Neonics_contact_test, t_unit != "U/g_protein")
temporary_data$c_mean = temporary_data$c_mean/1000
temporary_data$c_SD = temporary_data$c_SD/1000
temporary_data$c_SE = temporary_data$c_SE/1000
temporary_data$t_mean = temporary_data$t_mean/1000
temporary_data$t_SD = temporary_data$t_SD/1000
temporary_data$t_SE = temporary_data$t_SE/1000
temporary_data[temporary_data == "U/g_protein"] = "U/mg_protein"
Neonics_contact_test = rbind(temporary_data, temporary_data2)

# Transform "Lipid Peroxidation" in "Malondialdehyde Concentration"
Neonics_contact_test$endpoint = as.character(Neonics_contact_test$endpoint)
Neonics_contact_test[Neonics_contact_test == "Lipid Peroxidation"] = "Malondialdehyde Concentration"

# Change labels from TOMEX table

Neonics_contact_test[Neonics_contact_test == "Reactive Oxygen Species Production"] = "ROS"
Neonics_contact_test[Neonics_contact_test == "Catalase Activity"] = "CAT"
Neonics_contact_test[Neonics_contact_test == "Glutathione S-Transferase Activity"] = "GST"
Neonics_contact_test[Neonics_contact_test == "Carboxylesterase activity"] = "CarE"
Neonics_contact_test[Neonics_contact_test == "Glutathione Reductase Activity"] = "GR"
Neonics_contact_test[Neonics_contact_test == "Glutathione Concentration"] = "GSH"
Neonics_contact_test[Neonics_contact_test == "Malondialdehyde Concentration"] = "MDA"

Neonics_contact_test$endpoint = as.factor(Neonics_contact_test$endpoint)

# Delete temporary data
rm(temporary_data, temporary_data2)

# Overview
table(Neonics_contact_test$t_unit,Neonics_contact_test$endpoint)


# transform to correct data type ----

log_means = escalc(measure = "ROM", m1i = c_mean, m2i = t_mean, 
                   sd1i = c_SD, sd2i = t_SD,
                   n1i = c_n_replicates, n2i = t_n_replicates,
                   data = Neonics_contact_test)

# Positive yi value means that the endpoint is higher in treatment compared to control 
log_means$yi = log_means$yi * (-1) 

log_means$endpoint = factor(log_means$endpoint , levels = c("MDA",
                                                            "GSH",
                                                            "GR",
                                                            "CarE",
                                                            "GST",
                                                            "CAT",
                                                            "ROS"))


# analysis with all endpoints ----
individual_level = 1:nrow(log_means)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = log_means)
res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
             angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = "left", legend.pos = "none") +
  scale_fill_manual(values = c("gold3","slateblue3","purple3","firebrick4","mediumpurple3","deepskyblue3","gray47")) +
  scale_colour_manual(values = c("gold3","slateblue3","purple3","firebrick4","mediumpurple3","deepskyblue3","gray47"))

data.frame(log_means %>% 
             group_by(endpoint) %>% 
             summarize(mean_sd = sd(yi)))
