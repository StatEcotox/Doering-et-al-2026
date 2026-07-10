
# loading packages ----

library("metafor")


# load Neonic data ----

new_Neo_data = read.csv2("Data/new_Neonics_data.csv", stringsAsFactors = T, fileEncoding = "latin1")


# calculate log transformed ratio of means ----

# Calculate yi (effect sizes) and vi (sampling variances)
log_means = escalc(measure = "ROM", m1i = c_mean, m2i = t_mean, 
                   sd1i = c_SD, sd2i = t_SD,
                   n1i = c_n_replicates, n2i = t_n_replicates,
                   data = new_Neo_data)

# Exclude 1 data point with vi = 0
# log_means['vi'][log_means['vi'] == 0] = 0.0000001
log_means = subset(log_means, vi != 0)

# Positive yi value means that the endpoint is higher in treatment compared to control 
log_means$yi = log_means$yi * (-1) 


# export effect sizes ----

write.csv2(log_means, "Data/Neonics_log_means.csv", row.names=FALSE)








