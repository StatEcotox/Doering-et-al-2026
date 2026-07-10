
# loading packages ----

library("metafor")


# load MP data ----

new_NMP_data = read.csv2("Data/new_NMP_data.csv", stringsAsFactors = T, fileEncoding = "latin1")


# calculate log transformed ratio of means ----

# calculate yi (effect sizes) and vi (sampling variances)
log_means = escalc(measure = "ROM", m1i = c_mean, m2i = t_mean, 
                   sd1i = c_SD, sd2i = t_SD,
                   n1i = c_n_replicates, n2i = t_n_replicates,
                   data = new_NMP_data)

# Exclude 20 data points with vi = 0 and 9 data points with "NA"
# log_means['vi'][log_means['vi'] == 0] = 0.0000001
log_means = subset(log_means, vi != 0) 

# positive yi value means that the endpoint is higher in treatment compared to control 
log_means$yi = log_means$yi * (-1) 


# export effect sizes ----

write.csv2(log_means, "Data/NMP_log_means.csv", row.names = FALSE)




