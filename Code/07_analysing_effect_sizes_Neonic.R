
# loading packages ----

library("metafor")
library("orchaRd")
library("ggplot2")
library("dplyr")


# load MP data ----

log_means = read.csv2("Data/Neonics_log_means.csv", stringsAsFactors = T, fileEncoding="latin1")

log_means$endpoint = factor(log_means$endpoint , levels = c("PC",
                                                            "8-OHdG",
                                                            "MDA",
                                                            "GSH",
                                                            "GST mRNA expression",
                                                            "CAT mRNA expression", 
                                                            "SOD mRNA expression", 
                                                            "GR",
                                                            "CarE",
                                                            "POD",
                                                            "GST",
                                                            "SOD",
                                                            "CAT",
                                                            "OH",
                                                            "ROS"))

# subset data in all endpoints ----

ros_content = subset(log_means, endpoint == "ROS")
hydroxyl_levels = subset(log_means, endpoint == "OH")
cat_activity = subset(log_means, endpoint == "CAT")
sod_activity = subset(log_means, endpoint == "SOD")
gst_activity = subset(log_means, endpoint == "GST")
pod_activity = subset(log_means, endpoint == "POD")
care_activity = subset(log_means, endpoint == "CarE")
gr_activity = subset(log_means, endpoint == "GR") 
sod_mRNA_expression= subset(log_means, endpoint == "SOD mRNA expression")
cat_mRNA_expression = subset(log_means, endpoint == "CAT mRNA expression")
gst_mRNA_expression = subset(log_means, endpoint == "GST mRNA expression")
gsh_content = subset(log_means, endpoint == "GSH")
mda_content = subset(log_means, endpoint == "MDA")
OhdG_content = subset(log_means, endpoint == "8-OHdG")
pc_concentration = subset(log_means, endpoint == "PC")


# analysis with all endpoints ----
individual_level = 1:nrow(log_means)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = log_means)
res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
             angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = -2.7, legend.pos = "none") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1), labels = c(-3,-2,-1,0,1), expand = c(0,0), limits = c(-3.13,1.8)) +
  scale_fill_manual(values = c("lavenderblush4","wheat4","gold3","slateblue3","mediumpurple4","deepskyblue4","orange4","purple3","firebrick4","cyan4","mediumpurple3","orange3","deepskyblue3","khaki4","gray47")) +
  scale_colour_manual(values = c("lavenderblush4","wheat4","gold3","slateblue3","mediumpurple4","deepskyblue4","orange4","purple3","firebrick4","cyan4","mediumpurple3","orange3","deepskyblue3","khaki4","gray47"))


data.frame(log_means %>% 
             group_by(endpoint) %>% 
             summarize(mean_sd = sd(yi)))


# calculate within-particleID variances ----

log_means_comparisons = subset(log_means, !is.na(neonic) &
                                 !is.na(combined_concentration) &
                                 !is.na(exposure_duration_effectsize_days))

# determining which predictor is needed for model ----
individual_level = 1:nrow(log_means_comparisons)
res_all_predictors = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+neonic+log_combined_concentration+exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")

res_all_without_neonics = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+log_combined_concentration+exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")
anova(res_all_predictors, res_all_without_neonics) # full model not significantly better

res_all_without_concentration = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+neonic+exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")
anova(res_all_predictors, res_all_without_concentration) # full model significantly better

res_all_without_time = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+neonic+log_combined_concentration-1, data = log_means_comparisons, method = "ML")
anova(res_all_predictors, res_all_without_time) # full model not significantly better


# analysis of interactions ----

log_means_comparisons = subset(log_means, !is.na(neonic) &
                                 !is.na(combined_concentration) &
                                 !is.na(exposure_duration_effectsize_days))

individual_level = 1:nrow(log_means_comparisons)
res_all_interaction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + neonic + log_combined_concentration + exposure_duration_effectsize_days + endpoint:neonic + endpoint:log_combined_concentration + endpoint:exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")

res_all_noNeonicInteraction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + neonic + log_combined_concentration + exposure_duration_effectsize_days + endpoint:log_combined_concentration + endpoint:exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")
anova(res_all_interaction, res_all_noNeonicInteraction) 

res_all_noConcentrationInteraction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + neonic + log_combined_concentration + exposure_duration_effectsize_days + endpoint:neonic + endpoint:exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")
anova(res_all_interaction, res_all_noConcentrationInteraction) 

res_all_noTimeInteraction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + neonic + log_combined_concentration + exposure_duration_effectsize_days + endpoint:neonic + endpoint:log_combined_concentration-1, data = log_means_comparisons, method = "ML")
anova(res_all_interaction, res_all_noTimeInteraction) 


# how much residual heterogeneity is explained by the inclusion of experimental design and pollutant properties ----
log_means_comparisons = subset(log_means, !is.na(neonic) &
                                 !is.na(combined_concentration) &
                                 !is.na(exposure_duration_effectsize_days))

individual_level = 1:nrow(log_means_comparisons)
res_only_endpoint = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint, data = log_means_comparisons)
res_all_predictors = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + neonic + log_combined_concentration + exposure_duration_effectsize_days + endpoint:neonic + endpoint:log_combined_concentration + endpoint:exposure_duration_effectsize_days, data = log_means_comparisons)

(sum(res_only_endpoint$sigma2) - sum(res_all_predictors$sigma2)) / sum(res_only_endpoint$sigma2)*100
# REML:21.18939 %
# ML:37.8704 %

W = diag(1/res_only_endpoint$vi)
X = model.matrix(res_only_endpoint)
P = W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_only_endpoint$sigma2) / (sum(res_only_endpoint$sigma2) + (res_only_endpoint$k-res_only_endpoint$p)/sum(diag(P)))
# 99.6518 % of the total unaccounted variance comprises residual heterogeneity

W = diag(1/res_all_predictors$vi)
X = model.matrix(res_all_predictors)
P = W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_all_predictors$sigma2) / (sum(res_all_predictors$sigma2) + (res_all_predictors$k-res_all_predictors$p)/sum(diag(P)))
# 99.4876 % of the total unaccounted variance comprises residual heterogeneity




# supplement ----

## analysis of endpoints over time ----

layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10, 11, 12,
                13, 14, 15), nrow = 5, byrow = TRUE))

# ros content
individual_level = 1:nrow(ros_content)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = ros_content)
res_time
xs = seq(-1,44,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# OH content
individual_level = 1:nrow(hydroxyl_levels)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = hydroxyl_levels)
res_time
xs = seq(-1,58,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# cat activity
individual_level = 1:nrow(cat_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = cat_activity)
res_time
xs = seq(-2,58,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# sod activity
individual_level = 1:nrow(sod_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = sod_activity)
res_time
xs = seq(-2,59,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# gst activity
individual_level = 1:nrow(gst_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = gst_activity)
res_time
xs = seq(-3,59,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# pod activity
individual_level = 1:nrow(pod_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = pod_activity)
res_time
xs = seq(-1,29,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# care activity
individual_level = 1:nrow(care_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = care_activity)
res_time
xs = seq(-1,29,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# gr activity
individual_level = 1:nrow(gr_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = gr_activity)
res_time
xs = seq(-1,30,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# sod mRNA expression
individual_level = 1:nrow(sod_mRNA_expression)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = sod_mRNA_expression)
res_time
xs = seq(-1,44,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# cat mRNA expression
individual_level = 1:nrow(cat_mRNA_expression)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = cat_mRNA_expression)
res_time
xs = seq(-1,44,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# gst mRNA expression
individual_level = 1:nrow(gst_mRNA_expression)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = gst_mRNA_expression)
res_time
xs = seq(1,22,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# gsh concentration
individual_level = 1:nrow(gsh_content)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = gsh_content)
res_time
xs = seq(-1,47,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# mda content
individual_level = 1:nrow(mda_content)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = mda_content)
res_time
xs = seq(-2,59,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# 8-OhdG content
individual_level = 1:nrow(OhdG_content)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = OhdG_content)
res_time
xs = seq(-1,58,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# protein carbonyl concentration
individual_level = 1:nrow(pc_concentration)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = pc_concentration)
res_time
xs = seq(-1,44,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)


## analysis of endpoints at different neonic concentrations ----

layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10, 11, 12,
                13, 14, 15), nrow = 5, byrow = TRUE))
# ros content
individual_level = 1:nrow(ros_content)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = ros_content)
res_conc
xs = seq(-14,-3,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01)), labels = c(0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# OH content
individual_level = 1:nrow(hydroxyl_levels)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = hydroxyl_levels)
res_conc
xs = seq(-12,-9, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# cat activity
individual_level = 1:nrow(cat_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = cat_activity)
res_conc
xs = seq(-17,-4,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# sod activity
individual_level = 1:nrow(sod_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = sod_activity)
res_conc
xs = seq(-14,-4, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# gst activity
individual_level = 1:nrow(gst_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = gst_activity)
res_conc
xs = seq(-17,-3, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# pod activity
individual_level = 1:nrow(pod_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = pod_activity)
res_conc
xs = seq(-12,-5, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# care activity
individual_level = 1:nrow(care_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = care_activity)
res_conc
xs = seq(-15,-3, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# gr activity
individual_level = 1:nrow(gr_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = gr_activity)
res_conc
xs = seq(-12,-3, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# sod mRNA expression
individual_level = 1:nrow(sod_mRNA_expression)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = sod_mRNA_expression)
res_conc
xs = seq(-12,-8, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# cat mRNA expression
individual_level = 1:nrow(cat_mRNA_expression)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = cat_mRNA_expression)
res_conc
xs = seq(-12,-8, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# gst mRNA expression
individual_level = 1:nrow(gst_mRNA_expression)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = gst_mRNA_expression)
res_conc
xs = seq(-14,-10, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# gsh concentration
individual_level = 1:nrow(gsh_content)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = gsh_content)
res_conc
xs = seq(-17,-3, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# mda content
individual_level = 1:nrow(mda_content)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = mda_content)
res_conc
xs = seq(-14,-5, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# 8-OhdG content
individual_level = 1:nrow(OhdG_content)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = OhdG_content)
res_conc
xs = seq(-12,-9, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)

# protein carbonyl concentration
individual_level = 1:nrow(pc_concentration)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = pc_concentration)
res_conc
xs = seq(-12,-8, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "neonic concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01)), labels = c(0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01), cex.axis = 1.5)


## analysis of neonic effects ----

# Only Clothianidin and Imidacloprid

data.frame(sort(table(log_means$neonic), decreasing = TRUE))

clothianidin = subset(log_means, neonic == "Clothianidin")
individual_level = 1:nrow(clothianidin)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = clothianidin)
res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
             angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = -1.5, legend.pos = "none") +
  scale_y_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), labels = c(-1.5,-1,-0.5,0,0.5,1,1.5), expand = c(0,0), limits = c(-1.77,1.55)) +
  scale_fill_manual(values = c("lavenderblush4","gold3","cyan4","orange3","mediumpurple3","deepskyblue3","gray47")) +
  scale_colour_manual(values = c("lavenderblush4","gold3","cyan4","orange3","mediumpurple3","deepskyblue3","gray47"))


imidacloprid = subset(log_means, neonic == "Imidacloprid")
imidacloprid$endpoint = as.character(imidacloprid$endpoint)
imidacloprid[imidacloprid == "SOD mRNA expression"] = "SOD1"
imidacloprid[imidacloprid == "CAT mRNA expression"] = "CAT1"
imidacloprid[imidacloprid == "GST mRNA expression"] = "GST1"
imidacloprid$endpoint = as.factor(imidacloprid$endpoint)
imidacloprid$endpoint = factor(imidacloprid$endpoint , levels = c("MDA",
                                                                  "GST1",
                                                                  "CAT1",
                                                                  "SOD1",
                                                                  "CarE",
                                                                  "POD",
                                                                  "GST",
                                                                  "SOD",
                                                                  "CAT",
                                                                  "ROS"))

individual_level = 1:nrow(imidacloprid)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = imidacloprid)
res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
             angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = -1.5, legend.pos = "none") +
  scale_y_continuous(breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5), labels = c(-1.5,-1,-0.5,0,0.5,1,1.5), expand = c(0,0), limits = c(-1.77,1.55)) +
  scale_fill_manual(values = c("gold3","mediumpurple4","deepskyblue4","orange4","firebrick4","cyan4","mediumpurple3","orange3","deepskyblue3","gray47")) +
  scale_colour_manual(values = c("gold3","mediumpurple4","deepskyblue4","orange4","firebrick4","cyan4","mediumpurple3","orange3","deepskyblue3","gray47"))
