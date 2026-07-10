
# loading packages ----

library("metafor")
library("orchaRd")
library("ggplot2")
library("dplyr")


# load NMP data ----

log_means = read.csv2("Data/NMP_log_means.csv", stringsAsFactors = T, fileEncoding = "latin1")

log_means$endpoint = factor(log_means$endpoint , levels = c("8-OHdG",
                                                            "MDA",
                                                            "GSH",
                                                            "TAC",
                                                            "GST mRNA expression",
                                                            "CAT mRNA expression",
                                                            "SOD mRNA expression",
                                                            "POD",
                                                            "GST",
                                                            "SOD",
                                                            "CAT",
                                                            "ROS"))


# subset data in all endpoints (ROS markers) ----

ros_content = subset(log_means, endpoint == "ROS")
cat_activity = subset(log_means, endpoint == "CAT")
sod_activity = subset(log_means, endpoint == "SOD")
gst_activity = subset(log_means, endpoint == "GST")
pod_activity = subset(log_means, endpoint == "POD")
sod_mRNA_expression= subset(log_means, endpoint == "SOD mRNA expression")
cat_mRNA_expression = subset(log_means, endpoint == "CAT mRNA expression")
gst_mRNA_expression = subset(log_means, endpoint == "GST mRNA expression")
TAC_activity = subset(log_means, endpoint == "TAC")
gsh_content = subset(log_means, endpoint == "GSH")
mda_content = subset(log_means, endpoint == "MDA")
OhdG_content = subset(log_means, endpoint == "8-OHdG")



# analysis with all endpoints ----

individual_level = 1:nrow(log_means)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = log_means)
res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
             angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = -2.6, legend.pos = "none") +
  scale_y_continuous(breaks = c(-3,-2,-1,0,1), labels = c(-3,-2,-1,0,1), expand = c(0,0), limits = c(-3.13,1.8)) +
  scale_fill_manual(values = c("wheat4","gold3","slateblue3","coral4","mediumpurple4","deepskyblue4","orange4","cyan4","mediumpurple3","orange3","deepskyblue3","gray47")) +
  scale_colour_manual(values = c("wheat4","gold3","slateblue3","coral4","mediumpurple4","deepskyblue4","orange4","cyan4","mediumpurple3","orange3","deepskyblue3","gray47")) 

data.frame(log_means %>% 
             group_by(endpoint) %>% 
             summarize(mean_sd = sd(yi)))


# calculate within-particleID variances ----

log_means_comparisons = subset(log_means, !is.na(combined_concentration) &
                                 !is.na(shape) &
                                 !is.na(polymer_type) &
                                 !is.na(combined_size_mm) &
                                 !is.na(exposure_duration_effectsize_days))


# determining which predictor is needed for model ----

individual_level = 1:nrow(log_means_comparisons)
res_all_predictors = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+log_combined_concentration+shape+polymer_type+combined_size_mm+exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")

res_all_without_time = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+log_combined_concentration+shape+polymer_type+combined_size_mm-1, data = log_means_comparisons, method = "ML")
anova(res_all_predictors, res_all_without_time) 

res_all_without_concentration = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+shape+polymer_type+combined_size_mm+exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")
anova(res_all_predictors, res_all_without_concentration) 

res_all_without_shape = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+log_combined_concentration+polymer_type+combined_size_mm+exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")
anova(res_all_predictors, res_all_without_shape) 

res_all_without_polymer = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+log_combined_concentration+shape+combined_size_mm+exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")
anova(res_all_predictors, res_all_without_polymer) 

res_all_without_size = rma.mv(yi, vi, random = ~1|DOI/individual_level,mods = ~endpoint+log_combined_concentration+shape+polymer_type+exposure_duration_effectsize_days-1, data = log_means_comparisons, method = "ML")
anova(res_all_predictors, res_all_without_size)


# analysis of interactions ----

log_means_comparisons = subset(log_means, !is.na(shape) &
                                 !is.na(polymer_type) &
                                 !is.na(combined_size_mm))

individual_level = 1:nrow(log_means_comparisons)
res_all_interaction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + log_combined_concentration + exposure_duration_effectsize_days + polymer_type + shape + combined_size_mm + endpoint:log_combined_concentration + endpoint:exposure_duration_effectsize_days + endpoint:polymer_type + endpoint:shape + endpoint:combined_size_mm-1, data = log_means_comparisons, method = "ML")

res_all_noConcentrationInteraction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + log_combined_concentration + exposure_duration_effectsize_days + polymer_type + shape + combined_size_mm + endpoint:exposure_duration_effectsize_days + endpoint:polymer_type + endpoint:shape + endpoint:combined_size_mm-1, data = log_means_comparisons, method = "ML")
anova(res_all_interaction, res_all_noConcentrationInteraction) 

res_all_noTimeInteraction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + log_combined_concentration + exposure_duration_effectsize_days + polymer_type + shape + combined_size_mm + endpoint:log_combined_concentration + endpoint:polymer_type + endpoint:shape + endpoint:combined_size_mm-1, data = log_means_comparisons, method = "ML")
anova(res_all_interaction, res_all_noTimeInteraction) 

res_all_noPolymerInteraction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + log_combined_concentration + exposure_duration_effectsize_days + polymer_type + shape + combined_size_mm + endpoint:log_combined_concentration + endpoint:exposure_duration_effectsize_days + endpoint:shape + endpoint:combined_size_mm-1, data = log_means_comparisons, method = "ML")
anova(res_all_interaction, res_all_noPolymerInteraction) 

res_all_noShapeInteraction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + log_combined_concentration + exposure_duration_effectsize_days + polymer_type + shape + combined_size_mm + endpoint:log_combined_concentration + endpoint:exposure_duration_effectsize_days + endpoint:polymer_type + endpoint:combined_size_mm-1, data = log_means_comparisons, method = "ML")
anova(res_all_interaction, res_all_noShapeInteraction)

res_all_noSizeInteraction = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + log_combined_concentration + exposure_duration_effectsize_days + polymer_type + shape + combined_size_mm + endpoint:log_combined_concentration + endpoint:exposure_duration_effectsize_days + endpoint:polymer_type + endpoint:shape-1, data = log_means_comparisons, method = "ML")
anova(res_all_interaction, res_all_noSizeInteraction)


# how much residual heterogeneity is explained by the inclusion of experimental design and pollutant properties ----

log_means_comparisons = subset(log_means, !is.na(shape) &
                                 !is.na(polymer_type) &
                                 !is.na(combined_size_mm))

individual_level = 1:nrow(log_means_comparisons)
res_only_endpoint = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = log_means_comparisons)
res_all_predictors = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint + log_combined_concentration + exposure_duration_effectsize_days + polymer_type + shape + combined_size_mm + endpoint:log_combined_concentration + endpoint:exposure_duration_effectsize_days + endpoint:polymer_type + endpoint:shape + endpoint:combined_size_mm-1, data = log_means_comparisons)

(sum(res_only_endpoint$sigma2) - sum(res_all_predictors$sigma2)) / sum(res_only_endpoint$sigma2)*100
# REML: 46.71053 %
# ML: 46.85243 %

W = diag(1/res_only_endpoint$vi)
X = model.matrix(res_only_endpoint)
P = W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_only_endpoint$sigma2) / (sum(res_only_endpoint$sigma2) + (res_only_endpoint$k-res_only_endpoint$p)/sum(diag(P)))
# 99.91417 % of the total unaccounted variance comprises residual heterogeneity

W = diag(1/res_all_predictors$vi)
X = model.matrix(res_all_predictors)
P = W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(res_all_predictors$sigma2) / (sum(res_all_predictors$sigma2) + (res_all_predictors$k-res_all_predictors$p)/sum(diag(P)))
# 99.80663 % of the total unaccounted variance comprises residual heterogeneity




# supplement ----

## analysis of endpoints over time ----

layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10, 11, 12), nrow = 4, byrow = TRUE))

# ros content
individual_level = 1:nrow(ros_content)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = ros_content)
res_time
xs = seq(5,44,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# cat activity
individual_level = 1:nrow(cat_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = cat_activity)
res_time
xs = seq(5,58,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# sod activity
individual_level = 1:nrow(sod_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = sod_activity)
res_time
xs = seq(5,58,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# gst activity
individual_level = 1:nrow(gst_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = gst_activity)
res_time
xs = seq(5,30,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# pod activity
individual_level = 1:nrow(pod_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = pod_activity)
res_time
xs = seq(5,30,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# sod mRNA expression activity
individual_level = 1:nrow(sod_mRNA_expression)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = sod_mRNA_expression)
res_time
xs = seq(5,44,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# cat mRNA expression
individual_level = 1:nrow(cat_mRNA_expression)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = cat_mRNA_expression)
res_time
xs = seq(5,30,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# gst mRNA expression
individual_level = 1:nrow(gst_mRNA_expression)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = gst_mRNA_expression)
res_time
xs = seq(12,30,length = 1000)
tmp = predict(res_time, newmods = cbind(xs))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# TAC activity
individual_level = 1:nrow(TAC_activity)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = TAC_activity)
res_time
xs = seq(5,30,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# gsh content
individual_level = 1:nrow(gsh_content)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = gsh_content)
res_time
xs = seq(5,30,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# mda content
individual_level = 1:nrow(mda_content)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = mda_content)
res_time
xs = seq(5,44,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)

# 8-OHdG content
individual_level = 1:nrow(OhdG_content)
res_time = rma.mv(yi, vi, mods = ~exposure_duration_effectsize_days + I(exposure_duration_effectsize_days^2), random = ~1|DOI/individual_level, data = OhdG_content)
res_time
xs = seq(5,30,length = 1000)
tmp = predict(res_time, newmods = cbind(xs, xs^2))
regplot(res_time, mod = "exposure_duration_effectsize_days", pred = tmp, xlab = "Days", ylab = "logROM", xvals = xs, cex.lab = 1.5, cex.axis = 1.5)


## analysis of endpoints at different MP concentrations ----

layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10, 11, 12), nrow = 4, byrow = TRUE))

# ros content
individual_level = 1:nrow(ros_content)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = ros_content)
res_conc
xs = seq(-6,3,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# cat activity
individual_level = 1:nrow(cat_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = cat_activity)
res_conc
xs = seq(-6,4, length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# sod activity
individual_level = 1:nrow(sod_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = sod_activity)
res_conc
xs = seq(-13,4,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# gst activity
individual_level = 1:nrow(gst_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = gst_activity)
res_conc
xs = seq(-6,4,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# pod activity
individual_level = 1:nrow(pod_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = pod_activity)
res_conc
xs = seq(-5,4,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# sod mRNA expression
individual_level = 1:nrow(sod_mRNA_expression)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = sod_mRNA_expression)
res_conc
xs = seq(-4,1,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# cat mRNA expression
individual_level = 1:nrow(cat_mRNA_expression)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = cat_mRNA_expression)
res_conc
xs = seq(-4,1,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# gst mRNA expression
individual_level = 1:nrow(gst_mRNA_expression)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = gst_mRNA_expression)
res_conc
xs = seq(-3.1,-1.5,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.05, 0.1, 0.2)), labels = c(0.05, 0.1, 0.2), cex.axis = 1.5)

# tac
individual_level = 1:nrow(TAC_activity)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = TAC_activity)
res_conc
xs = seq(-4,-1,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.05, 0.1, 0.2)), labels = c(0.05, 0.1, 0.2), cex.axis = 1.5)

# gsh content
individual_level = 1:nrow(gsh_content)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = gsh_content)
res_conc
xs = seq(-13,3,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# mda content
individual_level = 1:nrow(mda_content)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = mda_content)
res_conc
xs = seq(-6,4,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)

# 8-OHdG content
individual_level = 1:nrow(OhdG_content)
res_conc = rma.mv(yi, vi, mods = ~log_combined_concentration + I(log_combined_concentration^2), random = ~1|DOI/individual_level, data = OhdG_content)
res_conc
xs = seq(-4,1,length = 1000)
tmp = predict(res_conc, newmods = cbind(xs, xs^2))
regplot(res_conc, mod = "log_combined_concentration", pred = tmp, xvals = xs, xaxt = "n", xlab = "NMP concentration (% w/w)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)
axis(side = 1, at = log(c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10)), labels = c(0.00001, 0.0001, 0.001, 0.01, 0.1, 1, 10), cex.axis = 1.5)


## analysis of endpoints at different MP sizes ----

layout(matrix(c(1, 2, 3,
                4, 5, 6,
                7, 8, 9,
                10, 11, 12), nrow = 4, byrow = TRUE))

# ros content
individual_level = 1:nrow(ros_content)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = ros_content)
res_size
xs = seq(0,23,length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# cat activity
individual_level = 1:nrow(cat_activity)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = cat_activity)
res_size
xs = seq(0,22, length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# sod activity
individual_level = 1:nrow(sod_activity)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = sod_activity)
res_size
xs = seq(0,22, length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# gst activity
individual_level = 1:nrow(gst_activity)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = gst_activity)
res_size
xs = seq(0,21,length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# pod activity
individual_level = 1:nrow(pod_activity)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = pod_activity)
res_size
xs = seq(0,19,length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# sod mRNA expression
individual_level = 1:nrow(sod_mRNA_expression)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = sod_mRNA_expression)
res_size
xs = seq(2,14,length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# cat mRNA expression
individual_level = 1:nrow(cat_mRNA_expression)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = cat_mRNA_expression)
res_size
xs = seq(2,14,length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# gst mRNA expression
# individual_level = 1:nrow(gst_mRNA_expression)
# res_conc = rma.mv(yi, vi, mods = ~combined_size_mm, random = ~1|DOI/individual_level, data = gst_mRNA_expression)
# res_conc
# regplot(res_conc, xlab = "MP size (mm)", ylab = "log transformed ratio of means")

# TAC activity
# individual_level = 1:nrow(TAC_activity)
# res_conc = rma.mv(yi, vi, mods = ~combined_size_mm, random = ~1|DOI/individual_level, data = TAC_activity)
# res_conc
# regplot(res_conc, xlab = "MP size (mm)", ylab = "log transformed ratio of means")

# gsh content
individual_level = 1:nrow(gsh_content)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = gsh_content)
res_size
xs = seq(1,23,length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# mda content
individual_level = 1:nrow(mda_content)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = mda_content)
res_size
xs = seq(0,24,length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)

# 8-OHdG content
individual_level = 1:nrow(OhdG_content)
res_size = rma.mv(yi, vi, mods = ~combined_size_mm + I(combined_size_mm^2), random = ~1|DOI/individual_level, data = OhdG_content)
res_size
xs = seq(2,21,length = 1000)
tmp = predict(res_size, newmods = cbind(xs, xs^2))
regplot(res_size, mod = "combined_size_mm", pred = tmp, xvals = xs, xlab = "NMP size (mm)", ylab = "logROM", cex.lab = 1.5, cex.axis = 1.5)


## impact of MP polymer on endpoints ----

log_means_for_plots = log_means
log_means_for_plots$endpoint = as.character(log_means_for_plots$endpoint)
log_means_for_plots[log_means_for_plots == "ROS"] = "a"
log_means_for_plots[log_means_for_plots == "CAT"] = "b"
log_means_for_plots[log_means_for_plots == "SOD"] = "c"
log_means_for_plots[log_means_for_plots == "GST"] = "d"
log_means_for_plots[log_means_for_plots == "POD"] = "e"
log_means_for_plots[log_means_for_plots == "SOD mRNA expression"] = "f"
log_means_for_plots[log_means_for_plots == "CAT mRNA expression"] = "g"
log_means_for_plots[log_means_for_plots == "GST mRNA expression"] = "h"
log_means_for_plots[log_means_for_plots == "TAC"] = "i"
log_means_for_plots[log_means_for_plots == "GSH"] = "j"
log_means_for_plots[log_means_for_plots == "MDA"] = "k"
log_means_for_plots[log_means_for_plots == "8-OHdG"] = "l"

log_means_for_plots$endpoint = as.factor(log_means_for_plots$endpoint)

sort(table(log_means_for_plots$polymer_type))

# ps
ps = subset(log_means_for_plots, polymer_type == "PS")
ps$endpoint = factor(ps$endpoint , levels = c("l",
                                                            "k",
                                                            "j",
                                                            "i",
                                                            "h",
                                                            "g",
                                                            "f",
                                                            "e",
                                                            "d",
                                                            "c",
                                                            "b",
                                                            "a"))

individual_level = 1:nrow(ps)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = ps)

res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
             angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = -1.05, legend.pos = "none") +
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1,1.5), labels = c(-1,-0.5,0,0.5,1,1.5), expand = c(0,0), limits = c(-1.28,1.6)) +
  scale_fill_manual(values = c("wheat4","gold3","slateblue3","coral4","mediumpurple4","deepskyblue4","orange4","cyan4","mediumpurple3","orange3","deepskyblue3","gray47")) +
  scale_colour_manual(values = c("wheat4","gold3","slateblue3","coral4","mediumpurple4","deepskyblue4","orange4","cyan4","mediumpurple3","orange3","deepskyblue3","gray47"))

# ldpe
ldpe = subset(log_means_for_plots, polymer_type == "LDPE")

ldpe$endpoint = factor(ldpe$endpoint , levels = c("l",
                                              "k",
                                              "j",
                                              "i",
                                              "h",
                                              "g",
                                              "f",
                                              "e",
                                              "d",
                                              "c",
                                              "b",
                                              "a"))
individual_level = 1:nrow(ldpe)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = ldpe)
res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
                         angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = -1.05, legend.pos = "none") +
  scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1,1.5), labels = c(-1,-0.5,0,0.5,1,1.5), expand = c(0,0), limits = c(-1.28,1.6)) +
  scale_fill_manual(values = c("wheat4","gold3","cyan4","mediumpurple3","orange3","deepskyblue3","gray47")) +
  scale_colour_manual(values = c("wheat4","gold3","cyan4","mediumpurple3","orange3","deepskyblue3","gray47"))



## impact of MP shape on endpoints ----

# fragment
fragment = subset(log_means_for_plots, shape == "fragment")
fragment$endpoint = factor(fragment$endpoint , levels = c("l",
                                                  "k",
                                                  "j",
                                                  "i",
                                                  "h",
                                                  "g",
                                                  "f",
                                                  "e",
                                                  "d",
                                                  "c",
                                                  "b",
                                                  "a"))

individual_level = 1:nrow(fragment)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = fragment)
res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
                             angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = -0.79, legend.pos = "none") +
  scale_y_continuous(breaks = c(-0.5,0,0.5,1), labels = c(-0.5,0,0.5,1), expand = c(0,0), limits = c(-0.94,0.87)) +
  scale_fill_manual(values = c("wheat4","gold3","slateblue3","coral4","mediumpurple4","deepskyblue4","orange4","mediumpurple3","orange3","deepskyblue3","gray47")) +
  scale_colour_manual(values = c("wheat4","gold3","slateblue3","coral4","mediumpurple4","deepskyblue4","orange4","mediumpurple3","orange3","deepskyblue3","gray47"))

# bead
bead = subset(log_means_for_plots, shape == "bead")
bead$endpoint = factor(bead$endpoint , levels = c("l",
                                                          "k", 
                                                          "j", 
                                                          "i", 
                                                          "h",
                                                          "g",
                                                          "f",
                                                          "e",
                                                          "d",
                                                          "c",
                                                          "b",
                                                          "a"))
individual_level = 1:nrow(bead)
res_all = rma.mv(yi, vi, random = ~1|DOI/individual_level, mods = ~endpoint-1, data = bead)
res_all

orchard_plot(res_all, mod = "endpoint", group = "DOI", xlab = "log transformed ratio of means",
                         angle = 0, trunk.size = 1.1, alpha = 0.4, branch.size = 1.5, k.pos = -0.79, legend.pos = "none") +
  scale_y_continuous(breaks = c(-0.5,0,0.5,1), labels = c(-0.5,0,0.5,1), expand = c(0,0), limits = c(-0.94,0.87)) +
  scale_fill_manual(values = c("gold3","cyan4","mediumpurple3","orange3","deepskyblue3","gray47")) +
  scale_colour_manual(values = c("gold3","cyan4","mediumpurple3","orange3","deepskyblue3","gray47"))

mean(fragment$combined_size_mm)
sd(fragment$combined_size_mm)
mean(bead$combined_size_mm)
sd(bead$combined_size_mm)



