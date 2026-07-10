
# Load Data ----
sessionInfo()
all_NMP_data = read.csv2("Data/Doeringetal_2026_NMP.csv", stringsAsFactors = T, fileEncoding = "latin1")


# Remove columns that are not used in the analysis  ---- 

# Delete all columns with only one value (in most cases "NA", "yes", or "no")
new_NMP_data = Filter(function(x)(length(unique(x))>1), all_NMP_data)


# Transform factors into numerics/characters ----
new_NMP_data$c_unit = as.character(new_NMP_data$c_unit)
new_NMP_data$t_unit = as.character(new_NMP_data$t_unit)
new_NMP_data$endpoint = as.character(new_NMP_data$endpoint)
new_NMP_data$c_n_replicates = as.integer(new_NMP_data$c_n_replicates)
new_NMP_data$t_n_replicates = as.integer(new_NMP_data$t_n_replicates)
new_NMP_data$c_mean = as.numeric(levels(new_NMP_data$c_mean))[new_NMP_data$c_mean]
new_NMP_data$t_mean = as.numeric(levels(new_NMP_data$t_mean))[new_NMP_data$t_mean]
new_NMP_data$c_SD = as.numeric(levels(new_NMP_data$c_SD))[new_NMP_data$c_SD]
new_NMP_data$t_SD = as.numeric(levels(new_NMP_data$t_SD))[new_NMP_data$t_SD]
new_NMP_data$c_SE = as.numeric(levels(new_NMP_data$c_SE))[new_NMP_data$c_SE]
new_NMP_data$t_SE = as.numeric(levels(new_NMP_data$t_SE))[new_NMP_data$t_SE]
new_NMP_data$nominal_dose_altern_type = as.numeric(levels(new_NMP_data$nominal_dose_altern_type))[new_NMP_data$nominal_dose_altern_type]
new_NMP_data$nominal_dose_mass = as.numeric(levels(new_NMP_data$nominal_dose_mass))[new_NMP_data$nominal_dose_mass]
new_NMP_data$nominal_dose_mass_units = as.character(levels(new_NMP_data$nominal_dose_mass_units))[new_NMP_data$nominal_dose_mass_units]
new_NMP_data$size_length_mm_nominal_mean = as.numeric(new_NMP_data$size_length_mm_nominal_mean)
new_NMP_data$size_length_mm_measured_mean = as.numeric(new_NMP_data$size_length_mm_measured_mean)

# Unify units ----

# Transform "µmol/l_homogenate" in "nmol/ml_homogenate"
new_NMP_data[new_NMP_data == "µmol/l_homogenate"] = "nmol/ml_homogenate"

# Transform "µmol/g_protein" in "nmol/mg_protein"
new_NMP_data[new_NMP_data == "µmol/g_protein"] = "nmol/mg_protein"

# Transform "mmol/mg_protein" in "nmol/mg_protein"
temporary_data = subset(new_NMP_data, t_unit == "mmol/mg_protein")
temporary_data2 = subset(new_NMP_data, t_unit != "mmol/mg_protein")
temporary_data$c_mean = temporary_data$c_mean*1000000
temporary_data$c_SD = temporary_data$c_SD*1000000
temporary_data$c_SE = temporary_data$c_SE*1000000
temporary_data$t_mean = temporary_data$t_mean*1000000
temporary_data$t_SD = temporary_data$t_SD*1000000
temporary_data$t_SE = temporary_data$t_SE*1000000
temporary_data[temporary_data == "mmol/mg_protein"] = "nmol/mg_protein"
new_NMP_data = rbind(temporary_data, temporary_data2)

# Transform "fluorescence/µg_protein" in "fluorescence/mg_protein"
temporary_data = subset(new_NMP_data, t_unit == "fluorescence/µg_protein")
temporary_data2 = subset(new_NMP_data, t_unit != "fluorescence/µg_protein")
temporary_data$c_mean = temporary_data$c_mean*1000
temporary_data$c_SD = temporary_data$c_SD*1000
temporary_data$c_SE = temporary_data$c_SE*1000
temporary_data$t_mean = temporary_data$t_mean*1000
temporary_data$t_SD = temporary_data$t_SD*1000
temporary_data$t_SE = temporary_data$t_SE*1000
temporary_data[temporary_data == "fluorescence/µg_protein"] = "fluorescence/mg_protein"
new_NMP_data = rbind(temporary_data, temporary_data2)

# Transform "mg/g_protein" in "ng/mg_protein"
temporary_data = subset(new_NMP_data, t_unit == "mg/g_protein")
temporary_data2 = subset(new_NMP_data, t_unit != "mg/g_protein")
temporary_data$c_mean = temporary_data$c_mean*1000
temporary_data$c_SD = temporary_data$c_SD*1000
temporary_data$c_SE = temporary_data$c_SE*1000
temporary_data$t_mean = temporary_data$t_mean*1000
temporary_data$t_SD = temporary_data$t_SD*1000
temporary_data$t_SE = temporary_data$t_SE*1000
temporary_data[temporary_data == "mg/g_protein"] = "ng/mg_protein"
new_NMP_data = rbind(temporary_data, temporary_data2)

# Transform "U/g_protein" in "U/mg_protein"
temporary_data = subset(new_NMP_data, t_unit == "U/g_protein")
temporary_data2 = subset(new_NMP_data, t_unit != "U/g_protein")
temporary_data$c_mean = temporary_data$c_mean/1000
temporary_data$c_SD = temporary_data$c_SD/1000
temporary_data$c_SE = temporary_data$c_SE/1000
temporary_data$t_mean = temporary_data$t_mean/1000
temporary_data$t_SD = temporary_data$t_SD/1000
temporary_data$t_SE = temporary_data$t_SE/1000
temporary_data[temporary_data == "U/g_protein"] = "U/mg_protein"
new_NMP_data = rbind(temporary_data, temporary_data2)

# Transform "nmol/min/mg_protein" in "U/mg_protein"
temporary_data = subset(new_NMP_data, t_unit == "nmol/min/mg_protein")
temporary_data2 = subset(new_NMP_data, t_unit != "nmol/min/mg_protein")
temporary_data$c_mean = temporary_data$c_mean/1000
temporary_data$c_SD = temporary_data$c_SD/1000
temporary_data$c_SE = temporary_data$c_SE/1000
temporary_data$t_mean = temporary_data$t_mean/1000
temporary_data$t_SD = temporary_data$t_SD/1000
temporary_data$t_SE = temporary_data$t_SE/1000
temporary_data[temporary_data == "nmol/min/mg_protein"] = "U/mg_protein"
new_NMP_data = rbind(temporary_data, temporary_data2)

# Transform "percentage_to_control" in "fold_change"
temporary_data = subset(new_NMP_data, t_unit == "percentage_to_control")
temporary_data2 = subset(new_NMP_data, t_unit != "percentage_to_control")
temporary_data$c_mean = temporary_data$c_mean/100
temporary_data$c_SD = temporary_data$c_SD/100
temporary_data$c_SE = temporary_data$c_SE/100
temporary_data$t_mean = temporary_data$t_mean/100
temporary_data$t_SD = temporary_data$t_SD/100
temporary_data$t_SE = temporary_data$t_SE/100
temporary_data[temporary_data == "percentage_to_control"] = "fold_change"
new_NMP_data = rbind(temporary_data, temporary_data2)

# Transform "change_ratio" in "fold_change"
temporary_data = subset(new_NMP_data, t_unit == "change_ratio")
temporary_data2 = subset(new_NMP_data, t_unit != "change_ratio")
temporary_data$c_mean = temporary_data$c_mean+1
temporary_data$c_SD = temporary_data$c_SD+1
temporary_data$c_SE = temporary_data$c_SE+1
temporary_data$t_mean = temporary_data$t_mean+1
temporary_data$t_SD = temporary_data$t_SD+1
temporary_data$t_SE = temporary_data$t_SE+1
temporary_data[temporary_data == "change_ratio"] = "fold_change"
new_NMP_data = rbind(temporary_data, temporary_data2)

# Transform "ng/l_homogenate" of 8-OHdG (283.24 g/mol) in "nmol/ml_homogenate"
temporary_data = subset(new_NMP_data, endpoint == "8-OHdG Concentration")
temporary_data2 = subset(new_NMP_data, endpoint != "8-OHdG Concentration")
temporary_data3 = subset(temporary_data, t_unit == "ng/l_homogenate") 
temporary_data4 = subset(temporary_data, t_unit != "ng/l_homogenate")
temporary_data3$c_mean = (temporary_data3$c_mean/(283.24*1000))
temporary_data3$c_SD = (temporary_data3$c_SD/(283.24*1000))
temporary_data3$c_SE = (temporary_data3$c_SE/(283.24*1000))
temporary_data3$t_mean = (temporary_data3$t_mean/(283.24*1000))
temporary_data3$t_SD = (temporary_data3$t_SD/(283.24*1000))
temporary_data3$t_SE = (temporary_data3$t_SE/(283.24*1000))
temporary_data3[temporary_data3 == "ng/l_homogenate"] = "nmol/ml_homogenate"
new_NMP_data = rbind(temporary_data2, temporary_data3, temporary_data4)

# Transform "ng/mg_protein" of 8-OHdG (283.24 g/mol) in "nmol/mg_protein"
temporary_data = subset(new_NMP_data, endpoint == "8-OHdG Concentration")
temporary_data2 = subset(new_NMP_data, endpoint != "8-OHdG Concentration")
temporary_data3 = subset(temporary_data, t_unit == "ng/mg_protein") 
temporary_data4 = subset(temporary_data, t_unit != "ng/mg_protein")
temporary_data3$c_mean = (temporary_data3$c_mean/(283.24))
temporary_data3$c_SD = (temporary_data3$c_SD/(283.24))
temporary_data3$c_SE = (temporary_data3$c_SE/(283.24))
temporary_data3$t_mean = (temporary_data3$t_mean/(283.24))
temporary_data3$t_SD = (temporary_data3$t_SD/(283.24))
temporary_data3$t_SE = (temporary_data3$t_SE/(283.24))
temporary_data3[temporary_data3 == "ng/mg_protein"] = "nmol/mg_protein"
new_NMP_data = rbind(temporary_data2, temporary_data3, temporary_data4)

# Transform "ng/mg_protein" of Glutathione (307.33 g/mol) in "nmol/mg_protein"
temporary_data = subset(new_NMP_data, endpoint == "Glutathione Concentration")
temporary_data2 = subset(new_NMP_data, endpoint != "Glutathione Concentration")
temporary_data3 = subset(temporary_data, t_unit == "ng/mg_protein") 
temporary_data4 = subset(temporary_data, t_unit != "ng/mg_protein")
temporary_data3$c_mean = (temporary_data3$c_mean/(307.33))
temporary_data3$c_SD = (temporary_data3$c_SD/(307.33))
temporary_data3$c_SE = (temporary_data3$c_SE/(307.33))
temporary_data3$t_mean = (temporary_data3$t_mean/(307.33))
temporary_data3$t_SD = (temporary_data3$t_SD/(307.33))
temporary_data3$t_SE = (temporary_data3$t_SE/(307.33))
temporary_data3[temporary_data3 == "ng/mg_protein"] = "nmol/mg_protein"
new_NMP_data = rbind(temporary_data2, temporary_data3, temporary_data4)

# Transform "ng/mg_protein" of Glutathione (307.33 g/mol) in "nmol/mg_protein"
temporary_data = subset(new_NMP_data, endpoint == "Glutathione Concentration")
temporary_data2 = subset(new_NMP_data, endpoint != "Glutathione Concentration")
temporary_data3 = subset(temporary_data, t_unit == "ng/mg_protein") 
temporary_data4 = subset(temporary_data, t_unit != "ng/mg_protein")
temporary_data3$c_mean = (temporary_data3$c_mean/(307.33))
temporary_data3$c_SD = (temporary_data3$c_SD/(307.33))
temporary_data3$c_SE = (temporary_data3$c_SE/(307.33))
temporary_data3$t_mean = (temporary_data3$t_mean/(307.33))
temporary_data3$t_SD = (temporary_data3$t_SD/(307.33))
temporary_data3$t_SE = (temporary_data3$t_SE/(307.33))
temporary_data3[temporary_data3 == "ng/mg_protein"] = "nmol/mg_protein"
new_NMP_data = rbind(temporary_data2, temporary_data3, temporary_data4)

# Transform "Thiobarbituric Acid Reactive Species" in "Malondialdehyde Concentration"
new_NMP_data[new_NMP_data == "Thiobarbituric Acid Reactive Species"] = "Malondialdehyde Concentration"

# Combine measured and nominal particle concentration in one column
temporary_data = subset(new_NMP_data, nominal_dose_mass != "NA" &  is.na(new_NMP_data$nominal_dose_altern_type))
temporary_data2 = subset(new_NMP_data, is.na(new_NMP_data$nominal_dose_mass) & nominal_dose_altern_type != "NA")
temporary_data3 = subset(new_NMP_data, nominal_dose_mass != "NA" & nominal_dose_altern_type != "NA")
temporary_data4 = subset(new_NMP_data, is.na(new_NMP_data$nominal_dose_mass) & is.na(new_NMP_data$nominal_dose_altern_type))
temporary_data$combined_concentration = temporary_data$nominal_dose_mass
temporary_data$combined_concentration_unit = temporary_data$nominal_dose_mass_units
temporary_data$combined_concentration = temporary_data$combined_concentration/10000
temporary_data$combined_concentration_unit = "% w w-1"
temporary_data2$combined_concentration = temporary_data2$nominal_dose_altern_type
temporary_data2$combined_concentration_unit = temporary_data2$nominal_dose_altern_type_units
temporary_data3$combined_concentration = temporary_data3$nominal_dose_altern_type
temporary_data3$combined_concentration_unit = temporary_data3$nominal_dose_altern_type_units
temporary_data4$combined_concentration = temporary_data4$nominal_dose_altern_type
temporary_data4$combined_concentration_unit = temporary_data4$nominal_dose_altern_type_units
new_NMP_data = rbind(temporary_data, temporary_data2, temporary_data3, temporary_data4)
new_NMP_data$log_combined_concentration = log(new_NMP_data$combined_concentration)

# Combine measured and nominal particle size in one column
temporary_data = subset(new_NMP_data, size_length_mm_nominal_mean != "NA" &  is.na(new_NMP_data$size_length_mm_measured_mean))
temporary_data2 = subset(new_NMP_data, is.na(new_NMP_data$size_length_mm_nominal_mean) & size_length_mm_measured_mean != "NA")
temporary_data3 = subset(new_NMP_data, size_length_mm_nominal_mean != "NA" & size_length_mm_measured_mean != "NA")
temporary_data4 = subset(new_NMP_data, is.na(new_NMP_data$size_length_mm_nominal_mean) & is.na(new_NMP_data$size_length_mm_measured_mean))
temporary_data$combined_size_mm = temporary_data$size_length_mm_nominal_mean
temporary_data2$combined_size_mm = temporary_data2$size_length_mm_measured_mean
temporary_data3$combined_size_mm = temporary_data3$size_length_mm_nominal_mean
temporary_data4$combined_size_mm = temporary_data4$size_length_mm_nominal_mean
new_NMP_data = rbind(temporary_data, temporary_data2, temporary_data3, temporary_data4)

# Abbreviate ROS marker
new_NMP_data[new_NMP_data == "Reactive Oxygen Species Production"] = "ROS"
new_NMP_data[new_NMP_data == "Malondialdehyde Concentration"] = "MDA"
new_NMP_data[new_NMP_data == "Catalase Activity"] = "CAT"
new_NMP_data[new_NMP_data == "Superoxide Dismutase Activity"] = "SOD"
new_NMP_data[new_NMP_data == "Glutathione S-Transferase Activity"] = "GST"
new_NMP_data[new_NMP_data == "Peroxidase Activity"] = "POD"
new_NMP_data[new_NMP_data == "8-OHdG Concentration"] = "8-OHdG"
new_NMP_data[new_NMP_data == "Total antioxidant capacity (T-AOC)"] = "TAC"
new_NMP_data[new_NMP_data == "sod mRNA expression"] = "SOD mRNA expression"
new_NMP_data[new_NMP_data == "Catalase mRNA expression"] = "CAT mRNA expression"
new_NMP_data[new_NMP_data == "Glutathione Concentration"] = "GSH"
new_NMP_data[new_NMP_data == "Glutathione S-Transferase mRNA expression"] = "GST mRNA expression"

# Delete temporary data
rm(temporary_data, temporary_data2, temporary_data3, temporary_data4)

# Overview
table(new_NMP_data$t_unit, new_NMP_data$endpoint)

# Export prepared data ----
write.csv2(new_NMP_data, "Data/new_NMP_data.csv", row.names = FALSE)


