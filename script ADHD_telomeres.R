# ADHD vs telomere length

library(TwoSampleMR)
install.packages("Cairo")
install.packages("markdown")
install.packages("TwoSampleMR")
install.packages("remotes")
remotes::install_github("MRCIEU/TwoSampleMR")

# Two sample Mendelian Randomization 

## Exposure: Migraine Outcome:ADHD

#Get instruments 

exposure_dat <- read_exposure_data(
  filename = "adhd_demontis_2022_clean_sig_clumped.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  eaf_col = "EAF",
  pval_col = "P",
  units_col = "",
  gene_col = "",
  samplesize_col = "N",)

exposure_dat$exposure <- "adhd_2022"
exposure_dat$id.exposure <- "adhd_2022"

#write.csv(telomere_Codd_2021, file="telomere_Codd_2021.csv")

#Outcome
outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,
                                 filename = "telomere_Codd_2021.csv",
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col = "SE",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 eaf_col = "FREQ",
                                 pval_col = "P",
                                 units_col = "",
                                 gene_col = "",
                                 samplesize_col = "N",
)

outcome_dat$outcome <- "telomeres_2021"
outcome_dat$id.outcome <- "telomeres_2021"


write.csv(outcome_dat, file="outcome_dat_ADHD2022_telomeres.csv", row.names=F)


# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, file="dat_ADHD2022_telomeres.csv", row.names=T)


# Perform MR
res <- mr(dat)
View(res)

write.csv(res, file="res_ADHD2022_telomeres.csv")



## Directionality test

directionality_test <- directionality_test(dat)

write.csv(directionality_test, file="directionality_test_adhd_telomeres.csv", row.names=F)


## Heterogeneity test

heterogeneity_test <- mr_heterogeneity(dat)

## i²
isquared <- (((heterogeneity_test$Q)-(heterogeneity_test$Q_df))/heterogeneity_test$Q)
heterogeneity_test$isquared <- isquared

write.csv(heterogeneity_test, file="heterogeneity_test_adhd_telomeres.csv")

##mean F
#Rename required columns
## F and meanF

add_metadata()

dat$F <- (((dat$effective_n.exposure-27)-1)/27) * (dat$rsq.exposure/(1-dat$rsq.exposure)) 
F <- dat$F
mF  = mean(F)
exposure_withF <- dat
write.csv(dat, file="exposure_withrsq.csv")

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger


b_exp <- dat$beta.exposure
b_out <- dat$beta.outcome
se_exp <- dat$se.exposure
se_out <- dat$se.outcome


## MR-RAPS

#install.packages("mr.raps")
library(mr.raps)


#mr_raps_basic <- mr.raps(b_exp, b_out, se_exp, se_out, over.dispersion = TRUE,
mr_raps_overdispersed_robust <- mr.raps.overdispersed.robust(b_exp, b_out, se_exp, se_out, diagnosis = TRUE)

#cat(capture.output(print(mr_raps_basic), file="mr_raps_basic_ACR_ADHD.txt"))
cat(capture.output(print(mr_raps_overdispersed_robust), file="mr_raps_overdispersed_robust_ADHD_telomeres.txt"))


dat$samplesize.outcome <- "472174"
r.exposure <- get_r_from_lor(
  exposure_dat$beta.exposure,
  exposure_dat$eaf.exposure,
  adhd_demontis_2022_clean_sig_clumped$Nca,
  adhd_demontis_2022_clean_sig_clumped$Nco,
  0.02,
  model = "logit",
  correction = FALSE
)

exposure_dat$r.exposure <- r.exposure

r.outcome <- get_r_from_pn(outcome_dat$pval.outcome, outcome_dat$N)


#Report
report_ADHD_2022_telomeres <- mr_report(dat)


# Contamination mixture method

install.packages("MendelianRandomization")
library(MendelianRandomization)


MRInput_ADHD_telomeres <- mr_input(bx = dat$beta.exposure,
                                  bxse = dat$se.exposure,
                                  by = dat$beta.outcome,
                                  byse = dat$se.outcome)


mr_conmix_ADHD_telomeres <- mr_conmix(MRInput_ADHD_telomeres, psi = 0, CIMin = NA, CIMax = NA, CIStep = 0.01, alpha = 0.05)

# Exportando o resultado do mr_conmix

# Assuming your S4 object is named 'your_s4_object'
object_slots <- slotNames(mr_conmix_ADHD_telomeres)

# Create a list to store slot values
slot_values_list <- list()

# Iterate over the slots and extract values
for (slot_name in object_slots) {
  slot_values_list[[slot_name]] <- slot(mr_conmix_ADHD_telomeres, slot_name)
}

# The resulting 'slot_values_list' will contain the values of each slot

cat(capture.output(print(slot_values_list), file="mr_conmix_ADHD_telomeres.txt"))


# Steiger filtering

dat <- steiger_filtering(dat)

dat <- dat[-c(1, 3, 15), ]  # Exclui as linhas 1, 3 e 15

# Perform MR
res <- mr(dat)
View(res)

write.csv(res, file="res_migraine_ADHD2022.csv")

