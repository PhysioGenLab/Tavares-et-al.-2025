# Telomere Length on ADHD

library(TwoSampleMR)
install.packages("TwoSampleMR")

# Two sample Mendelian Randomization 

## Exposure:Telomere Outcome:ADHD

#Get instruments 

exposure_dat <- read_exposure_data(
  filename = "telomere_Codd_2021_sig_clumped.csv",
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
  samplesize_col = "N",)

exposure_dat$exposure <- "telomeres_2021"
exposure_dat$id.exposure <- "telomeres_2021"



#Outcome
outcome_dat <- read_outcome_data(snps = exposure_dat$SNP,
                                 filename = "adhd_demontis_2022_clean.csv",
                                 sep = ",",
                                 snp_col = "SNP",
                                 beta_col = "BETA",
                                 se_col = "beta_SE",
                                 effect_allele_col = "A1",
                                 other_allele_col = "A2",
                                 eaf_col = "EAF",
                                 pval_col = "P",
                                 units_col = "",
                                 gene_col = "",
                                 samplesize_col = "N",
)

outcome_dat$outcome <- "ADHD2022"
outcome_dat$id.outcome <- "ADHD2022"

write.csv(outcome_dat, file="outcome_dat_telomeres_ADHD2022.csv", row.names=F)


# Harmonise the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)
write.csv(dat, file="dat_telomeres_ADHD2022.csv", row.names=T)


# Perform MR
res <- mr(dat)
View(res)

write.csv(res, file="res_telomeres_ADHD2022.csv")



## Directionality test

directionality_test <- directionality_test(dat)

write.csv(directionality_test, file="directionality_test_telomeres_adhd.csv", row.names=F)


## Heterogeneity test

heterogeneity_test <- mr_heterogeneity(dat)

## i²
isquared <- (((heterogeneity_test$Q)-(heterogeneity_test$Q_df))/heterogeneity_test$Q)
heterogeneity_test$isquared <- isquared

write.csv(heterogeneity_test, file="heterogeneity_test_telomeres_adhd.csv")

##mean F
#Rename required columns
exposure_dat$BetaXG<-exposure_dat$beta.exposure
exposure_dat$seBetaXG<-exposure_dat$se.exposure
BetaXG   = exposure_dat$BetaXG
seBetaXG = exposure_dat$seBetaXG 

# Calculate F statistics
# and I-squared statistics
# to measure Instrument 
# strength for MR-Egger

BXG=abs(BetaXG)
F   = BXG^2/seBetaXG^2
mF  = mean(F)

exposure_dat$F <- F
write.csv(exposure_dat, file="telomere_exposure.csv")


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
cat(capture.output(print(mr_raps_overdispersed_robust), file="mr_raps_overdispersed_robust_telomere_ADHD.txt"))


#Report
report_telomeres_ADHD_2022 <- mr_report(dat)


# Contamination mixture method

#install.packages("MendelianRandomization")
library(MendelianRandomization)


MRInput_telomeres_ADHD <- mr_input(bx = dat$beta.exposure,
                                  bxse = dat$se.exposure,
                                  by = dat$beta.outcome,
                                  byse = dat$se.outcome)


mr_conmix_telomeres_ADHD <- mr_conmix(MRInput_telomeres_ADHD, psi = 0, CIMin = NA, CIMax = NA, CIStep = 0.01, alpha = 0.05)

# Exportando o resultado do mr_conmix

# Assuming your S4 object is named 'your_s4_object'
object_slots <- slotNames(mr_conmix_telomeres_ADHD)

# Create a list to store slot values
slot_values_list <- list()

# Iterate over the slots and extract values
for (slot_name in object_slots) {
  slot_values_list[[slot_name]] <- slot(mr_conmix_telomeres_ADHD, slot_name)
}

# The resulting 'slot_values_list' will contain the values of each slot

cat(capture.output(print(slot_values_list), file="mr_conmix_telomeres_ADHD.txt"))


# Steiger filtering

steiger_filtering(dat)

dat <- dat[-c(1, 3, 15), ]  # Exclui as linhas 1, 3 e 15

# Perform MR
res <- mr(dat)
View(res)

write.csv(res, file="res_telomeres_ADHD2022.csv")

