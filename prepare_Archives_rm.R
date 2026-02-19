library(TwoSampleMR)

# prepairing archives

#ADHD Demontis 2023

adhd_demontis_2022_clean_sig <- adhd_demontis_2022_clean [which(adhd_demontis_2022_clean$P<= 0.00000005),]
write.csv(adhd_demontis_2022_clean_sig, file="adhd_demontis_2022_clean_sig.csv", row.names=F)

#perform LD clump
adhd_demontis_2022_clean_sig_clumped <- clump_data(adhd_demontis_2022_clean_sig)
write.csv(adhd_demontis_2022_clean_sig_clumped, file="adhd_demontis_2022_clean_sig_clumped.csv", row.names=F)

#Telomere Length Codd 2021
#telomere_Codd_2021.csv

library(readr)
telomere_Codd_2021 <- read_csv("telomere_Codd_2021.csv")
View(telomere_Codd_2021)

telomere_Codd_2021_sig <- telomere_Codd_2021 [which(telomere_Codd_2021$P<= 0.00000005),]
write.csv(telomere_Codd_2021_sig, file="telomere_Codd_2021_sig.csv", row.names=F)


#perform LD clump
telomere_Codd_2021_sig_clumped <- clump_data(telomere_Codd_2021_sig)
write.csv(telomere_Codd_2021_sig_clumped, file="telomere_Codd_2021_sig_clumped.csv", row.names=F)
