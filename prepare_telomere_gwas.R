# Rename column B to NewName
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "variant_id"] <- "SNP"
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "p_value"] <- "P"
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "chromosome"] <- "CHR"
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "base_pair_location"] <- "BP"
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "effec_allele"] <- "A1"
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "other_allele"] <- "A2"
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "effect_allele_frequency"] <- "FREQ"
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "beta"] <- "BETA"
colnames(telomere_Codd_2021)[colnames(telomere_Codd_2021) == "standard_error"] <- "SE"

# salvar arquivo

write.csv(telomere_Codd_2021, file="telomere_Codd_2021.csv")
