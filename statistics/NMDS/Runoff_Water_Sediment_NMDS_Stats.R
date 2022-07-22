library(phyloseq)
library(tidyverse)
library(vegan)

set.seed(8675309)

phy <- readRDS("Armstrong_Rainfall_Phy.RDS")
permutations <- 9999

runoff_water_sediment_manure <- subset_samples(phy, matrix == "manure" | matrix == "runoff_water" | matrix == "runoff_sediment")

all_data <- data.frame(sample_data(runoff_water_sediment_manure)) %>%
  mutate(Matrix = case_when(matrix == "runoff_water" ~ "Runoff Water",
                            matrix == "runoff_sediment" ~ "Runoff Sediment",
                            matrix == "manure" ~ "Manure")) %>%
  mutate(Treatment = case_when(treatment == "ACM" ~ "No Strip + Manure",
                               treatment == "ACS" ~ "Strip + No Manure",
                               treatment == "ACSM" ~ "Strip + Manure",
                               treatment == "AM" ~ "Manure"))

sample_names(runoff_water_sediment_manure)

rownames(all_data) <- all_data$sample_name

sample_data(runoff_water_sediment_manure) <- all_data

all_dist <- phyloseq::distance(runoff_water_sediment_manure, method = "bray")

adonis(all_dist ~ Matrix, data = as(sample_data(runoff_water_sediment_manure), "data.frame"), permutations = permutations)

adonis(all_dist ~ Treatment, data = as(sample_data(runoff_water_sediment_manure), "data.frame"), permutations = permutations)

#Subset matrices
runoff_water <- subset_samples(runoff_water_sediment_manure, Matrix == "Runoff Water")
runoff_sediment <- subset_samples(runoff_water_sediment_manure, Matrix == "Runoff Sediment")

##WATER COMPARISONS
pairwise_p_water <- numeric()

#Strip + No Manure and Strip + Manure Water
ACS_ACSM_water <- subset_samples(runoff_water, treatment == "ACS" | treatment == "ACSM")

ACS_ACSM_water_dist <- phyloseq::distance(ACS_ACSM_water, method = "bray")

ACS_ACSM_water_test <- adonis(ACS_ACSM_water_dist ~ Treatment, data = as(sample_data(ACS_ACSM_water), "data.frame"), permutations = permutations)

pairwise_p_water["ACS_ACSM"] <- ACS_ACSM_water_test[["aov.tab"]][["Pr(>F)"]][1]

#Strip + Manure and No Strip + Manure Water
ACSM_ACM_water <- subset_samples(runoff_water, treatment == "ACSM" | treatment == "ACM")

ACSM_ACM_water_dist <- phyloseq::distance(ACSM_ACM_water, method = "bray")

ACSM_ACM_water_test <- adonis(ACSM_ACM_water_dist ~ Treatment, data = as(sample_data(ACSM_ACM_water), "data.frame"), permutations = permutations)

pairwise_p_water["ACSM_ACM"] <- ACSM_ACM_water_test[["aov.tab"]][["Pr(>F)"]][1]

#No Strip + Manure and Strip + No Manure Water
ACM_ACS_water <- subset_samples(runoff_water, treatment == "ACM" | treatment == "ACS")

ACM_ACS_water_dist <- phyloseq::distance(ACM_ACS_water, method = "bray")

ACM_ACS_water_test <- adonis(ACM_ACS_water_dist ~ Treatment, data = as(sample_data(ACM_ACS_water), "data.frame"), permutations = permutations)

pairwise_p_water["ACM_ACS"] <- ACM_ACS_water_test[["aov.tab"]][["Pr(>F)"]][1]

pairwise_p_water_adjust <- p.adjust(pairwise_p_water)

##SEDIMENT COMPARISONS
pairwise_p_sediment <- numeric()

#Strip + No Manure and Strip + Manure Sediment
ACS_ACSM_sediment <- subset_samples(runoff_sediment, treatment == "ACS" | treatment == "ACSM")

ACS_ACSM_sediment_dist <- phyloseq::distance(ACS_ACSM_sediment, method = "bray")

ACS_ACSM_sediment_test <- adonis(ACS_ACSM_sediment_dist ~ Treatment, data = as(sample_data(ACS_ACSM_sediment), "data.frame"), permutations = permutations)

pairwise_p_sediment["ACS_ACMS"] <- ACS_ACSM_sediment_test[["aov.tab"]][["Pr(>F)"]][1]

#Strip + Manure and No Strip + Manure Sediment
ACSM_ACM_sediment <- subset_samples(runoff_sediment, treatment == "ACSM" | treatment == "ACM")

ACSM_ACM_sediment_dist <- phyloseq::distance(ACSM_ACM_sediment, method = "bray")

ACSM_ACM_sediment_test <- adonis(ACSM_ACM_sediment_dist ~ Treatment, data = as(sample_data(ACSM_ACM_sediment), "data.frame"), permutations = permutations)

pairwise_p_sediment["ACSM_ACM"] <- ACSM_ACM_sediment_test[["aov.tab"]][["Pr(>F)"]][1]

#No Strip + Manure and Strip + No Manure Sediment
ACM_ACS_sediment <- subset_samples(runoff_sediment, treatment == "ACM" | treatment == "ACS")

ACM_ACS_sediment_dist <- phyloseq::distance(ACM_ACS_sediment, method = "bray")

ACM_ACS_sediment_test <- adonis(ACM_ACS_sediment_dist ~ Treatment, data = as(sample_data(ACM_ACS_sediment), "data.frame"), permutations = permutations)

pairwise_p_sediment["ACM_ACS"] <- ACM_ACS_sediment_test[["aov.tab"]][["Pr(>F)"]][1]

pairwise_p_sediment_adjust <- p.adjust(pairwise_p_sediment)



#No Strip + Manure Sediment and Manure
sediment_manure <- subset_samples(runoff_water_sediment_manure, Matrix == "Runoff Sediment" | Matrix == "Manure")
ACM_sediment_manure <- subset_samples(sediment_manure, treatment == "ACM" | treatment == "AM")

ACM_sediment_manure_dist <- phyloseq::distance(ACM_sediment_manure, method = "bray")

ACM_sediment_manure_test <- adonis(ACM_sediment_manure_dist ~ Treatment, data = as(sample_data(ACM_sediment_manure), "data.frame"), permutations = permutations)

pairwise_p_sediment["ACM_Manure"] <- ACM_sediment_manure_test[["aov.tab"]][["Pr(>F)"]][1]

pairwise_p_sediment_adjust <- p.adjust(pairwise_p_sediment)
