library("readxl")
library("tidyverse")
library("writexl")

#Manure Processing
raw_data_manure <- read_excel("Armstrong_Manure_qPCR_Raw.xlsx", sheet = "Raw_data")
meta_data_manure <- read_excel("Armstrong_Manure_qPCR_Raw.xlsx", sheet = "Meta_data")

sample_data_manure <- raw_data_manure %>% mutate(Replicate = paste(Sample, Primer, sep = ":")) %>%
  filter(!is.na(Ct))

two_tech_reps_manure <- sample_data_manure %>% group_by(Replicate) %>%
  tally() %>%
  filter(n > 1)

final_data_manure <- sample_data_manure %>% filter(Replicate %in% two_tech_reps_manure$Replicate) %>%
  group_by(Replicate) %>%
  summarise(AVG_Ct = mean(Ct), .groups = "drop") %>%
  separate(Replicate, c("Sample", "Primer"), ":", remove = TRUE) %>%
  pivot_wider(names_from = Primer, values_from = AVG_Ct) %>%
  replace(is.na(.), 0)

Delta_Ct_manure <- final_data_manure %>%
  mutate_at(vars(-matches("Sample"), -matches("16S_Eub_338_518")), list(~ . - `16S_Eub_338_518`)) %>%
  select(-2) %>%
  replace(. < 0, 0) %>%
  mutate_at(vars(-matches("Sample")), list(~ ifelse(. > 0, 2^-., .))) %>%
  merge(., meta_data_manure, by = "Sample") %>%
  rename(aaC3 = `aac3-Via`, tet33 = tet33_241_388, intl1 = intI1F165_clinical, aph3 = aphA3,
         `erm35` = `erm(35)`, `ermB` = `erm(B)`, `ermF` = `erm(F)`, `ermQ` = `erm(Q)`,
         `tet36` = `tet(36)`, `sul1` = `sul1 NEW`)

data.long_manure <- Delta_Ct_manure %>%
  pivot_longer(., cols = -c(Sample, Plot, Time, Treatment, Matrix),
               names_to = "Primer", values_to = "Relative_Abundance") %>%
  mutate(Gene_Class = case_when(Primer %in% c("aaC3", "aadD", "aadE", "aph3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm35") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfonamide",
                                Primer %in% c("intl1", "intl2", "intl3", "IS1247", "IS6100", "tnpA2",
                                              "tnpA5", "tnpA6") ~ "MGE",
                                TRUE ~ "Tetracycline"))

#Water Processing
raw_data_water <- read_excel("Armstrong_Water_Biomark_Raw.xlsx", sheet = "Raw_data")
meta_data_water <- read_excel("Armstrong_Water_Biomark_Raw.xlsx", sheet = "Meta_data")

sample_data_water <- raw_data_water %>% filter(Type == "Unknown") %>%
  filter(Primer != "16S_155y_1392") %>%
  mutate(Replicate = paste(Sample, Primer, sep = ":")) %>%
  mutate(Ct = replace(Ct, Ct > 27, NA)) %>%
  filter(!is.na(Ct))

two_tech_reps_water <- sample_data_water %>% group_by(Replicate) %>%
  tally() %>%
  filter(n > 1)

final_data_water <- sample_data_water %>% filter(Replicate %in% two_tech_reps_water$Replicate) %>%
  group_by(Replicate) %>%
  summarise(AVG_Ct = mean(Ct), .groups = "drop") %>%
  separate(Replicate, c("Sample", "Primer"), ":", remove = TRUE) %>%
  pivot_wider(names_from = Primer, values_from = AVG_Ct) %>%
  replace(is.na(.), 0)

Delta_Ct_water <- final_data_water %>%
  mutate_at(vars(-matches("Sample"), -matches("16S_Eub_338_518")), list(~ . - `16S_Eub_338_518`)) %>%
  select(-2) %>%
  replace(. < 0, 0) %>%
  mutate_at(vars(-matches("Sample")), list(~ ifelse(. > 0, 2^-., .))) %>%
  merge(., meta_data_water, by = "Sample") %>%
  rename(aaC3 = `aac3-Via`, tet33 = tet33_241_388, intl1 = intI1F165_clinical, aph3 = aphA3,
         `erm35` = `erm(35)`, `ermB` = `erm(B)`, `ermF` = `erm(F)`, `ermQ` = `erm(Q)`,
         `tet36` = `tet(36)`, `sul1` = `sul1 NEW`)

data.long_water <- Delta_Ct_water %>%
  pivot_longer(., cols = -c(Sample, Plot, Time, Treatment, Matrix),
               names_to = "Primer", values_to = "Relative_Abundance") %>%
  mutate(Gene_Class = case_when(Primer %in% c("aaC3", "aadD", "aadE", "aph3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm35") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfonamide",
                                Primer %in% c("intl1", "intl2", "intl3", "IS1247", "IS6100", "tnpA2",
                                              "tnpA5", "tnpA6") ~ "MGE",
                                TRUE ~ "Tetracycline"))

#Sediment Processing
raw_data_sediment <- read_excel("Armstrong_Water_Sediment_Biomark_Raw.xlsx", sheet = "Raw_data")
meta_data_sediment <- read_excel("Armstrong_Water_Sediment_Biomark_Raw.xlsx", sheet = "Meta_data")
null_primers_sediment <- read_excel("Armstrong_Water_Sediment_Biomark_Raw.xlsx", sheet = "Null_Primers")
plot6_3_sediment <- read_excel("Armstrong_Water_Sediment_Biomark_Raw.xlsx", sheet = "Plot6_3")

sample_data_sediment <- raw_data_sediment %>% filter(Type == "Unknown") %>%
  filter(Primer != "16S_155y_1392") %>%
  mutate(Replicate = paste(Sample, Primer, sep = ":")) %>%
  mutate(Ct = replace(Ct, Ct > 27, NA)) %>%
  filter(!is.na(Ct))

two_tech_reps_sediment <- sample_data_sediment %>% group_by(Replicate) %>%
  tally() %>%
  filter(n > 1)

final_data_sediment <- sample_data_sediment %>% filter(Replicate %in% two_tech_reps_sediment$Replicate) %>%
  group_by(Replicate) %>%
  summarise(AVG_Ct = mean(Ct), .groups = "drop") %>%
  separate(Replicate, c("Sample", "Primer"), ":", remove = TRUE) %>%
  pivot_wider(names_from = Primer, values_from = AVG_Ct) %>%
  replace(is.na(.), 0)

Delta_Ct_sediment <- final_data_sediment %>%
  mutate_at(vars(-matches("Sample"), -matches("16S_Eub_338_518")), list(~ . - `16S_Eub_338_518`)) %>%
  select(-2) %>%
  replace(. < 0, 0) %>%
  mutate_at(vars(-matches("Sample")), list(~ ifelse(. > 0, 2^-., .))) %>%
  merge(., null_primers_sediment, by = "Sample") %>%
  rbind(., plot6_3_sediment) %>%
  merge(., meta_data_sediment, by = "Sample") %>%
  rename(aaC3 = `aac3-Via`, tet33 = tet33_241_388, intl1 = intI1F165_clinical, aph3 = aphA3,
         `erm35` = `erm(35)`, `ermB` = `erm(B)`, `ermF` = `erm(F)`, `ermQ` = `erm(Q)`,
         `tet36` = `tet(36)`, `sul1` = `sul1 NEW`)

data.long_sediment <- Delta_Ct_sediment %>%
  pivot_longer(., cols = -c(Sample, Plot, Time, Treatment, Matrix),
               names_to = "Primer", values_to = "Relative_Abundance") %>%
  mutate(Gene_Class = case_when(Primer %in% c("aaC3", "aadD", "aadE", "aph3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm35") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfonamide",
                                Primer %in% c("intl1", "intl2", "intl3", "IS1247", "IS6100", "tnpA2",
                                              "tnpA5", "tnpA6") ~ "MGE",
                                TRUE ~ "Tetracycline"))

manure_water <- rbind(data.long_manure, data.long_water)

manure_water_sediment <- rbind(manure_water, data.long_sediment)

write_xlsx(manure_water_sediment, "manure_water_sediment_combined.xlsx")
