library("readxl")
library("tidyverse")
library("writexl")

#Manure Processing
raw_data_manure <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Raw_Data") %>%
  filter(Matrix == "Manure")

meta_data_manure <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Meta_Data") %>%
  filter(Matrix == "Manure")

sample_data_manure <- raw_data_manure %>% mutate(Replicate = paste(Sample, Primer, sep = ":")) %>%
  filter(!is.na(Ct))

final_data_manure <- sample_data_manure %>%
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
  merge(., meta_data_manure, by = "Sample")

data.long_manure <- Delta_Ct_manure %>%
  pivot_longer(., cols = -c(Sample, Plot, Time, Treatment, Matrix),
               names_to = "Primer", values_to = "Relative_Abundance") %>%
  mutate(Gene_Class = case_when(Primer %in% c("aadD", "aadE", "aphA3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm(35)") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfonamide",
                                Primer %in% c("intl1(clinical)", "intl2", "IS1247", "IS6100", "tnpA-02",
                                              "tnpA-05", "tnpA-06") ~ "MGE",
                                Primer %in% c("tet(36)", "tet(44)", "tet(L)", "tet(M)", "tet(O)", "tet(Q)", "tet(T)", "tet(W)", "tet(X)") ~ "Tetracycline"))

#Water Processing
raw_data_water <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Raw_Data") %>%
  filter(Matrix == "Water")

meta_data_water <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Meta_Data") %>%
  filter(Matrix == "Water")

null_primers_water <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Water_Null_Primers")

sample_data_water <- raw_data_water %>%
  mutate(Replicate = paste(Sample, Primer, sep = ":")) %>%
  mutate(Ct = replace(Ct, Ct > 27, NA)) %>%
  filter(!is.na(Ct))

three_reps_water <- sample_data_water %>% group_by(Replicate) %>%
  tally() %>%
  filter(n > 2)

final_data_water <- sample_data_water %>% filter(Replicate %in% three_reps_water$Replicate) %>%
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
  merge(., null_primers_water, by = "Sample") %>%
  merge(., meta_data_water, by = "Sample")

data.long_water <- Delta_Ct_water %>%
  pivot_longer(., cols = -c(Sample, Plot, Time, Treatment, Matrix),
               names_to = "Primer", values_to = "Relative_Abundance") %>%
  mutate(Gene_Class = case_when(Primer %in% c("aadD", "aadE", "aphA3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm(35)") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfonamide",
                                Primer %in% c("intl1(clinical)", "intl2", "IS1247", "IS6100", "tnpA-02",
                                              "tnpA-05", "tnpA-06") ~ "MGE",
                                Primer %in% c("tet(36)", "tet(44)", "tet(L)", "tet(M)", "tet(O)", "tet(Q)", "tet(T)", "tet(W)", "tet(X)") ~ "Tetracycline"))

#Sediment Processing
raw_data_sediment <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Raw_Data") %>%
  filter(Matrix == "Sediment")

meta_data_sediment <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Meta_Data") %>%
  filter(Matrix == "Sediment")

null_primers_sediment <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Sediment_Null_Primers")

null_plots_sediment <- read_excel("Armstrong_qPCR_raw.xlsx", sheet = "Sediment_Null_Plots")

sample_data_sediment <- raw_data_sediment %>%
  mutate(Replicate = paste(Sample, Primer, sep = ":")) %>%
  mutate(Ct = replace(Ct, Ct > 27, NA)) %>%
  filter(!is.na(Ct))

three_reps_sediment <- sample_data_sediment %>% group_by(Replicate) %>%
  tally() %>%
  filter(n > 2)

final_data_sediment <- sample_data_sediment %>% filter(Replicate %in% three_reps_sediment$Replicate) %>%
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
  rbind(., null_plots_sediment) %>%
  merge(., meta_data_sediment, by = "Sample")

data.long_sediment <- Delta_Ct_sediment %>%
  pivot_longer(., cols = -c(Sample, Plot, Time, Treatment, Matrix),
               names_to = "Primer", values_to = "Relative_Abundance") %>%
  mutate(Gene_Class = case_when(Primer %in% c("aadD", "aadE", "aphA3", "sat4") ~ "Aminoglycoside",
                                Primer %in% c("ermB", "ermF", "ermQ", "ermT", "erm(35)") ~ "MLSB",
                                Primer %in% c("sul1", "sul2") ~ "Sulfonamide",
                                Primer %in% c("intl1(clinical)", "intl2", "IS1247", "IS6100", "tnpA-02",
                                              "tnpA-05", "tnpA-06") ~ "MGE",
                                Primer %in% c("tet(36)", "tet(44)", "tet(L)", "tet(M)", "tet(O)", "tet(Q)", "tet(T)", "tet(W)", "tet(X)") ~ "Tetracycline"))

manure_water <- rbind(data.long_manure, data.long_water)

manure_water_sediment <- rbind(manure_water, data.long_sediment)

write_xlsx(manure_water_sediment, "Armstrong_qPCR_processed.xlsx")
