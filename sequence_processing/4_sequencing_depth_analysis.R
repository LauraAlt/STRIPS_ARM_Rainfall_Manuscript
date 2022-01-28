library("phyloseq")
library("data.table")
library("writexl")

phy = readRDS("Armstrong_Rainfall_Phy.RDS")
soil = subset_samples(phy, matrix == "soil")
manure = subset_samples(phy, matrix == "manure")
runoff_water = subset_samples(phy, matrix == "runoff_water")
runoff_sediment = subset_samples(phy, matrix == "runoff_sediment")

phy
# phyloseq-class experiment-level object
# otu_table()   OTU Table:         [ 60924 taxa and 827 samples ]
# sample_data() Sample Data:       [ 827 samples by 13 sample variables ]
# tax_table()   Taxonomy Table:    [ 60924 taxa by 6 taxonomic ranks ]
# refseq()      DNAStringSet:      [ 60924 reference sequences ]

sum(sample_sums(phy))
# 12984501

mean(sample_sums(phy))
# 15700.73
ntaxa(phy)
# 60924
min(sample_sums(phy))
# 804
max(sample_sums(phy))
# 45537

phy_soil = prune_taxa(taxa_sums(soil) > 0, soil)
any(taxa_sums(phy_soil) == 0)
mean(sample_sums(phy_soil))
# 16079.06
ntaxa(phy_soil)
# 54880
min(sample_sums(phy_soil))
# 1864
max(sample_sums(phy_soil))
# 45537

phy_manure = prune_taxa(taxa_sums(manure) > 0, manure)
any(taxa_sums(phy_manure) == 0)
mean(sample_sums(phy_manure))
# 23666.36
ntaxa(phy_manure)
# 2159
min(sample_sums(phy_manure))
# 17148
max(sample_sums(phy_manure))
# 29105

phy_runoff_water = prune_taxa(taxa_sums(runoff_water) > 0, runoff_water)
any(taxa_sums(phy_runoff_water) == 0)
mean(sample_sums(phy_runoff_water))
# 10601.04
ntaxa(phy_runoff_water)
# 3074
min(sample_sums(phy_runoff_water))
# 804
max(sample_sums(phy_runoff_water))
# 20631

phy_runoff_sediment = prune_taxa(taxa_sums(runoff_sediment) > 0, runoff_sediment)
any(taxa_sums(phy_runoff_sediment) == 0)
mean(sample_sums(phy_runoff_sediment))
# 10182.08
ntaxa(phy_runoff_sediment)
# 8895
min(sample_sums(phy_runoff_sediment))
# 2914
max(sample_sums(phy_runoff_sediment))
# 30541

unique(tax_table(phy)[, "Phylum"])
table(tax_table(phy)[, "Phylum"], exclude = NULL)

sdt = data.table(as(sample_data(phy), "data.frame"),
                 TotalReads = sample_sums(phy), keep.rownames = TRUE)

setnames(sdt, "rn", "SampleID")

write_xlsx(x = sdt, 
           path = "Armstrong_Rainfall_Samples_Sequence_Depth.xlsx")
