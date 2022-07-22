library(phyloseq)

## Making Phyloseq Objects

# Sequencing Run 1 (USDA)
seq_table_1 <- readRDS("sequence_data/sequence_files/STRIPS_Run_1/seq_table.RDS")
tax_table_1 <- readRDS("sequence_data/sequence_files/STRIPS_Run_1/tax_table.RDS")

phy_1 <- phyloseq(otu_table(seq_table_1, taxa_are_rows = TRUE),
                    tax_table(tax_table_1))

# Sequencing Run 2 (USDA)
seq_table_2 <- readRDS("sequence_data/sequence_files/STRIPS_Run_2/seq_table.RDS")
tax_table_2 <- readRDS("sequence_data/sequence_files/STRIPS_Run_2/tax_table.RDS")

phy_2 <- phyloseq(otu_table(seq_table_2, taxa_are_rows = TRUE),
                    tax_table(tax_table_2))

# Sequencing Run 3 (Argonne)
seq_table_3 <- readRDS("sequence_data/sequence_files/STRIPS_Run_3/seq_table.RDS")
tax_table_3 <- readRDS("sequence_data/sequence_files/STRIPS_Run_3/tax_table.RDS")

phy_3 <- phyloseq(otu_table(seq_table_3, taxa_are_rows = TRUE),
                  tax_table(tax_table_3))

## Merging Phyloseq Objects and removing excess sample labels

merged_phy <- merge_phyloseq(phy_1, phy_2, phy_3)

sample_names(merged_phy) <- gsub('_L001_R1_001', '', sample_names(merged_phy))
sample_names(merged_phy) <- gsub('(.*)_\\w+', '\\1', sample_names(merged_phy))
sample_names(merged_phy)

## Removing Mock Communities

merged_phy <- prune_samples(sample_names(merged_phy) != "mock_01" &
                             sample_names(merged_phy) != "mock_02" &
                             sample_names(merged_phy) != "mock_03" &
                             sample_names(merged_phy) != "mock_04" &
                             sample_names(merged_phy) != "mock_05" &
                             sample_names(merged_phy) != "mock_06" &
                             sample_names(merged_phy) != "mock_07" &
                             sample_names(merged_phy) != "mock_08" &
                             sample_names(merged_phy) != "mock_09" &
                             sample_names(merged_phy) != "mock_10",
                           merged_phy)

## Listing Remaining Samples

remaining_sample_names <- data.frame(unique(sample_names(merged_phy)))

## Listing Removed Samples

total_sample_names <- read.csv(file = "sequence_data/Armstrong_Sequences_Metadata.csv")

removed_samples <- data.frame(setdiff(total_sample_names$sample_name,
                           remaining_sample_names$unique.sample_names.merged_phy..))

removed_samples_list <- setdiff(total_sample_names$sample_name,
                                      remaining_sample_names$unique.sample_names.merged_phy..)

## Create New Metadata

meta_data <- dplyr::filter(total_sample_names, !sample_name %in% removed_samples_list)

## Assign ASV

dna <- Biostrings::DNAStringSet(taxa_names(merged_phy))
names(dna) <- taxa_names(merged_phy)
merged_phy <- merge_phyloseq(merged_phy, dna)
taxa_names(merged_phy) <- paste0("ASV", seq(ntaxa(merged_phy)))

## Add Metadata

sam.info <- sample_data(meta_data)
rownames(sam.info) <- sam.info$sample_name

merged_phy_meta <- merge_phyloseq(merged_phy, sam.info)
head(sample_data(merged_phy_meta))

## Subset to Bacteria

merged_phy_meta_bacteria = subset_taxa(merged_phy_meta, Kingdom == "Bacteria")
table(tax_table(merged_phy_meta_bacteria)[, "Kingdom"], exclude = NULL)

#Remove empty OTUs
any(taxa_sums(merged_phy_meta_bacteria) == 0)

sum(taxa_sums(merged_phy_meta_bacteria) == 0)

merged_phy_meta_bacteria0 = prune_taxa(taxa_sums(merged_phy_meta_bacteria) > 0, merged_phy_meta_bacteria)


## Save Phyloseq Object

saveRDS(merged_phy_meta_bacteria0, file = "Armstrong_Rainfall_Phy.RDS")

