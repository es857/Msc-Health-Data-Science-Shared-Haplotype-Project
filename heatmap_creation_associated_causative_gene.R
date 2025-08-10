library(dplyr)
library(ggplot2)
library(data.table)
library(ggplot2)
library(reshape2)
library(dplyr)
library(stringr)

colnames(geno) <- c("Position", "Sample", "Genotype")
# 1. Prepare gene_data (ensure chromosome_name is character)
gene_data <- gene_data %>%
  mutate(chromosome_name = as.character(chromosome_name))

gene_data_chr19 <- gene_data %>%
  filter(hgnc_symbol == "Chro_19")

gene_data_other <- gene_data %>%
  filter(hgnc_symbol != "Chro_19")

# 2. Process geno to extract Chr/Pos (handles X, Y, M)
geno <- geno %>%
  mutate(
    Chr = str_extract(Position, "^chr([0-9XYM]+):") %>% 
      str_remove("^chr") %>% 
      str_remove(":"),
    Pos = as.numeric(str_remove(Position, "^chr[0-9XYM]+:"))
  )

# 3. Filter for target gene regions and label variants
buffer_bp <- 10000
geno_other <- geno %>%
  inner_join(gene_data_other, by = c("Chr" = "chromosome_name")) %>%
  filter(Pos >= start_position - buffer_bp & Pos <= end_position + buffer_bp) %>%
  mutate(
    GenotypeCode = case_when(
      Genotype %in% c("0/0") ~ "Hom_REF",
      Genotype %in% c("0/1", "1/0") ~ "Het",
      Genotype %in% c("1/1") ~ "Hom_ALT",
      TRUE ~ "Missing"
    ),
    # Add gene labels to positions for the x-axis
    Pos_Label = paste0(hgnc_symbol, ":", Pos)  # e.g., "FICD:48140234"
  ) %>%
  arrange(Pos)  # Sort by genomic position

# 4. Plot (all genes combined)
ggplot(geno_other, aes(x = Pos_Label, y = Sample, fill = GenotypeCode)) +
  geom_tile() +
  scale_fill_manual(
    values = c("Het" = "yellow", "Hom_ALT" = "red", 
               "Hom_REF" = "green", "Missing" = "gray")
  ) +
  labs(
    title = "Genotype Heatmap: All genes except Chromosome 19 Haplotype",
    x = "Genomic Position (Gene:BasePair)",
    y = "Sample"
  ) +
  scale_x_discrete(breaks = geno_other$Pos_Label[seq(1, nrow(geno_other), by = 100)]) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8)
  )



geno_chr19 <- geno %>%
  inner_join(gene_data_chr19, by = c("Chr" = "chromosome_name")) %>%
  filter(Pos >= start_position & Pos <= end_position) %>%
  mutate(
    GenotypeCode = case_when(
      Genotype %in% c("0/0") ~ "Hom_REF",
      Genotype %in% c("0/1", "1/0") ~ "Het",
      Genotype %in% c("1/1") ~ "Hom_ALT",
      TRUE ~ "Missing"
    ),
    Pos_Label = paste0(hgnc_symbol, ":", Pos)
  ) %>%
  arrange(Pos)

# Define custom order of samples (top to bottom in the heatmap)
custom_order <- c("WG0091", "WG0094", "WG0347", "WG0600", "WG1537",
                  "WG0225", "WG0872", "WG0631", "WG0153", "WG0154",
                  "WG0161", "WG1726", "WG0867", "WG1322", "WG0596",
                  "WG0512", "WG1068", "WG0364", "WG0361", "WG0513",
                  "WG0367", "WG1094", "WG0878", "WG0158", "WG1078",
                  "WG0363", "WG0117", "WG1758", "WG0718", "WG0366")

# Apply to both geno_other and geno_chr19
geno_other$Sample <- factor(geno_other$Sample, levels = custom_order)
geno_chr19$Sample <- factor(geno_chr19$Sample, levels = custom_order)


ggplot(geno_chr19, aes(x = Pos_Label, y = Sample, fill = GenotypeCode)) +
  geom_tile() +
  scale_fill_manual(
    values = c("Het" = "yellow", "Hom_ALT" = "red", 
               "Hom_REF" = "green", "Missing" = "gray")
  ) +
  labs(
    title = "Genotype Heatmap: Chromosome 19 Haplotype Gene Only",
    x = "Genomic Position (Gene:BasePair)",
    y = "Sample"
  ) +
  scale_x_discrete(breaks = geno_chr19$Pos_Label[seq(1, nrow(geno_chr19), by = 750)]) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
    axis.text.y = element_text(size = 8)
  )
