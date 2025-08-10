library(ggplot2)
library(dplyr)

# Minimum length filter in base pairs
min_length <- 550

# Clean WG0718 and WG1094 data frames (remove NA and empty)
WG0718_clean <- WG0718 %>%
  filter_all(all_vars(!is.na(.) & . != "")) %>%
  setNames(c("chr", "start", "end"))

WG1094_clean <- WG1094 %>%
  filter_all(all_vars(!is.na(.) & . != "")) %>%
  setNames(c("chr", "start", "end"))

# Combine both datasets with sample labels
combined <- bind_rows(
  WG0718_clean %>% mutate(sample = "WG0718"),
  WG1094_clean %>% mutate(sample = "WG1094")
) %>%
  # Filter segments by minimum length
  filter((end - start) >= min_length) %>%
  mutate(
    y = as.numeric(factor(sample))  # y position: WG0718=1, WG1094=2
  )

# NARS2 gene region with 50kb flanking buffer
nars2_gene <- data.frame(
  gene = "NARS2",
  chr = "11",
  start = 78435620,
  end = 78575194,
  buffer = 50000
) %>%
  mutate(
    plot_start = start - buffer,
    plot_end = end + buffer
  )

# Plot segments as rectangles with minimal clutter
ggplot() +
  geom_rect(data = nars2_gene,
            aes(xmin = plot_start, xmax = plot_end, ymin = 0.5, ymax = 2.5),
            fill = "lightblue", alpha = 0.15) +
  geom_rect(data = nars2_gene,
            aes(xmin = start, xmax = end, ymin = 1.75, ymax = 2.25),
            fill = NA, color = "dodgerblue3", size = 0.7) +
  geom_rect(data = combined,
            aes(xmin = start, xmax = end, ymin = y - 0.2, ymax = y + 0.2, fill = sample),
            color = "black", linewidth = 0.15, show.legend = FALSE) +
  scale_y_continuous(breaks = 1:2, labels = c("WG0718", "WG1094"), limits = c(0.5, 2.5)) +
  scale_fill_manual(values = c("WG0718" = "#FDB863", "WG1094" = "#648FFF")) +
  labs(x = "Genomic Position (chr11)",
       y = "",
       title = paste0("Haplotype Regions Near NARS2 (Segments >= ", min_length/1000, " kb)")) +
  theme_minimal(base_size = 14) +
  theme(axis.text.y = element_text(size = 13, face = "bold"),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

