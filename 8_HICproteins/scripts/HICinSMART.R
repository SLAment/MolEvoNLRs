#!/usr/bin/env Rscript

### HICinSMART: Proportion of HIC proteins with different repeat types in the SMART database
#############################################################################
# ===========================================================================
# Sandra Lorena Ament Velasquez
# 2025/03/04
# ++++++++++++++++++++++++++++++++++++++++++++++

# ============================
# Load the necessary libraries
# ============================
library(dplyr)
library(ggplot2)
library(ggtext) # for element_markdown()
library(epitools) # for oddsratio()
library(ggrepel) # geom_text_repel()
library(cowplot) # plot_grid()
# ============================

# ============================
# Reading the data
# ============================
## Snakemake
treks <- read.table(snakemake@input$treks, header = F, quote = '', sep = "\t")
annotation <- read.table(snakemake@input$annotation, header = T, quote = '', sep = "\t")
PSIMmin <- snakemake@params$psimmin

# ============================
# Fix the data frame
# ============================
names(treks) <- c('Repeat', "Taxon", 'Group', 'Phylum', 'SeqID', 'Genus', 'Species', 'Repeat_len', 'Number_repeats', 'Psim', 'Range')

treks <- merge(treks, annotation %>% dplyr::select(!Group), by = c("SeqID", "Repeat")) %>% filter(Genus != "Quercus") # Quercus often has contaminations from endophytes

# Replace "OtherFungi" with "Other Fungi" in Group labels
treks <- treks %>%
  mutate(Group = recode(Group, "OtherFungi" = "Other Fungi")) %>% 
  mutate(Group = recode(Group, "OtherEukaryotes" = "Other Eukaryotes"))

# Create a factor for Group ordered by Taxon
treks <- treks %>%
  arrange(Taxon, Group) %>%
  mutate(Group = factor(Group, levels = unique(Group)))

# Create a factor for Repeat ordered to have LRR at the end
treks <- treks %>% mutate(Repeat = factor(Repeat, levels = c("ANK", "TPR", "WD40", "LRR") ))

# Define custom colors for each Taxon
taxon_colors <- c(
  "Animals" = "red",
  "Archaea" = "gold2",
  "Viridiplantae" = "forestgreen",
  "Fungi" = "darkorange",
  "Bacteria" = "steelblue",
  "Other Eukaryotes" = "brown"  # Add any other taxonomic groups here
)

# HIC bar fill colors
hic_colors <- c("with_HIC" = "black", "without_HIC" = "azure3")

# For statistics later
# Create a binary HIC variable
treks <- treks %>%
  mutate(HIC = Number_repeats > 0 & Psim >= PSIMmin)

# Simplify ATPase to "ATPase" vs "None"
treks <- treks %>%
  mutate(ATPase_class = ifelse(ATPase == "None", "None", "With ATPase"))

# Exclude proteins with an extra repeat domain, since those might be false
# positives for their respective main repeats Some proteins have repeated
# entries but the second repeat was not properly annotated, so it looks as if
# it had a single repeat type.
noextras <- treks %>% filter(ExtraRepeats == "None") %>% .$SeqID %>% unique()
extras <- treks %>% filter(ExtraRepeats != "None") %>% .$SeqID %>% unique() 
trickyones <- intersect(noextras,extras)

length(extras) + length(trickyones) # 815
length(unique(treks$SeqID)) # 160477

treks <- treks %>% filter(!SeqID %in% c(extras, trickyones))

# ============================
# ---- How many HIC proteins are there per taxonomic group? ----
# ============================

# Summarize into HIC and non-HIC counts
group_summary_long <- treks %>%
  group_by(Repeat, Taxon, Group, ATPase) %>%
  summarise(
    with_HIC = sum(Number_repeats > 0 & Psim >= PSIMmin),
    without_HIC = n() - with_HIC,
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(with_HIC, without_HIC), names_to = "HIC_status", values_to = "count")

# Reorder the Group by Taxon
group_summary_long <- group_summary_long %>%
  arrange(Taxon, Group) %>%
  mutate(Group = factor(Group, levels = unique(Group)))

# Create colored x-axis labels
group_labels <- treks %>%
  distinct(Group, Taxon) %>%
  mutate(label = paste0("<span style='color:", taxon_colors[Taxon], "'>", Group, "</span>"))

# Merge colored labels back in
group_summary_long <- group_summary_long %>%
  left_join(group_labels, by =  c("Group", "Taxon")) %>%
  mutate(label = factor(label, levels = unique(label)))

# Compute total proteins for label display
label_data <- group_summary_long %>%
  group_by(Repeat, Taxon, Group, ATPase, label) %>%
  summarise(total = sum(count), .groups = "drop")

### Perspective: Proportion of HIC proteins amongst all proteins for each repeat type and ATPase
(hicprop <- ggplot(group_summary_long, 
                   aes(x = label, y = count, fill = HIC_status)) +
    geom_bar(stat = "identity") +
    facet_grid(ATPase ~ Repeat, scales = "free") +
    scale_fill_manual(values = hic_colors, labels = c("With HIC", "Without HIC")) +
    ggrepel::geom_text_repel(
      data = label_data,
      aes(x = label, y = total, label = total),
      size = 2.3,
      inherit.aes = FALSE,
      direction = "y",          # only repel vertically (keeps things tidy)
      nudge_y = 30,           # push slightly upward
      # segment.color = NA        # remove connecting lines
      segment.color = "grey30",   # color of connecting lines
      segment.size = 0.3          # thinner lines (default is 0.5)
    ) +
    labs(
      y = "Number of proteins",
      fill = "HIC status"
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_markdown(angle = 45, hjust = 1),
      axis.title.x = element_blank(),
      plot.margin = margin(5, 30, 5, 5) # give extra space if labels overflow
    )
)

# prevent clipping of repelled labels
hicprop <- hicprop + coord_cartesian(clip = "off")

ggsave(plot = hicprop, 
       filename = snakemake@output$hicprop,
       width = 9.5, height = 7)

# ============================
### Perspective: How frequent are HIC repeats among WD, ANK and TPR repeats, where do they occur?
# ============================

# Summarize into HIC and non-HIC counts, without caring about the ATPase
group_summary_long_noATPase <- treks %>%
  group_by(Repeat, Taxon, Group) %>%
  summarise(
    with_HIC = sum(Number_repeats > 0 & Psim >= PSIMmin),
    without_HIC = n() - with_HIC,
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = c(with_HIC, without_HIC), names_to = "HIC_status", values_to = "count")

# Compute proportions of HIC proteins
label_data_pct <- group_summary_long_noATPase %>%
  tidyr::pivot_wider(names_from = HIC_status, values_from = count) %>%
  mutate(
    total = with_HIC + without_HIC,
    proportion_HIC = with_HIC / total,
    perc = sprintf("%.1f%%", proportion_HIC * 100)
  ) 
 
## Just the proportions of HIC so the coverage per group is not distracting
(totalprots <- ggplot(label_data_pct, aes(x = Group, y = proportion_HIC, fill = Taxon)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ Repeat, scales = "free") +
  scale_fill_manual(values = taxon_colors) +
  theme_classic() +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        strip.text = element_text(size = 10),
        legend.position = "none") +
  ylab("Percentage of proteins\nwith HIC") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  ggrepel::geom_text_repel(
    data = label_data_pct,
    aes(x = Group, y = proportion_HIC, label = perc),
    size = 2.5,
    nudge_y = 0.01,
    direction = "y",
    segment.size = 0.2) )

# ============================
### Perspective: What is the fraction of HIC ANK, TPR, WD repeats associated to ATPase domains? 
# ============================

# Filter for HIC proteins
hic_proteins <- treks %>%
  filter(Number_repeats > 0, Psim >= PSIMmin)

# Summarize
hic_atpase_summary <- hic_proteins %>%
  group_by(Repeat, Taxon, Group, ATPase) %>%
  summarise(count = n()) %>% #, .groups = "drop") %>%
  group_by(Repeat, Taxon, Group) %>%
  mutate(
    total_HIC = sum(count),
    proportion = count / total_HIC
  ) %>%
  ungroup()

# Replace the naming to make it less confusing
hic_atpase_summary <- hic_atpase_summary %>%
  mutate(ATPase = if_else(ATPase == "ATPase", "Other ATPases", ATPase))

HICprots <- ggplot(hic_atpase_summary %>% filter(ATPase != "None"), 
                   aes(x = Group, y = proportion, fill = Taxon)) +
  geom_bar(stat = "identity") +
  ggrepel::geom_text_repel(
    aes(label = sprintf("%.1f%%", proportion * 100)),
    size = 2.3,
    direction = "y",
    nudge_y = 0.01,  # vertical spacing
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  facet_grid(ATPase ~ Repeat) +
  scale_fill_manual(values = taxon_colors) +
  theme_bw() +
  scale_y_continuous(
    limits = c(0, 0.6),
    breaks = c(0, 0.25, 0.5),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ylab("Percentage of HIC proteins\nwith ATPase domains") +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.x= element_blank(),
        strip.text = element_text(size = 10),
        strip.background = element_rect(fill = "white"),
        panel.grid.major.y = element_blank(),
        legend.position="bottom",
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow = 1))

(allvshic <- cowplot::plot_grid(totalprots, HICprots, 
                                nrow=2, 
                                rel_heights = c(1,2),
                                labels=c('a', 'b'), 
                                align = "hv", 
                                axis = "lr") )

ggsave(plot = allvshic, 
       filename = snakemake@output$allvshic,
       width = 9.5, height = 8.7)

# ============================
# ---- Just the totals of HIC proteins with ATPase domain ----
# ============================

hic_atpase_nlrs <- hic_proteins %>% filter(HIC == TRUE) %>%
  group_by(Taxon, Group, ATPase_class) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Taxon, Group) %>%
  mutate(
    total_ATPase = sum(count),
    fracATPase = count / total_ATPase
  ) %>%
  ungroup()

( totalsATPaseHIC <- ggplot(hic_atpase_nlrs %>% filter(ATPase_class == "With ATPase"), 
       aes(x = Group, y = fracATPase, fill = Taxon)) +
  geom_bar(stat = "identity") +
  ggrepel::geom_text_repel(
    aes(label = sprintf("%.1f%%", fracATPase * 100)),
    size = 2.3,
    direction = "y",
    nudge_y = 0.01,  # vertical spacing
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  # facet_grid(ATPase ~ Repeat) +
  scale_fill_manual(values = taxon_colors) +
  theme_classic() +
  scale_y_continuous(
    # limits = c(0, 0.5),
    # breaks = c(0, 0.25, 0.5),
    labels = scales::percent_format(accuracy = 1)
  ) +
  ylab("Proportion of HIC proteins\nwith ATPase domains") +
  theme(axis.text.x = element_markdown(angle = 45, hjust = 1),
        axis.title.x= element_blank(),
        strip.text = element_text(size = 10),
        legend.position="bottom",
        legend.box = "horizontal") +
  guides(fill = guide_legend(nrow = 1)) )

ggsave(plot = totalsATPaseHIC, 
       filename = snakemake@output$totals,
       width = 8.5, height = 6)

# ============================
# ---- Statistics: Is HIC over-represented in NLRs? ----
# ============================

# Create contingency table
hic_table <- table(treks$HIC, treks$ATPase_class)
hic_table

hic_table[2,1]/sum(hic_table[,1])*100 # % HIC proteins without an ATPase domain
hic_table[2,2]/sum(hic_table[,2])*100 # % HIC proteins with an ATPase domain

# Chi-square test 
chisq.test(hic_table)
fisher.test(hic_table)


epitools::oddsratio(table(treks$ATPase_class, treks$HIC))
# Presence of an ATPase domain was strongly associated with the occurrence of
# HIC. Among proteins without an ATPase domain, 5.59% contained HIC, whereas
# 22.04% of proteins with an ATPase domain did. This corresponds to an odds ratio
# of 4.78 (95% CI: 4.51–5.06, Fisher’s exact test p < 2.2e-16).
