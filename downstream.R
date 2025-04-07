# Processing of random barcode count tables for COBALTse

# Libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(viridis)
library(readr)

# Outdir
outdir <- '/nemo/stp/babs/working/bootj/projects/swantonc/eva.gongross/eg879'

# Read in file
dat <-
  read_csv(
    paste0(outdir, '/allcounts.csv'),
    col_names = c(
      'rand_bc',
      'count',
      'sgRNA',
      'sample_id',
      'illumina_id',
      'lane_id',
      'run_id'
    )
  )

# Calculate diversity index for each sample
# Add ID col and calculate diversity index 
dat_score <- dat %>%
  group_by(sample_id) %>%
  mutate(Sum = sum(count)) %>%
  mutate(Pi = count/Sum) %>%
  mutate(logPi = log(Pi)) %>%
  mutate(Pi_x_logPi = Pi*logPi) %>%
  summarise(diversity_index = sum(Pi_x_logPi)*-1)

# Annotate 
meta <- read_xlsx('/nemo/stp/babs/working/bootj/github/cobaltSeq/DN20060_all_samples.xlsx',
                  skip = 1)
colnames(meta)[1] <- colnames(dat_score)[1]
meta$`Sample Treatment`[meta$`Sample Treatment` == 'A549_Invitro'] <- 'A549_invitro'
dat_score_anno <- merge(dat_score, meta, by = 'sample_id')

# Specify levels of treatment
dat_score_anno$`Sample Treatment` <- factor(
  dat_score_anno$`Sample Treatment`,
  levels = c(
    "H358_T0",
    "H358_Lung",
    "A549_T0",
    "A549_Lung",
    "A549_invitro",
    "H441_T0",
    "H441_Lung",
    "H2122_T0",
    "H2122_Lung",
    "IV, A549",
    "IV, H358",
    "Sub-cut, A549",
    "T0_cells, pre-inf_mix, A549",
    "T0_cells, post-inf_mix, A549",
    "T0_cells, pre-inf_mix, H358",
    "T0_cells, post-inf_mix, H358"
  )
)

# Visualise some results
p <- ggplot(data = dat_score_anno,
       aes(x = `Sample Treatment`, y = diversity_index, fill = `Sample Treatment`)) +
  geom_boxplot(outliers = F) +
  scale_fill_viridis(discrete = TRUE, alpha = 0.6) +
  geom_jitter(color = "black",
              size = 1,
              alpha = 0.9,
              position = position_jitter(width=0.2, height=0.1)) +
  theme_bw() +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  ),
  legend.position = 'none')
ggsave(plot = p,
       filename = paste0(outdir, '/testfig1.png'),
       height = 5,
       width = 5,
       dpi = 300)
  
# Find most abundant sgRNAs (regardless of clones) within each sample
# Sum counts for sgRNA
sgRNA_counts <- dat %>%
  group_by(sample_id, sgRNA) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  mutate(logTotalCount = log10(total_count))

# Annotate
sgRNA_counts_anno <- merge(sgRNA_counts, meta, by = 'sample_id')

# Summary plot
hm <- ggplot(sgRNA_counts_anno, aes(x = sgRNA, y = sample_id, fill = logTotalCount)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 4),
        axis.text.y = element_text(size = 6)) +
  facet_grid(~`Sample Treatment`)

ggsave(plot = hm,
       filename = paste0(outdir, '/hm1.png'),
       height = 25,
       width = 40,
       dpi = 300)


