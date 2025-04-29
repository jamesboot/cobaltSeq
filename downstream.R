# Processing of random barcode count tables for COBALTse

# Libraries
library(dplyr)
library(readxl)
library(ggplot2)
library(viridis)
library(readr)
library(networkD3)
library(tidyr)
library(webshot)
library(ComplexHeatmap)
library(tibble)

# Outdir
outdir <- '/nemo/stp/babs/working/bootj/projects/swantonc/eva.gongross/eg879'
setwd(outdir)

# Read in file
dat <-
  read_csv(
    paste0(outdir, '/allcounts.csv')
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

saveRDS(dat_score_anno, 'diversity_score_df.RDS')

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
gRNA_counts <- dat %>%
  group_by(sample_id, gRNA) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  mutate(logTotalCount = log10(total_count))

# Annotate
gRNA_counts_anno <- merge(gRNA_counts, meta, by = 'sample_id')

# Heatmap 

# Pivot wider
hmdf <- gRNA_counts_anno %>%
  select(sample_id, gRNA, logTotalCount) %>%
  pivot_wider(names_from = gRNA, values_from = logTotalCount) %>%
  column_to_rownames('sample_id')

# Convert to matrix
mat <- as.matrix(hmdf)
rownames(meta) <- meta$sample_id

# Prep colours for annotation
colrs1 <- c('#E6194B', '#3CB44B', '#4363D8', '#F58231', '#911EB4', '#42D4F4', '#BFEF45','#F032E6', '#469990','#FFE119','#FABED4','#9AACEF', '#A52A2A', '#808080','#FFD8B1', '#AA6E28')
names(colrs1) <- c(
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

colrs2 <- magma(7)
names(colrs2) <- unique(meta[rownames(mat),]$`Sample Genotype`)

# Create annotation
ha <- HeatmapAnnotation(df = data.frame(Treatment = meta[rownames(mat),]$`Sample Treatment`,
                                        Genotype = meta[rownames(mat),]$`Sample Genotype`),
                        col = list(Treatment = colrs1,
                                   Genotype = colrs2),
                        which = 'row')

# Heatmap
hm <- Heatmap(
  mat,
  na_col = "grey",
  col = viridis(100),
  row_split = meta[rownames(mat), ]$`Sample Treatment`,
  right_annotation = ha,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 6),
  show_row_names = F,
  row_gap = unit(1, 'mm'),
  column_names_gp = gpar(fontsize = 8)
)

# Save 
png('heatmap1.png',
    width = 15,
    height = 15,
    units = 'in',
    res = 300)
draw(hm)
dev.off()

# Try plotting percentage
gRNA_counts <- dat %>%
  group_by(sample_id, gRNA) %>%
  summarise(total_count = sum(count)) %>%
  arrange(desc(total_count)) %>%
  group_by(sample_id) %>%
  mutate(sample_counts = sum(total_count)) %>%
  mutate(percentage_counts = log10(total_count/sample_counts))
  
# Annotate
gRNA_counts_anno <- merge(gRNA_counts, meta, by = 'sample_id')

# Save for markdown report
saveRDS(gRNA_counts_anno, 'heatmapDF.RDS')

# Heatmap 

# Pivot wider
hmdf <- gRNA_counts_anno %>%
  select(sample_id, gRNA, percentage_counts) %>%
  pivot_wider(names_from = gRNA, values_from = percentage_counts) %>%
  column_to_rownames('sample_id')

# Convert to matrix
mat <- as.matrix(hmdf)
rownames(meta) <- meta$sample_id

# Prep colours for annotation
colrs1 <- c('#E6194B', '#3CB44B', '#4363D8', '#F58231', '#911EB4', '#42D4F4', '#BFEF45','#F032E6', '#469990','#FFE119','#FABED4','#9AACEF', '#A52A2A', '#808080','#FFD8B1', '#AA6E28')
names(colrs1) <- c(
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

colrs2 <- magma(7)
names(colrs2) <- unique(meta[rownames(mat),]$`Sample Genotype`)

# Create annotation
ha <- HeatmapAnnotation(df = data.frame(Treatment = meta[rownames(mat),]$`Sample Treatment`,
                                        Genotype = meta[rownames(mat),]$`Sample Genotype`),
                        col = list(Treatment = colrs1,
                                   Genotype = colrs2),
                        which = 'row')

# Heatmap
hm <- Heatmap(
  mat,
  na_col = "grey",
  col = viridis(100),
  row_split = meta[rownames(mat), ]$`Sample Treatment`,
  right_annotation = ha,
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 6),
  show_row_names = F,
  row_gap = unit(1, 'mm'),
  column_names_gp = gpar(fontsize = 8),
  width = unit(10, "in"), 
  height = unit(10, "in")
)

saveRDS(hm, 'heatmapObj.RDS')



# Save 
png('heatmap2.png',
    width = 15,
    height = 15,
    units = 'in',
    res = 300)
draw(hm)
dev.off()



# PART 2 ----
# Assess overall QC of reads - what happens to them all?
# Sankey plot to visualise

# Load gRNA counts CSV
gRNA_qc <- read.csv(paste0(outdir, '/all_gRNA_counts.csv'))

# Create a table of total read count per sample
total_reads <- gRNA_qc %>%
  select(sample_id, total_reads_R1) %>%
  group_by(sample_id) %>%
  distinct() %>%
  summarise(
    total_reads = sum(total_reads_R1)
  )

# Create a table of total number of each class per sample
gRNA_reads <- gRNA_qc %>%
  select(sample_id, valid_gRNA:rand_bc_too_short) %>%
  group_by(sample_id) %>%
  summarise(
    total_valid_gRNA_reads = sum(valid_gRNA),
    total_valid_rand_bc = sum(valid_rand_bc),
    total_rand_bc_seq_error = sum(rand_bc_seq_error),
    total_rand_bc_too_short = sum(rand_bc_too_short)
  )

# Merge above tables and get number of invalid gRNA reads
read_report <- inner_join(total_reads, gRNA_reads, by = 'sample_id') 
read_report <- as.data.frame(read_report)

read_report$total_reads <- as.numeric(as.character(read_report$total_reads))
read_report$total_valid_gRNA_reads <- as.numeric(as.character(read_report$total_valid_gRNA_reads))

read_report <- read_report %>%
  mutate(invalid_gRNA = total_reads - total_valid_gRNA_reads)

# Make a plot for all samples
all_report <- read_report %>%
  summarise(across(total_reads:invalid_gRNA, sum))

# Create object for plotting
sankeyList <- list(nodes = data.frame(
  name = c(
    'total_reads',
    'valid_gRNA',
    'invalid_gRNA',
    "valid_rand_bc",
    "rand_bc_seq_error",
    "rand_bc_too_short"
  )
),
links = data.frame(
  source = c(0, 0, 1, 1, 1),
  target = c(1, 2, 3, 4, 5),
  value = c(
    all_report$total_valid_gRNA_reads,
    all_report$invalid_gRNA,
    all_report$total_valid_rand_bc,
    all_report$total_rand_bc_seq_error,
    all_report$total_rand_bc_too_short
  )
))

# Plot 
sankeyNetwork(Links = sankeyList$links,
              Nodes = sankeyList$nodes,
              Source = 'source',
              Target = 'target',
              Value = 'value',
              NodeID = 'name',
              fontSize = 16,
              fontFamily = 'Arial')

# Make plot per sample
# Create directory for images
dir.create(paste0(outdir, '/figures'))
for (sample in unique(total_reads$sample_id)) {
  # Filter tables to sample of interest
  soi <- sample
  message(paste('Starting sample:', soi))
  sample_read_report <- read_report %>%
    filter(sample_id == soi)
  
  # Create object for plotting
  sankeyList <- list(nodes = data.frame(
    name = c(
      'total_reads',
      'valid_gRNA',
      'invalid_gRNA',
      "valid_rand_bc",
      "rand_bc_seq_error",
      "rand_bc_too_short"
    )
  ),
  links = data.frame(
    source = c(0, 0, 1, 1, 1),
    target = c(1, 2, 3, 4, 5),
    value = c(
      sample_read_report$total_valid_gRNA_reads,
      sample_read_report$invalid_gRNA,
      sample_read_report$total_valid_rand_bc,
      sample_read_report$total_rand_bc_seq_error,
      sample_read_report$total_rand_bc_too_short
    )
  ))
  
  # Plot
  sn <- sankeyNetwork(
    Links = sankeyList$links,
    Nodes = sankeyList$nodes,
    Source = 'source',
    Target = 'target',
    Value = 'value',
    NodeID = 'name',
    fontSize = 16,
    fontFamily = 'Arial'
  )

  # Save
  saveNetwork(sn, paste0(outdir, '/figures/', soi, ".html"))
}

# # Also do a plot that shows all the gRNAs
# sample_gRNA_qc <- gRNA_qc %>%
#   filter(sample_id == soi) %>%
#   group_by(sample_id, barcode_name) %>%
#   summarise(
#     valid_gRNA = sum(valid_gRNA),
#     valid_rand_bc = sum(valid_rand_bc),
#     rand_bc_seq_error = sum(rand_bc_seq_error),
#     rand_bc_too_short = sum(rand_bc_too_short)
#   )
# 
# # Pivot wider
# valid_gRNAs <- sample_gRNA_qc %>%
#   pivot_wider(id_cols = sample_id, names_from = barcode_name, values_from = valid_gRNA)
# 
# valid_rand_bc <- sample_gRNA_qc %>%
#   pivot_wider(id_cols = sample_id, names_from = barcode_name, values_from = valid_rand_bc)
# 
# rand_bc_seq_error <- sample_gRNA_qc %>%
#   pivot_wider(id_cols = sample_id, names_from = barcode_name, values_from = rand_bc_seq_error)
# 
# rand_bc_too_short <- sample_gRNA_qc %>%
#   pivot_wider(id_cols = sample_id, names_from = barcode_name, values_from = rand_bc_too_short)
# 
# # Create object for plotting
# n_gRNA <- length(colnames(valid_gRNAs)[-1])
# sankeyList <- list(nodes = data.frame(
#   name = c(
#     'total_reads',
#     colnames(valid_gRNAs)[-1],
#     'invalid_gRNA',
#     "valid_rand_bc",
#     "rand_bc_seq_error",
#     "rand_bc_too_short"
#   )
# ),
# links = data.frame(
#   source = c(0, rep(0, n_gRNA), 1:n_gRNA, 1:n_gRNA, 1:n_gRNA),
#   target = c(2, 1:n_gRNA, rep(3, n_gRNA), rep(4, n_gRNA), rep(5, n_gRNA)),
#   value = c(
#     sample_read_report$invalid_gRNA,
#     as.numeric(as.vector(valid_gRNAs[1, -1])),
#     as.numeric(as.vector(valid_rand_bc[1, -1])),
#     as.numeric(as.vector(rand_bc_seq_error[1, -1])),
#     as.numeric(as.vector(rand_bc_too_short[1, -1]))
#   )
# ))
# 
# # Plot 
# sankeyNetwork(Links = sankeyList$links,
#               Nodes = sankeyList$nodes,
#               Source = 'source',
#               Target = 'target',
#               Value = 'value',
#               NodeID = 'name')
