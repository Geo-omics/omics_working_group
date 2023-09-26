# reference sequence import

ref_path <- "mapping/reads/GLAMR_sxtAll_bwa/database/GLAMRsxtAll.fa"

ref_seq <- Biostrings::readDNAStringSet(ref_path) %>%
    as.character() %>%
    data.frame(ref_base = ., seqnames = names(.)) %>%
    mutate(seq_length = nchar(ref_base)) %>%
    separate_longer_position(ref_base, width = 1) %>%
    group_by(seqnames) %>%
    mutate(pos = row_number())

# An example BAM file path
bam_path = "/geomicro/data2/pdenuyl2/neurotoxin_thesis/FINAL1/mapping/reads/GLAMR_sxtAll_mm2/output/bam/samp_471_GLAMRsxtAll_mapped.bam"

# function to generate bam stats - reference cover & % id
bam_stats <- function(bam_path){
 
  bam <- Rsamtools::BamFile(file = bam_path,
                            index = paste0(bam_path,".bai"))

  out <- Rsamtools::pileup(bam,pileupParam = PileupParam(distinguish_strands = FALSE)) %>%
    full_join(., ref_seq, by = c("seqnames", "pos")) %>%  # join pileup file and alignment reference
    group_by(seqnames, pos, seq_length) %>%
    arrange(seqnames, pos, desc(count)) %>%
      mutate(depth = sum(count),                          # alignment depth at each position
             rel_abund_base = count/depth) %>%            # relative abundance of base read at position compared to all reads at position
    group_by(seqnames) %>%
      mutate(mean_depth = mean(depth),                    # mean depth of sequence
             depth_var = var(depth),
             depth_sd = sd(depth),
             bam_path = bam_path,
             sample_id = bam_path %>% str_remove(".*bam/") %>% str_remove("_GLAMRsxtAll.*")) %>%
    filter(nucleotide == ref_base) %>%                    # only look at bases that match reference
    group_by(seqnames) %>%
      mutate(percent_cover = (sum(count>0) / seq_length) * 100,      # number of bases matching reference divided by reference sequence length
      percent_id_ref_contig = (sum(count) / sum(depth)) * 100) %>%    #all reads matching ref / total reads for sequence          
    ungroup()
}

# create list of bam paths to process
bam_paths <- system("ls mapping/reads/GLAMR_sxtAll_mm2/output/bam/*.bam", intern = TRUE) %>%
  tibble(bam_path = .)

# run for bam stats!
#GLAMR_bam_stats_mm2 <- map_df(bam_paths$bam_path, ~ bam_stats(.x), .progress = TRUE) #RUN LAST - (date: Sept. 20, 2023)
#write_csv(GLAMR_bam_stats_mm2, file = "mapping/reads/GLAMR_sxtAll_mm2/output/GLAMR_bam_stats_mm2.csv")
GLAMR_bam_stats_mm2 <- read_csv(file = "mapping/reads/GLAMR_sxtAll_mm2/output/GLAMR_bam_stats_mm2.csv", col_names = TRUE)

# data table of samples with >50% cover for all three references
GLAMR_bam_stats_60_cov_mm2 <- GLAMR_bam_stats_mm2 %>%
    group_by(sample_id, seqnames) %>%
      #reframe(                                          # skim down data table to just data for each reference sequence
      #  sample_id = unique(sample_id),
      #  percent_cover = unique(percent_cover),
      #  percent_id_ref_contig = unique(percent_id_ref_contig)) %>%
  #select columns
      select(seqnames, sample_id, percent_cover, percent_id_ref_contig) %>% #check order of columns
      distinct() %>%
    group_by(sample_id) %>%
      mutate(gene_present = percent_cover > 60,                       # does reference sequence have > 60% coverage?
             all_present = all(gene_present)) %>%                     # do all three references have > 60% coverage?
    subset(all_present == TRUE)                                       # keep only samples that have >60% coverage for all three references

View(GLAMR_bam_stats_60_cov_mm2)

# data table of samples with >50% cover for all three references, but coverage from references is combined (one combined reference)
GLAMR_bam_stats_60_cov_comb_ref <- GLAMR_bam_stats_mm2 %>%
    subset(sample_id %in% unique(GLAMR_bam_stats_60_cov_mm2$sample_id)) %>%    # subset only samples that have >50% cover for all three references
    group_by(sample_id) %>%
      mutate(seq_length_comb_ref = sum(unique(seq_length)),                        # sum unique reference lengths (each unique sequence) to get length of combined reference
             percent_cover_comb_ref = (sum(count>0) / seq_length_comb_ref) * 100,                      # number of bases matching combined reference divided by combined reference sequence length
             percent_id_comb_ref = (sum(count) / sum(depth)) * 100) %>%            # all reads matching combined ref / total reads for combined ref sequence
      reframe(                                                                     # skim down data table to just data for each sample
        sample_id = unique(sample_id),
        percent_cover_comb_ref = unique(percent_cover_comb_ref),
        percent_id_comb_ref  = unique( percent_id_comb_ref))
 
View(GLAMR_bam_stats_60_cov_comb_ref)
#Summary: 22 samples, all >95% combined ref ID - September 20, 2023