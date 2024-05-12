library(miplicorn)
library(tidyverse)
library(here)
library(flextable)
library(binom)
extract_sampleID <- function(x) {
  # split all three by delimiter
  splt <- str_split_fixed(x,'-',4)
  # make a dataframe
  splt_ID <- as.data.frame(splt) %>%
    mutate(ID = if_else(
      V1 == 'UNC',.[[2]],
      .[[1]]))
  return(splt_ID$ID)
}

# Make a list of amino acids
amino_acids <- vector(mode="list", length=18)
names(amino_acids) <- c("Val", "Ala", "His","Leu","Arg","Pro","Lys","Thr","Cys","Met","Tyr","Ile","Asn","Glu","Ser","Gly","Asp","Phe")
amino_acids[[1]] <- "V"; amino_acids[[2]] <- "A"; amino_acids[[3]] <- "H"; amino_acids[[4]] <- "L"; amino_acids[[5]] <- "R"; amino_acids[[6]] <- "P"; amino_acids[[7]] <- "K"; amino_acids[[8]] <- "T"; amino_acids[[9]] <- "C"; amino_acids[[10]] <- "M"; amino_acids[[11]] <- "Y"; amino_acids[[12]] <- "I"; amino_acids[[13]] <- "N"; amino_acids[[14]] <- "E"; amino_acids[[15]] <- "S"; amino_acids[[16]] <- "G"; amino_acids[[17]] <- "D"; amino_acids[[18]] <- "F"

here::i_am("R/miplicorn/01_prevalence_table.R")

genotype_df <- vroom::vroom(
  here("data/tables/genotypes_AA_table.csv"),
  col_names = FALSE,
  show_col_types = FALSE
) %>%
  # Take the transpose of our matrix, making rows columns and columns rows.
  # This will allows us to keep all the data in our .csv file.
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(-rowname) %>%
  tidyr::pivot_wider(
    names_from = rowname,
    values_from = value
  ) %>%
  # Assign the column names of our tibble and clean them up
  dplyr::select(-name) %>%
  janitor::row_to_names(1) %>%
  # Convert our data to a long format
  tidyr::pivot_longer(
    cols = -c(1:6),
    names_to = "sample",
    values_to = "value"
  ) %>%
  janitor::clean_names() %>%
  dplyr::relocate(sample) %>%
  dplyr::rename('genotype_vals' := value)

genotype_df$sample <- extract_sampleID(genotype_df$sample)

metadata<-read_tsv(here("data/metadata/TAN_metadata.tsv"))

genotype_df <- dplyr::left_join(genotype_df, metadata, by = c("sample" = "sample_ID"))

genotype_df <- genotype_df %>%
  mutate(region = if_else(str_detect(Sample_Set, "ZAN") == TRUE,"ZAN","TAN"))
genotype_df$genotype_vals <- as.numeric(genotype_df$genotype_vals)

for(i in names(amino_acids)) {
  genotype_df$mutation_name <- str_replace_all(genotype_df$mutation_name,i,amino_acids[[i]])
}

to_filter <- c("YB5004","YB5005","YB5006","YB5007","YB5042","YB5181","YB5187","YB5196","YB5198","YB5209","YB5252","YB5257","YB5316","YB5451")
genotype_df <- genotype_df %>%
  filter(!sample %in% to_filter)

write_csv(genotype_df,"results/genotype_df.csv")

# get_wilson_CI <- function(p_hat,n,alpha = 0.05) {
#   #-----------------------------------------------------------------------------
#   # Compute the Wilson (aka Score) confidence interval for a popn. proportion
#   # \text{CI lower bound} = \frac{p + \frac{z^2}{2n} - z \sqrt{\frac{p(1 - p)}{n} + \frac{z^2}{4n^2}}}{1 + \frac{z^2}{n}}
#   # \text{CI upper bound} = \frac{p + \frac{z^2}{2n} + z \sqrt{\frac{p(1 - p)}{n} + \frac{z^2}{4n^2}}}{1 + \frac{z^2}{n}}
#   #-----------------------------------------------------------------------------
#   # p_hat calculated prevalence
#   # n sample size
#   # alpha    1 - (confidence level)
#   #-----------------------------------------------------------------------------
#   z <- qnorm(1 - alpha / 2)
#   omega <- n / (n + z^2)
#   A <- p_hat + z^2 / (2 * n)
#   B <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
#   lower <- omega * (A - B)
#   upper <- omega * (A + B)
#
#   # Ensure CI bounds are within [0, 1]
#   lower <- max(lower, 0)
#   upper <- min(upper, 1)
#
#   CI <- c('lower' = lower, 'upper' = upper)
#   return(CI)
# }

prev_table_TAN <- genotype_df %>%
  filter(region == "TAN") %>%
  group_by(mutation_name) %>%
  summarise(total_ALT_M = sum(genotype_vals >= 1, na.rm = TRUE),
            total_REF_M = sum(genotype_vals == 0, na.rm = TRUE),
            N_M = total_ALT_M + total_REF_M,
            prevalence_M = total_ALT_M/N_M) %>%
  mutate(
    CI = pmap(list(total_ALT_M, N_M), ~ binom.confint(x = ..1, n = ..2, methods = "exact")),
    ci.lower_M = map_dbl(CI, 'lower'),
    ci.upper_M = map_dbl(CI, 'upper')
  ) %>%
  select(-CI)

prev_table_ZAN <- genotype_df %>%
  filter(region == "ZAN") %>%
  group_by(mutation_name) %>%
  summarise(total_ALT_Z = sum(genotype_vals >= 1, na.rm = TRUE),
            total_REF_Z = sum(genotype_vals == 0, na.rm = TRUE),
            N_Z = total_ALT_Z + total_REF_Z,
            prevalence_Z = total_ALT_Z/N_Z) %>%
  mutate(
    CI = pmap(list(total_ALT_Z, N_Z), ~ binom.confint(x = ..1, n = ..2, methods = "exact")),
    ci.lower_Z = map_dbl(CI, 'lower'),
    ci.upper_Z = map_dbl(CI, 'upper')
  ) %>%
  select(-CI)

prev_main_zan_tbl <- left_join(prev_table_ZAN,
                               prev_table_TAN,
                               by = "mutation_name") %>%
  mutate_if(is.numeric, round,3) %>%
  unite(CI_Z,ci.lower_Z,ci.upper_Z,remove = T, sep = "-") %>%
  unite(CI_M,ci.lower_M,ci.upper_M,remove = T, sep = "-")

mutations_for_tbl <- c("crt-K76T","crt-M74I","crt-N75E", "dhfr-ts-S108N","dhfr-ts-N51I","dhfr-ts-C59R","dhps-A437G","dhfr-ts-I164L","dhfr-ts-S108T","dhps-A581G", "dhps-K540E", "k13-K189N", "k13-K189T", "mdr1-N86Y", "mdr1-D1246Y", "mdr1-Y184F","mdr2-I492V","dhfr-ts-A16V")

prev_table_filtered <- prev_main_zan_tbl %>%
  filter(mutation_name %in% mutations_for_tbl)

# ggplot(prev_table_filtered, aes(x = mutation_name, y = prevalence, fill = region)) +
#   geom_bar(stat = 'identity',position=position_dodge(), width =.4)+
#   coord_cartesian(ylim=c(0,1)) +
#   theme_classic() +
#   labs(y= "Prevalence",
#        #title = "Prevalence Based On Region"
#   ) +
#   guides(fill=guide_legend(title="Region")) +
#   geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.5,
#                 position=position_dodge(.4)) +
#   scale_fill_manual(values = c("#2538be","#Bb141d")) +
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.title.x = element_blank(),
#         axis.text = element_text(size = 25),
#         axis.title = element_text(size = 35),
#         legend.text = element_text(size = 25),
#         legend.title = element_text(size = 25))

border_style = officer::fp_border(color="black", width=3)
table_sub <-
  prev_table_filtered %>%
  select(-c("total_ALT_Z","total_REF_Z","total_ALT_M","total_REF_M")) %>%
  select(mutation_name,prevalence_Z,CI_Z,N_Z,prevalence_M,CI_M,N_M) %>%
  mutate(mutation_name = if_else(
    str_detect(mutation_name,"dhfr") == TRUE,
    paste("dhfr-",str_extract(mutation_name, "[^-]+$"),sep = ""),
    mutation_name)) %>%
  # order alphabetically
  arrange(stringr::str_extract(mutation_name, ".+?(?=-)"),
          # then order by the numeric in the amino acid sequence
          as.integer(stringr::str_extract(gsub("^.*-", "", mutation_name), "\\d+"))) %>%
  mutate(mutation_name = paste("Pf",mutation_name,sep = "")) %>%
  # mutate(
  #   prevalence_Z = case_when(
  #     prevalence_Z < 0.005 ~ 'ND†',
  #     TRUE ~ as.character(prevalence_Z)
  #   ),
  #   prevalence_M = case_when(
  #     prevalence_M < 0.005 ~ 'ND†',
  #     TRUE ~ as.character(prevalence_M)
  #   ),
  #   CI_Z = case_when(
  #     prevalence_Z == 'ND†' ~ '-',
  #     TRUE ~ as.character(CI_Z)
  #   ),
  #   CI_M = case_when(
  #     prevalence_M == 'ND†' ~ '-',
  #     TRUE ~ as.character(CI_M)
  #   )
  #   ) %>%
  flextable::flextable() %>%
  flextable::autofit() %>%
  flextable::width(j=1, width = 2) %>%
  flextable::width(j=c(2,4), width = 1) %>%
  flextable::width(j=c(3,5), width = 1.5) %>%
  flextable::add_header_row(
    top = TRUE,
    values = c("Mutation",
               "Zanzibar",
               "",
               "",
               "Mainland",
               "",
               "")) %>%
  flextable::set_header_labels(
    mutation_name = "",
    prevalence_Z = "Mutant Allele Prevalence‡",
    CI_Z = "CI‡",
    N_Z = "# Genotyped Samples§",
    prevalence_M = "Mutant Allele Prevalence‡",
    CI_M = "CI‡",
    N_M = "# Genotyped Samples§") %>%
  # merge mutation name header
  flextable::merge_at(i = 1:2, j = 1, part = "header") %>%
  # merge Zanzibar part
  flextable::merge_at(i = 1, j = 2:4, part = "header") %>%
  # merge Mainland part
  flextable::merge_at(i = 1, j = 5:7, part = "header") %>%
  # Remove all existing borders
  flextable::border_remove() %>%
  # add  box theme
  flextable::theme_box() %>%
  # change borders
  flextable::border(border = border_style) %>%
  flextable::border(border = border_style, part = "header") %>%
  flextable::align(align = "center", j = c(2:7), part = "all")
flextable::save_as_image(table_sub, path = "results/selected_mts_prev.pdf",webshot = "webshot2")
flextable::save_as_docx(table_sub, path = "results/selected_mts_prev.docx")
