# Clean datasets for analysis

source(here::here("code/00_config.R"))


# ----- Load data -----
wwscan_2025 <- read_csv(here("data/original/WWSCAN_AllPlants_SDR_2025.csv"))

casper_ra_files <- dir_ls(
  here("data/original/casper_export_20260604"),
  recurse = TRUE,
  glob = "*relative_abundance.csv"
)

epic_sheets <- excel_sheets(here("data/original/all_epic.xlsx"))
epic_dfs <- epic_sheets %>% 
  set_names() %>% 
  map(~ read_xlsx(here("data/original/all_epic.xlsx"), sheet = .x))


# ----- Clean dfs -----
## WWSCAN -----
colnames(wwscan_2025)

wwscan_2025 <- wwscan_2025 %>%
  rename(
    HAdV_F_gc_g_dry_weight     = HAdV_F_gc_g_dry_weigh,
    HAdV_F_gc_g_dry_weight_LCI = HAdV_F_gc_g_dry_weigh_LCI,
    HAdV_F_gc_g_dry_weight_UCI = HAdV_F_gc_g_dry_weigh_UCI,
    C_auris_gc_g_dry_weight_UCI = `C_auris_gc_g_dry_weight_ UCI`,
    InfA_H1_gc_g_dry_weight = `InfA_H1_gc_g_dry_weigh`
  )

## CASPER -----
# Col names are consistent; add variable identifying location associated with data
casper_combined <- casper_ra_files %>%
  map(~ read_csv(.x) %>%
        mutate(location = basename(dirname(.x)))) %>%
  list_rbind()

# Add a col for state name/abbr and lab
casper_combined <- casper_combined %>%
  mutate(
    state_abbr = str_extract(location, "(?<=, ).{2}"),
    state_name = state.name[match(state_abbr, state.abb)],
    lab_id = ifelse(str_detect(location, "SB"),"SB","MU")
  )

# Merge population estimates from Table 1 of CASPER preprint
casper_combined <- casper_combined %>%
  mutate(
    population = case_when(
      str_detect(location, "Boston DITP North, MA") ~ 1100000,
      str_detect(location, "Boston DITP South, MA") ~ 1400000,
      str_detect(location, "NYC") ~ NA,
      location == "Chicago (CHI-A), IL" ~ 1300000,
      location == "Chicago (CHI-B), IL" ~ 1100000,
      location == "Chicago (CHI-C), IL" ~ NA,
      location == "Chicago (CHI-D1), IL" ~ 1100000,
      location == "Chicago (CHI-D2), IL" ~ 1100000,
      location == "Chicago (CHI-D3), IL" ~ 2200000,
      location == "Ottumwa WPCF, IA" ~ 26000,
      location == "Columbia WWTP, MO" ~ 120000,
      location == "Kansas City Blue River WWTP, MO" ~ 490000,
      location == "Kansas City Westside WWTP, MO" ~ 61000,
      location == "Milan WWTF, MO" ~ 1800,
      location == "Monett WWTP, MO" ~ 9100,
      str_detect(location, "Miami-Dade CDWWTP, FL") ~ 900000,
      location == "Central Oklahoma (OK-A), OK" ~ 110000,
      location == "Central Oklahoma (OK-B), OK" ~ 440000,
      location == "Ontario, CA" ~ 890000,
      location == "Palo Alto RWQCP, CA" ~ 240000,
      location == "Riverside WQCP, CA" ~ 350000,
      location == "Sacramento, CA" ~ 1600000,
      str_detect(location, "Southern California, CA") ~ NA,
      location == "Boise RWF, ID" ~ 270000,
      T ~ NA # any sites since preprint posting are NA 
    )
  )

# Convert clade counts to wide format
casper_clade_counts_wide <- casper_combined %>%
  select(
    location,
    state_abbr,
    state_name,
    lab_id,
    population,
    date,
    label,
    total_read_pairs,
    n_reads_profiled,
    pmmov_counts,
    pmmov_counts_scaled,
    clade_counts
  ) %>%
  pivot_wider(
    .,
    names_from = label,
    values_from = clade_counts,
    names_glue = "{label}_{.value}",
  )

colnames(casper_clade_counts_wide) <- gsub(" ", "_", colnames(casper_clade_counts_wide))
colnames(casper_clade_counts_wide) <- gsub("-", "_", colnames(casper_clade_counts_wide))
colnames(casper_clade_counts_wide) <- gsub("\\(", "", colnames(casper_clade_counts_wide))
colnames(casper_clade_counts_wide) <- gsub("\\)", "", colnames(casper_clade_counts_wide))

# Aggregate DITP N and S
# (keep MU and SB separate)
DITP_MU_data <- casper_clade_counts_wide %>% 
  filter(location == "Boston DITP North, MA (MU)" |
           location == "Boston DITP South, MA (MU)") %>% 
  group_by(date) %>% 
  summarise(n_samples = n())
# ^N and S have same sampling dates for MU

DITP_SB_data <- casper_clade_counts_wide %>% 
  filter(location == "Boston DITP North, MA (SB)" |
           location == "Boston DITP South, MA (SB)") %>% 
  group_by(date) %>% 
  summarise(n_samples = n())
# ^N and S have same sampling dates for SB

DITP_data <- casper_clade_counts_wide %>%
  filter(grepl("DITP", location)) %>%
  group_by(state_abbr, state_name, lab_id, date) %>%
  summarise(across(
    ends_with("_counts") |
      ends_with("_pairs") |
      ends_with("_profiled"),
    ~ sum(.x, na.rm = T) # sum N + S values
  ), .groups = "drop") %>%
  mutate(location = case_when(
    lab_id == "MU" ~ "Boston DITP ALL, MA (MU)",
    lab_id == "SB" ~ "Boston DITP ALL, MA (SB)",
    T ~ ""
  ),
  population = 1100000 + 1400000, # sum of N and S populations
  pmmov_counts_scaled = pmmov_counts * total_read_pairs / n_reads_profiled) %>%
  select(location,state_abbr,state_name,lab_id,population,date,total_read_pairs,
         n_reads_profiled,pmmov_counts,pmmov_counts_scaled,everything()) # order columns

# Append to df
casper_clade_counts_wide <- rbind(casper_clade_counts_wide,DITP_data)

## Epic ----
# Convert df to long format
epic_dfs_long <- map(epic_dfs,
                     ~ pivot_longer(
                       .x,
                       cols = c(`Sep 7 – Sep 13, 2024`:`Mar 28 – Mar 29, 2026`),
                       names_to = "week",
                       values_to = "count",
                       values_transform = list(count = as.character)
                     ) %>% 
                       mutate(week_end = as.Date(sub(".*–\\s*", "", week), format = "%b %d, %Y"),
                              count = as.numeric(count)) %>% 
                       filter(week_end!=as.Date("2026-03-29"))) # omit the 2-d week

saveRDS(epic_dfs_long, file = here("data/clean/all_epic_long.rds"))


# ----- PMMoV normalization -----
## WWSCAN -----
summary(wwscan_2025$PMMoV_gc_g_dry_weight)
colnames(wwscan_2025)

wwscan_2025 <- wwscan_2025 %>%
  mutate(across(C_auris_gc_g_dry_weight:SC2_N_gc_g_dry_weight, 
                ~ .x == 0, .names = "{.col}_is_nd"))

if (impute_nondetects) {
  wwscan_2025 <- wwscan_2025 %>%
    mutate(across(C_auris_gc_g_dry_weight:SC2_N_gc_g_dry_weight, 
                  ~ ifelse(.x == 0, 1000/2, .x)))
}

wwscan_2025_PMMoV_norm <- wwscan_2025 %>%
  filter(PMMoV_gc_g_dry_weight > 0 &
           !is.na(PMMoV_gc_g_dry_weight)) %>%
  select(!ends_with(c("_UCI", "_LCI"))) %>%
  relocate(City:Year, BCoV_Recovery, PMMoV_gc_g_dry_weight, .before = 1) %>%
  mutate(
    across(
      C_auris_gc_g_dry_weight:SC2_N_gc_g_dry_weight,
      ~ .x / PMMoV_gc_g_dry_weight,
      .names = "{.col}_PMMoV_norm"
    )
  )

# Aggregate Lander St and West Boise in ID
lander_st_data <- wwscan_2025_PMMoV_norm %>% 
  filter(Plant == "Lander Street, Boise, ID")

west_boise_data <- wwscan_2025_PMMoV_norm %>% 
  filter(Plant == "West Boise, ID") 

# ^ virtually the same sampling dates so use pop-weighted average by date
boise_data <- wwscan_2025_PMMoV_norm %>%
  filter(Plant == "Lander Street, Boise, ID" | Plant == "West Boise, ID") %>%
  group_by(Month, Day, Year) %>%
  summarise(
    across(
      ends_with("BCoV_Recovery") |
        ends_with("_gc_g_dry_weight") |
        ends_with("_PMMoV_norm"),
      ~ weighted.mean(.x, w = Population_Served, na.rm = T)
    ),
    across(ends_with("_is_nd"), ~ all(.x, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    City = "Boise",
    State = "Idaho",
    State_Abbr = "ID",
    Plant = "All Boise",
    Site_Name = "Lander St and West Boise",
    Population_Served = 108556 + 186901
  ) %>% 
  select(City,State,State_Abbr,Plant,Site_Name,Population_Served,
         everything()) # order columns

# Append to df
wwscan_2025_PMMoV_norm <- rbind(wwscan_2025_PMMoV_norm,boise_data)

write_csv(wwscan_2025_PMMoV_norm,
          here(paste0("data/clean/wwscan_2025_PMMoV_norm", nd_suffix, ".csv")))

## CASPER -----
summary(casper_clade_counts_wide$pmmov_counts_scaled)
colnames(casper_clade_counts_wide)

casper_clade_counts_wide <- casper_clade_counts_wide %>%
  mutate(across(ends_with("_clade_counts"), ~ .x == 0, .names = "{.col}_is_nd"))

if (impute_nondetects) {
  casper_clade_counts_wide <- casper_clade_counts_wide %>%
    mutate(across(ends_with("_clade_counts"), ~ ifelse(.x == 0, 1/2, .x)))
}

casper_clade_counts_wide <- casper_clade_counts_wide %>%
  filter(pmmov_counts_scaled > 0) %>%
  mutate(
    across(
      ends_with("_clade_counts"),
      ~ .x / pmmov_counts_scaled,
      .names = "{.col}_PMMoV_scaled_norm"
    )
  )

sort(unique(casper_combined$label)) # note: Human_mastadenovirus_F = HAdV F (includes F40/41)

write_csv(
  casper_clade_counts_wide,
  here(paste0("data/clean/casper_clade_counts_combined_PMMoV_scaled_norm", nd_suffix, ".csv"))
)
