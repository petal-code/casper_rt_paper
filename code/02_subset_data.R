# Subset datasets for analysis

source(here::here("code/00_config.R"))


# ---- Load data ----
epic_dfs <- read_rds(here("data/clean/all_epic_long.rds"))
wwscan_wwtp_data <- read_csv(here(paste0("data/clean/wwscan_2025_PMMoV_norm", nd_suffix, ".csv")))
casper_wwtp_data <- read_csv(here(paste0("data/clean/casper_clade_counts_combined_PMMoV_scaled_norm", nd_suffix, ".csv")))

# Convert dates to Date format
wwscan_wwtp_data$date <- make_date(wwscan_wwtp_data$Year, wwscan_wwtp_data$Month, wwscan_wwtp_data$Day)
casper_wwtp_data$date <- as.Date(casper_wwtp_data$date)


# ---- Subset state-level clinical data (Epic) ----
epic_states_include <- c("CA","MA","ID","IA")

epic_state_name_map <- c(CA = "California", MA = "Massachusetts", ID = "Idaho", IA = "Iowa")

# Maps virus key -> source df, filter column, filter value
epic_virus_spec <- list(
  SC2  = list(state_df = "sc2_diagnosis_state",        county_df = "sc2_diagnosis_county",
              col = "Diagnosis (Encounter)",       state_val = "COVID_derived",          county_val = "COVID_derived"),
  IAV  = list(state_df = "iav_lab_no_ind_state",        county_df = "iav_lab_no_ind_county",
              col = "Influenza A Lab Results",    state_val = "Positive_derived",          county_val = "Positive_derived"),
  IBV  = list(state_df = "ibv_lab_no_ind_state",        county_df = "ibv_lab_no_ind_county",
              col = "Influenza B Lab Results",    state_val = "Positive_derived",          county_val = "Positive_derived"),
  RSV  = list(state_df = "rsv_diagnosis_state",  county_df = "rsv_diagnosis_county",
              col = "Diagnosis (Encounter)",      state_val = "RSV_derived",       county_val = "RSV_derived"),
  noro = list(state_df = "noro_diagnosis_state", county_df = "noro_diagnosis_county",
              col = "Diagnosis (Encounter)",      state_val = "Norovirus_derived", county_val = "Norovirus_derived"),
  rota = list(state_df = "rota_diagnosis_state", county_df = "rota_diagnosis_county",
              col = "Diagnosis (Encounter)",      state_val = "Rotavirus_derived", county_val = "Rotavirus_derived")
)

epic_viruses_include <- names(epic_virus_spec)

epic_state_dfs <- setNames(
  lapply(epic_states_include, function(x) setNames(vector("list", length(epic_viruses_include)), epic_viruses_include)),
  epic_states_include
)

for (s in epic_states_include) {
  for (v in epic_viruses_include) {
    spec <- epic_virus_spec[[v]]
    epic_state_dfs[[s]][[v]] <- epic_dfs[[spec$state_df]] %>%
      filter(`State of Residence` == epic_state_name_map[s],
             .data[[spec$col]] == spec$state_val) %>%
      select(week_end, count) %>%
      rename(date = week_end, confirm = count)
  }
}

saveRDS(epic_state_dfs, file = here("data/clean/subsetted/epic_state_viruses.rds"))


# ---- Subset local-level data ----
## Epic ----
epic_counties_include <- c("Riverside_CA", "Suffolk_MA", "Ada_ID")

# use same viruses as state

epic_county_string_map <- c(Riverside_CA = "RIVERSIDE, CA", Suffolk_MA = "SUFFOLK, MA", Ada_ID = "ADA, ID")

epic_county_dfs <- setNames(
  lapply(epic_counties_include, function(x) setNames(vector("list", length(epic_viruses_include)), epic_viruses_include)),
  epic_counties_include
)

for (cty in epic_counties_include) {
  for (v in epic_viruses_include) {
    spec <- epic_virus_spec[[v]]
    epic_county_dfs[[cty]][[v]] <- epic_dfs[[spec$county_df]] %>%
      filter(`County of Residence` == epic_county_string_map[cty],
             .data[[spec$col]] == spec$county_val) %>%
      select(week_end, count) %>%
      rename(date = week_end, confirm = count)
  }
}

saveRDS(epic_county_dfs, file = here("data/clean/subsetted/epic_county_viruses.rds"))

## WWSCAN ----
wwscan_wwtp_data <- wwscan_wwtp_data %>%
  mutate(Plant = gsub(", ","_",Plant),
         Plant = gsub(" ","_",Plant)) %>%
  select(City:Year, date, ends_with("PMMoV_norm"), ends_with("_is_nd"))

wwscan_wwtps_include <- c("All_Boise", 
                          "Boston_MA",
                          "Ottumwa_IA",
                          "Riverside_CA")

wwscan_viruses_include <- c("SC2_N","Influenza_A","InfA_H1","InfA_H3","InfA_H5","Influenza_B",
                            "RSV","HMPV_4","EVD68","Parvo_B19","Noro_G2",
                            "Rotavirus","HAdV_F","HAV")

wwscan_wwtp_dfs <- setNames(
  lapply(wwscan_wwtps_include, function(x){
    setNames(vector("list", length(wwscan_viruses_include)), wwscan_viruses_include)
  }),
  wwscan_wwtps_include
)

for(w in 1:length(wwscan_wwtps_include)){
  for(v in 1:length(wwscan_viruses_include)){
    tmp <- wwscan_wwtp_data %>%
      filter(Plant == wwscan_wwtps_include[w]) %>%
      select(date,
             matches(wwscan_viruses_include[v]) & ends_with("_norm"),
             matches(wwscan_viruses_include[v]) & ends_with("_is_nd")) %>%
      rename(confirm = 2, is_nd = 3) %>%
      filter(!is.na(confirm))
    wwscan_wwtp_dfs[[wwscan_wwtps_include[w]]][[wwscan_viruses_include[v]]] <- tmp
  }
}

saveRDS(wwscan_wwtp_dfs, file = here(paste0("data/clean/subsetted/wwscan_wwtps_viruses", nd_suffix, ".rds")))

## CASPER ----
casper_wwtp_data <- casper_wwtp_data %>%
  mutate(location = gsub(", ","_",location),
         location = gsub(" ","_",location),
         location = gsub("-","_",location),
         location = gsub("\\(","_",location),
         location = gsub("\\)","_",location)) %>%
  select(location:date, ends_with("PMMoV_scaled_norm"), ends_with("_is_nd"))

# Only include MU labs for this analysis
casper_wwtp_data <- casper_wwtp_data %>% 
  filter(lab_id == "MU")

casper_wwtps_include <- c("Boise_RWF_ID", 
                          "Boston_DITP_ALL_MA__MU_",
                          "Ottumwa_WPCF_IA",
                          "Riverside_WQCP_CA")

casper_viruses_include <- c("SARS_CoV_2","Influenza_A_clade_counts","Influenza_A_H1N1","Influenza_A_H3N1","Influenza_A_H5N1",
                            "Influenza_B","RSV_clade_counts","Metapneumovirus","Enterovirus_D68",
                            "Parvovirus_B19","Norovirus_genogroup_II","Rotavirus_clade_counts","Human_mastadenovirus_F","Hepatitis_A","Hepatitis_B","Hepatitis_C",
                            "Aichivirus","Astrovirus_clade_counts","Sapovirus","Human_bocavirus","RSV_A","RSV_B","Rotavirus_A",
                            "Rotavirus_B","Rotavirus_C_clade_counts","Human_adenovirus_F40","Human_adenovirus_F41","Influenza_C","Parainfluenza_virus_1",
                            "Parainfluenza_virus_2","Parainfluenza_virus_3","Parainfluenza_virus_4","Norovirus_genogroup_I_",# underscore needed so GII not captured again
                            "Norovirus_genogroup_IV","Rhinovirus_A","Rhinovirus_B","Rhinovirus_C_clade_counts","Human_coronavirus_229E","Human_coronavirus_HKU1",
                            "Human_coronavirus_NL63","Human_coronavirus_OC43","Enterovirus_A","Enterovirus_B","Enterovirus_C_clade_counts","Human_mastadenovirus_A",
                            "Human_mastadenovirus_B","Human_mastadenovirus_C_clade_counts","Human_mastadenovirus_D","Norwalk","Enterovirus_clade_counts"
                            )

casper_wwtp_dfs <- setNames(
  lapply(casper_wwtps_include, function(x){
    setNames(vector("list", length(casper_viruses_include)), casper_viruses_include)
  }),
  casper_wwtps_include
)

for(w in 1:length(casper_wwtps_include)){
  for(v in 1:length(casper_viruses_include)){
    tmp <- casper_wwtp_data %>%
      filter(location == casper_wwtps_include[w]) %>%
      select(date,
             matches(casper_viruses_include[v]) & ends_with("_norm"),
             matches(casper_viruses_include[v]) & ends_with("_is_nd")) %>%
      rename(confirm = 2, is_nd = 3) %>%
      filter(!is.na(confirm))
    casper_wwtp_dfs[[casper_wwtps_include[w]]][[casper_viruses_include[v]]] <- tmp
  }
}

saveRDS(casper_wwtp_dfs, file = here(paste0("data/clean/subsetted/casper_wwtps_viruses", nd_suffix, ".rds")))
