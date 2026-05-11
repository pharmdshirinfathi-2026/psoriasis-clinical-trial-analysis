# 1. Libraries
# Run this script from the repository/project folder that contains meta_data_psoriasis.csv.
# Example: Rscript Clinical_trial_analysis.R
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(purrr)
library(tibble)
library(ggplot2)
library(lubridate)
library(grid)
# 2. Input and output folders
input_file <- "meta_data_psoriasis.csv"
output_dir <- paste0("meta_data_psoriasis_", Sys.Date())
table_dir  <- file.path(output_dir, "tables")
qc_dir     <- file.path(output_dir, "QC_manual_classification")
figure_dir <- file.path(output_dir, "figures")
dir.create(output_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)
dir.create(qc_dir, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)
# 3. Read input file
ctg_raw <- read_csv(input_file, show_col_types = FALSE)
required_cols <- c(
  "NCT Number",
  "Study Title",
  "Study URL",
  "Brief Summary",
  "Interventions",
  "Primary Outcome Measures",
  "Phases",
  "Study Design",
  "Completion Date"
)
missing_cols <- setdiff(required_cols, names(ctg_raw))
if (length(missing_cols) > 0) {
  stop(
    "These required columns are missing from the input file: ",
    paste(missing_cols, collapse = ", ")
  )
}
# 4. Helper functions
normalize_text <- function(x) {
  x <- ifelse(is.na(x), "", x)
  x |>
    str_to_lower() |>
    str_replace_all("[\u2010\u2011\u2012\u2013\u2014\u2212]", "-") |>
    str_replace_all("[_/]", "-") |>
    str_replace_all("[^a-z0-9\\-\\s/\\+\\(\\)]", " ") |>
    str_squish()
}
paste_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}
extract_year_from_completion <- function(x) {
  x <- as.character(x)
  year1 <- str_extract(x, "^\\d{4}")
  year2 <- str_extract(x, "\\d{4}$")
  out <- ifelse(!is.na(year1), year1, year2)
  as.integer(out)
}
# 5. Drug dictionary with numeric class codes
# class_code:0 = non-drug / non-relevant and 1 = small molecule and 2 = biologic and 3 = conventional systemic / csDMARD / conventional topical
drug_dictionary <- tribble(
  ~standard_name, ~class_code, ~class_label, ~target_group, ~aliases,
# 0. Non-drug / non-relevant
  "placebo", 0, "Non-drug / non-relevant", "Comparator",
  "placebo|placebo tablet|placebo injection|placebo comparator|aba placebo|ust placebo|secukinumab placebo|placebo methotrexate|placebo followed by secukinumab|placebo to adalimumab|placebo to etanercept|placebo to methotrexate|placebo to upadacitinib|placebo for risankizumab|placebo matching bms[-\\s]?945429|placebo controlled period",
  "vehicle", 0, "Non-drug / non-relevant", "Vehicle/control",
  "vehicle|vehicle cream|vehicle ointment",
  "standard of care", 0, "Non-drug / non-relevant", "Standard care/control",
  "standard of care",
  "lab draws and imaging", 0, "Non-drug / non-relevant", "Procedure",
  "lab draws and imaging|lab draw|imaging",
  "autoinjector", 0, "Non-drug / non-relevant", "Device",
  "autoinjector|autoinjector a|device",
  "other", 0, "Non-drug / non-relevant", "Other",
  "^other$",
  "phototherapy", 0, "Non-drug / non-relevant", "Phototherapy",
  "phototherapy|uvb|uva|nb[-\\s]?uvb|puva|excimer|laser",
  "metformin", 0, "Non-drug / non-relevant", "Non-psoriasis metabolic drug",
  "metformin",
  "fumaric acid esters", 0, "Non-drug / non-relevant", "Excluded/non-main-analysis drug",
  "fumaric acid|fumaric acid esters|dimethyl fumarate|fumarate",
  "naltrexone", 0, "Non-drug / non-relevant", "Non-main-analysis drug",
  "naltrexone",
  "herpes zoster vaccine", 0, "Non-drug / non-relevant", "Vaccine",
  "herpes zoster|hz vaccine|zoster vaccine",
  # 2. Biologics: TNF-alpha inhibitors
  "adalimumab", 2, "Biologic", "TNF-alpha",
  "adalimumab|humira|abrilada|amjevita|cyltezo|hadlima|hulio|hyrimoz|idacio|yuflyma|yusimry|adalimumab \\(ada\\)",
  "infliximab", 2, "Biologic", "TNF-alpha",
  "infliximab|remicade|inflectra|renflexis|avsola|ixifi|zymfentra|ta[-\\s]?650|ifx",
  "etanercept", 2, "Biologic", "TNF-alpha",
  "etanercept|enbrel|erelzi|chs[-\\s]?0214|active comparator etanercept|etanercept / autoinjector a",
  "certolizumab pegol", 2, "Biologic", "TNF-alpha",
  "certolizumab|certolizumab pegol|cimzia|cdp[-\\s]?870|czp|czp 200 mg q2w|czp 400 mg q4w",
  "golimumab", 2, "Biologic", "TNF-alpha",
  "golimumab|simponi|simponi aria|cnto[-\\s]?148|placebo golimumab",
  # 2. Biologics: IL-17 axis
  "secukinumab", 2, "Biologic", "IL-17A",
  "secukinumab|cosentyx|ain[-\\s]?457|ain457|ain457a|secukinumab 75 mg|secukinumab 150 mg|secukinumab 300 mg|secukinumab low dose|secukinumab high dose",
  "ixekizumab", 2, "Biologic", "IL-17A",
  "ixekizumab|taltz|ly[-\\s]?2439821|ixekizumab 80 mg|80 mg ixekizumab|ixekizumab 80 mg q2w|ixekizumab 80 mg q4w",
  "brodalumab", 2, "Biologic", "IL-17RA",
  "brodalumab|siliq|lumicef|amg[-\\s]?827|140 mg brodalumab|210 mg brodalumab",
  "bimekizumab", 2, "Biologic", "IL-17A/F",
  "bimekizumab|bimzelx|ucb[-\\s]?4940",
  "sonelokimab", 2, "Biologic", "IL-17A/F",
  "sonelokimab|m[-\\s]?1095",
  "netakimab", 2, "Biologic", "IL-17A",
  "netakimab|bcd[-\\s]?085|efleira",
  "izokibep / ABY-035", 2, "Biologic", "IL-17A",
  "izokibep|aby[-\\s]?035|aby035",
  "remtolumab / ABT-122", 2, "Biologic", "TNF-alpha/IL-17A",
  "remtolumab|abt[-\\s]?122|abt122",
  # 2. Biologics: IL-12/23 and IL-23
  "ustekinumab", 2, "Biologic", "IL-12/23p40",
  "ustekinumab|stelara|cnto[-\\s]?1275|cnto 1275 63 mg|ustekinumab 45 mg|ustekinumab 90 mg",
  "guselkumab", 2, "Biologic", "IL-23p19",
  "guselkumab|tremfya|cnto[-\\s]?1959|cnto 1959 50 mg|ctno 1959 100 mg|cn_to_1959",
  "tildrakizumab", 2, "Biologic", "IL-23p19",
  "tildrakizumab|ilumya|ilumetri|mk[-\\s]?3222|sunpg18_07|tildrakizumab 100 mg|tildrakizumab 200 mg",
  "risankizumab", 2, "Biologic", "IL-23p19",
  "risankizumab|skyrizi|bi[-\\s]?655066|abbv[-\\s]?066|placebo for risankizumab",
  "mirikizumab", 2, "Biologic", "IL-23p19",
  "mirikizumab|omvoh|ly[-\\s]?3074828",
  "BMS-945429", 2, "Biologic", "IL-23/immune biologic",
  "bms[-\\s]?945429|placebo matching bms[-\\s]?945429",
  # 2. Biologics: IL-36 and other biologic targets
  "spesolimab", 2, "Biologic", "IL-36R",
  "spesolimab|spesolima|spevigo|bi[-\\s]?655130",
  "imsidolimab", 2, "Biologic", "IL-36R",
  "imsidolimab|anb[-\\s]?019",
  "efalizumab", 2, "Biologic", "CD11a/LFA-1",
  "efalizumab|raptiva|hu1124",
  "alefacept", 2, "Biologic", "CD2/LFA-3",
  "alefacept|amevive",
  "siplizumab", 2, "Biologic", "CD2",
  "siplizumab|medi[-\\s]?507",
  "galiximab", 2, "Biologic", "CD80",
  "galiximab|idec[-\\s]?114",
  "abatacept", 2, "Biologic", "CD80/CD86",
  "abatacept|orencia|ctla4[-\\s]?ig|aba",
  "cetuximab", 2, "Biologic", "EGFR",
  "cetuximab|erbitux",
  # 1. Small molecules
  "apremilast", 1, "Small molecule", "PDE4",
  "apremilast|otezla|cc[-\\s]?10004|apremilast 20mg|apremilast 30mg|apremilast 30 mg|placebo \\+ 20 mg apremilast|placebo \\+ 30 mg apremilast",
  "roflumilast", 1, "Small molecule", "PDE4",
  "roflumilast|zoryve|daliresp|apta[-\\s]?2217|b9302[-\\s]?107|by[-\\s]?217",
  "crisaborole", 1, "Small molecule", "PDE4",
  "crisaborole|eucrisa|an[-\\s]?2728",
  "Hemay005", 1, "Small molecule", "PDE4",
  "mufemilast|hemay[-\\s]?005|hemay005",
  "tofacitinib", 1, "Small molecule", "JAK",
  "tofacitinib|xeljanz|cp[-\\s]?690[\\s,]?550|tofacitinib 5 mg bid|tofacitinib 10 mg bid",
  "upadacitinib", 1, "Small molecule", "JAK1",
  "upadacitinib|upadacitini|rinvoq|abt[-\\s]?494|placebo to upadacitinib",
  "filgotinib", 1, "Small molecule", "JAK1",
  "filgotinib|jyseleca|glpg[-\\s]?0634|gs[-\\s]?6034",
  "baricitinib", 1, "Small molecule", "JAK1/2",
  "baricitinib|olumiant|ly[-\\s]?3009104",
  "peficitinib", 1, "Small molecule", "JAK",
  "peficitinib|smyraf|asp[-\\s]?015k",
  "ruxolitinib", 1, "Small molecule", "JAK1/2",
  "ruxolitinib|jakafi|jakavi|opzelura|incb[-\\s]?018424",
  "abrocitinib", 1, "Small molecule", "JAK1",
  "abrocitinib|cibinqo|pf[-\\s]?04965842",
  "itacitinib", 1, "Small molecule", "JAK1",
  "itacitinib|incb[-\\s]?039110|ndi[-\\s]?034858",
  "solcitinib", 1, "Small molecule", "JAK1",
  "solcitinib|gsk[-\\s]?2586184",
  "brepocitinib", 1, "Small molecule", "TYK2/JAK1",
  "brepocitinib|pf[-\\s]?06700841",
  "deucravacitinib", 1, "Small molecule", "TYK2",
  "deucravacitinib|sotyktu|bms[-\\s]?986165",
  "BMS-986202", 1, "Small molecule", "TYK2",
  "bms[-\\s]?986202",
  "zasocitinib", 1, "Small molecule", "TYK2",
  "zasocitinib|tak[-\\s]?279",
  "tapinarof", 1, "Small molecule", "AhR",
  "tapinarof|vtama|gsk[-\\s]?2894512|wbi[-\\s]?1001",
  "ponesimod", 1, "Small molecule", "S1PR1",
  "ponesimod|ponvory|act[-\\s]?128800|rg[-\\s]?3477",
  "piclidenoson", 1, "Small molecule", "A3AR",
  "piclidenoson|cf[-\\s]?101|ib[-\\s]?meca|ibmeca",
  "belumosudil", 1, "Small molecule", "ROCK2",
  "belumosudil|rezurock|kd[-\\s]?025",
  "vimirogant", 1, "Small molecule", "RORgammat",
  "vimirogant|vimirogant 43742|vtp[-\\s]?43742",
  "PF-06763809", 1, "Small molecule", "RORgammat",
  "pf[-\\s]?06763809",
  "S18-000003", 1, "Small molecule", "RORgammat",
  "s18[-\\s]?000003",
  "GSK2981278", 1, "Small molecule", "RORgammat",
  "gsk[-\\s]?2981278|gsk 2981278",
  "icotrokinra", 1, "Small molecule", "IL-23 receptor",
  "icotrokinra|icotyde|jnj[-\\s]?2113|pn[-\\s]?235",
  # 3. Conventional systemic / csDMARD / conventional topical
  "methotrexate", 3, "Conventional systemic / csDMARD / conventional topical", "Antimetabolite/csDMARD",
  "methotrexate|mtx|methotrexate \\(mtx\\)|otrexup|rasuvo|trexall|placebo methotrexate|placebo to methotrexate",
  "cyclosporine", 3, "Conventional systemic / csDMARD / conventional topical", "Calcineurin inhibitor",
  "cyclosporine|cyclosporins|ciclosporin|ciclosporine|neoral|sandimmune|gengraf",
  "acitretin", 3, "Conventional systemic / csDMARD / conventional topical", "Retinoid",
  "acitretin|soriatane|neotigason",
  "leflunomide", 3, "Conventional systemic / csDMARD / conventional topical", "csDMARD",
  "leflunomide|arava",
  "sulfasalazine", 3, "Conventional systemic / csDMARD / conventional topical", "csDMARD",
  "sulfasalazine|azulfidine",
  "hydroxychloroquine", 3, "Conventional systemic / csDMARD / conventional topical", "Conventional systemic antirheumatic",
  "hydroxychloroquine|plaquenil",
  "clobetasol propionate", 3, "Conventional systemic / csDMARD / conventional topical", "Topical corticosteroid",
  "clobetasol|clobetasol propionate|temovate|clobex|impoyz",
  "halobetasol propionate", 3, "Conventional systemic / csDMARD / conventional topical", "Topical corticosteroid",
  "halobetasol|halobetasol propionate|ultravate|lexette|bryhali",
  "betamethasone", 3, "Conventional systemic / csDMARD / conventional topical", "Topical corticosteroid",
  "betamethasone|diprolene|luxiq|enstilar",
  "methylprednisolone", 3, "Conventional systemic / csDMARD / conventional topical", "Systemic corticosteroid",
  "methylprednisolone",
  "triamcinolone", 3, "Conventional systemic / csDMARD / conventional topical", "Corticosteroid",
  "triamcinolone",
  "systemic corticosteroids", 3, "Conventional systemic / csDMARD / conventional topical", "Systemic corticosteroid",
  "systemic corticosteroids|systemic corticosteroids/methotrexate",
  "calcipotriol / calcipotriene", 3, "Conventional systemic / csDMARD / conventional topical", "Vitamin D analog",
  "calcipotriol|calcipotriene|dovonex|sorilux",
  "calcitriol", 3, "Conventional systemic / csDMARD / conventional topical", "Vitamin D analog",
  "calcitriol|vectical",
  "tazarotene", 3, "Conventional systemic / csDMARD / conventional topical", "Topical retinoid",
  "tazarotene|tazorac|arazlo|fabior",
  "retinoids", 3, "Conventional systemic / csDMARD / conventional topical", "Retinoid class",
  "retinoids",
  "tacrolimus", 3, "Conventional systemic / csDMARD / conventional topical", "Topical calcineurin inhibitor",
  "tacrolimus|protopic",
  "pimecrolimus", 3, "Conventional systemic / csDMARD / conventional topical", "Topical calcineurin inhibitor",
  "pimecrolimus|elidel",
  "anthralin", 3, "Conventional systemic / csDMARD / conventional topical", "Topical anti-proliferative",
  "anthralin|dithranol"
)
# 6. Parse intervention column into intervention-level table
intervention_long <- ctg_raw |>
  mutate(row_id = row_number()) |>
  select(
    row_id,
    `NCT Number`,
    `Study Title`,
    `Study URL`,
    Interventions,
    Phases,
    `Study Design`,
    `Completion Date`
  ) |>
  separate_rows(Interventions, sep = "\\|") |>
  mutate(
    intervention_raw = str_squish(Interventions),
    intervention_type = if_else(
      str_detect(intervention_raw, ":"),
      str_squish(str_extract(intervention_raw, "^[^:]+")),
      NA_character_
    ),
    intervention_name_original = if_else(
      str_detect(intervention_raw, ":"),
      str_squish(str_replace(intervention_raw, "^[^:]+:", "")),
      intervention_raw
    ),
    intervention_name_clean = normalize_text(intervention_name_original)
  ) |>
  filter(intervention_name_original != "") |>
  select(-Interventions)
# 7. Classification function
classify_intervention <- function(x) {
  x_norm <- normalize_text(x)
  hits <- drug_dictionary |>
    rowwise() |>
    mutate(
      matched = str_detect(
        x_norm,
        regex(paste0("\\b(", aliases, ")\\b"), ignore_case = TRUE)
      )
    ) |>
    ungroup() |>
    filter(matched)
  if (nrow(hits) == 0) {
    return(tibble(
      standard_name = NA_character_,
      class_code = NA_integer_,
      class_label = "UNMATCHED",
      target_group = NA_character_,
      match_status = "unmatched_dictionary"
    ))
  }
  priority_order <- c(2, 1, 3, 0)
  selected_code <- priority_order[priority_order %in% hits$class_code][1]
  selected_label <- hits |>
    filter(class_code == selected_code) |>
    pull(class_label) |>
    unique() |>
    first()
  selected_targets <- hits |>
    filter(class_code == selected_code) |>
    pull(target_group) |>
    unique() |>
    paste(collapse = "; ")
  tibble(
    standard_name = paste_unique(hits$standard_name),
    class_code = selected_code,
    class_label = selected_label,
    target_group = selected_targets,
    match_status = "matched_dictionary"
  )
}
classification_long <- intervention_long |>
  rowwise() |>
  mutate(class_info = list(classify_intervention(intervention_name_original))) |>
  unnest(class_info) |>
  ungroup()
# 8. Unmatched file for manual classification
unmatched_for_manual <- classification_long |>
  filter(match_status == "unmatched_dictionary") |>
  group_by(
    intervention_name_original,
    intervention_name_clean,
    intervention_type
  ) |>
  summarise(
    n_occurrences = n(),
    n_unique_trials = n_distinct(`NCT Number`),
    example_nct_ids = paste(head(unique(`NCT Number`), 8), collapse = "; "),
    example_titles = paste(head(unique(`Study Title`), 3), collapse = " || "),
    suggested_manual_code = NA_integer_,
    suggested_standard_name = NA_character_,
    manual_note = NA_character_,
    .groups = "drop"
  ) |>
  arrange(desc(n_occurrences), intervention_name_original)
write_csv(
  unmatched_for_manual,
  file.path(qc_dir, "UNMATCHED_interventions_for_manual_coding.csv")
)
# 9. Optional manual override file
# suggested_manual_code:
# 0 = non-drug / non-relevant
# 1 = small molecule
# 2 = biologic
# 3 = conventional systemic / csDMARD / conventional topical
manual_template_file <- file.path(qc_dir, "manual_intervention_classification_TEMPLATE.csv")
if (!file.exists(manual_template_file)) {
  write_csv(unmatched_for_manual, manual_template_file)
}
manual_completed_file <- file.path(qc_dir, "manual_intervention_classification_completed.csv")
if (file.exists(manual_completed_file)) {
  manual_completed <- read_csv(manual_completed_file, show_col_types = FALSE) |>
    mutate(
      intervention_name_clean = normalize_text(intervention_name_clean),
      suggested_manual_code = as.integer(suggested_manual_code)
    ) |>
    select(
      intervention_name_clean,
      suggested_manual_code,
      suggested_standard_name,
      manual_note
    )
  classification_long <- classification_long |>
    left_join(manual_completed, by = "intervention_name_clean") |>
    mutate(
      class_code = if_else(
        match_status == "unmatched_dictionary" & !is.na(suggested_manual_code),
        suggested_manual_code,
        class_code
      ),
      standard_name = if_else(
        match_status == "unmatched_dictionary" &
          !is.na(suggested_standard_name) &
          suggested_standard_name != "",
        suggested_standard_name,
        standard_name
      ),
      class_label = case_when(
        class_code == 0 ~ "Non-drug / non-relevant",
        class_code == 1 ~ "Small molecule",
        class_code == 2 ~ "Biologic",
        class_code == 3 ~ "Conventional systemic / csDMARD / conventional topical",
        TRUE ~ class_label
      ),
      match_status = if_else(
        match_status == "unmatched_dictionary" & !is.na(suggested_manual_code),
        "manual_override",
        match_status
      )
    )
}
# 10. Save intervention-level classified table
write_csv(
  classification_long,
  file.path(table_dir, "01_intervention_level_classified.csv")
)
matched_summary <- classification_long |>
  count(
    standard_name,
    class_code,
    class_label,
    target_group,
    match_status,
    sort = TRUE,
    name = "n_occurrences"
  )
write_csv(
  matched_summary,
  file.path(table_dir, "02_intervention_matched_summary.csv")
)
# 11. Trial-level classification
trial_classification <- classification_long |>
  group_by(
    row_id,
    `NCT Number`,
    `Study Title`,
    `Study URL`,
    Phases,
    `Study Design`,
    `Completion Date`
  ) |>
  summarise(
    all_interventions = paste_unique(intervention_name_original),
    matched_standard_names = paste_unique(standard_name),
    all_class_codes = paste_unique(as.character(class_code)),
    has_code_0 = any(class_code == 0, na.rm = TRUE),
    has_code_1 = any(class_code == 1, na.rm = TRUE),
    has_code_2 = any(class_code == 2, na.rm = TRUE),
    has_code_3 = any(class_code == 3, na.rm = TRUE),
    has_unmatched = any(is.na(class_code), na.rm = TRUE),
    n_intervention_items = n(),
    n_unmatched_items = sum(is.na(class_code)),
    .groups = "drop"
  ) |>
  mutate(
    trial_class_code = case_when(
      has_code_2 & !has_code_1 & !has_code_3 ~ "2",
      has_code_1 & !has_code_2 & !has_code_3 ~ "1",
      has_code_3 & !has_code_1 & !has_code_2 ~ "3",
      has_code_2 & (has_code_1 | has_code_3) ~ "mixed_biologic_plus_nonbiologic",
      has_code_1 & has_code_3 & !has_code_2 ~ "mixed_small_plus_conventional",
      !has_code_1 & !has_code_2 & !has_code_3 & has_code_0 ~ "0",
      TRUE ~ "unmatched_or_unclear"
    ),
    trial_class_label = case_when(
      trial_class_code == "2" ~ "Biologic only",
      trial_class_code == "1" ~ "Small molecule only",
      trial_class_code == "3" ~ "Conventional systemic/topical only",
      trial_class_code == "mixed_biologic_plus_nonbiologic" ~
        "Combination: biologic + small molecule or conventional",
      trial_class_code == "mixed_small_plus_conventional" ~
        "Combination: small molecule + conventional",
      trial_class_code == "0" ~ "Non-drug / non-relevant only",
      TRUE ~ "Unmatched / unclear"
    ),
    completion_year = extract_year_from_completion(`Completion Date`),
    phase_clean = case_when(
      str_detect(Phases, "PHASE1") & str_detect(Phases, "PHASE2") ~ "Phase 1/2",
      str_detect(Phases, "PHASE2") & str_detect(Phases, "PHASE3") ~ "Phase 2/3",
      str_detect(Phases, "PHASE2") ~ "Phase 2",
      str_detect(Phases, "PHASE3") ~ "Phase 3",
      str_detect(Phases, "PHASE4") ~ "Phase 4",
      TRUE ~ Phases
    ),
    phase_clean = factor(
      phase_clean,
      levels = c("Phase 1/2", "Phase 2", "Phase 2/3", "Phase 3", "Phase 4")
    )
  )
write_csv(
  trial_classification,
  file.path(table_dir, "03_trial_level_classified.csv")
)
# 12. Summary tables
summary_trial_class <- trial_classification |>
  count(trial_class_label, name = "n_trials") |>
  arrange(desc(n_trials))
summary_phase_by_class <- trial_classification |>
  count(phase_clean, trial_class_label, name = "n_trials") |>
  arrange(phase_clean, desc(n_trials))
summary_year_by_class <- trial_classification |>
  filter(!is.na(completion_year)) |>
  count(completion_year, trial_class_label, name = "n_trials") |>
  arrange(completion_year, trial_class_label)
summary_intervention_code <- classification_long |>
  count(class_code, class_label, name = "n_intervention_records") |>
  arrange(class_code)
summary_unmatched_by_trial <- trial_classification |>
  filter(has_unmatched | n_unmatched_items > 0) |>
  arrange(desc(n_unmatched_items))
write_csv(summary_trial_class, file.path(table_dir, "summary_01_trial_class_distribution.csv"))
write_csv(summary_phase_by_class, file.path(table_dir, "summary_02_phase_by_trial_class.csv"))
write_csv(summary_year_by_class, file.path(table_dir, "summary_03_completion_year_by_trial_class.csv"))
write_csv(summary_intervention_code, file.path(table_dir, "summary_04_intervention_code_distribution.csv"))
write_csv(summary_unmatched_by_trial, file.path(table_dir, "summary_05_trials_with_unmatched_items.csv"))
# 13. Figures
pub_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "bottom"
  )
class_colors <- c(
  "Biologic only" = "#D55E00",
  "Small molecule only" = "#0072B2",
  "Conventional systemic/topical only" = "#CC79A7",
  "Combination: biologic + small molecule or conventional" = "#009E73",
  "Combination: small molecule + conventional" = "#56B4E9",
  "Non-drug / non-relevant only" = "#999999",
  "Unmatched / unclear" = "#666666"
)
save_plot <- function(p, name, width = 8, height = 5) {
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(figure_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 600)
}
# Figure 1. Trial-class distribution
p1 <- ggplot(
  summary_trial_class,
  aes(x = reorder(trial_class_label, -n_trials), y = n_trials, fill = trial_class_label)
) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = n_trials),
    hjust = -0.15,
    size = 4,
    fontface = "bold"
  ) +
  coord_flip() +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Distribution of completed phase II-IV psoriasis/PsA trials by therapeutic class",
    x = NULL,
    y = "Number of trials",
    fill = "Trial class"
  ) +
  pub_theme +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      title.theme = element_text(face = "bold", size = 13),
      label.theme = element_text(face = "bold", size = 11)
    )
  ) +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12)
  )
print(p1)
save_plot(p1, "Figure_1_trial_class_distribution", width = 12, height = 5.5)
# Figure 2. Phase-by-class distribution
p2 <- ggplot(
  summary_phase_by_class,
  aes(x = phase_clean, y = n_trials, fill = trial_class_label)
) +
  geom_col(position = "stack", color = "white", linewidth = 0.2) +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  labs(
    title = "Therapeutic-class distribution across clinical trial phases",
    x = "Clinical trial phase",
    y = "Number of trials",
    fill = "Trial class"
  ) +
  pub_theme +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      title.theme = element_text(face = "bold", size = 13),
      label.theme = element_text(face = "bold", size = 11)
    )
  ) +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(face = "bold", size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.x = unit(0.35, "cm")
  )
print(p2)
save_plot(p2, "Figure_2_phase_by_trial_class", width = 14, height = 5.5)
# Figure 3. Completion year by class
p3 <- ggplot(
  summary_year_by_class,
  aes(x = completion_year, y = n_trials, fill = trial_class_label)
) +
  geom_col(color = "white", linewidth = 0.15) +
  scale_fill_manual(values = class_colors, drop = FALSE) +
  scale_x_continuous(
    breaks = seq(
      min(summary_year_by_class$completion_year, na.rm = TRUE),
      max(summary_year_by_class$completion_year, na.rm = TRUE),
      by = 2
    )
  ) +
  labs(
    title = "Temporal distribution of completed psoriasis/PsA trials by therapeutic class",
    x = "Completion year",
    y = "Number of trials",
    fill = "Trial class"
  ) +
  pub_theme +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      title.theme = element_text(face = "bold", size = 13),
      label.theme = element_text(face = "bold", size = 11)
    )
  ) +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 11, angle = 45, hjust = 1),
    axis.text.y  = element_text(face = "bold", size = 12),
    
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(face = "bold", size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.x = unit(0.35, "cm")
  )
print(p3)
save_plot(p3, "Figure_3_completion_year_by_trial_class", width = 11, height = 5.5)
# Figure 4. Intervention-code distribution
summary_intervention_code <- summary_intervention_code |>
  mutate(
    class_label_plot = case_when(
      class_code == 0 ~ "Non-drug / non-relevant",
      class_code == 1 ~ "Small molecule",
      class_code == 2 ~ "Biologic",
      class_code == 3 ~ "Conventional systemic/topical",
      TRUE ~ class_label
    )
  )
p4 <- ggplot(
  summary_intervention_code,
  aes(x = factor(class_code), y = n_intervention_records, fill = class_label_plot)
) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  geom_text(
    aes(label = n_intervention_records),
    vjust = -0.25,
    size = 4,
    fontface = "bold"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Intervention-level classification based on numeric therapeutic codes",
    x = "Class code",
    y = "Number of intervention records",
    fill = "Class"
  ) +
  pub_theme +
  guides(
    fill = guide_legend(
      nrow = 2,
      byrow = TRUE,
      title.theme = element_text(face = "bold", size = 13),
      label.theme = element_text(face = "bold", size = 11)
    )
  ) +
  theme(
    axis.title.x = element_text(face = "bold", size = 14),
    axis.title.y = element_text(face = "bold", size = 14),
    axis.text.x  = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold", size = 13),
    legend.text  = element_text(face = "bold", size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.spacing.x = unit(0.35, "cm")
  )
print(p4)
save_plot(p4, "Figure_4_intervention_code_distribution", width = 9.5, height = 5)
# 14. Console report
cat("\n=====================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("=====================================================\n\n")
cat("Input file:", input_file, "\n")
cat("Total trials in CSV:", nrow(ctg_raw), "\n")
cat("Total intervention records:", nrow(classification_long), "\n")
cat("Unmatched intervention names:", nrow(unmatched_for_manual), "\n\n")
cat("Trial-level class distribution:\n")
print(summary_trial_class)
cat("\nIntervention-level code distribution:\n")
print(summary_intervention_code)
cat("\nFiles saved in:\n")
cat(output_dir, "\n\n")
cat("Most important manual QC file:\n")
cat(file.path(qc_dir, "UNMATCHED_interventions_for_manual_coding.csv"), "\n\n")
