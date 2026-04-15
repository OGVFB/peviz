# R/globals.R

# Suppress R CMD check notes for variables used in dplyr NSE (Non-Standard Evaluation)
utils::globalVariables(c(
  "AA", "Protein", "Type", "begin", "begin_num", "category", "clinical",
  "description", "end", "end_num", "evidences", "fasta", "group", "index",
  "is_aa", "is_pathogenic", "mutatedType", "name", "pos_label", "position",
  "ref_pos", "type", "seq", "true_seq", "count", "names_formatted", "id",
  "wildType", "consequenceType", "clinicalSignificances", "Position", "character"
))
