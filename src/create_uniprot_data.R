library(dplyr)
system('wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz')
uniprotDB <- Biostrings::readAAStringSet('uniprot_sprot.fasta.gz')
#extract protein names
proteins <- uniprotDB@ranges %>% data.frame() %>% as_tibble(rownames = 'index') %>%
  mutate(org = stringr::str_extract(names, 'OS.*OX') %>% gsub('OS\\=| OX','',.),
         gene = gsub('\\w+\\|\\w+\\|','', names) %>% gsub('_.*','',.)) %>%
  rowwise() %>%
  mutate(id = stringr::str_split(names, '\\|| ')[[1]][2]) %>%
  ungroup()
save(proteins, uniprotDB, file = 'inst/app/data/peviz_uniprot_data.Rdata', compress = FALSE)
