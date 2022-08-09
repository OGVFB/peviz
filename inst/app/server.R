library(shiny)
library(msa)
library(dplyr)
library(httr)
library(jsonlite)
library(xml2)
library(ggplot2)
# input in uniprot db
#uniprotDB <- Biostrings::readAAStringSet('./data/uniprot_sprot.fasta.gz')
# extract protein names
# proteins <- uniprotDB@ranges %>% data.frame() %>% as_tibble(rownames = 'index') %>%
#   mutate(org = stringr::str_extract(names, 'OS.*OX') %>% gsub('OS\\=| OX','',.),
#          gene = gsub('\\w+\\|\\w+\\|','', names) %>% gsub(' .*','',.)) %>%
#   rowwise() %>%
#   mutate(id = stringr::str_split(names, '\\|| ')[[1]][2]) %>%
#   ungroup()
# save(proteins, uniprotDB, file = 'inst/app/data/peviz_uniprot_data.Rdata', compress = FALSE)
load('./data/peviz_uniprot_data.Rdata')


server <- function(input, output, session) {
  observe({
    query <- parseQueryString(session$clientData$url_search)

    if (is.null(query[['uniprot']])){
      updateSelectizeInput(session, 'uniprot',
                           choices = proteins$names,
                           options = list(placeholder = 'Type to search'),
                           selected = '',
                           server = TRUE)
    }
    observeEvent(input$uniprot, {
      updateSelectizeInput(session, 'primary_protein',
                           choices = input$uniprot,
                           options = list(placeholder = 'Optional'),
                           selected = '',
                           server = TRUE)
    })
  })



  draw_msa <- eventReactive(input$Draw, {
    # select proteins
    indi <- proteins %>%  filter(names %in% input$uniprot) %>% pull(index)
    msa_align <- msa(uniprotDB[as.integer(indi)])

    # build color
    color_pick <- ggmsa:::scheme_AA$Chemistry_AA
    names(color_pick) <- row.names(ggmsa:::scheme_AA)

    # extract msa order
    order <- unmasked(msa_align) %>% data.frame() %>% row.names()
    org_order <- stringr::str_extract(order, 'OS.*OX') %>% gsub('OS\\=| OX','',.)
    gene_order <- stringr::str_split(order, '\\|| ') %>% purrr::map(3) %>% unlist()
    id_order <- stringr::str_split(order, '\\|| ')%>% purrr::map(2) %>% unlist()
    order_full <- paste(org_order, gene_order, id_order, sep = ' | ')

    # if uniprot annotations desired
    if (input$primary_protein != ''){
      cat(input$primary_protein)
      selected_protein <- stringr::str_split(input$primary_protein,
                                             '\\|')[[1]][2]

      requestURL <- paste0("https://www.ebi.ac.uk/proteins/api/features?offset=0&size=-1&accession=", selected_protein)
      r <- GET(requestURL, accept("application/json"))

      #stop_for_status(r)

      json <- toJSON(content(r))
      uniprot_domains <- fromJSON(json)$features[[1]]  %>% filter(category == 'DOMAINS_AND_SITES')
      # expand sequence range
      expand_seq <- uniprot_domains %>%
        dplyr::select(evidences, begin, end) %>%
        transmute(evidences, seq = purrr::map2(begin, end, seq, by = 1)) %>%
        tidyr::unnest(cols = seq) %>%
        distinct() %>%
        data.frame() %>%
        left_join(uniprot_domains %>% dplyr::select(type, description, evidences), by = 'evidences')
      # create new data frame with UniProt annotation
      select_protein_aa <- ggmsa::tidy_msa(msa_align@unmasked ) %>% filter(name == input$primary_protein)
      select_protein_aa <- select_protein_aa %>% mutate(seq = row_number()) %>% left_join(expand_seq, by = 'seq') %>%
        filter(evidences != 'NULL') %>%
        mutate(name = paste0(type, ': ', description)) %>%
        dplyr::select(name, position) %>%
        mutate(character = '^')
      plot_data <- ggmsa::tidy_msa(msa_align@unmasked ) %>%
        mutate(Type = 'Amino Acid') %>%
        bind_rows(select_protein_aa %>% mutate(Type = 'UniProt Annotation')) %>%
        mutate(group =
                 substr(formatC(position, width = 5, format = "d", flag = "0"), 1, 3),
               Position =
                 substr(formatC(position, width = 5, format = "d", flag = "0"), 4, 5)) %>%
        left_join(., proteins[indi,], by = c('name' = 'names')) %>%
        mutate(Protein = paste(org, gene, id, sep = ' | '),
               Protein = factor(Protein, levels = rev(order_full))) %>%
        rename(AA = character)
      plot_data$group <- factor(plot_data$group)
      levels(plot_data$group) = paste0(plot_data$group, " (Hundreds Position)") %>% unique()
      plot <- plot_data %>%
        ggplot(aes(x=Position,y=Protein, label = AA, fill = AA)) +
        ggforce::facet_col(~ group + Type, scales = 'free_y', space = 'free')
    } else {
      plot_data <- ggmsa::tidy_msa(msa_align@unmasked ) %>%
        mutate(group =
                 substr(formatC(position, width = 5, format = "d", flag = "0"), 1, 3),
               Position =
                 substr(formatC(position, width = 5, format = "d", flag = "0"), 4, 5)) %>%
        left_join(., proteins[indi,], by = c('name' = 'names')) %>%
        mutate(Protein = paste(org, gene, id, sep = ' | '),
               Protein = factor(Protein, levels = rev(order_full))) %>%
        rename(AA = character)
      plot_data$group <- factor(plot_data$group)
      levels(plot_data$group) = paste0(plot_data$group, " (Hundreds Position)") %>% unique()
      plot <- plot_data %>%
        ggplot(aes(x=Position,y=Protein, label = AA, fill = AA)) +
        facet_wrap(~group, scales = 'free_y', ncol = 1)
    }


    #output <- list
    #output$msa <- msa_align
    #output$plot <- plot
    #output
    plot +
      geom_raster() +
      geom_text() +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())  +
      scale_fill_manual(values = color_pick) +
      scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
      xlab("Position (tens)")
  })

  draw_tree <- eventReactive(input$Draw, {
    indi <- proteins %>%  filter(names %in% input$uniprot) %>% pull(index)
    msa_align <- msa(uniprotDB[as.integer(indi)])
    options(ignore.negative.edge=TRUE)
    msa_align2 <- msaConvert(msa_align, type="seqinr::alignment")
    msa_align2$nam <- stringr::str_extract(msa_align2$nam, 'OS.*OX') %>% gsub('OS\\=| OX','',.)
    d <- seqinr::dist.alignment(msa_align2, "identity")
    Tree <- ape::nj(d)
    ggtree::ggtree(Tree) + ggtree::geom_tiplab()
  })

  output$msa <- renderPlot({
    draw_msa()
  })

  output$tree <- renderPlot({
    draw_tree()
  })

}
