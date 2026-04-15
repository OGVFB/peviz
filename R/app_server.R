#' The application server-side
#'
#' @import shiny
#' @import dplyr
#' @import ggplot2
#' @import msa
#' @importFrom DT renderDT datatable renderDataTable
#' @importFrom data.table fread
#' @importFrom httr GET stop_for_status content accept
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom Biostrings readAAStringSet AAStringSet unmasked
#' @importFrom stringr str_split
#' @importFrom purrr map2
#' @importFrom tidyr unnest
#' @importFrom ggmsa tidy_msa
#' @importFrom ggforce facet_col
#' @importFrom ape nj
#' @importFrom seqinr dist.alignment write.fasta
#' @importFrom ggtree ggtree geom_tiplab
#' @importFrom utils URLencode
#' @importFrom stats na.omit
#' @noRd
app_server <- function(input, output, session) {

  options(shiny.maxRequestSize = 5000 * 1024^2)

  # --- REACTIVE VALUES FOR API SEARCH ----
  api_db_val <- reactiveVal(list(uniprotDB = NULL, proteins = NULL))
  api_search_results <- reactiveVal(NULL)

  # --- UNIPROT SEARCH LOGIC ----
  observeEvent(input$search_api, {
    req(input$api_query)

    query <- URLencode(input$api_query)
    url <- paste0("https://rest.uniprot.org/uniprotkb/search?query=", query,
                  "&format=tsv&fields=accession,id,gene_names,protein_name,organism_name,length&size=500")

    showNotification("Searching UniProt...", id = "uniprot_search", duration = NULL)

    tryCatch({
      df <- data.table::fread(url, stringsAsFactors = FALSE)
      api_search_results(df)
      removeNotification("uniprot_search")
    }, error = function(e) {
      removeNotification("uniprot_search")
      showNotification(paste("Search failed:", e$message), type = "error")
      api_search_results(NULL)
    })
  })

  output$api_results_table <- DT::renderDT({
    req(api_search_results())
    DT::datatable(api_search_results(), selection = 'multiple',
                  options = list(pageLength = 5, scrollX = TRUE))
  })

  # --- FETCH SELECTED SEQUENCES FROM API ----
  observeEvent(input$fetch_api_seqs, {
    req(api_search_results())
    selected_rows <- input$api_results_table_rows_selected

    if (length(selected_rows) == 0) {
      showNotification("Please select at least one row from the table.", type = "warning")
      return()
    }

    acc_col <- if ("Entry" %in% names(api_search_results())) "Entry" else names(api_search_results())[1]
    selected_accessions <- api_search_results()[[acc_col]][selected_rows]

    showNotification("Downloading FASTA sequences...", id = "fetch_seqs", duration = NULL)

    tryCatch({
      query_string <- paste0("accession:", paste(selected_accessions, collapse = "+OR+accession:"))
      fasta_url <- paste0("https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=", query_string)

      r <- httr::GET(fasta_url)
      httr::stop_for_status(r)
      fasta_text <- httr::content(r, as = "text", encoding = "UTF-8")

      tmp_fasta <- tempfile(fileext = ".fasta")
      writeLines(fasta_text, tmp_fasta)

      uniprotDB <- Biostrings::readAAStringSet(tmp_fasta)

      proteins <- uniprotDB@ranges %>% data.frame() %>% as_tibble(rownames = 'index') %>%
        mutate(names_formatted = parse_uniprot_names(names)) %>%
        rowwise() %>%
        mutate(id = stringr::str_split(names, '\\|| ')[[1]][2]) %>%
        ungroup()

      api_db_val(list(uniprotDB = uniprotDB, proteins = proteins))
      updateSelectizeInput(session, 'uniprot', choices = proteins$names, selected = proteins$names, server = TRUE)

      removeNotification("fetch_seqs")
      showNotification("Sequences loaded successfully!", type = "message")

    }, error = function(e) {
      removeNotification("fetch_seqs")
      showNotification(paste("Failed to fetch FASTA:", e$message), type = "error")
    })
  })

  # --- DATABASE LOAD ----
  db <- reactive({
    req(input$Database)
    if (input$Database == 'Swiss-Prot'){
      # Update path to look in the package installation directory
      data_path <- system.file("extdata", "peviz_uniprot_data.Rdata", package = "peviz")
      load(data_path)
      return(list(uniprotDB = uniprotDB, proteins = proteins))
    } else if (input$Database == 'Local UniProt Fasta') {
      req(input$local_fasta)
      datapath <- input$local_fasta$datapath
      if (tools::file_ext(datapath) == 'gz') {
        system(paste('gunzip', datapath))
        datapath <- gsub('.gz', '', datapath)
      }
      uniprotDB <- Biostrings::readAAStringSet(datapath)

      proteins <- uniprotDB@ranges %>% data.frame() %>% as_tibble(rownames = 'index') %>%
        mutate(names_formatted = parse_uniprot_names(names)) %>%
        rowwise() %>%
        mutate(id = stringr::str_split(names, '\\|| ')[[1]][2]) %>%
        ungroup()

      return(list(uniprotDB = uniprotDB, proteins = proteins))
    } else if (input$Database == 'UniProt API Search') {
      req(api_db_val()$uniprotDB)
      return(api_db_val())
    }
  })

  # --- DYNAMIC UI UPDATES ----
  observe({
    query <- parseQueryString(session$clientData$url_search)
    if (is.null(query[['uniprot']])){
      req(input$Database)
      if (input$Database != 'UniProt API Search') {
        updateSelectizeInput(session, 'uniprot', choices = db()$proteins$names, server = TRUE)
      }
    }
    if (is.null(query[['color']])){
      scheme_AA <- utils::getFromNamespace("scheme_AA", "ggmsa")
      valid_colors <- colnames(scheme_AA)[colnames(scheme_AA) != 'CN6']
      updateSelectizeInput(session, 'color', choices = valid_colors, selected = 'Chemistry_AA')
    }
  })

  observeEvent(input$uniprot, {
    updateSelectizeInput(session, 'primary_protein', choices = input$uniprot, server = TRUE)

    ref_selection <- if(is.null(input$reference_seq) || input$reference_seq == "") input$uniprot[1] else input$reference_seq
    updateSelectizeInput(session, 'reference_seq', choices = input$uniprot, selected = ref_selection, server = TRUE)
  })

  # ---- UNIPROT API FETCHING ----
  uniprot_features <- reactive({
    if (is.null(input$primary_protein) || input$primary_protein == '') return(NULL)

    selected_protein <- stringr::str_split(input$primary_protein, '\\|')[[1]][2]

    # 1. Fetch Standard Features
    req_url <- paste0("https://www.ebi.ac.uk/proteins/api/features?offset=0&size=-1&accession=", selected_protein)
    r_feat <- httr::GET(req_url, httr::accept("application/json"))

    if (r_feat$status_code != 200) return(NULL)
    feat_data <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r_feat)))$features[[1]]

    if (!is.null(feat_data) && nrow(feat_data) > 0 && "category" %in% names(feat_data)) {
      features <- feat_data %>% filter(toupper(type) != 'VARIANT', toupper(category) != 'VARIANTS')

      if ("category" %in% names(features)) features$category <- clean_list_col(features$category)
      if ("type" %in% names(features)) features$type <- clean_list_col(features$type)
      if ("description" %in% names(features)) features$description <- clean_list_col(features$description)
    } else {
      features <- data.frame()
    }

    # 2. Fetch Detailed Variations
    variants <- data.frame()

    if (input$fetch_variants) {
      var_url <- paste0("https://www.ebi.ac.uk/proteins/api/variation?offset=0&size=-1&accession=", selected_protein)
      r_var <- httr::GET(var_url, httr::accept("application/json"))

      if (r_var$status_code == 200) {
        var_data <- jsonlite::fromJSON(jsonlite::toJSON(httr::content(r_var)))$features[[1]]
        if (!is.null(var_data) && nrow(var_data) > 0) {
          variants <- var_data

          if ("consequenceType" %in% names(variants)) {
            variants$description <- clean_list_col(variants$consequenceType)
          } else {
            variants$description <- "Variant"
          }

          if ("clinicalSignificances" %in% names(variants)) {
            variants$clinical <- clean_list_col(variants$clinicalSignificances, extract_col = "type")
          } else {
            variants$clinical <- "Unknown"
          }

          valid_cols <- intersect(names(variants), c("begin", "end", "wildType", "mutatedType", "description", "clinical"))
          variants <- variants[, valid_cols, drop = FALSE]
          variants$category <- "VARIANTS"
          variants$type <- "VARIANT"
          variants$clinical[is.na(variants$clinical)] <- 'unknown'
        }
      }
    }

    # 3. Combine
    combined <- bind_rows(features, variants)
    return(combined)
  })

  observeEvent(uniprot_features(), {
    req(features <- uniprot_features())
    valid_categories <- na.omit(unique(features$category))
    updateSelectizeInput(session, 'uniprot_category', choices = valid_categories, selected = 'DOMAINS_AND_SITES')
  })

  observeEvent(c(uniprot_features(), input$uniprot_category), {
    req(features <- uniprot_features())
    if (length(input$uniprot_category) > 0) {
      features <- features %>% filter(category %in% input$uniprot_category)
    }
    valid_types <- na.omit(unique(features$type))
    updateSelectizeInput(session, 'uniprot_type', choices = valid_types, selected = character(0))
  })

  # ---- MULTIPLE SEQUENCE ALIGNMENT ----
  draw_msa <- eventReactive(input$Draw, {
    req(input$uniprot)
    indi <- db()$proteins %>% filter(names %in% input$uniprot) %>% pull(index) %>% as.integer()

    scheme_AA <- utils::getFromNamespace("scheme_AA", "ggmsa")
    color_pick <- scheme_AA[, input$color]
    names(color_pick) <- row.names(scheme_AA)
    is_single_seq <- length(indi) == 1

    if (is_single_seq) {
      seqs <- db()$uniprotDB[indi]
      order_full <- parse_uniprot_names(names(seqs))
      tidy_msa_df <- ggmsa::tidy_msa(seqs)
      fasta_out <- tibble(names = names(seqs), fasta = as.character(seqs))
    } else {
      msa_align <- msa(db()$uniprotDB[indi], method = input$method)
      order_full <- c('Consensus', parse_uniprot_names(names(Biostrings::unmasked(msa_align))))
      consensus <- Biostrings::AAStringSet(x = gsub("?", "X", msaConsensusSequence(msa_align), fixed = TRUE))
      names(consensus) <- 'Consensus'
      all_seqs <- c(msa_align@unmasked, consensus)
      tidy_msa_df <- ggmsa::tidy_msa(all_seqs)
      fasta_out <- tibble(names = names(all_seqs), fasta = as.character(all_seqs))
    }

    ref_name <- if (is.null(input$reference_seq) || input$reference_seq == "" || !(input$reference_seq %in% tidy_msa_df$name)) tidy_msa_df$name[1] else input$reference_seq

    ref_mapping <- tidy_msa_df %>%
      filter(name == ref_name) %>%
      arrange(position) %>%
      mutate(is_aa = character != "-", ref_pos = cumsum(is_aa)) %>%
      select(position, ref_pos, is_aa)

    tidy_msa_df <- tidy_msa_df %>% left_join(ref_mapping, by = "position")
    order_full <- c(' Sequence Numbering', order_full)

    # --- APPLY ANNOTATIONS & VARIANTS ----
    features <- uniprot_features()
    has_valid_annotations <- FALSE

    select_protein_aa <- data.frame()
    variant_counts <- data.frame()

    if (!is.null(features)) {
      if (length(input$uniprot_category) > 0) features <- features %>% filter(category %in% input$uniprot_category)
      if (length(input$uniprot_type) > 0) features <- features %>% filter(type %in% input$uniprot_type)

      variant_features <- features %>% filter(toupper(type) == "VARIANT")
      standard_features <- features %>% filter(toupper(type) != "VARIANT")

      if (nrow(standard_features) > 0) {
        expand_seq <- standard_features %>%
          mutate(
            begin_num = suppressWarnings(as.numeric(begin)),
            end_num = suppressWarnings(as.numeric(end))
          ) %>%
          filter(!is.na(begin_num) & !is.na(end_num)) %>%
          select(evidences, begin_num, end_num, type, description) %>%
          mutate(seq = purrr::map2(begin_num, end_num, ~ seq(.x, .y, by = 1))) %>%
          tidyr::unnest(cols = seq) %>%
          distinct()

        select_protein_aa <- tidy_msa_df %>%
          filter(name == input$primary_protein) %>%
          arrange(position) %>%
          mutate(is_aa = character != "-", true_seq = cumsum(is_aa)) %>%
          filter(is_aa) %>%
          left_join(expand_seq, by = c('true_seq' = 'seq'), relationship = "many-to-many") %>%
          filter(!is.na(evidences), evidences != 'NULL') %>%
          mutate(Protein = paste0(type, ': ', description), character = '^') %>%
          select(Protein, position, character)
      }

      if (nrow(variant_features) > 0) {
        if (!"clinical" %in% names(variant_features)) {
          variant_features$clinical <- "Unknown"
        }

        clean_variants <- variant_features %>%
          mutate(begin_num = suppressWarnings(as.numeric(begin))) %>%
          filter(!is.na(begin_num)) %>%
          select(begin_num, mutatedType, description, clinical) %>%
          distinct() %>%
          mutate(is_pathogenic = grepl("pathogenic", tolower(clinical)))

        variant_counts <- tidy_msa_df %>%
          filter(name == input$primary_protein) %>%
          arrange(position) %>%
          mutate(is_aa = character != "-", true_seq = cumsum(is_aa)) %>%
          filter(is_aa) %>%
          left_join(clean_variants, by = c('true_seq' = 'begin_num'), relationship = "many-to-many") %>%
          filter(!is.na(description)) %>%
          group_by(position, is_pathogenic) %>%
          summarise(count = n(), .groups = 'drop') %>%
          mutate(
            Protein = if_else(is_pathogenic, "Variants (Pathogenic)", "Variants (Other)"),
            character = as.character(count)
          ) %>%
          select(Protein, position, character)
      }

      has_valid_annotations <- nrow(select_protein_aa) > 0 || nrow(variant_counts) > 0
    }

    plot_data <- format_positions(tidy_msa_df) %>%
      mutate(Type = 'Amino Acid') %>%
      left_join(db()$proteins[indi, ], by = c('name' = 'names')) %>%
      mutate(Protein = parse_uniprot_names(name)) %>%
      rename(AA = character)

    if (nrow(select_protein_aa) > 0) {
      order_full <- c(order_full, unique(select_protein_aa$Protein))
      annot_data <- format_positions(select_protein_aa %>% mutate(Type = 'UniProt Annotation')) %>% rename(AA = character)
      plot_data <- bind_rows(plot_data, annot_data)
    }

    if (nrow(variant_counts) > 0) {
      var_levels <- c("Variants (Pathogenic)", "Variants (Other)")
      var_levels_present <- intersect(var_levels, unique(variant_counts$Protein))
      order_full <- c(order_full, var_levels_present)

      var_data <- format_positions(variant_counts %>% mutate(Type = 'UniProt Annotation')) %>% rename(AA = character)
      plot_data <- bind_rows(plot_data, var_data)
    }

    pos_track <- plot_data %>%
      filter(Type == 'Amino Acid') %>%
      select(position, group, Position, ref_pos, is_aa) %>%
      distinct() %>%
      mutate(
        Protein = ' Sequence Numbering', Type = 'Amino Acid', AA = NA,
        pos_label = if_else(!is.na(ref_pos) & ref_pos %% 10 == 0 & is_aa, as.character(ref_pos), "")
      )

    plot_data <- bind_rows(plot_data %>% mutate(pos_label = NA), pos_track) %>%
      mutate(Protein = factor(Protein, levels = rev(order_full)))

    # --- PLOTTING ----
    multiplier <- if(is_single_seq) 80 else if(has_valid_annotations) 35 else 25
    height <- nrow(unique(select(plot_data, group, Protein))) * multiplier

    base_plot <- ggplot(plot_data, aes(x = Position, y = Protein))
    if (has_valid_annotations) {
      base_plot <- base_plot + ggforce::facet_col(~ group + Type, scales = 'free_y', space = 'free')
    } else {
      base_plot <- base_plot + facet_wrap(~ group, scales = 'free_y', ncol = 1)
    }

    ref_split <- stringr::str_split(ref_name, "\\|")[[1]]
    clean_ref <- if(length(ref_split) >= 3) ref_split[3] else ref_name
    axis_label <- paste("Reference numbering:", clean_ref)

    plot <- base_plot +
      geom_raster(aes(fill = AA), data = filter(plot_data, Protein != ' Sequence Numbering')) +
      geom_text(aes(label = AA), data = filter(plot_data, Protein != ' Sequence Numbering')) +
      geom_text(aes(label = pos_label), data = filter(plot_data, Protein == ' Sequence Numbering'), size = 3.5, fontface = "bold") +
      scale_fill_manual(values = color_pick, na.value = "transparent", na.translate = FALSE) +
      labs(x = axis_label) +
      theme_minimal() +
      theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
            strip.background = element_blank(), strip.text = element_blank())

    list(plot = plot, height = height, fasta = fasta_out)
  })

  # --- PHYLOGENETIC TREE ----
  draw_tree <- eventReactive(input$Draw, {
    indi <- db()$proteins %>% filter(names %in% input$uniprot) %>% pull(index) %>% as.integer()
    if (length(indi) < 3) {
      return(ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Select at least 3 sequences to draw a phylogenetic tree.", color = "grey50"))
    }

    msa_align <- msaConvert(msa(db()$uniprotDB[indi], method = input$method), type="seqinr::alignment")
    msa_align$nam <- parse_uniprot_names(msa_align$nam)

    Tree <- ape::nj(seqinr::dist.alignment(msa_align, "identity"))
    ggtree::ggtree(Tree) + ggtree::geom_tiplab() + coord_cartesian(clip="off") + theme(plot.margin = margin(1,8,1,1, "cm"))
  })

  # --- OUTPUTS ----
  output$msa <- renderPlot({ draw_msa()$plot }, height = eventReactive(input$Draw, { draw_msa()$height }))
  output$tree <- renderPlot({ draw_tree() })
  output$fasta_table <- DT::renderDataTable({
    draw_msa()$fasta %>% mutate(fasta = gsub('-', '', fasta)) %>% DT::datatable()
  })

  # --- DOWNLOADS ----
  output$fasta_aligned_download <- downloadHandler(
    filename = function() paste0('data-', Sys.Date(), '_aligned.fasta'),
    content = function(con) {
      fa <- draw_msa()$fasta %>% filter(names != 'Consensus')
      seqinr::write.fasta(as.list(fa$fasta), fa$names, file.out = con)
    }
  )

  output$fasta_download <- downloadHandler(
    filename = function() paste0('data-', Sys.Date(), '.fasta'),
    content = function(con) {
      fa <- draw_msa()$fasta %>% filter(names != 'Consensus')
      seqinr::write.fasta(as.list(gsub('-', '', fa$fasta)), fa$names, file.out = con)
    }
  )
}
