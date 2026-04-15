#' Run the Shiny Application
#'
#' @param ... arguments to pass to shinyApp()
#' @export
run_app <- function(...) {
  shiny::shinyApp(ui = app_ui, server = app_server, ...)
}
