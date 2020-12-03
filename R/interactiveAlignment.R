#' Run interactive alignment app
#'
#' This app allows to interactively explore the alignment of two snowprofiles, which are either
#' given as input to this function, or are uploaded to the app interactively as caaml files.
#' Example profiles are also provided in the app.
#'
#' @import shiny
#' @import grid
#' @import sarp.snowprofile
#'
#' @param query an optional query snowprofile
#' @param ref an optional reference snowprofile
#'
#' @return An interactive session will be started
#'
#' @examples
#' if (FALSE){  # this example won't be started in package tests.
#'
#' ## start app and choose profiles from within the app:
#' interactiveAlignment()
#'
#' ## start app with package internal profile data (from `sarp.snowprofile`):
#' interactiveAlignment(query = SPpairs$A_modeled, ref = SPpairs$A_manual)
#'
#' }
#'
#' @author fherla
#'
#' @export
interactiveAlignment <- function(query = NaN, ref = NaN) {


  ## --- handle input profiles ----
  if (!all(is.na(query)) && !all(is.na(ref))) {  # start app with input profiles
    values <- reactiveValues(providedInput = TRUE, query = query, ref = ref)
  } else {
    values <- reactiveValues(providedInput = FALSE)
  }
  ## --- UI ----
  ui <- fluidPage(
    withMathJax(),
    tags$head(tags$style(type="text/css", ".container-fluid {  max-width: 1300px; }"),
              tags$head(tags$script('
                                var dimension = [0, 0];
                                $(document).on("shiny:connected", function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                                $(window).resize(function(e) {
                                    dimension[0] = window.innerWidth;
                                    dimension[1] = window.innerHeight;
                                    Shiny.onInputChange("dimension", dimension);
                                });
                            '))
              ),
    headerPanel("DTW snow profile alignment"),
    fluidRow(
      column(4, offset = 0,
             ## Source - Select whether to browse
             conditionalPanel(condition = "!output.onInput",
                              selectInput(inputId = "source", label = "Source",
                                           choices = c("1) Variable local warping", "2) Pre/post storm", "3) Human vs modeled",
                                                       "Browse for your own CAAML files",
                                                       "Browse for your own PRF file"),
                                           selected = "1) Variable local warping"),
                              ## Action button
                              actionButton("do", "Load", style = "background-color: silver;  color:navy"),
                              br(),
                              hr()
                             ),
             ## Reverse checkbox
             checkboxInput("reverse", "Switch query and reference", value = FALSE),
             ## Scale profiles to equal heights before alignment
             checkboxInput("initialScaling", "Scale profiles to equal height before aligning", value = TRUE),
             ## resampling rate:
             checkboxInput("resampling", "Resample profiles onto regular grid", value = TRUE),
             conditionalPanel(condition = "input.resampling",
                              sliderInput("resamplingRate", "sampling rate (cm)", min = 0.25, max = 5, value = 0.5, step = 0.25)),
             ## use weighted grain Similarity matrix:
             checkboxInput("layerWeighting", "Apply a layer weighting scheme to the grain similarity matrix for preferential alignment; unweighted matrix is still used for similarity measure Phi", value = TRUE)
             ),


      column(8,
             ## Browse CAAML or PRF files:
             conditionalPanel(condition = "(!output.onInput && input.source == 'Browse for your own CAAML files')",
                              div(style = "display: inline-block; vertical-align: top;",
                                fileInput("queryFile", h5("Query"), accept = c(".caaml"))),
                              div(style = "display: inline-block; vertical-align: top; height: 1px; width: 25px"),
                              div(style = "display: inline-block; vertical-align: top;",
                                  fileInput("refFile", h5("Reference"), accept = c(".caaml")))
                              ),
             conditionalPanel(condition = "(!output.onInput && input.source == 'Browse for your own PRF file')",
                              div(style = "display: inline-block; vertical-align: top;",
                                  fileInput("prfFile", h5("PRF file"), accept = c(".prf"))),
                              helpText("Note, the PRF file needs to contain at least two profiles, the first and second of which will be the query and reference profiles.")
                              )
      ),
      column(12, hr())
    ),

    sidebarPanel(
      # tags$head(tags$style(type="text/css", "select { min-width: 300px; max-width: 300px;}"),
      #           tags$style(type="text/css", ".span4 { min-width: 320px; max-width: 320px;}"),
      #           tags$style(type="text/css", ".well { min-width: 300px; max-width: 300px;}")
      #           ),
      width = 3,
      checkboxInput("openEnd", "Open End alignment", value = TRUE),
      conditionalPanel(condition = "input.openEnd",
                       checkboxInput("checkGlobal", "Check global alignment", value = TRUE)),
      radioButtons("alignDir", "Direction of alignment",
                   choices = list("Bottom-up (BU)", "Top-down (TD)", "BU/TD"),
                   selected = "BU/TD"),
      checkboxInput("ddate", "Add deposition date info", value = FALSE),
      conditionalPanel(condition = "input.ddate",
                       sliderInput("ddateNorm",
                                   label = "Date normalization factor (unit: days)",
                                   min = 1, max = 20, value = 3)
      ),

      ## Set weights
      h4("Weights"),
      sliderInput(inputId = "weightSlider",
                  label = " grain type | hardness",
                           # div(style="width: 300px;",
                           #    div(style='float:left;', 'grain type'),
                           #    div(style='margin-left: 65%;', 'hardness')),
                  min = 0, max = 1, value = 0.6),

      conditionalPanel(
        condition = "input.ddate",
        sliderInput(inputId = "weightSlider2",
                    ## how to top-align following labels??
                    ## tried vertical-align:top without success
                    label = "{grain type & hardness} | ddate",
                            # div(style='width:300px;',
                            #     div(style='float:left;', 'grain type'),
                            #     div(style='margin-left:40%;', 'hardness'),
                            #     div(style='float:right;', 'ddate')
                            # ),
                    min = 0, max = 1, value = 0.7)
      ),
      ## Set warping window
      h4("Warping window"),
      sliderInput(inputId = "wsize", "percentage of layers/height",
                  min = 0, max = 1, value = 0.3),
      # conditionalPanel(condition = "input.ddate",
      #                  sliderInput(inputId = "dwsize", "number of days",
      #                              min = 5, max = 40, value = 40)),

      ## Step pattern
      h4("Local slope constraint"),
      helpText("symmetricP1 limits layer stretching and compressing to double/half the original thickness"),
      radioButtons(inputId = "stepPattern_m", label = "",
                  choices = list("unconstrained", "symmetricP1", "symmetricP2"), selected = "symmetricP1",
                  inline = FALSE)
    ),

    mainPanel(
      tags$style(HTML("
    .tabbable > .nav > li > a                  {background-color: silver;  color:navy}
    .tabbable > .nav > li[class=active]    > a {background-color: white; color:black}
    ")),
      tabsetPanel(
        tabPanel("Profile Alignment",
                 br(),
                 fixedRow(column(width = 5, HTML("Normalized DTW <b>distance</b>: ")),
                          column(width = 2, strong(textOutput("normDist")))),
                 fixedRow(column(width = 5, HTML("<b>similarity</b> measure Phi: ")),
                          column(width = 2, strong(textOutput("simSP")))),
                 br(),
                 fixedRow(column(width = 7, checkboxInput("verboseSim", "Print detailed similarity to console", value = FALSE),
                                 textOutput("verboseSim"))),
                 conditionalPanel(condition = 'input.ddate',
                                  br(),
                                  fixedRow(column(width = 6, checkboxInput("labelDdate",
                                                                           "Label deposition date", value = FALSE)))),
                 plotOutput("alignmentPlot", height = "650px")
        ),
        tabPanel("Cost Density & Warping Path",
                 br(),
                 fixedRow(column(width = 5, HTML("Normalized DTW <b>distance</b>: ")),
                          column(width = 2, strong(textOutput("normDistII")))),
                 fixedRow(column(width = 5, HTML("<b>similarity</b> measure Phi: ")),
                          column(width = 2, strong(textOutput("simSPII")))),
                 br(),
                 fixedRow(column(width = 3, offset = 3, radioButtons("labelHeight", "Units",
                                                         choices = list("Layer #" = FALSE, "Height (cm)" = TRUE),
                                                         selected = FALSE, inline = TRUE)),
                          column(width = 3, radioButtons("localCost", "Cost",
                                                         choices = list("Global" = FALSE, "Local" = TRUE),
                                                         selected = TRUE, inline = TRUE))),
                 plotOutput("costDensity",  height = "650px")
        )
      )
    )
  )
  ## --- server function ----
  server <- function(input, output, session) {

   ## ---- initialize reactive profile data ----
    ## dependend on action button or provided input
    isolate({
      if (values$providedInput) profiles <- reactiveValues(query = values$query, ref = values$ref)
      else profiles <- reactiveValues(query = NULL, ref = NULL)
    })

    ## ---- load profiles upon action button ----
    ## subroutine
    get_profiles <- eventReactive(input$do, {
      SPpairs <- sarp.snowprofile::SPpairs
      if (input$source == "3) Human vs modeled") {
        query <- SPpairs$A_modeled
        ref <- SPpairs$A_manual
      } else if (input$source == "2) Pre/post storm") {
        query <- SPpairs$C_day1
        ref <- SPpairs$C_day2
      } else if (input$source == "1) Variable local warping") {
        query <- SPpairs$D_generalAlignment1
        ref <- SPpairs$D_generalAlignment2
      } else if (input$source == "Browse for your own PRF file") {
        prfRead <- snowprofilePrf(input$prfFile$datapath)
        query <- prfRead[[1]]
        ref <- prfRead[[2]]
      } else {
        req(input$queryFile, input$refFile)
        query <- snowprofileCaaml(input$queryFile$datapath)
        ref <- snowprofileCaaml(input$refFile$datapath)
      }
      list(ref = ref, query = query)
    })
    ## update profiles with subroutine
    ## i.e. workaround to be able to start off with other values written into profiles!
    observeEvent(get_profiles(), {
      p <- get_profiles()
      profiles$query = p$query
      profiles$ref = p$ref
    })

    ## ---- store and update dims and weights ----
    properties <- reactiveValues()
    MINdwsize <- reactiveVal(0)
    observeEvent(
      c(input$weightSlider, input$weightSlider2, input$ddate),
      priority = 3,
      {
        if (input$ddate) {
          ## check if ddate info available in profiles:
          if (!"ddate" %in% names(isolate(profiles$query$layers)) |
              !"ddate" %in% names(isolate(profiles$ref$layers))) {
            showModal(modalDialog(
              title = "Oooops!",
              "At least one of your profiles doesn't have any deposition date information.",
              easyClose = TRUE
            ))
            ## change input$ddate to wrong and continue with subsequent if clause:
            updateCheckboxInput(session, "ddate", value = FALSE)
          } else {
            ## ddate info is available:
            properties$dims = c("gtype", "hardness", "ddate")
            properties$weights = c(gType = input$weightSlider2 * input$weightSlider,
                                   hardness = input$weightSlider2 * (1 - input$weightSlider),
                                   dDate = 1-input$weightSlider2)
            # ## calculate min dwsize which still yields a warping path:
            # ## i.e. find minimum dwsize for every layer and then take maximum of that vector
            # mv <- sapply(profiles$query$layers$ddate, function(x, y) min(abs(x - y), na.rm = TRUE),
            #              y = profiles$ref$layers$ddate)
            # MINdwsize(max(mv)+1)
          }
        }
        if(!input$ddate) {
          properties$dims = c("gtype", "hardness")
          properties$weights = c(gType = input$weightSlider,
                                 hardness = (1-input$weightSlider))
        }
      }
    )

    ## ---- store and UPDATE alignment ----
    ## update align whenever changes happen:
    align <- reactiveValues()
    observeEvent(

      ## UPDATE whenever detect changes in:
      c(profiles$ref, profiles$query, properties$weights, input$reverse, input$ddateNorm, input$wsize, input$dwsize,
        input$layerWeighting, input$initialScaling, input$openEnd, input$checkGlobal, input$stepPattern_m,
        input$resampling, input$resamplingRate, input$alignDir), {

        ## calculate profile alignment:
        ## requires profiles to continue
        req(profiles$query, profiles$ref)
        ## 'reverse' UI tick box
        if (input$reverse) profiles <- list(ref = profiles$query, query = profiles$ref)

        ## rescaling and resampling need to go hand in hand:
        if (input$resampling & !input$initialScaling) {
          updateCheckboxInput(session, "resampling", value = FALSE)
          showModal(modalDialog(
            title = "Sorry!",
            HTML(paste0("Resampling is not supported for profiles of different total snow height. Turning off resampling..")),
            easyClose = TRUE
          ))
          rrate <- NA
        } else {
          ## resampling rate:
          if (input$resampling) rrate <- input$resamplingRate
          else rrate <- NA
        }

        ## don't allow to go below MINdwsize:
        if (FALSE) {
          ## uncomment if date warping window dwsize is wanted:
          # (input$dwsize < MINdwsize()) {
          #   updateSliderInput(session, "dwsize", value =  MINdwsize())
          #   showModal(modalDialog(
          #     title = "Oha!",
          #     HTML(paste0("No warping path available, if you go below ",  MINdwsize(), " days of date window size!<br>
          #            Check out the cost density to better understand why.")),
          #     easyClose = TRUE
          #   ))
        } else {
          ## dwsize is large enough:

          ## choose grain similarity matrix:
          grainDist <- sim2dist(grainSimilarity_align(FALSE))
          if (input$layerWeighting) layerWeights <- layerWeightingMat(FALSE)
          else layerWeights <- NA

          ## choose step pattern:
          if (input$stepPattern_m == "unconstrained") stepPat <- symmetric2
          else if (input$stepPattern_m == "symmetricP1") stepPat <- symmetricP1
          else if (input$stepPattern_m == "symmetricP2") stepPat <- symmetricP2

          ## alignment direction:
          if (input$alignDir == "Bottom-up (BU)") {
            BU <- TRUE
            TD <- FALSE
          } else if (input$alignDir == "Top-down (TD)") {
            BU <- FALSE
            TD <- TRUE
          } else if (input$alignDir == "BU/TD") {
            BU <- TRUE
            TD <- TRUE
          }

          ## call to dtw:
          align$alignment <- tryCatch({
            dtwSP(query = mergeIdentLayers(profiles$query, properties = properties$dims),
                  ref = mergeIdentLayers(profiles$ref, properties = properties$dims),
                  grain_type_distMat = grainDist,
                  prefLayerWeights = layerWeights,
                  dims = properties$dims, weights = properties$weights,
                  resamplingRate = rrate,
                  rescale2refHS = input$initialScaling,
                  ddateNorm = input$ddateNorm,
                  windowFunction = warpWindowSP,
                  window.size = input$wsize, ddate.window.size = input$dwsize,
                  step.pattern = stepPat,
                  open.end = input$openEnd, checkGlobalAlignment = input$checkGlobal,
                  keep.internals = TRUE, bottom.up = BU, top.down = TD)
          }, error = function(err) {
            WSIZE <- input$wsize
            alignmentLoop <- function() {
              catch <- tryCatch({
                dtwSP(query = mergeIdentLayers(profiles$query, properties = properties$dims),
                      ref = mergeIdentLayers(profiles$ref, properties = properties$dims),
                      grain_type_distMat = grainDist,
                      prefLayerWeights = layerWeights,
                      dims = properties$dims, weights = properties$weights,
                      resamplingRate = rrate,
                      rescale2refHS = input$initialScaling,
                      ddateNorm = input$ddateNorm,
                      windowFunction = warpWindowSP,
                      window.size = WSIZE, ddate.window.size = input$dwsize,
                      step.pattern = stepPat,
                      open.end = input$openEnd, checkGlobalAlignment = input$checkGlobal,
                      keep.internals = TRUE, bottom.up = BU, top.down = TD)
              }, error = function(err) {
                return(NA)
              })
            }
            align$alignment <- alignmentLoop()
            while (all(is.na(align$alignment))) {
              if (WSIZE >= 1) break()
              WSIZE = WSIZE + 0.01
              align$alignment <- alignmentLoop()
            }
            if (WSIZE < 1) {
              updateSliderInput(session, "wsize", value = WSIZE)
              showModal(modalDialog(
                title = "Oha!",
                HTML(paste0("No warping path available, if you go below a window size of ", WSIZE, "!<br>
                     Check out the cost density to better understand why. <br> Switch between the axis units 'Layer #' and 'Height (cm)' if necessary.")),
                easyClose = TRUE
              ))
            } else {
              showModal(modalDialog(
                title = "Sorry!",
                HTML(paste0("Something has gone quite wrong and no alignment could be calculated. <br>
                            Can you reverse what you just did?")),
                easyClose = TRUE
              ))
            }

            return(align$alignment)
          })  # close tryCatch call to assign align$alignment
        }
    })

    ## ---- render alignment figure ----
    output$alignmentPlot <- renderPlot({
      ## requires alignment to continue:
      req(align$alignment)
      ## plot:
      label.ddate <- ifelse((input$labelDdate & input$ddate), TRUE, FALSE)
      plotSPalignment(query = NA, ref = NA, dtwAlignment = align$alignment, label.ddate = label.ddate,
                      keep.alignment = TRUE)
      ## write text into figure:
      weightText <- paste(sapply(seq(length(properties$dims)),
                                 function(x) paste(names(properties$weights)[x], "=", properties$weights[x])),
                          collapse = " | ")
      grid::grid.text(paste("weights:", weightText),
                      x = 0.45, y = 0.93,
                      gp = grid::gpar(fontsize=12, col="grey"))
    })

    ## ---- render density figure ----
    output$costDensity <- renderPlot({
      ## requires alignment to continue:
      req(align$alignment)
      ## plot:
      align$alignment$localCostMatrix[is.na(align$alignment$costMatrix)] <- NA
      plotCostDensitySP(align$alignment, localCost = as.logical(input$localCost),
                        labelHeight = as.logical(input$labelHeight), marginalPros = TRUE)
    })

    ## ---- render text ----
    ## note: next line's hack is necessary to use that identical output in two different tabs!
    output$normDist <- output$normDistII <-  renderText({
      ## requires alignment to continue:
      req(align$alignment)
      ## text:
      paste(formatC(align$alignment$normalizedDistance, format = "f", digits = 3))
    })

    output$simSP <- output$simSPII <- renderText({
      ## requires alignment to continue:
      req(align$alignment)
      req(align$alignment$queryWarped)
      ## text:
      paste(formatC(simSP(align$alignment$reference, align$alignment$queryWarped,
                              gtype_distMat = sim2dist(grainSimilarity_evaluate(triag = FALSE))),
                    format = "f", digits = 3))
    })

    output$verboseSim <- renderText({
      ## requires alignment
      req(align$alignment)
      req(align$alignment$queryWarped)
      ## verbose text:
      if (input$verboseSim) {
        tmp <- simSP(align$alignment$reference, align$alignment$queryWarped,
                  gtype_distMat = sim2dist(grainSimilarity_evaluate(triag = FALSE)),
                  verbose = TRUE)
        "(printed outside app)"
      }
    })

    ## ---- communication between UI and server ----
    ## needed in order to make action button disappear when input profiles are provided
    output$onInput <- reactive({
      if (values$providedInput) TRUE
      else FALSE
    })
    outputOptions(output, 'onInput', suspendWhenHidden = FALSE)
  }

  ## --- run app ----
  shinyApp(ui = ui, server = server)

}
