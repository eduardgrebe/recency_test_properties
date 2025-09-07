# Copyright 2022 UNAIDS, World Health Organization, and individual contributors.
# Author: Eduard Grebe
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.  This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.  You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

library(shiny)
library(tidyverse)
library(arrow)
library(foreach)
library(parallel)
library(doParallel)
library(future)
library(promises)
library(ipc)

source("R/utility_functions.R")
source("R/mdri_adjustment.R")

# Load the dataset once across sessions
mdridat <- arrow::read_parquet("data/mdri_results.parquet")
frrdat <- arrow::read_parquet("data/PRt_fits_alltimes.parquet")

# Set up cluster for parallel computation
#cluster <- parallel::makeCluster(parallel::detectCores(), outfile="")
#doParallel::registerDoParallel(cluster)
plan(multisession, workers = availableCores()-2)


server <- function(input, output, session) {

  observeEvent(input$help_rita, {
    showModal(modalDialog(
      p("In order to estimate the mean duration of recent infection (MDRI), it is necessary to specify the recent infection testing algorithm (RITA). The tool looks up calibration data for the specified RITA and then applies adjustments for the subtype mix in the population of interest, and for early diagnosis and treatment, both of which can impact MDRI."),
      tags$b("Time cutoff:"),
      p("A cut-off time T must be chosen for a RITA, in terms of which the properties of the RITA, including MDRI, are defined."),
      tags$b("Primary assay:"),
      p("This is the primary recency assay, usually an immunoassay, by which infections are classified as 'recent' or 'longstanding'."),
      tags$b("Primary assay and viral load thresholds:"),
      p("Different thresholds are associated with different MDRIs. Calibration data is used to look up the MDRI associated with a particular threshold combination. Even without prior diagnosis screening or ARV detection, the distribution of times from infection to viral suppression must be adjusted for when viral load is included in the RITA, since calibration data is obtained from untreated individuals."),
      tags$b("Prior diagnosis:"),
      p("Please indicate whether the RITA will classify HIV-infected individuals who previously received an HIV diagnosis as longstanding. MDRI calibration data do not account for this additional 'route' out of the recent state, and if individuals with previous diagnoses are classified as longstanding, the distribution of times from infection to diagnosis must be provided to adjust for it. Shorter times to diagnosis are associated with a reduction in the effective MDRI."),
      tags$b("ARV detection:"),
      p("Please indicate whether the RITA will test for the presence of ARVs. MDRI calibration data do not account for the use of ARV detection in RITAs, and if individuals on treatment are classified as longstanding, the distribution of times from infection to treatment must be provided to adjust for it. Shorter times to treatment initiation are associated with a reduction in the effective MDRI."),
      tags$b("Screening sensitivity:"),
      p("The unadjusted MDRIs stored in this tool are expressed relative to an HIV case definition of 'detectable by a viral load assay with a detection threshold of 1 RNA copy/mL. Any practical survey must therefore use an MDRI adjusted for the window period of the least sensitive screening assay on which positivity is required, typically a serological assay. Please select the screening assay class or provide a relative window period in days."),

      title = "Specifying the RITA",
      footer = modalButton("Dismiss"),
      easyClose = TRUE,
      size = "l"
    ))
  })

  observeEvent(input$help_subtype, {
    showModal(modalDialog(
      p("MDRI has been demonstrated to vary according to HIV subtype. While subtype confirmation is not logistically feasible in most settings, external data can be used to specify the subtype mix and the population where a survey will be conducted. Use this section to specify the approximate proportions of the HIV-infected population infected with each HIV-1 subtype. Unfortunately calibration data are only available for the subtypes listed."),

      title = "Adjusting for subtype distribution",
      footer = modalButton("Dismiss"),
      easyClose = TRUE,
      size = "l"
    ))
  })

  observeEvent(input$help_ttd, {
    showModal(modalDialog(
      p(""),

      title = "Adjusting for the distribution of times from infection to diagnosis or treatment",
      footer = modalButton("Dismiss"),
      easyClose = TRUE,
      size = "l"
    ))
  })

  observeEvent(input$help_frr, {
    showModal(modalDialog(
      p("FRR is highly dependent on treatment coverage. When neither viral load nor ARV exposure testing is included in a RITA, the FRR in treated invidivuals can be extremely high (above 50%). It is not recommended to use a RITA without viral load or ARV testing and this tool does not support estimating the context-specific FRR for such a RITA. Further note that the FRR is estimated in untreated individuals (weighted for the distribution of times from infection to diagnosis), and assumed to be zero in treated individuals. In the CEPHIA datasets underlying this tool all treated individuals are virally suppressed. In reality some individuals are treated but not virally suppressed, but the FRR in this group is not known."),
      p(),
      p("Please note that it is up to the user to ensure that the time-to-diagnosis/treatment weighting function is consistent with the specified treatment coverage. For example, the tool will accept inputs specifying the treatment coverage as 20% and a median time to diagnosis of 1 year, but this combination is exceedingly unlikely in reality. If both the treatment coverage and time-to-diagnosis parameters specified are consistent with empiric data or epidemiological modelling (e.g., Shiny90), the inputs are likely to be consistent."),

      title = "Treatment coverage for FRR estimation",
      footer = modalButton("Dismiss"),
      easyClose = TRUE,
      size = "l"
    ))
  })

  v <- reactiveValues(
    adjusted_mdri_ci = NULL,
    fully_adjusted_mdri_ci = NULL,
    frr_untreated_ci = NULL,
    frr_weighted_ci = NULL
  )

  observe({
    set_to_no <- input$vl_threshold == 0 & input$diagnosed_screened_out == "no" & input$arv_detection == "no"
    if(set_to_no) {
      updateRadioButtons(
        session,
        "adjust_ttd",
        selected = "no"
      )
    }
  })

  observe({
    median <- input$median_ttd
    updateSliderInput(
      session,
      "quantile_0_25",
      value = 0.75 * median,
      min = 0.1,
      max = median,
      step = 0.05
    )
  })

  observe({
    total_nonC <- input$subtype_a_perc + input$subtype_d_perc + input$subtype_b_perc
    updateSliderInput(
      session,
      "subtype_c_perc",
      value = min(input$subtype_c_perc, 100 - total_nonC),
      min = 0,
      max = 100 - total_nonC,
      step = 1
    )
    total_nonA <- input$subtype_c_perc + input$subtype_d_perc + input$subtype_b_perc
    updateSliderInput(
      session,
      "subtype_a_perc",
      value = min(input$subtype_a_perc, 100 - total_nonA),
      min = 0,
      max = 100 - total_nonA,
      step = 1
    )
    total_nonD <- input$subtype_c_perc + input$subtype_a_perc + input$subtype_b_perc
    updateSliderInput(
      session,
      "subtype_d_perc",
      value = min(input$subtype_d_perc, 100 - total_nonD),
      min = 0,
      max = 100 - total_nonD,
      step = 1
    )
    total_nonB <- input$subtype_c_perc + input$subtype_a_perc + input$subtype_d_perc
    updateSliderInput(
      session,
      "subtype_b_perc",
      value = min(input$subtype_b_perc, 100 - total_nonB),
      min = 0,
      max = 100 - total_nonB,
      step = 1
    )
  })



  # reactive(
  #   if ("Other" %in% input$subtypes) {
  #     showNotification("Calibration data only available for subtypes A, B, C and D.", type = "warning", duration = 10)
  #   }
  # )


  wid <- NULL
  observe({
    #browser()
    subtype_total <- input$subtype_a_perc + input$subtype_b_perc + input$subtype_c_perc + input$subtype_d_perc
    if (input$weight_by_subtype == "yes" & subtype_total != 100 & is.null(wid)) {
      wid <<- showNotification(
        "Subtype percentages do not add up to 100%",
        type = "warning",
        duration = 0,
        closeButton = FALSE
      )
    }
    if ((input$weight_by_subtype == "no" | subtype_total == 100) & !is.null(wid)) {
      removeNotification(wid)
      wid <<- NULL
    }
  })

  # Show warning if subtype proportions do not add up
  # wid <- NULL
  # observeEvent(input$subtypes, {
  #   if ("Other" %in% input$subtypes & is.null(wid)) {
  #     wid <<- showNotification(
  #       "Calibration data only available for subtypes A, B, C and D.",
  #       type = "error",
  #       duration = 0,
  #       closeButton = FALSE
  #       )
  #   }
  #   if (!("Other" %in% input$subtypes) & !is.null(wid)) {
  #     removeNotification(wid)
  #     wid <<- NULL
  #   }
  # })

  get_threshold <- reactive({
    "LAg-Sedia"
    "LAg-Maxim"
    "BioRadAvidity-CDC"
    "BED"
    "Asante Electronic"
    threshold <- ifelse(
      input$assay == "LAg-Sedia",
      input$assay_threshold_lags,
      ifelse(
        input$assay == "LAg-Maxim",
        input$assay_threshold_lagm,
        ifelse(
          input$assay == "BioRadAvidity-CDC",
          input$assay_threshold_bra,
          ifelse(
            input$assay == "BED",
            input$assay_threshold_bed,
            ifelse(
              input$assay == "Asante Electronic",
              input$assay_threshold_asante,
              NULL
            )
          )
        )
      )
    )
    return(threshold)
  })

  filter_assay <- reactive({
    #browser()

    assaydat <- filter_mdri_data_assay(
      df = mdridat,
      assay = input$assay,
      threshold = get_threshold(),
      viral_load_threshold = as.numeric(input$vl_threshold),
      time_cutoff = as.numeric(input$bigT)
    )
    return(assaydat)
  })

  filter_assay_frr <- reactive({
    assaydat <- filter_frr_data_assay(
      df = frrdat,
      assay = input$assay,
      threshold = get_threshold(),
      viral_load_threshold = as.numeric(input$vl_threshold),
      time_cutoff = as.numeric(input$bigT)
    )
    return(assaydat)
  })

  # observe({
  #   mdridat_assay <- filter_mdri_data_assay(
  #     df = mdridat,
  #     assay = input$assay,
  #     threshold = input$assay_threshold,
  #     viral_load_threshold = as.numeric(input$vl_threshold),
  #     time_cutoff = as.numeric(input$bigT)
  #   )
  #   makeReactiveBinding("mdridat_assay")
  # })

  filter_subtype <- reactive({
    if (input$weight_by_subtype == "yes") {
      subtype_vector <- c()
      if (input$subtype_a_perc > 0) subtype_vector <- append(subtype_vector, "A1")
      if (input$subtype_b_perc > 0) subtype_vector <- append(subtype_vector, "B")
      if (input$subtype_c_perc > 0) subtype_vector <- append(subtype_vector, "C")
      if (input$subtype_d_perc > 0) subtype_vector <- append(subtype_vector, "D")
      if (is.null(subtype_vector)) {
        subtypedat <- filter_mdri_data_subtype(
          df = filter_assay(),
          subtypes = "All"
        )
      } else {
        subtypedat <- filter_mdri_data_subtype(
          df = filter_assay(),
          subtypes = subtype_vector
        )
      }
    } else if (input$weight_by_subtype == "no") {
      subtypedat <- filter_mdri_data_subtype(
        df = filter_assay(),
        subtypes = "All"
      )
    }
    return(subtypedat)
  })

  # observe({
  #   if (input$weight_by_subtype == "yes") {
  #     subtype_vector <- c()
  #     if (input$subtype_a_perc > 0) subtype_vector <- append(subtype_vector, "A1")
  #     if (input$subtype_b_perc > 0) subtype_vector <- append(subtype_vector, "B")
  #     if (input$subtype_c_perc > 0) subtype_vector <- append(subtype_vector, "C")
  #     if (input$subtype_d_perc > 0) subtype_vector <- append(subtype_vector, "D")
  #     if (is.null(subtype_vector)) {
  #       mdridat_subtype <- filter_mdri_data_subtype(
  #         df = mdridat_assay,
  #         subtype == "All"
  #       )
  #     } else {
  #       mdridat_subtype <- filter_mdri_data_subtype(
  #         df = mdridat_assay,
  #         subtype == subtype_vector
  #       )
  #     }
  #   } else if (input$weight_by_subtype == "no") {
  #       mdridat_subtype <- filter_mdri_data_subtype(
  #         df = mdridat_assay,
  #         subtype == "All"
  #       )
  #   }
  #   makeReactiveBinding("mdridat_subtype")
  # })

  output$mdri_table_assay <- renderTable({
    mdridat_assay <- filter_assay()
    return(mdridat_assay)
  })

  output$mdri_table_subtype <- renderTable({
    mdridat_subtype <- filter_subtype()
    return(mdridat_subtype)
  })

  # output$final_mdri < renderText({
  #   if (input$weight_by_subtype == "yes") {
  #     return("Must weight")
  #   } else {
  #     return("Must not weight")
  #   }
  # })

  # output$distPlot <- renderPlot({
  #   hist(rnorm(input$obs), col = "darkgray", border = "white")
  # })
  output$textOut <- renderText({
    paste0(
      paste0("Subtype C: ", input$subtype_c_perc, "%;  "),
      paste0("Subtype A: ", input$subtype_a_perc, "%;  "),
      paste0("Subtype D: ", input$subtype_d_perc, "%;  "),
      paste0("Subtype B: ", input$subtype_b_perc, "%; T/F: "),
      input$weight_by_subtype
    )
  })

  weight_by_subtype <- reactive({
    if (input$weight_by_subtype == "yes") {
      #browser()
      mdris <- filter_assay() %>%
        filter(subtype != "All") %>%
        arrange(subtype)
      weighted_mdri <- weighted.mean(
        x = mdris$mdriPE,
        w = c(
          input$subtype_a_perc/100,
          input$subtype_b_perc/100,
          input$subtype_c_perc/100,
          input$subtype_d_perc/100
        )
      )
      sigmas <- mdris$mdriSE
      W <- c(
        input$subtype_a_perc/100,
        input$subtype_b_perc/100,
        input$subtype_c_perc/100,
        input$subtype_d_perc/100
      )
      sigma_weighted_mdri <- sqrt(sum(sigmas^2 * W^2))
      weighted_mdri_lb <- weighted_mdri - 1.96 * sigma_weighted_mdri
      weighted_mdri_ub <- weighted_mdri + 1.96 * sigma_weighted_mdri
    } else if (input$weight_by_subtype == "no") {
      mdris <- filter_assay() %>%
        filter(subtype == "All")
      weighted_mdri <- mdris$mdriPE
      weighted_mdri_lb <- mdris$mdriLB
      weighted_mdri_ub <- mdris$mdriUB
      sigma_weighted_mdri <- mdris$mdriSE
    }
    return(c(weighted_mdri, weighted_mdri_lb, weighted_mdri_ub, sigma_weighted_mdri))
  })

  find_scale_shape <- reactive({
    #browser()
    if (input$shape_method == "shape_parameter") {
      shape_param <- input$weibull_shape
      scale_param <- find_scale(
        target_median = input$median_ttd * 365.25,
        beta = shape_param,
        max_t = 20000
      )
      delta_q25 <- NA
      delta_q50 <- abs( weibull_quantile(q = 0.5, alpha = scale_param, beta = shape_param)/365.25 - input$median_ttd)
    } else if (input$shape_method == "quantiles") {
      beta_start <- 2.5
      scale_param <- find_scale(
        target_median = input$median_ttd * 365.25,
        beta = beta_start,
        max_t = 20000
      )
      shape_param <- find_shape_quantile(
        target = input$quantile_0_25 * 365.25,
        quantile = 0.25,
        alpha = scale_param,
        max_t = 20000
      )
      delta_q25 <- abs( weibull_quantile(q = 0.25, alpha = scale_param, beta = shape_param)/365.25 - input$quantile_0_25)
      delta_q50 <- abs( weibull_quantile(q = 0.5, alpha = scale_param, beta = shape_param)/365.25 - input$median_ttd)
      while (delta_q25 > input$fitting_tolerance | delta_q50 > input$fitting_tolerance) {
        scale_param <- find_scale(
          target_median = input$median_ttd * 365.25,
          beta = shape_param,
          max_t = 20000
        )
        shape_param <- find_shape_quantile(
          target = input$quantile_0_25 * 365.25,
          quantile = 0.25,
          alpha = scale_param,
          max_t = 20000
        )
        delta_q25 <- abs( weibull_quantile(q = 0.25, alpha = scale_param, beta = shape_param)/365.25 - input$quantile_0_25)
        delta_q50 <- abs( weibull_quantile(q = 0.5, alpha = scale_param, beta = shape_param)/365.25 - input$median_ttd)
      }
    }
    return(c(scale_param, shape_param, delta_q50, delta_q25))
  })

  output$weighting_function <- renderPlot({
    params <- find_scale_shape()
    scale_param = params[1]
    shape_param = params[2]
    #browser()
    # scale <- find_scale(
    #   target_median = input$median_ttd*365.25,
    #   beta = input$weibull_shape,
    #   max_t = 20000
    # )
    plot <- plot_weibull(
      from = 0,
      to = input$median_ttd * 365.25 * 3,
      alpha = scale_param,
      beta = shape_param,
      vertline_T = as.numeric(input$bigT),
      vertline_median = input$median_ttd * 365.25
    )
    return(plot)
  })

  output$weighting_function_frr <- renderPlot({
    params <- find_scale_shape()
    scale_param = params[1]
    shape_param = params[2]
    #browser()
    # scale <- find_scale(
    #   target_median = input$median_ttd*365.25,
    #   beta = input$weibull_shape,
    #   max_t = 20000
    # )
    plot <- plot_weibull(
      from = 0,
      to = input$median_ttd * 365.25 * 3,
      alpha = scale_param,
      beta = shape_param,
      vertline_T = as.numeric(input$bigT),
      vertline_median = input$median_ttd * 365.25
    )
    return(plot)
  })


  output$quantiles <- renderText({
    params <- find_scale_shape()
    scale_param = params[1]
    shape_param = params[2]
    delta_q50 = params[3]
    delta_q25 = params[4]
    return(
      paste0(
        "25%: ",
        round(weibull_quantile(q = 0.25, alpha = scale_param, beta = shape_param, max_t = 20000) / 365.25, 2),
        "yrs; 50%: ",
        round(weibull_quantile(q = 0.5, alpha = scale_param, beta = shape_param, max_t = 20000) / 365.25, 2),
        "yrs; 75%: ",
        round(weibull_quantile(q = 0.75, alpha = scale_param, beta = shape_param, max_t = 20000) / 365.25, 2),
        "yrs. Delta (median): ",
        round(delta_q50, 4),
        "yrs. Delta (25%): ",
        round(delta_q25, 4)
      )
    )
  })


  output$quantile_0_75 <- renderText({
    params <- find_scale_shape()
    scale_param = params[1]
    shape_param = params[2]
    value <- round(weibull_quantile(q = 0.75, alpha = scale_param, beta = shape_param, max_t = 20000) / 365.25, 1)
    return(value)
  })


  weight_by_subtype_and_ttd <- reactive({
    #browser()
    params <- find_scale_shape()
    scale_param = params[1]
    shape_param = params[2]

    if (input$weight_by_subtype == "yes") {

      paramsets <- filter_assay() %>%
        filter(subtype != "All") %>%
        arrange(subtype)
      ttd_weighted_mdris <- numeric()

      for (i in 1:nrow(paramsets)) {
        params <- c(
          paramsets$beta0[i],
          paramsets$beta1[i],
          paramsets$beta2[i],
          paramsets$beta3[i]
        )
        ttd_weighted_mdri <- cubature::adaptIntegrate(
          PprimeR,
          lowerLimit = 0,
          upperLimit = as.numeric(input$bigT),
          pr_parms = params,
          alpha = scale_param,
          beta = shape_param
        )$integral
        ttd_weighted_mdris[i] <- ttd_weighted_mdri
      }
      weighted_mdri <- weighted.mean(
        x = ttd_weighted_mdris,
        w = c(
          input$subtype_a_perc / 100,
          input$subtype_b_perc / 100,
          input$subtype_c_perc / 100,
          input$subtype_d_perc / 100
        )
      )
    } else if (input$weight_by_subtype == "no") {
      mdris <- filter_assay() %>%
        filter(subtype == "All")
      params <- c(
        mdris$beta0[1],
        mdris$beta1[1],
        mdris$beta2[1],
        mdris$beta3[1]
      )
      # scale <- find_scale(
      #   target_median = input$median_ttd * 365.25,
      #   beta = input$weibull_shape,
      #   max_t = 20000
      # )
      weighted_mdri <- cubature::adaptIntegrate(
        PprimeR,
        lowerLimit = 0,
        upperLimit = as.numeric(input$bigT),
        pr_parms = params,
        alpha = scale_param,
        beta = shape_param
      )$integral
    }
    return(weighted_mdri)
  })

  weight_by_ttd_frr <- reactive({
    #browser()
    params <- find_scale_shape()
    scale_param = params[1]
    shape_param = params[2]

      frrdat <- filter_assay_frr()
      params <- c(
        frrdat$beta0[1],
        frrdat$beta1[1],
        frrdat$beta2[1],
        frrdat$beta3[1]
      )

      frr_untreated_ttd_weighted <- frr_untreated_weighted(
        pr_parms = params,
        bigT = as.numeric(input$bigT),
        alpha = scale_param,
        beta = shape_param
      )

    return(frr_untreated_ttd_weighted)
  })


  # outout$quantile_0_75 <- renderText({
  #   scale <- find_scale(
  #     target_median = input$median_ttd*365.25,
  #     beta = input$weibull_shape,
  #     max_t = 20000
  #     )
  #   shape <- input$weibull_shape
  #   quantile_0_75 <- round(weibull_quantile(q = 0.75, alpha = scale, beta = shape, max_t = 20000)/365.25, 1)
  #   return(quantile_0_75)
  # })

  output$weighted_mdri <- renderText({
    values <- weight_by_subtype()
    #c(weighted_mdri, weighted_mdri_lb, weighted_mdri_ub, sigma_weighted_mdri)
    paste0(
      round(values[1], 1),
      " (95% CI: ",
      round(values[2], 1),
      ",",
      round(values[3], 1),
      "; RSE: ",
      round(values[4]/values[1]*100, 1),
      "%)"
    )
  })

  output$screening_wp_echo <- renderText({
    if (input$screening_adjustment == "window_period") {
      return(as.numeric(input$screening_wp_wp))
    } else if (input$screening_adjustment == "screening_assay") {
      return(as.numeric(input$screening_wp_assay))
    }

  })

  output$adjusted_mdri <- renderText({
    pe = round(weight_by_subtype_and_ttd(), 1)
    #se_ci = ttd_weighted_ci()
    if (is.null(v$adjusted_mdri_ci)) {
      returnval = paste0(
        round(pe, 1),
        " (95% CI not calculated)"
      )
    } else {
      returnval = paste0(
        round(pe, 1),
        " (95% CI: ",
        round(v$adjusted_mdri_ci[1], 1),
        ",",
        round(v$adjusted_mdri_ci[2], 1),
        "; RSE: ",
        round(v$adjusted_mdri_ci[3]/pe*100, 1),
        "%)"
      )
    }
    return(returnval)
  })

  output$final_mdri <- renderText({
    adjusted_mdri_pe <- weight_by_subtype_and_ttd()
    if (input$screening_adjustment == "window_period") {
      screening_adjustment_value <- as.numeric(input$screening_wp_wp)
    } else if (input$screening_adjustment == "screening_assay") {
      screening_adjustment_value <- as.numeric(input$screening_wp_assay)
    }
    pe <- adjusted_mdri_pe - screening_adjustment_value
    if (is.null(v$fully_adjusted_mdri_ci)) {
      returnval = paste0(
        round(pe, 1),
        " (95% CI not calculated)"
      )
    } else {
      returnval = paste0(
        round(pe, 1),
        " (95% CI: ",
        round(v$fully_adjusted_mdri_ci[1], 1),
        ",",
        round(v$fully_adjusted_mdri_ci[2], 1),
        "; RSE: ",
        round(v$fully_adjusted_mdri_ci[3]/pe*100, 1),
        "%)"
      )
    }
    return(returnval)
  })

  get_bs_params <- reactive({
    runs <- filter_subtype()$run
    bs_params <- read_bs_params_dataset(runs, n = as.integer(input$n_bootstraps))
    return(bs_params)
  })

  get_bs_params_frr <- reactive({
    runs <- filter_assay_frr()$run
    bs_params <- read_bs_params_frr_dataset(runs, as.integer(input$n_bootstraps))
    return(bs_params)
  })

  #ttd_weighted_ci <- eventReactive(input$compute_ci_mdri_button, {
  observeEvent(input$compute_ci_mdri_button, {
    #browser()
    # showModal(modalDialog(
    #   title = "Confidence intervals",
    #   "Computing confidence intervals may take several minutes",
    #   easyClose = FALSE
    # ))

    weighting_params <- find_scale_shape()
    scale_param = weighting_params[1]
    shape_param = weighting_params[2]

    if (input$weight_by_subtype == "yes") {
      subtypedat <- filter_subtype()
      run_subtype <- tibble(
        subtype = subtypedat$subtype,
        run = subtypedat$run,
        w_s = case_when(
          subtype == "A1" ~ input$subtype_a_perc / 100,
          subtype == "B" ~ input$subtype_b_perc / 100,
          subtype == "C" ~ input$subtype_c_perc / 100,
          subtype == "D" ~ input$subtype_d_perc / 100
        )
      )
      paramsets <- get_bs_params() %>%
        left_join(run_subtype) %>%
        group_by(run) %>%
        mutate(bsno = 1:n()) %>%
        ungroup()

      n <- nrow(paramsets)
      withProgress(
        message = 'Weighting by time to diagnosis',
        detail = 'Calculation may take several minutes', value = 0, {
          # cluster <- parallel::makeCluster(parallel::detectCores(), outfile="")
          # doParallel::registerDoParallel(cluster)
          # subtype_ttd_weighted_mdris <- foreach::foreach(i = 1:100,
          #                                                .combine = c,
          #                                                .inorder = FALSE,
          #                                                .export = c("session",
          #                                                            "find_scale")
          # ) %dopar%
          #   {
          #     # shiny::incProgress(
          #     #   amount = 1/n,
          #     #   session = session
          #     #   )
          #     print(i)
          #     params <- c(
          #       paramsets$beta0[i],
          #       paramsets$beta1[i],
          #       paramsets$beta2[i],
          #       paramsets$beta3[i]
          #     )
          #     scale <- find_scale(
          #       target_median = input$median_ttd * 365.25,
          #       beta = input$weibull_shape,
          #       max_t = 20000
          #     )
          #     ttd_weighted_mdri <- cubature::adaptIntegrate(
          #       PprimeR,
          #       lowerLimit = 0,
          #       upperLimit = as.numeric(input$bigT),
          #       pr_parms = params,
          #       alpha = scale,
          #       beta = input$weibull_shape
          #     )$integral
          #   return(ttd_weighted_mdri)
          #   }

          #!!! Try parallel again
          # ttd_weighted_mdris <- foreach::foreach(
          #   i = 1:100,
          #   .combine = c,
          #   .inorder = FALSE,
          #   .export = c("session", "find_scale")
          # ) %dopar% {
          #   # shiny::incProgress(
          #   #       amount = 1/n,
          #   #       session = session
          #   #       )
          #   params <- c(
          #     paramsets$beta0[i],
          #     paramsets$beta1[i],
          #     paramsets$beta2[i],
          #     paramsets$beta3[i]
          #   )
          #   scale <- find_scale(
          #     target_median = input$median_ttd * 365.25,
          #     beta = input$weibull_shape,
          #     max_t = 20000
          #   )
          #   ttd_weighted_mdri <- cubature::adaptIntegrate(
          #     PprimeR,
          #     lowerLimit = 0,
          #     upperLimit = as.numeric(input$bigT),
          #     pr_parms = params,
          #     alpha = scale_param,
          #     beta = shape_param
          #   )$integral
          #   return(ttd_weighted_mdri)
          # }




          # #!!! Try using future
          # browser()
          # progress = AsyncProgress$new(message="Complex analysis")
          # f <- future({
          #   ttd_weighted_mdris <- numeric()
          #   for (i in 1:n) {
          #     progress$inc(1/n)
          #     incProgress(1/n)
          #       params <- c(
          #         paramsets$beta0[i],
          #         paramsets$beta1[i],
          #         paramsets$beta2[i],
          #         paramsets$beta3[i]
          #       )
          #       scale <- find_scale(
          #         target_median = input$median_ttd * 365.25,
          #         beta = input$weibull_shape,
          #         max_t = 20000
          #       )
          #       ttd_weighted_mdri <- cubature::adaptIntegrate(
          #         PprimeR,
          #         lowerLimit = 0,
          #         upperLimit = as.numeric(input$bigT),
          #         pr_parms = params,
          #         alpha = scale_param,
          #         beta = shape_param
          #       )$integral
          #       ttd_weighted_mdris[i] <- ttd_weighted_mdri
          #   }
          #   progress$close()
          #   return(ttd_weighted_mdris)
          # })


          # for (i in 1:n) {
          #   incProgress(1/n)
          #   params <- c(
          #     paramsets$beta0[i],
          #     paramsets$beta1[i],
          #     paramsets$beta2[i],
          #     paramsets$beta3[i]
          #   )
          #   scale <- find_scale(
          #     target_median = input$median_ttd * 365.25,
          #     beta = input$weibull_shape,
          #     max_t = 20000
          #   )
          #   ttd_weighted_mdri <- cubature::adaptIntegrate(
          #     PprimeR,
          #     lowerLimit = 0,
          #     upperLimit = as.numeric(input$bigT),
          #     pr_parms = params,
          #     alpha = scale_param,
          #     beta = shape_param
          #   )$integral
          #   ttd_weighted_mdris[i] <- ttd_weighted_mdri
          # }





          #!! Compute sequentially !!
          ttd_weighted_mdris <- numeric()
          for (i in 1:n) {
            incProgress(1/n)
            params <- c(
              paramsets$beta0[i],
              paramsets$beta1[i],
              paramsets$beta2[i],
              paramsets$beta3[i]
            )
            scale <- find_scale(
              target_median = input$median_ttd * 365.25,
              beta = input$weibull_shape,
              max_t = 20000
            )
            ttd_weighted_mdri <- cubature::adaptIntegrate(
              PprimeR,
              lowerLimit = 0,
              upperLimit = as.numeric(input$bigT),
              pr_parms = params,
              alpha = scale_param,
              beta = shape_param
            )$integral
            ttd_weighted_mdris[i] <- ttd_weighted_mdri
          }
        })
      paramsets$ttd_weighted_mdri <- ttd_weighted_mdris
      subtype_ttd_weighted_mdris = paramsets %>%
        group_by(bsno) %>%
        summarise(
          subtype_ttd_weighted_mdri = weighted.mean(
            x = ttd_weighted_mdri,
            w = w_s
          )
        ) %>%
        pull(subtype_ttd_weighted_mdri)

      weighted_mdri_lb <- quantile(subtype_ttd_weighted_mdris, 0.025)
      weighted_mdri_ub <- quantile(subtype_ttd_weighted_mdris, 0.975)
      sigma_weighted_mdri <- sd(subtype_ttd_weighted_mdris)
      v$adjusted_mdri_ci <- c(weighted_mdri_lb, weighted_mdri_ub, sigma_weighted_mdri)

      # HERE HERE
      if (input$screening_adjustment == "window_period") {
        screening_adjustment_value <- as.numeric(input$screening_wp_wp)
      } else if (input$screening_adjustment == "screening_assay") {
        screening_adjustment_value <- as.numeric(input$screening_wp_assay)
      }
      subtype_ttd_weighted_mdris_scrn_adj = subtype_ttd_weighted_mdris - screening_adjustment_value
      weighted_mdri_lb_scrn_adj <- quantile(subtype_ttd_weighted_mdris_scrn_adj, 0.025)
      weighted_mdri_ub_scrn_adj <- quantile(subtype_ttd_weighted_mdris_scrn_adj, 0.975)
      sigma_weighted_mdri_scrn_adj <- sd(subtype_ttd_weighted_mdris_scrn_adj)
      v$fully_adjusted_mdri_ci <- c(weighted_mdri_lb_scrn_adj, weighted_mdri_ub_scrn_adj, sigma_weighted_mdri_scrn_adj)
      #return(c(sigma_weighted_mdri, weighted_mdri_lb, weighted_mdri_ub))
    } else if (input$weight_by_subtype == "no") {

    }

  })


  output$runs <- renderText({
    runs <- filter_subtype()$run
    return(runs)
  })

  # output$ci_mdris <- renderText({
  #   ci <- ttd_weighted_ci()
  #   paste0("95% CI:" )
  # })

  # output$bstab <- renderTable({
  #   return(get_bs_params())
  # })

  output$threshold_output <- renderText({
    return(get_threshold())
  })

  output$active_tab <- renderText({
    return(input$output_tab)
  })


  output$frr_untreated <- reactive({
    #browser()
    pe = weight_by_ttd_frr()
    if (is.null(v$frr_untreated_ci)) {
      returnval = paste0(
        round(pe*100, 2),
        "% (95% CI not calculated)"
      )
    } else {
      returnval = paste0(
        round(pe*100, 2),
        "% (95% CI: ",
        round(v$frr_untreated_ci[1]*100, 2),
        "%,",
        round(v$frr_untreated_ci[2]*100, 2),
        "%; RSE: ",
        round(v$frr_untreated_ci[3]/pe*100, 2),
        "%)"
      )
    }
    return(returnval)
  })

  output$frr_final <- reactive({
    #browser()
    c = as.numeric(input$rx_coverage)/100
    pe = (1 - c) * weight_by_ttd_frr()
    if (is.null(v$frr_weighted_ci)) {
      returnval = paste0(
        round(pe*100, 2),
        "% (95% CI not calculated)"
      )
    } else {
      returnval = paste0(
        round(pe*100, 2),
        "% (95% CI: ",
        round(v$frr_weighted_ci[1]*100, 2),
        "%,",
        round(v$frr_weighted_ci[2]*100, 2),
        "%; RSE: ",
        round(v$frr_weighted_ci[3]/pe*100, 2),
        "%)"
      )
    }
    return(returnval)
  })

  observeEvent(input$compute_ci_frr_button, {
    #browser()
    weighting_params <- find_scale_shape()
    scale_param = weighting_params[1]
    shape_param = weighting_params[2]

    paramsets <- get_bs_params_frr()

    # subtypedat <- filter_subtype()
    # run_subtype <- tibble(
    #   subtype = subtypedat$subtype,
    #   run = subtypedat$run,
    #   w_s = case_when(
    #     subtype == "A1" ~ input$subtype_a_perc / 100,
    #     subtype == "B" ~ input$subtype_b_perc / 100,
    #     subtype == "C" ~ input$subtype_c_perc / 100,
    #     subtype == "D" ~ input$subtype_d_perc / 100
    #   )
    # )
    # paramsets <- get_bs_params() %>%
    #   left_join(run_subtype) %>%
    #   group_by(run) %>%
    #   mutate(bsno = 1:n()) %>%
    #   ungroup()

    n <- nrow(paramsets)
    withProgress(
      message = 'Weighting by time to diagnosis',
      detail = 'Calculation may take several minutes', value = 0, {
        #!! Compute sequentially !!
        ttd_weighted_untreated_frrs <- numeric()
        for (i in 1:n) {
          incProgress(1/n)
          params <- c(
            paramsets$beta0[i],
            paramsets$beta1[i],
            paramsets$beta2[i],
            paramsets$beta3[i]
          )
          frr_untreated_ttd_weighted <- frr_untreated_weighted(
            pr_parms = params,
            bigT = as.numeric(input$bigT),
            alpha = scale_param,
            beta = shape_param
          )
          ttd_weighted_untreated_frrs[i] <- frr_untreated_ttd_weighted
        }
      })

    frr_untreated_lb <- quantile(ttd_weighted_untreated_frrs, 0.025)
    frr_untreated_ub <- quantile(ttd_weighted_untreated_frrs, 0.975)
    sigma_untreated_frr <- sd(ttd_weighted_untreated_frrs)
    v$frr_untreated_ci <- c(frr_untreated_lb, frr_untreated_ub, sigma_untreated_frr)

    c = as.numeric(input$rx_coverage)/100
    c_weighted_untreated_frrs <- (1-c) * ttd_weighted_untreated_frrs
    frr_overall_lb <- quantile(c_weighted_untreated_frrs, 0.025)
    frr_overall_ub <- quantile(c_weighted_untreated_frrs, 0.975)
    sigma_frr_overall <- sd(c_weighted_untreated_frrs)
    v$frr_weighted_ci <- c(frr_overall_lb, frr_overall_ub, sigma_frr_overall)

  })

}



ui <- fluidPage(
  fluidRow(
    titlePanel(
      "Context-specific recency test property estimation tool",
      windowTitle = "MDRI & FRR estimation"
    )
  ),
  fluidRow(
    tabsetPanel(
      type = "tabs",
      tabPanel(
        "Estimate MDRI and FRR",
    column(
      12,
      fluidRow(
        column(
          3,
          wellPanel(
            h4("Recent infection testing algorithm"),
            selectInput(
              "bigT",
              "Time cutoff (T)",
              choices = c(
                "1 year" = 365.25,
                "18 months" = 548,
                "2 years" = 730.5
              ),
              selected = 730.5
            ),
            selectInput(
              "assay",
              "Primary assay",
              choices = c(
                "Sedia Limiting Antigen Avidity EIA" = "LAg-Sedia",
                "Maxim Limiting Antigen Avidity EIA" = "LAg-Maxim",
                "BioRad Avidity (CDC protocol)" = "BioRadAvidity-CDC",
                "Sedia BED" = "BED",
                "Sedia Asanté (electronic reader)" = "Asante Electronic"
              ),
              selected = "LAg-Sedia", # c("Sedia Limiting Antigen Avidity EIA" = "LAg-Sedia"),
              multiple = FALSE
            ),
            conditionalPanel(
              condition = "input.assay == 'LAg-Maxim'",
              sliderInput(
                "assay_threshold_lagm",
                "Assay threshold",
                min = 1.0,
                max = 2.5,
                value = 1.5,
                step = 0.25
              )
              # selectInput(
              #   "assay_threshold",
              #   "Assay threshold",
              #   choices = c(
              #     "≤1.0" = 1,
              #     "≤1.25" = 1.25,
              #     "≤1.5" = 1.5,
              #     "≤1.75" = 1.75,
              #     "≤2.0" = 2,
              #     "≤2.25" = 2.25,
              #     "≤2.5" = 2.5,
              #     "Other" = 9999
              #   ),
              #   selected = 1.5
              # )
            ),
            conditionalPanel(
              condition = "input.assay == 'LAg-Sedia'",
              sliderInput(
                "assay_threshold_lags",
                "Assay threshold",
                min = 1.0,
                max = 2.5,
                value = 1.5,
                step = 0.25
              )
              # selectInput(
              #   "assay_threshold",
              #   "Assay threshold",
              #   choices = c(
              #     "≤1.0" = 1,
              #     "≤1.25" = 1.25,
              #     "≤1.5" = 1.5,
              #     "≤1.75" = 1.75,
              #     "≤2.0" = 2,
              #     "≤2.25" = 2.25,
              #     "≤2.5" = 2.5,
              #     "Other" = 9999
              #   ),
              #   selected = 1.5
              # )
            ),
            conditionalPanel(
              condition = "input.assay == 'BioRadAvidity-CDC'",
              sliderInput(
                "assay_threshold_bra",
                "Assay threshold",
                min = 10.0,
                max = 80.0,
                value = 40.0,
                step = 10.0
              )
              # selectInput(
              #   "assay_threshold",
              #   "Assay threshold",
              #   choices = c(
              #     "≤10" = 10,
              #     "≤20" = 20,
              #     "≤30" = 30,
              #     "≤40" = 40,
              #     "≤50" = 50,
              #     "≤60" = 60,
              #     "≤70" = 70,
              #     "≤80" = 80,
              #     "Other" = 9999
              #   ),
              #   selected = 40
              # )
            ),
            conditionalPanel(
              condition = "input.assay == 'BED'",
              sliderInput(
                "assay_threshold_bed",
                "Assay threshold",
                min = 0.4,
                max = 1.8,
                value = 0.8,
                step = 0.2
              )
              # selectInput(
              #   "assay_threshold",
              #   "Assay threshold",
              #   choices = c(
              #     "≤0.4" = 0.4,
              #     "≤0.6" = 0.6,
              #     "≤0.8" = 0.8,
              #     "≤1.0" = 1.0,
              #     "≤1.2" = 1.2,
              #     "≤1.4" = 1.4,
              #     "≤1.6" = 1.6,
              #     "≤1.8" = 1.8,
              #     "Other" = 9999
              #   ),
              #   selected = 0.8
              # )
            ),
            conditionalPanel(
              condition = "input.assay == 'Asante Electronic'",
              sliderInput(
                "assay_threshold_asante",
                "Assay threshold",
                min = 2.00,
                max = 4.00,
                value = 3.00,
                step = 0.25
              )
              # selectInput(
              #   "assay_threshold",
              #   "Assay threshold",
              #   choices = c(
              #     "≤2.00" = 2.00,
              #     "≤2.25" = 2.25,
              #     "≤2.50" = 2.50,
              #     "≤2.75" = 2.75,
              #     "≤3.00" = 3.00,
              #     "≤3.25" = 3.25,
              #     "≤3.50" = 3.50,
              #     "≤3.75" = 3.75,
              #     "≤4.00" = 4.00,
              #     "Other" = 9999
              #   ),
              #   selected = 3.0
              # )
            ),
            selectInput(
              "vl_threshold",
              "Viral load threshold",
              choices = c(
                "no VL threshold" = 0,
                ">75 c/mL" = 75,
                ">200 c/mL" = 200,
                ">400 c/mL" = 400,
                ">1000 c/mL" = 1000
              ),
              selected = 75
            ),
            radioButtons(
              "diagnosed_screened_out",
              "Does the RITA include prior HIV diagnosis?",
              choices = c(
                "Yes" = "yes",
                "No" = "no"
              ),
              selected = "yes"
            ),
            radioButtons(
              "arv_detection",
              "Does the RITA include ARV testing?",
              choices = c(
                "Yes" = "yes",
                "No" = "no"
              ),
              selected = "yes"
            ),
            radioButtons(
              "screening_adjustment",
              "Adjust for screening sensitivity",
              choices = c(
                "Select screening assay" = "screening_assay",
                "Provide window period" = "window_period"
              ),
              selected = "screening_assay"
            ),
            conditionalPanel(
              condition = "input.screening_adjustment == 'screening_assay'",
              selectInput(
                "screening_wp_assay",
                "Least sensitive assay required to define positivity",
                choices = c(
                  "Ab rapid test" = 22,
                  "Ab/Ag rapid test" = 13,
                  "Third gen lab-based serology" = 17,
                  "Fourth gen lab-based serology" = 16,
                  "Geenius" = 29,
                  "Western blot" = 30,
                  "NAT 20-30c/mL" = 4,
                  "NAT 100c/mL" = 6
                ),
                selected = 22
              )
            ),
            conditionalPanel(
              condition = "input.screening_adjustment == 'window_period'",
              sliderInput(
                "screening_wp_wp",
                label = "Window period relative to 1c/mL NAT",
                min = 1,
                max = 35,
                value = 22,
                step = 1
              )
            ),
            div(
              actionButton("help_rita", "Help"),
              align = "center"
            )


            # ,
            # radioButtons(
            #   "second_assay_yn",
            #   "Does the RITA include a second recency assay?",
            #   choices = c(
            #     "No" = "no",
            #     "Yes" = "yes"
            #              ),
            #   selected = "no"
            #   ),
            # conditionalPanel(
            #   condition = "input.second_assay_yn == 'yes'",
            #   selectInput(
            #     "secondary_assay",
            #     "Second assay",
            #     choices = c(
            #       "Viral load" = "viral_load"
            #     ),
            #     selected = "viral_load"
            #   )
            # )
          ),
          wellPanel(
            h4("Manage computation"),
            selectInput(
              "n_bootstraps",
              "Bootstrap iterations",
              choices = c(
                "100" = 100,
                "1000" = 1000,
                "5000" = 5000,
                "10000" = 10000
              ),
              selected = 1000
            ),
            em("Select a small number for testing and at least 5,000 for final estimates.")
          )
        ),
        column(
          3,
          conditionalPanel(
            condition = "input.output_tab == 'Weighted MDRI' | input.output_tab == 'MDRI data'",
            wellPanel(
              h4("Subtype distribution"),
              # This creates a character vector
              # selectInput(
              #   "subtypes",
              #   "Subtypes present in population",
              #   choices = c(
              #     "Subtype C" = "C",z
              #     "Subtype A" = "A1",
              #     "Subtype D" = "D",
              #     "Subtype B" = "B",
              #     "Other" = "Other"
              #   ),
              #   selected = NULL,
              #   multiple = TRUE
              # ),
              radioButtons(
                "weight_by_subtype",
                "Adjust for subtype distribution",
                choices = c(
                  "Yes" = "yes",
                  "No" = "no"
                ),
                selected = "yes"
              ),
              conditionalPanel(
                condition = "input.weight_by_subtype == 'yes'",
                sliderInput(
                  "subtype_c_perc",
                  "Percentage subtype C",
                  min = 0,
                  max = 100, #maxC,
                  value = 5,
                  step = 1
                ),
                sliderInput(
                  "subtype_a_perc",
                  "Percentage subtype A",
                  min = 0,
                  max = 100, #maxA,
                  value = 50,
                  step = 1
                ),
                sliderInput(
                  "subtype_d_perc",
                  "Percentage subtype D",
                  min = 0,
                  max = 100, #maxD,
                  value = 45,
                  step = 1
                ),
                sliderInput(
                  "subtype_b_perc",
                  "Percentage subtype B",
                  min = 0,
                  max = 100, #maxB,
                  value = 0,
                  step = 1
                )
              ),
              div(
                actionButton("help_subtype", "Help"),
                align = "center"
              )
              # radioButtons(
              #   "weight_by_sex",
              #   "Adjust for sex distribution",
              #   choices = c(
              #     "Yes" = "yes",
              #     "No" = "no"
              #   ),
              #   selected = "no"
              # ),
              # conditionalPanel(
              #   condition = "input.weight_by_sex == 'yes'",
              #   sliderInput(
              #     "female_perc",
              #     "Percentage female",
              #     min = 0,
              #     max = 100,
              #     value = 50,
              #     step = 1
              #   ),
              #   em("Sex adjustment is not functional at present", style = "color:red")
              # )
            )
          ),
          conditionalPanel(
            condition = "input.output_tab == 'Weighted FRR'",
            wellPanel(
              h4("FRR parameters"),
              sliderInput(
                "rx_coverage",
                "Treatment coverage (%)",
                min = 0,
                max = 100,
                value = 80,
                step = 1
              ),
              div(
                actionButton("help_frr", "Help"),
                align = "center"
              )
            )
          ),
          wellPanel(
            h4("Adjustment for early diagnosis/treatment"),
            conditionalPanel(
              condition = "input.vl_threshold == 0 & input.diagnosed_screened_out == 'no' & input.arv_detection == 'no'",
              em("No adjustment for early diagnosis or treatment is necessary. RITAs without viral load are not recommended, since very high false recent ratios would be expected in individuals on antiretroviral treatment.", style = "color:red")
            ),
            conditionalPanel(
              condition = "(input.vl_threshold > 0 | input.diagnosed_screened_out == 'yes' | input.arv_detection == 'yes') & input.adjust_ttd == 'no'",
              em("Adjustment for early diagnosis or treatment is recommended.", style = "color:red")
            ),
            conditionalPanel(
              condition = "input.arv_detection == 'yes' & input.diagnosed_screened_out == 'no'",
              em("For the specified RITA, the distribution of times from infection to treatment initiation is most relevant.")
            ),
            conditionalPanel(
              condition = "input.vl_threshold > 0 & input.arv_detection == 'no' & input.diagnosed_screened_out == 'no'",
              em("For the specified RITA, the distribution of times from infection to viral suppression is most relevant.")
            ),
            conditionalPanel(
              condition = "input.diagnosed_screened_out == 'yes'",
              em("For the specified RITA, the distribution of times from infection to diagnosis is most relevant.")
            ),
            radioButtons(
              "adjust_ttd",
              "Adjust for early diagnosis or treatment",
              choices = c(
                "Yes" = "yes",
                "No" = "no"
              ),
              selected = "yes"
            ),
            conditionalPanel(
              condition = "input.adjust_ttd == 'yes'",
              sliderInput(
                "median_ttd",
                "Best estimate of median time to diagnosis/treatment (years)",
                min = 0.25,
                max = 10,
                value = 2,
                step = 0.05
              ),
              radioButtons(
                "shape_method",
                "Specify shape by",
                choices = c(
                  "Weibull shape parameter" = "shape_parameter",
                  "Quantiles" = "quantiles"
                ),
                selected = "shape_parameter"
              ),
              conditionalPanel(
                condition = "input.shape_method == 'shape_parameter'",
                sliderInput(
                  "weibull_shape",
                  "Weibull shape parameter",
                  min = 1,
                  max = 10,
                  value = 2.5,
                  step = 0.25
                )
              ),
              conditionalPanel(
                condition = "input.shape_method == 'quantiles'",
                sliderInput(
                  "quantile_0_25",
                  "Time when 25% of individuals are diagnosed (years)",
                  min = 0.1,
                  max = 2.0,
                  value = 1,
                  step = 0.05
                ),
                numericInput(
                  "fitting_tolerance",
                  "Fitting tolerance (advanced)",
                  value = 0.025,
                  min = 0.01,
                  max = 0.1,
                  step = 0.005
                ),
                tags$b("Time when 75% of individuals are diagnosed (years)"),
                textOutput("quantile_0_75"),
                # sliderInput(
                #   "quantile_0_75",
                #   "75th percentile for time to diagnosis",
                #   min = 2.1,
                #   max = 1.9,
                #   value = 1,
                #   step = 0.1
                # ),
                #em("Quantile-based estimation not active", style = "color:red")
              )
            ),
            div(
              actionButton("help_ttd", "Help"),
              align = "center"
            )
          )


          # ** Quantiles options:
          # .25 and 0.75 (and show 0.5)
          # pick 0.25 default midpoint between 0 and 0.5
          # constrain 0.75 to possible values that a parameter could be found for

          # **!! Add more instructional text in each well relating to the inputs
          #

          # ,
          # wellPanel(
          #   h4("Placeholder!"),
          #   sliderInput("obs", "Number of observations:", min = 10, max = 500, value = 100)
          # )
        ),
        column(
          6,
          tabsetPanel(
            id = "output_tab",
            type = "tabs",
            tabPanel(
              "Weighted MDRI",
              wellPanel(
                tags$b("Subtype-weighted MDRI:"),
                textOutput("weighted_mdri"),
                p(), p(),
                tags$b("Time to diagnosis-adjusted MDRI:"),
                textOutput("adjusted_mdri"),
                p(), p(),
                tags$b("Screening assay adjustment:"),
                textOutput("screening_wp_echo"),
                p(), p(),
                tags$b("Fully adjusted MDRI:"),
                textOutput("final_mdri"),
                div(
                  actionButton("compute_ci_mdri_button", "(Re)calculate 95% CI"),
                  align = "center"
                ),
                p(), p(), p(),
                tags$b("Weighting function (survival in undiagnosed/untreated state):"),
                p(),
                plotOutput("weighting_function"),
                # div(
                #   em("Note: Only 1,000 bootstrap iterations used to reduce compute times before production deployment.", style = "color:red"),
                #   align = "center"
                # )
              )
            ),
            tabPanel(
              "MDRI data",
              wellPanel(
                h4("Subtype-specific MDRIs for the selected RITA:"),
                tableOutput("mdri_table_assay"),
                p(),
                h4("Subtype-specific MDRIs to be used:"),
                tableOutput("mdri_table_subtype")
              )
            ),
            tabPanel(
              "Weighted FRR",
              wellPanel(
                # Put FRR results here
                tags$b("FRR in treated invidivuals:"),
                p("0.0%"),
                tags$b("TTD-weighted FRR in untreated invidivuals (infected for > T):"),
                textOutput("frr_untreated"),
                p(), p(),
                tags$b("TTD and Rx weighted FRR:"),
                textOutput("frr_final"),
                div(
                  actionButton("compute_ci_frr_button", "(Re)calculate 95% CI"),
                  align = "center"
                ),

                p(), p(), p(),
                tags$b("Weighting function (survival in undiagnosed/untreated state):"),
                p(),
                plotOutput("weighting_function_frr")
              )
            ) #,
            # tabPanel(
            #   "TTD weighting function",
            #   wellPanel(
            #     h4("Weighting function (survival in undiagnosed state)"),
            #     plotOutput("weighting_function"),
            #     tags$b("Quantiles:"), textOutput("quantiles")
            #   )
            # ),
            # tabPanel(
            #   "Debugging output",
            #   wellPanel(
            #     h4("Subtype weighting setup:"),
            #     textOutput("textOut"),
            #     h4("Assay threshold:"),
            #     textOutput("threshold_output"),
            #     h4("Runs:"),
            #     textOutput("runs"),
            #     h4("Output tab"),
            #     textOutput("active_tab")
            #     #tableOutput("bstab")
            #   )
            # )
          )
        )
      )#,
    )

      ),
    tabPanel(
      "Documentation",
      wellPanel(
        p("This tool uses public datasets from the Consortium for the Evaluation
                  and Performance of HIV Incidence Assays (CEPHIA) to compute
                  context-specific MDRI and FRR estimates, for specified recent infection testing
                  algorithms (RITAs), adjusted for subtype distribution and early
                  diagnosis and treatment."),
        p("Click the help button in each input parameter section for further information.")
      )
    )
    )
    # Close tabsetpanel here
  ),
  fluidRow(
    wellPanel(
      div(
        em("Copyright © UNAIDS and World Health Organisation 2022. Code available under the GPLv3. App version 0.98."),
        align = "center"
      ),
      p(),
      div(
        img(src = "https://www.unaids.org/sites/default/files/media/images/UNAIDS-SDG-2019_en.png", width="200px"),
        align = "center"
      ),
      br(),
      div(
        img(src = "https://www.who.int/images/default-source/fallback/header-logos/h-logo-blue.svg", width="200px"),
        align = "center"
      )

    )
  )
  #  ),
  #   tabPanel(
  #     "Compute new MDRI fits",
  #     fluidRow(
  #       wellPanel(
  #         p("No live calculation: submit a job and receive results by email")
  #       )
  #     ),
  #     fluidRow(
  #       wellPanel(
  #         p("Here UI for computing new fits")
  #       )
  #
  #     )
  #
  #   ),
  #   tabPanel(
  #     "Help",
  #     wellPanel(
  #       p("Help page goes here")
  #     )
  #
  #   )
  # )
  #)
)

shinyApp(ui = ui, server = server)

# Run with:
# runApp('~/dev/recency_apps/mdri_estimation/app')
