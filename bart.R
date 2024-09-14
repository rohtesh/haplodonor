
library(shiny)
options(mc.cores = 8)
library(BART3)

## load the BART fit
post <- readRDS('sparse20.rds')
class. <- class(post$tx.train)[1]
if(class. == 'matrix') names. <- dimnames(post$tx.train)[[2]]
else if(class. == 'data.frame') names. <- names(post$tx.train)

# Fit the Cox model with the updated factor levels and labels
## cox_model_standard <- coxph(Surv(tpfs, pfs) ~ disease_stage_rsm + 
##                               age_category + kps + recipientcmv + hctcigp + 
##                               dqb1_gvh + drb1_gvh + yeartx + dnrage + 
##                               b_leader + tce, data = data_clean)

b_leader_labels <- c('Match' = 0, 'Mismatch' = 1)
drb1_gvh_labels <- c('Match' = 0, 'Missing' = 0.5, 'Mismatch' = 1)
dqb1_gvh_labels <- c('Match' = 0, 'Missing' = 0.5, 'Mismatch' = 1)
disease_stage_rsm_labels <- c('CR1/Early' = 1, 'CR2+' = 2, 'Advanced/Active' = 4)
hctcigp_labels <- c('None' = 0, '1' = 1, '2' = 2, '3+' = 3)
kps_labels <- c('10:80' = 0, 'Missing' = 0.5, '90:100' = 1)
tce_labels <- c('No non-permissive mismatch' = 0, 'Missing' = 0.5, 
                'Non-permissive mismatch' = 1)
recipientcmv_labels <- c('Negative' = 0, 'Missing' = 0.5, 'Positive' = 1)

data_clean <- data.frame(age_category = 1:80)
data_clean$age_category <- factor(data_clean$age_category)

# --- Shiny App UI ---
ui <- fluidPage(
  titlePanel("PFS Probability"),
  sidebarLayout(
    sidebarPanel(
      selectInput("b_leader", "B Leader:", 
                  choices = setNames(names(b_leader_labels), b_leader_labels)),
      selectInput("drb1_gvh", "DRB1 GVH:", 
                  choices = setNames(names(drb1_gvh_labels), drb1_gvh_labels)),
      selectInput("dqb1_gvh", "DQB1 GVH:", 
                  choices = setNames(names(dqb1_gvh_labels), dqb1_gvh_labels)), 
      numericInput("dnrage", "Donor Age:", value = 30, min = 0),
      selectInput("age_category", "Age Category:", 
                  choices = levels(data_clean$age_category)),
      selectInput("disease_stage_rsm", "Disease Stage:", 
                  choices = setNames(names(disease_stage_rsm_labels), disease_stage_rsm_labels)),
      selectInput("hctcigp", "HCT-CI:", 
                  choices = setNames(names(hctcigp_labels), hctcigp_labels)),
      selectInput("kps", "KPS:", 
                  choices = setNames(names(kps_labels), kps_labels)),
      selectInput("recipientcmv", "Recipient CMV:", 
                  choices = setNames(names(recipientcmv_labels), recipientcmv_labels)),
      selectInput("tce", "DPB1:", 
                  choices = setNames(names(tce_labels), tce_labels)),
      actionButton("submit", "Calculate")
    ),
    mainPanel(
      plotOutput("survival_curve"),
      h4("Predicted PFS Probabilities:"),
      verbatimTextOutput("prob_12_months"),
      verbatimTextOutput("prob_24_months"),
      verbatimTextOutput("prob_36_months"),
      
      # Disclaimer
      p(strong("Disclaimer:"),
        "This application and the results provided are from preliminary analysis and for model building only. They are not intended for clinical decision-making or any other purpose. Please consult with a qualified healthcare professional for any health-related concerns.")
    )
  )
)

# --- Shiny App Server ---
server <- function(input, output) {
  output$survival_curve <- renderPlot({
    ## Create new_data inside reactive expression
    new_data <- data.frame(
      times = 0, ## place holder
      disease =  input$disease,
      sex = input$sex, ## 1:M, 2:F
      dnrage = input$dnrage,
      graftype = input$graftype,
      yeartx = 2016, ## 3-year follow-up is pre-pandemic 
      age = input$age_category, 
      numhires10 = input$numhires10,
      ethgp = input$ethgp,
      kps = input$kps, 
      hctcigp = input$hctcigp, 
      condint = input$condint,
      a_gvh = input$a_gvh,
      b_gvh = input$b_gvh,
      c_gvh = input$c_gvh,
      dqb1_gvh = input$dqb1_gvh, 
      drb1_gvh = input$drb1_gvh, 
      dx2tx_months = 7,  ## Set to the median
      b_leader = input$b_leader, 
      tce = input$tce, 
      dcmvpr = input$donorcmv, 
      rcmvpr = input$recipientcmv, 
      sibling = input$sibling,
      parent = input$parent,
      child = input$child,
      relative_miss = 0,
      stage = input$disease_stage_rsm, 
      tce_miss = input$tce_miss, 
      dsex = input$dsex 
    )
    
    # Input validation (optional but recommended):
    validate(
      need(nrow(new_data) > 0, "Please provide all input values."),
      need(!any(is.na(new_data)), "Input values cannot be missing.")
      need(ncol(new_data) == length(names.), "Number of variables mismatch.")
      need(all(names(new_data) == names.), "Variable order mismatch.")
    )
    
    ##Calculate survival probabilities
    new_data. <- NULL
    for(i in 1:post$K) new_data. <- rbind(new_data., new_data)
    new_data.$times <- post$times
    surv_obj <- predict(post, new_data.)

    ##surv_obj <- survfit(cox_model_standard, newdata = new_data)

    ##Generate survival curve
    plot(c(0, post$times), c(1, surv_obj$surv.test.mean), type = 's',
         ylim = 0:1, xlim = 0:60, xlab = "Time (Months)",
         ylab = "Survival Probability: S(t|x)")

    ## ggsurvplot(surv_obj,
    ##            data = data_clean, 
    ##            conf.int = TRUE,
    ##            xlab = "Time (Months)",
    ##            ylab = "Survival Probability",
    ##            title = "Predicted Survival Curve"
    ## )
  })
  
  # Calculate and display probabilities at specific time points
  observe({
    ## Create new_data inside reactive expression
    new_data <- data.frame(
      times = 0, ## place holder
      disease =  input$disease,
      sex = input$sex, ## 1:M, 2:F
      dnrage = input$dnrage,
      graftype = input$graftype,
      yeartx = 2016, ## 3-year follow-up is pre-pandemic 
      age = input$age_category, 
      numhires10 = input$numhires10,
      ethgp = input$ethgp,
      kps = input$kps, 
      hctcigp = input$hctcigp, 
      condint = input$condint,
      a_gvh = input$a_gvh,
      b_gvh = input$b_gvh,
      c_gvh = input$c_gvh,
      dqb1_gvh = input$dqb1_gvh, 
      drb1_gvh = input$drb1_gvh, 
      dx2tx_months = 7,  ## Set to the median
      b_leader = input$b_leader, 
      tce = input$tce, 
      dcmvpr = input$donorcmv, 
      rcmvpr = input$recipientcmv, 
      sibling = input$sibling,
      parent = input$parent,
      child = input$child,
      relative_miss = 0,
      stage = input$disease_stage_rsm, 
      tce_miss = input$tce_miss, 
      dsex = input$dsex 
    )    
    
    # Predicting survival probabilities
    times <- c(12, 24, 36) ## corresponds to year() 
    year <- c('1' = 10, '2' = 13, '3' = 17)

    new_data. <- NULL
    for(i in 1:post$K) new_data. <- rbind(new_data., new_data)
    new_data.$times <- post$times
    surv_obj <- predict(post, new_data.)

    ##surv_summary <- summary(survfit(cox_model_standard, newdata = new_data), times = times)

    output$prob_12_months <- renderPrint({
      cat("12 Months:", round(surv_obj$surv.test.mean[year[1]], 4),
          " (95% CI:", round(surv_obj$surv.test.lower[year[1]], 4), ":", 
          round(surv_obj$surv.test.upper[year[1]], 4), ")")
    })
    output$prob_24_months <- renderPrint({
      cat("24 Months:", round(surv_obj$surv.test.mean[year[2]], 4),
          " (95% CI:", round(surv_obj$surv.test.lower[year[2]], 4), ":", 
          round(surv_obj$surv.test.upper[year[2]], 4), ")")
    })
    output$prob_36_months <- renderPrint({
      cat("36 Months:", round(surv_obj$surv.test.mean[year[3]], 4),
          " (95% CI:", round(surv_obj$surv.test.lower[year[3]], 4), ":", 
          round(surv_obj$surv.test.upper[year[3]], 4), ")")
    })
  })
}

# --- Run the App ---
shinyApp(ui = ui, server = server)
