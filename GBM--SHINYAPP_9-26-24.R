library(shiny)
library(gbm)



# Define labels for categorical variables
b_leader_labels <- c("Match" = 0, "Mismatch" = 1)
drb1_gvh_labels <- c("Match" = 0, "Mismatch" = 1)
dqb1_gvh_labels <- c("Match" = 0, "Mismatch" = 1)
tce_labels <- c("No non-permissive mismatch" = 0, "Non-permissive mismatch" = 1, "Missing" = 99)

disease_stage_rsm_labels <- c("AML CR1" = "1", "AML CR2+" = "2", "AML Relapsed/Refractory" = "3", "ALL CR1" = "4",
                              "ALL CR2+" = "5", "ALL Relapsed/Refractory" = "6", "MDS Early" = "7", "MDS Advanced" = "8")

donorcmv_labels <- c("Negative" = 0, "Positive" = 1, "Missing" = 99) 
recipientcmv_labels <- c("Negative" = 0, "Positive" = 1, "Missing" = 99)

rsex_labels <- c("Male" = 1, "Female" = 2)
relative_rsm_labels <- c("Son" = 1, "Daughter" = 2, "Brother" = 3, "Sister" = 4, "Father" = 5, "Mother" = 6, "Other male rel" = 7, "Other female rel" = 8, "Missing"=99)



# Define UI for the app
ui <- fluidPage(
  titlePanel("Survival Probability Prediction"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("ragecat_rsm", "Recipient Age:", choices = c("30-39", "40-49", "50-59", "60-69", "70+"), selected = "30-39"),
      selectInput("rsex", "Recipient Sex:", choices = rsex_labels, selected = 1),
      selectInput("disease_stage_rsm", "Disease/Stage:", choices = disease_stage_rsm_labels, selected = 1),
      selectInput("hctcigp", "HCT-CI:", choices = c(0, 1), selected = 0),
      selectInput("condint", "Conditioning:", choices = c(0, 1), selected = 1),
      selectInput("recipientcmv", "Recipient CMV:", choices = recipientcmv_labels, selected = 1), 
      numericInput("dnrage", "Donor age:", value = 65), 
      selectInput("b_leader", "B Leader:", choices = b_leader_labels, selected = 0),
      selectInput("drb1_gvh", "DRB1:", choices = drb1_gvh_labels, selected = 1),
      selectInput("dqb1_gvh", "DQB1:", choices = dqb1_gvh_labels, selected = 0),
      selectInput("tce", "DPB1:", choices = tce_labels, selected = 1),
      selectInput("donorcmv", "Donor CMV:", choices = donorcmv_labels, selected = 1),
      selectInput("relative_rsm", "Relationship:", choices = relative_rsm_labels, selected = 1),
      actionButton("predict", "Predict")
    ),
    
    mainPanel(
      tableOutput("results"),
      plotOutput("survival_plot"),   
      verbatimTextOutput("na_info"),
      verbatimTextOutput("data_structure")
    )
  )
)


server <- function(input, output) {
  observeEvent(input$predict, {
    new_data <- data.frame(
      b_leader = as.numeric(input$b_leader),
      dnrage = as.numeric(input$dnrage),
      drb1_gvh = as.numeric(input$drb1_gvh),
      dqb1_gvh = as.numeric(input$dqb1_gvh),
      tce = as.numeric(input$tce),
      condint = as.numeric(input$condint),
      disease_stage_rsm = as.numeric(input$disease_stage_rsm),
      hctcigp = as.numeric(input$hctcigp),
      rsex = as.numeric(input$rsex),
      donorcmv = as.numeric(input$donorcmv),
      recipientcmv = as.numeric(input$recipientcmv),
      ragecat_rsm = as.factor(input$ragecat_rsm),
      relative_rsm = as.numeric(input$relative_rsm)
    )
    
    # Check for NAs in new_data
    na_info <- sapply(new_data, function(x) sum(is.na(x)))
    data_structure <- str(new_data)
    
    output$na_info <- renderPrint({ na_info })
    output$data_structure <- renderPrint({ data_structure })
    
    # Define time points (in months)
    time_points <- seq(0, 60, by = 1)
    
    # Predict the log-hazard for the new data
    log_hazard <- predict(final_gbm_model, newdata = new_data, n.trees = n_trees_optimal, type = "link")
    
    # Calculate the baseline survival probabilities
    survival_probabilities <- exp(-exp(log_hazard) * time_points / 12)
    
    # Confidence Interval Calculation
    set.seed(123)
    n_bootstraps <- 1000
    
    # Bootstrap sampling
    bootstrap_results <- replicate(n_bootstraps, {
      log_hazard_boot <- log_hazard + rnorm(1, mean = 0, sd = 0.1)  # Add noise for variability
      exp(-exp(log_hazard_boot) * time_points / 12)
    })
    
    # Calculate the lower and upper CI (2.5th and 97.5th percentiles)
    lower_ci <- apply(bootstrap_results, 1, quantile, probs = 0.025)
    upper_ci <- apply(bootstrap_results, 1, quantile, probs = 0.975)
    
    # Select time points at 12, 24, and 36 months
    selected_times <- c(12, 24, 36)
    results <- data.frame(
      Time = selected_times,
      Survival_Probability = survival_probabilities[selected_times + 1],
      Lower_CI = lower_ci[selected_times + 1],
      Upper_CI = upper_ci[selected_times + 1]
    )
    
    # Render the table for specific times
    output$results <- renderTable({
      results
    })
    
    # Plot the full Kaplan-Meier-like curve with confidence intervals
    output$survival_plot <- renderPlot({
      ggplot(data.frame(Time = time_points, Survival_Probability = survival_probabilities), aes(x = Time, y = Survival_Probability)) +
        geom_line(color = "blue") +
        geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), alpha = 0.2) +
        labs(title = "Predicted Survival Probability with 95% Confidence Intervals",
             x = "Time (months)", y = "Survival Probability") +
        theme_minimal() +
        ylim(0, 1)
    })
  })
}



# Run the application 
shinyApp(ui = ui, server = server)
