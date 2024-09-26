library(shiny)
library(gbm)



# Define labels for categorical variables
ragecat_rsm_levels <- c("<=18", "18-29", "30-39", "40-49", "50-59", "60+") 
hctcigp_levels <- c("0", "1", "2", "3", "4")

b_leader_labels <- c("Match" = 0, "Mismatch" = 1)
drb1_gvh_labels <- c("Match" = 0, "Mismatch" = 1)
dqb1_gvh_labels <- c("Match" = 0, "Mismatch" = 1)
tce_labels <- c("No non-permissive mismatch" = 0, "Non-permissive mismatch" = 1, "Missing" = 99)

disease_stage_rsm_labels <- c("AML CR1" = "1", "AML CR2+" = "2", "AML Relapsed/Refractory" = "3", 
                              "ALL CR1" = "4", "ALL CR2+" = "5", "ALL Relapsed/Refractory" = "6", 
                              "MDS Early" = "7", "MDS Advanced" = "8")

donorcmv_labels <- c("Negative" = 0, "Positive" = 1, "Missing" = 99) 
recipientcmv_labels <- c("Negative" = 0, "Positive" = 1, "Missing" = 99)

rsex_labels <- c("Male" = 1, "Female" = 2)
condint_labels <- c("MAC" = 1, "RIC/NMA" = 2)
relative_rsm_labels <- c("Son" = 1, "Daughter" = 2, "Brother" = 3, "Sister" = 4, 
                         "Father" = 5, "Mother" = 6, "Other male rel" = 7, 
                         "Other female rel" = 8, "Missing" = 99)



# Define UI for the app
ui <- fluidPage(
  titlePanel("Survival Probability Prediction"),
  
  sidebarLayout(
    sidebarPanel(
      # Recipient Information
      selectInput("ragecat_rsm", "Recipient Age:", choices = ragecat_rsm_levels, selected = "30-39"),
      selectInput("rsex", "Recipient Sex:", choices = rsex_labels, selected = 1),
      selectInput("disease_stage_rsm", "Disease/Stage:", choices = disease_stage_rsm_labels, selected = 1),
      selectInput("hctcigp", "HCT-CI:", choices = hctcigp_levels),
      selectInput("condint", "Conditioning:", choices = condint_labels, selected = 1), 
      selectInput("recipientcmv", "Recipient CMV:", choices = recipientcmv_labels, selected = 1), 
      
      # Donor 1 Information
      h4("Donor 1"),
      numericInput("dnrage1", "Donor age:", value = 65), 
      selectInput("b_leader1", "B Leader:", choices = b_leader_labels, selected = 0),
      selectInput("drb1_gvh1", "DRB1:", choices = drb1_gvh_labels, selected = 1),
      selectInput("dqb1_gvh1", "DQB1:", choices = dqb1_gvh_labels, selected = 0),
      selectInput("tce1", "DPB1:", choices = tce_labels, selected = 1),
      selectInput("donorcmv1", "Donor CMV:", choices = donorcmv_labels, selected = 1),
      selectInput("relative_rsm1", "Relationship:", choices = relative_rsm_labels, selected = 99),
      
      # Donor 2 Information
      h4("Donor 2"),
      numericInput("dnrage2", "Donor age:", value = 20), 
      selectInput("b_leader2", "B Leader:", choices = b_leader_labels, selected = 0),
      selectInput("drb1_gvh2", "DRB1:", choices = drb1_gvh_labels, selected = 1),
      selectInput("dqb1_gvh2", "DQB1:", choices = dqb1_gvh_labels, selected = 0),
      selectInput("tce2", "DPB1:", choices = tce_labels, selected = 1),
      selectInput("donorcmv2", "Donor CMV:", choices = donorcmv_labels, selected = 1),
      selectInput("relative_rsm2", "Relationship:", choices = relative_rsm_labels, selected = 99),
      
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
    # Function to create new data for a donor
    create_donor_data <- function(dnrage, b_leader, drb1_gvh, dqb1_gvh, tce, donorcmv, relative_rsm) {
      data.frame(
        b_leader = as.numeric(b_leader),
        dnrage = as.numeric(dnrage),
        drb1_gvh = as.numeric(drb1_gvh),
        dqb1_gvh = as.numeric(dqb1_gvh),
        tce = as.numeric(tce),
        condint = as.numeric(input$condint),
        disease_stage_rsm = as.numeric(input$disease_stage_rsm),
        hctcigp = as.numeric(input$hctcigp),
        rsex = as.numeric(input$rsex),
        donorcmv = as.numeric(donorcmv),
        recipientcmv = as.numeric(input$recipientcmv),
        ragecat_rsm = as.factor(input$ragecat_rsm),
        relative_rsm = as.numeric(relative_rsm)
      )
    }
    
    # Create data for Donor 1 
    new_data_donor1 <- create_donor_data(
      dnrage = input$dnrage1, 
      b_leader = input$b_leader1, 
      drb1_gvh = input$drb1_gvh1,
      dqb1_gvh = input$dqb1_gvh1,
      tce = input$tce1,
      donorcmv = input$donorcmv1,
      relative_rsm = input$relative_rsm1
    )
    
    # Create data for Donor 2
    new_data_donor2 <- create_donor_data(
      dnrage = input$dnrage2,
      b_leader = input$b_leader2,
      drb1_gvh = input$drb1_gvh2,
      dqb1_gvh = input$dqb1_gvh2,
      tce = input$tce2,
      donorcmv = input$donorcmv2,
      relative_rsm = input$relative_rsm2
    )
    
    # Define time points (in months)
    time_points <- seq(0, 60, by = 1)
    
    # Select time points at 12, 24, and 36 months
    selected_times <- c(12, 24, 36)
    
    # --- Calculate survival probabilities and CIs for Donor 1 ---
    # Predict the log-hazard 
    log_hazard_donor1 <- predict(final_gbm_model, newdata = new_data_donor1, n.trees = n_trees_optimal, type = "link")
    
    # Calculate baseline survival probabilities
    survival_probabilities_donor1 <- exp(-exp(log_hazard_donor1) * time_points / 12)
    
    # Bootstrap sampling for Donor 1
    set.seed(123)
    n_bootstraps <- 1000
    bootstrap_results_donor1 <- replicate(n_bootstraps, {
      log_hazard_boot <- log_hazard_donor1 + rnorm(1, mean = 0, sd = 0.1) 
      exp(-exp(log_hazard_boot) * time_points / 12)
    })
    
    # Calculate CI for Donor 1
    lower_ci_donor1 <- apply(bootstrap_results_donor1, 1, quantile, probs = 0.025)
    upper_ci_donor1 <- apply(bootstrap_results_donor1, 1, quantile, probs = 0.975)
    
    # --- Calculate survival probabilities and CIs for Donor 2 ---
    # Predict the log-hazard 
    log_hazard_donor2 <- predict(final_gbm_model, newdata = new_data_donor2, n.trees = n_trees_optimal, type = "link")
    
    # Calculate baseline survival probabilities
    survival_probabilities_donor2 <- exp(-exp(log_hazard_donor2) * time_points / 12)
    
    # Bootstrap sampling for Donor 2
    set.seed(123) 
    bootstrap_results_donor2 <- replicate(n_bootstraps, {
      log_hazard_boot <- log_hazard_donor2 + rnorm(1, mean = 0, sd = 0.1) 
      exp(-exp(log_hazard_boot) * time_points / 12)
    })
    
    # Calculate CI for Donor 2
    lower_ci_donor2 <- apply(bootstrap_results_donor2, 1, quantile, probs = 0.025)
    upper_ci_donor2 <- apply(bootstrap_results_donor2, 1, quantile, probs = 0.975)
    
    # Select time points at 12, 24, and 36 months
    selected_times <- c(12, 24, 36)
    
    # Create separate result tables for each donor
    results_donor1 <- data.frame(
      Time = selected_times,
      Survival_Probability = survival_probabilities_donor1[selected_times + 1],
      Lower_CI = lower_ci_donor1[selected_times + 1],
      Upper_CI = upper_ci_donor1[selected_times + 1]
    )
    
    results_donor2 <- data.frame(
      Time = selected_times,
      Survival_Probability = survival_probabilities_donor2[selected_times + 1],
      Lower_CI = lower_ci_donor2[selected_times + 1],
      Upper_CI = upper_ci_donor2[selected_times + 1]
    )
    
    # Combine results into a tabular format
    results <- data.frame(
      Time = paste(selected_times, "months"), 
      "Donor_1" = paste0(round(results_donor1$Survival_Probability, 2), 
                                      " (", 
                                      round(results_donor1$Lower_CI, 2), 
                                      "-", 
                                      round(results_donor1$Upper_CI, 2),
                                      ")"),
      "Donor_2" = paste0(round(results_donor2$Survival_Probability, 2), 
                                      " (", 
                                      round(results_donor2$Lower_CI, 2), 
                                      "-", 
                                      round(results_donor2$Upper_CI, 2),
                                      ")")
    )
    
    # Render the table
    output$results <- renderTable({ 
      results 
    })
    
    # Plot survival curves for both donors
    survival_data <- rbind(
      data.frame(Time = time_points, Survival_Probability = survival_probabilities_donor1, Donor = "Donor 1", Lower_CI = lower_ci_donor1, Upper_CI = upper_ci_donor1),
      data.frame(Time = time_points, Survival_Probability = survival_probabilities_donor2, Donor = "Donor 2", Lower_CI = lower_ci_donor2, Upper_CI = upper_ci_donor2)
    )
    
    output$survival_plot <- renderPlot({
      ggplot(survival_data, aes(x = Time, y = Survival_Probability, color = Donor)) +
        geom_line() +
        geom_ribbon(aes(ymin = Lower_CI, ymax = Upper_CI, fill = Donor), alpha = 0.05, colour = NA) + # Remove border
        labs(title = "Predicted Survival Probability with 95% Confidence Intervals",
             x = "Time (months)", y = "Survival Probability") +
        theme_minimal() +
        ylim(0, 1)
    })
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
