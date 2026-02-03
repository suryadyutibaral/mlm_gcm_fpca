#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(MASS)
library(plotly)

ui <- fluidPage(
  titlePanel("Covariance Topology Explorer"),

sidebarPanel(
  width = 3,
  tags$style("
    .control-label { font-size: 12px; }
    .shiny-input-container { margin-bottom: 6px; }
  "),

  sliderInput("n_subj", "Subjects",
              min = 50, max = 2000, value = 1000, step = 50),

  sliderInput("n_time", "Time points",
              min = 5, max = 500, value = 100, step = 5),

  numericInput("beta0", "β₀", value = 1),
  numericInput("beta1", "β₁", value = 2),
  numericInput("sigma_eps", "Noise SD", value = 0.5, min = 0),

  tags$details(
    tags$summary("Random effects"),
    numericInput("sigma_b0", "Var(b₀)", value = 1.0, min = 0),
    numericInput("sigma_b1", "Var(b₁)", value = 1.0, min = 0),
    numericInput("sigma_b01", "Cov(b₀,b₁)", value = 0)
  ),

  actionButton("simulate", "Update", class = "btn-primary btn-sm")
),

    mainPanel(
      plotlyOutput("cov_plot", height = "600px")
    )
  )


server <- function(input, output) {

  sim_data <- eventReactive(input$simulate, {

    set.seed(123)

    # Time grid
    time <- seq(0, 1, length.out = input$n_time)

    # Random effects covariance
    Sigma_b <- matrix(
      c(input$sigma_b0,
        input$sigma_b01,
        input$sigma_b01,
        input$sigma_b1),
      nrow = 2
    )

    # Random effects
    b <- MASS::mvrnorm(input$n_subj, mu = c(0, 0), Sigma = Sigma_b)

    # Generate curves
    X <- matrix(NA, input$n_subj, input$n_time)

    for (i in seq_len(input$n_subj)) {
      X[i, ] <-
        (input$beta0 + b[i, 1]) +
        (input$beta1 + b[i, 2]) * time +
        rnorm(input$n_time, sd = input$sigma_eps)
    }

    # Center
    mu_hat <- colMeans(X)
    X_centered <- sweep(X, 2, mu_hat)

    # Sample covariance
    G_hat <- t(X_centered) %*% X_centered / input$n_subj

    # Long format for plotting
    cov_df <- expand.grid(
      t1 = time,
      t2 = time
    )

    cov_df$cov <- as.vector(G_hat)

    cov_df
  })

  output$cov_plot <- renderPlotly({

    cov_df <- sim_data()

    cmax <- max(abs(cov_df$cov))

plot_ly(
  data = cov_df,
  x = ~t1, y = ~t2, z = ~cov,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 2,
    opacity = 0.6,
    color = ~cov,
    colorscale = "Viridis",
    cmin = min(cov_df$cov),
    cmax = max(cov_df$cov),
    showscale = TRUE,
    colorbar = list(title = "Covariance")
  )
) %>%
      layout(
        title = "Empirical Covariance Point Cloud",
        scene = list(
          xaxis = list(title = "Time t"),
          yaxis = list(title = "Time s"),
          zaxis = list(title = "Covariance")
        )
      )
  })
}

shinyApp(ui, server)
