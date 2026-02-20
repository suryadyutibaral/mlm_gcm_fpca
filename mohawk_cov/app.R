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
library(fdapace)
library(pracma)
ui <- fluidPage(
  
  titlePanel("Covariance Topology Explorer"),
  
  sidebarLayout(
    
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
      
      numericInput("beta", "β", value = 1),
      numericInput("sigma_eps", "Noise SD", value = 0.5, min = 0),
      
      textInput(
        "f_t",
        "f(time) =",
        value = "time"
      ),
      
      checkboxInput("Orthonormalize", "Orthonormalize", FALSE), 
      
      tags$details(
        tags$summary("Random effects"),
        numericInput("sigma_b", "Var(b)", value = 1.0, min = 0)
      ),
      
      actionButton("simulate", "Update", class = "btn-primary btn-sm")
    ),
    
    mainPanel(
      tabsetPanel(
    
        tabPanel(
          "Covariance Surface",
          plotlyOutput("cov_plot", height = "600px")
        ),
    
        tabPanel(
          "FPCA Eigenfunction",
          checkboxInput("compute_fpca", "Compute FPCA", value = FALSE),
          plotlyOutput("eig_plot", height = "600px")
        )
    
      )
    )

    
  )
)

server <- function(input, output) {

sim_data <- eventReactive(input$simulate, {

  set.seed(123)

  time <- seq(0, 1, length.out = input$n_time)

  f_time <- tryCatch({
    eval(parse(text = input$f_t), envir = list(time = time))
  }, error = function(e) {
    time
  })

  if (length(f_time) != input$n_time) {
    f_time <- time
  }
  
  if(input$Orthonormalize){
    f_time <- f_time/sqrt(trapz(time, f_time^2))
  } else{
    f_time <- f_time
  }

  b <- rnorm(input$n_subj, mean = 0, sd = input$sigma_b)

  X <- matrix(NA, input$n_subj, input$n_time)

  for (i in seq_len(input$n_subj)) {
    X[i, ] <-
      (input$beta + b[i]) * f_time +
      rnorm(input$n_time, sd = input$sigma_eps)
  }

  mu_hat <- colMeans(X)
  X_centered <- sweep(X, 2, mu_hat)

  G_hat <- t(X_centered) %*% X_centered / input$n_subj

  cov_df <- expand.grid(
    t1 = time,
    t2 = time
  )

  cov_df$cov <- as.vector(G_hat)

  # ----- FPCA -----
  Ly <- lapply(seq_len(nrow(X)), function(i) X[i, ])
  Lt <- lapply(seq_len(nrow(X)), function(i) time)
  
  fpca_fit <- NULL

  if (input$compute_fpca) {
    Ly <- lapply(seq_len(nrow(X)), function(i) X[i, ])
    Lt <- lapply(seq_len(nrow(X)), function(i) time)
  
    fpca_fit <- FPCA(Ly, Lt, optns = list(methodMuCovEst = "smooth", shrink =  TRUE, methodXi = "IN", error = TRUE, dataType = "DenseWithMV"))
  }

  list(
    cov_df = cov_df,
    time = time,
    fpca = fpca_fit,
    f_time = f_time
  )
})

output$eig_plot <- renderPlotly({

  sim <- sim_data()
  time <- sim$time
  f_time <- sim$f_time
  p <- plot_ly()
  
    p <- p %>%
    add_lines(
      x = time,
      y = f_time,
      name = "f(t)",
      line = list(dash = "dash")
    )
    
    if (is.null(sim$fpca)) {
    
    return(
      p %>%
        layout(
          title = "FPCA Not Computed",
          xaxis = list(title = "Time"),
          yaxis = list(title = "Value")
        )
    )
  }
    
  fpca_fit <- sim$fpca
  phi_mat <- fpca_fit$phi[, 1, drop = FALSE]
  
  p <- p %>%
      add_lines(
        x = time,
        y = phi_mat,
        name = paste0("Eigenfunction ", 1,
                      " (λ=", round(fpca_fit$lambda[1], 3), ")")
      )
  p %>%
    layout(
      title = "FPCA Eigenfunctions",
      xaxis = list(title = "Time"),
      yaxis = list(title = "Eigenfunction Value")
    )
})
output$cov_plot <- renderPlotly({ 
  cov_df <- sim_data() 
  cmax <- max(abs(cov_df$cov_df)) 
  plot_ly(data = cov_df$cov_df, 
          x = ~t1, y = ~t2, z = ~cov, 
          type = "scatter3d", mode = "markers", 
          marker = list( size = 2, opacity = 0.6, 
                         color = ~cov, colorscale = "Viridis", 
                         cmin = min(cov_df$cov), cmax = max(cov_df$cov), 
                         showscale = TRUE, 
                         colorbar = list(title = "Covariance") ) ) %>% 
    layout( title = "Empirical Covariance Point Cloud", 
            scene = list( xaxis = list(title = "Time t"), 
                          yaxis = list(title = "Time s"), 
                          zaxis = list(title = "Covariance") ) ) 
}) 
}

shinyApp(ui, server)
