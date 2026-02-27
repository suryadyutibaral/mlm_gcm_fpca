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
  
  withMathJax(),
  
  tags$div(
    style = "margin-bottom:20px; font-size:14px; line-height:1.5;",
    
    tags$h4("Welcome!"),
    
    tags$p("This app provides an interactive demonstration of how latent 
    functional structure generates covariance surfaces (covariance matrices 
    visualized over continuous time). The covariance surface serves as a 
    bridge between parametric longitudinal modeling approaches, such as 
    Structural Equation Models (SEM) and Multilevel Models (MLM), and the 
    nonparametric framework of Functional Principal Component Analysis (FPCA)."),
    
    tags$p("In growth curve settings, the time-dependent loadings in SEM and the 
    functional forms specified in MLM both describe systematic temporal change 
    processes. FPCA recovers this structure directly from the data. Its 
    eigenfunctions provide nonparametric approximations to the underlying
    functions of time that drive between-subject variability. Indeed, the 
    loadings in SEM and functional forms of the MLM are both nonparametrically 
    obtained from the functional principal components."),
    
    tags$p("By exploring how changes in the generating function and variance components 
    alter the covariance surface and its eigenfunctions, one can see how these 
    seemingly distinct analytical approaches are unified through the covariance. 
    This perspective highlights how parametric and nonparametric methods converge 
    in the analysis of longitudinal data pushing us further towards realizing 
    Nesselroade's aim of fully understanding both the 'warp and woof of the 
    developmental fabric'."),
    
    tags$h4("Method"),
    
    tags$p("Data are generated from the multilevel model:"),
    
    helpText("$$ y_i(t) = \\alpha_i f(t) + e_i(t) $$"),
    
    tags$p("where "),
    
    tags$ul(
      tags$li("f(t) is a user-defined function on the interval [0,1],"),
      tags$li("\\( \\alpha_i = \\beta + b_i \\) is a subject-specific coefficient,"),
      tags$li("\\( b_i \\sim N(0, \\sigma_b^2) \\) represents between-subject variability,"),
      tags$li("\\( e_i(t) \\sim N(0, \\sigma_\\varepsilon^2) \\) represents measurement noise.")
    ),
    
    tags$h4("Instructions"),
    
    tags$p("We have options to vary number of subjects and time points. As you know from Stats 101, more 
    people and time points lead to better estimation of the covariance matrix. That is true here too. However,
    the maximum possible values here are set to reasonable values to ensure the app does not crash. Further, 
    when computing FPCA please set smaller values for time points to enable faster computation."),
    
    tags$p("There are buttons for each element of the data generating model. Feel free to make changes and explore on 
    your own. I would encourage you to change f(t) and measurement error variance and see how that affects the 
    estimation of f(t). Here is an fun activity. Plot the covariance surface with Noise SD = .5, and then plot 
    the same with Noie SD = 1. Can you explain why this happens?"),
    
    tags$p("When changing the f(t), please use 'time' and not 't'. For example, if you want f(t) = 
    sin(t), you should write sin(time) in the 'f(time) =' box. Most of the functions you use will probably not be 
    normalized, so we normalize them for you. This is particularly important when you compare the
    estimated functions to the generating functions. Remember, FPCA returns only orthonormal functions."),
    
    tags$p("The app was coded to be used in a laptop or a desktop. You can use it on your mobile, but note that the covariance 
    surface may be difficult to 'explore.' The FPCA computation should work out better. You might want to hold you mobile
    in landscape mode for a better experience."),
    
    
    tags$p("Feel free to visualize the empirical covariance surface and examine how functional principal component analysis (FPCA)
    recovers the underlying generating function. Last bit of advice"),
    tags$h5("Be creative, and most importantly,..."),
    tags$h5("Have Fun!")
  ),
  
  sidebarLayout(
    
    sidebarPanel(
      width = 3,
      
      tags$style("
        .control-label { font-size: 12px; }
        .shiny-input-container { margin-bottom: 6px; }
      "),
      
      sliderInput("n_subj", "Subjects",
                  min = 50, max = 1000, value = 1000, step = 50),
      
      sliderInput("n_time", "Time points",
                  min = 5, max = 300, value = 100, step = 5),
      
      numericInput("beta", "β", value = 1),
      numericInput("sigma_eps", HTML("Noise SD (σ<sub>ε</sub>)"), value = 0.5, min = 0),
      
      textInput(
        "f_t",
        "f(time) =",
        value = "time"
      ),
      
      checkboxInput("Orthonormalize", "Normalize", TRUE), 
      
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
