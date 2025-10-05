library(ggplot2)
library(shiny)

# Regression function
m_fun <- function(x){
  res = 1/((x/3)+0.1)
  return(sin(res))
}

# Creating a dataframe of X and corresponding Y
generate_sample <- function(n, alpha, beta, sigma2 = 1){
  set.seed(123)
  x <- rbeta(n, alpha, beta)
  y <- m_fun(x) + rnorm(n, mean = 0, sd = sqrt(sigma2))
  data.frame(X = x, Y = y)
}
# Fit 4th degree polynomial (only if >=5 points)
fit_poly4 <- function(df) {
  if (nrow(df) < 5) {
    return(NULL)
  }
  lm(Y ~ poly(X, 4, raw = TRUE), data = df)
}

# Compute second derivative of fitted polynomial at x
poly4_second_derivative <- function(model, x) {
  coefs <- coef(model)
  b2 <- ifelse(length(coefs) >= 3, coefs[3], 0)
  b3 <- ifelse(length(coefs) >= 4, coefs[4], 0)
  b4 <- ifelse(length(coefs) >= 5, coefs[5], 0)
  2*b2 + 6*b3*x + 12*b4*x^2
}

# Estimate theta22 and sigma2 (using blocking)
estimate_thetasigma <- function(data, N) {
  n <- nrow(data)
  if (N > n) return(NULL) # making sure we don't have more blocks than observations
  if (N==1){
    indices = rep(1,n)
  }else{
    indices <- cut(seq_len(n), breaks = N, labels = FALSE)
  }
  if(any(table(indices) < 5)) return(NULL)

  theta22_vals <- rep(NA, n)
  resid_vals <- rep(NA, n)
  
  for (j in 1:N) {
    block_idx <- which(indices == j)
    block_data <- data[block_idx, , drop = FALSE]
    fit <- fit_poly4(block_data)
    
    if (!is.null(fit)) {
      m_dd <- poly4_second_derivative(fit, block_data$X)
      theta22_vals[block_idx] <- m_dd^2
      resid_vals[block_idx] <- (fit$residuals)^2
    }
  }
  
  theta22_hat <- mean(theta22_vals)  # disregarding the NA values when taking the mean
  if (is.nan(theta22_hat)) theta22_hat <- NA_real_
  sigma2_hat <- sum(resid_vals) / (n - 5*N)  # avoid division by 0
  RSS <- sum(resid_vals)
  if (is.na(theta22_hat) || is.na(sigma2_hat)) return(NULL)
  list(theta22 = theta22_hat, sigma2 = sigma2_hat, RSS = RSS)
}


# AMISE bandwidth formula
h_AMISE <- function(n, sigma2, theta22, support_length) {
  if (is.null(sigma2) || is.null(theta22) || is.na(sigma2) || is.na(theta22) ||
      sigma2 <= 0 || theta22 <= 0) return(NA_real_)
  n^(-1/5) * ((35 * sigma2 * support_length) / theta22)^(1/5)
}

# data: data.frame with X and Y
# N: candidate number of blocks
# Returns C_p(N)
compute_Cp <- function(data, N, RSS_Nmax) {
  n <- nrow(data)
  Nmax <- max(min(floor(n / 20), 5), 1)
  # Compute RSS(N)
  est <- estimate_thetasigma(data, N)
  if(!is.null(est)){
    RSS_N <- est$RSS
  } else{
    RSS_N <- NA_real_
  }
  # Check for invalid RSS
  if (is.na(RSS_N) || is.na(RSS_Nmax) || RSS_Nmax <= 0) return(NA_real_)
  
  Cp <- RSS_N / (RSS_Nmax / (n - 5*Nmax)) - (n - 10*N)
  Cp
}
# data: data.frame with X and Y
# N_candidates: vector of candidate N values
find_optimal_N <- function(data, N_candidates = 1:20) {
  n <- nrow(data)
  # Compute RSS(N_max) so we don't do it multiple times in the loop
  Nmax <- max(min(floor(n / 20), 5), 1)
  est_max <- estimate_thetasigma(data, Nmax)
  RSS_Nmax <- est_max$RSS
  
  Cp_values <- sapply(N_candidates, function(N) compute_Cp(data, N, RSS_Nmax))
  
  # remove NA values
  valid_idx <- which(!is.na(Cp_values))
  if (length(valid_idx) == 0) return(NA_integer_)
  
  N_candidates[valid_idx][which.min(Cp_values[valid_idx])]
}


# -----------------------------
# Shiny app
# -----------------------------

ui <- fluidPage(
  titlePanel("Local Linear Estimator: AMISE Bandwidth"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("n", "Sample size (n):", min = 100, max = 5000, value = 3000, step = 100),
      sliderInput("N", "Number of blocks (N):", min = 1, max = 10, value = 5),
      sliderInput("alpha", "Beta(α, β): α", min = 0.5, max = 5, value = 3, step = 0.1),
      sliderInput("beta", "Beta(α, β): β", min = 0.5, max = 5, value = 3, step = 0.1),
      sliderInput("sigma2", "σ² (variance of gaussian noise):", min = 0.2, max = 2, value = 1, step = 0.2),
    ),
    mainPanel(
      verbatimTextOutput("bandwidthOut"),
      plotOutput("samplePlot"),
      plotOutput("beta_hist", height = "300px"),
      plotOutput("heatmap_plot", height = "400px"),
      plotOutput("h_vs_n_plot", height = "300px"),
      plotOutput("h_vs_N_plot", height = "300px")
    )
  )
)

server <- function(input, output, session) {
  sim_data <- reactive({
    generate_sample(input$n, input$alpha, input$beta, input$sigma2)
  })
  
  estimates <- reactive({
    req(sim_data())
    estimate_thetasigma(sim_data(), input$N)
  })
  
  output$bandwidthOut <- renderPrint({
    est <- estimates()
    n <- input$n
    support_length <- 1  # Beta is supported on [0,1]
    h <- h_AMISE(n, est$sigma2, est$theta22, support_length)
    cat("For N:", input$N, "(Slider input) \n")
    cat("Estimated theta22:", est$theta22, "\n")
    cat("Estimated sigma²:", est$sigma2, "\n")
    cat("h_AMISE:", h, "\n")
  })
  
  output$samplePlot <- renderPlot({
    df <- sim_data()
    req(df)
    ggplot(df, aes(x = X, y = Y)) +
      geom_point(alpha = 0.5) +
      stat_function(fun = m_fun, color = "red", size = 1.2) +
      labs(title = "Simulated Sample with True Regression Function",
           subtitle = "Red curve = true m(x)") +
      theme_minimal()
  })
  output$beta_hist <- renderPlot({
    df <- sim_data()
    req(df)
    ggplot(df, aes(x = X)) +
      geom_histogram(aes(y = ..count..), bins = 30, fill = "skyblue", color = "white") +
      #stat_function(fun = function(x) dbeta(x, input$alpha, input$beta),
      #              color = "red", size = 1.2) +
      labs(title = paste0("Distribution of X ~ Beta(", input$alpha, ", ", input$beta, ")"),
           x = "X", y = "Number of observations") +
      theme_minimal()
  })
  output$heatmap_plot <- renderPlot({
    n <- input$n
    sigma2 <- input$sigma2
    
    # Define grid of alpha and beta values
    alpha_vals <- seq(0.5, 5, by = 0.5)
    beta_vals  <- seq(0.5, 5, by = 0.5)
    
    grid <- expand.grid(alpha = alpha_vals, beta = beta_vals)
    
    # For each (alpha, beta), generate sample and estimate h_AMISE
    h_vals <- sapply(1:nrow(grid), function(i) {
      df <- generate_sample(n, grid$alpha[i], grid$beta[i], sigma2)
      N <- find_optimal_N(df, N_candidates = 1:10)
      est <- estimate_thetasigma(df, N)
      h_AMISE(n, est$sigma2, est$theta22, support_length = 1)
    })
    
    grid$h <- h_vals
    
    ggplot(grid, aes(x = alpha, y = beta, fill = h)) +
      geom_tile() +
      scale_fill_viridis_c(option = "plasma") +
      labs(title = expression(paste("Heatmap of ", h[AMISE], " vs. alpha, beta (optimal N)")),
           x = "Alpha",
           y = "Beta",
           fill = "h_AMISE") +
      theme_minimal()
  })
  output$h_vs_n_plot <- renderPlot({
    n_vals <- seq(100, 5000, by = 50)
    h_vals <- numeric(length(n_vals))
    N_opt_vals <- numeric(length(n_vals))
    
    for (i in 1:length(n_vals)) {
      df <- generate_sample(n_vals[i], input$alpha, input$beta, input$sigma2)
      N_opt <- find_optimal_N(df, N_candidates = 1:10)
      if (!is.na(N_opt)){
      N_opt_vals[i] <- N_opt
      } else {
        N_opt_vals[i] <- 1
      }
      est <- estimate_thetasigma(df, N_opt)
      if (!is.null(est)) {
        h_vals[i] <- h_AMISE(n_vals[i], est$sigma2, est$theta22, support_length = 1)
      } else {
        h_vals[i] <- NA_real_
      }
    }
    
    df_plot <- data.frame(n = n_vals, h_AMISE = h_vals)
    
    ggplot(df_plot, aes(x = n, y = h_AMISE)) +
      geom_line(color = "blue", size = 1.2) +
      geom_point(color = "darkblue") +
      labs(title = expression(paste("h"[AMISE], " vs Sample Size n (Optimal N)")),
           x = "Sample size n", y = expression(h[AMISE])) +
      theme_minimal()
  })
  
  output$h_vs_N_plot <- renderPlot({
    N_vals <- 1:10
    h_vals <- numeric(length(N_vals))
    n <- input$n
    df <- generate_sample(n, input$alpha, input$beta, input$sigma2)
    
    for (i in 1:length(N_vals)) {
      est <- estimate_thetasigma(df, N_vals[i])
      if (!is.null(est)) {
        h_vals[i] <- h_AMISE(n, est$sigma2, est$theta22, support_length = 1)
      } else {
        h_vals[i] <- NA_real_
      }
    }
    
    df_plot <- data.frame(N = N_vals, h_AMISE = h_vals)
    df_plot <- df_plot[!is.na(df_plot$h_AMISE), ]
    
    ggplot(df_plot, aes(x = N, y = h_AMISE)) +
      geom_line(color = "purple", size = 1.2) +
      geom_point(color = "violetred") +
      scale_x_continuous(breaks = N_vals) +
      labs(title = expression(paste("h"[AMISE], " vs Number of Blocks N")),
           x = "Number of blocks N", y = expression(h[AMISE])) +
      theme_minimal()
  })

}

shinyApp(ui, server)
