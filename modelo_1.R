# primer modelo 1
# Primer modelo 
MCMC_1 <- function(y, B) {
      # tamaño y estadísticos
      n      <- length(y)
      sum_y  <- sum(y)
      mean_y <- mean(y)
      var_y  <- var(y)
      # hiperparametros
      mu0  <- 13.495 
      g20  <- 11.382
      nu0  <- 1
      s20  <- 1.182
      # almacenamiento
      THETA <- matrix(data = NA, nrow = B, ncol = 2)
      LL    <- matrix(data = NA, nrow = B, ncol = 1)
      # Valores iniciales
      # valor inicial: simular de la previa o estadísticos, ya que tiene sentido usarlos
      theta <- mean_y
      sig2  <- var_y
      # cadena
      set.seed(2023)
      for (b in 1:B) {
            # 2.1 actualizar el valor de theta
            g2n   <- 1/(1/g20 + n/sig2)      
            mun   <- (mu0/g20 + sum_y/sig2)*g2n
            theta <- rnorm(n = 1, mean = mun, sd = sqrt(g2n))
            # 2.2 actualizar el valor de sigma^2
            nun   <- nu0 + n
            s2n   <- (nu0*s20 + (n-1)*var_y + n*(mean_y - theta)^2)
            sig2 <- 1/rgamma(n = 1, shape = nun/2, rate = s2n/2)
            # 2.3 almacenar valores
            THETA[b,] <- c(theta,sig2)
            # log-verosimilitud
            LL[b] <- sum(dnorm(x = y, mean = THETA[b,1], sd = sqrt(THETA[b,2]), log = T))
      }
      # fin de la cadena
      # salida
      THETA <- as.data.frame(THETA[-c(1:1000),])
      LL    <- as.data.frame(LL[-c(1:1000),])
      colnames(THETA) <- c("theta", "sig2")
      colnames(LL) <- c("ll")
      return(list(THETA = THETA, LL = LL))
}
