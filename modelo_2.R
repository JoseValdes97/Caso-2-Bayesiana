# Segundo modelo 2 
# Segundo modelo 
MCMC_2 <- function(y, B, nj, yb, s2) {
      # tamaños
      n <- sum(nj)
      m <- length(nj)
      # hiperparametros
      mu0  <- 13.495 
      g20  <- 11.382
      eta0 <- 1  
      t20  <- 1.182
      nu0   <- 1
      s20  <- t20
      # valores iniciales
      theta <- yb
      sig2  <- mean(s2)
      mu    <- mean(theta)
      tau2  <- var(theta)
      # almacenamiento
      THETA <- matrix(data = NA, nrow = B, ncol = m+3)
      LL    <- matrix(data = NA, nrow = B, ncol = 1)
      # cadena
      set.seed(2023)
      for (b in 1:B) {
            # inicializar vtheta y actualizar theta
            vtheta <- 1/(1/tau2 + nj/sig2)
            theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + nj*yb/sig2), sd = sqrt(vtheta))
            # actualizar sigma^2
            sig2 <- 1/rgamma(n = 1, shape = 0.5*(nu0 + n), rate = 0.5*(nu0*s20 + sum((nj-1)*s2 + nj*(yb - theta)^2)))
            # inicializar vmu y actualizar mu
            vmu <- 1/(1/g20 + m/tau2)
            mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu)) 
            # actualizar tau^2
            tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0 + m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
            # almacenar valores
            THETA[b,] <- c(theta, sig2, mu, tau2)
            # log-verosimilitud
            LL[b] <- sum(dnorm(x = y, mean = rep(theta, nj), sd = sqrt(sig2), log = T))
      }
      # fin de la cadena
      # salida
      THETA <- as.data.frame(THETA[-c(1:1000),])
      LL    <- as.data.frame(LL[-c(1:1000),])
      colnames(THETA) <- c(paste0("theta",1:m), "sig2", "mu", "tau2")
      colnames(LL) <- c("ll")
      return(list(THETA = THETA, LL = LL))
}
