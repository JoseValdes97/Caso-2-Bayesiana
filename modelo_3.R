# Tercer modelo
MCMC_3 <- function(y,B, nj, yb, s2) {
      # tamaños
      n <- sum(nj)
      m <- length(nj)
      # hiperparametros
      mu0  <- 13.495 
      g20  <- 11.382
      eta0 <- 1  
      t20  <- 1.182
      nu   <- 1  
      al0  <- 1
      be0  <- 0.846 
      # valores iniciales
      theta <- yb
      sig2  <- s2  # sigma_j^2
      mu    <- mean(theta)
      tau2  <- var(theta)
      # iniciar sigma random round(rgamma(1,shape=al0/2,rate=be0/2),3)
      ups2  <- 1.33  # sigma^2
      # almacenamiento
      THETA <- matrix(data = NA, nrow = B, ncol = 2*m+4)
      LL    <- matrix(data = NA, nrow = B, ncol = 1)
      # cadena
      set.seed(2023)
      for (b in 1:B) {
            # actualizar theta
            vtheta <- 1/(1/tau2 + nj/sig2)
            theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + nj*yb/sig2), sd = sqrt(vtheta))
            # actualizar sigma_j^2
            sig2 <- 1/rgamma(n = m, shape = 0.5*(nu + nj), rate = 0.5*(nu*ups2 + (nj-1)*s2 + nj*(yb - theta)^2))
            # actualizar mu
            vmu <- 1/(1/g20 + m/tau2)
            mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu))
            # actualizar tau2
            tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0+m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
            # actualizar sigma^2
            ups2 <- rgamma(n = 1, shape = al0 + 0.5*m*nu, rate = be0 + 0.5*nu*sum(1/sig2))
            # almacenar
            THETA[b,] <- c(theta, sig2, mu, tau2, nu, ups2)
            # log-verosimilitud
            LL[b] <- sum(dnorm(x = y, mean = rep(theta, nj), sd = sqrt(rep(sig2, nj)), log = T))
      }
      # fin de la cadena
      # salida
      colnames(THETA) <- c(paste0("theta", 1:m), paste0("sig2", 1:m), "mu", "tau2", "nu", "ups2")
      colnames(LL) <- c("ll")
      THETA <- as.data.frame(THETA)
      LL    <- as.data.frame(LL)
      return(list(THETA = THETA, LL = LL))
}
