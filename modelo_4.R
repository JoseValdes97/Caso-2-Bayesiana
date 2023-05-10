# modelo 4
MCMC_4 <- function(y,B,nj) {
      # tamaños
      n <- sum(nj)
      # hiperparametros
      mu0  <- 13.495 
      g20  <- 11.382
      al0  <- 1
      be0  <- 0.846 
      kp   <- 3
      # valores iniciales
      theta <- mean(y)
      sig2  <- var(y)  # sigma_j^2
      # almacenamiento
      THETA <- matrix(data = NA, nrow = B, ncol = 2)
      LL    <- matrix(data = NA, nrow = B, ncol = 1)
      # cadena
      set.seed(2023)
      for (b in 1:B) {
            # actualizar varsigma^2
            varsig2 <- 1/rgamma(n = n, shape = 0.5*(kp + 1), rate = 0.5*(kp*sig2 + (y - theta)^2))
            # actualizar theta
            vtheta <- 1/(1/g20 + sum(1/varsig2))
            theta  <- rnorm(n = 1, mean = vtheta*(mu0/g20 + sum(y/varsig2)), sd = sqrt(vtheta))
            # actualizar sigma^2
            sig2 <- rgamma(n = 1, shape = 0.5*((kp*n)+ al0), rate = 0.5*((kp*sum(1/varsig2))+ be0))
            # almacenar
            THETA[b,] <- c(theta, sig2)
            # log-verosimilitud
            LL[b] <- sum(dt.scaled(x = y, df = kp,mean = theta, sd = sqrt(sig2),log = T))
      }
      # fin de la cadena
      # salida
      THETA <- as.data.frame(THETA[-c(1:1000),])
      LL    <- as.data.frame(LL[-c(1:1000),])
      colnames(THETA) <- c("theta", "sig2")
      colnames(LL) <- c("ll")
      return(list(THETA = THETA, LL = LL))
}
