# modelo 5
MCMC_5 <- function(y,B,nj,yb) {
      # tamaños
      n <- sum(nj)
      m <- length(nj)
      # hiperparametros
      mu0  <- 13.495 
      g20  <- 11.382
      eta0 <- 1
      t20  <- 1.182
      al0  <- 1
      be0  <- 0.846 
      kp   <- 3
      # valores iniciales
      theta <- yb
      mu    <- mean(theta)
      tau2  <- var(theta)
      set.seed(0510)
      sig2  <- rgamma(n=1,shape= al0/2,rate = be0/2)  
      # almacenamiento
      THETA <- matrix(data = NA, nrow = B, ncol = m+3)
      LL    <- matrix(data = NA, nrow = B, ncol = 1)
      # auxiliares sumas paralos parametros
      lims <- cumsum(nj)
      # cadena
      set.seed(2023)
      for (b in 1:B) {
            # actualizar varsigma^2
            varsig2 <- 1/rgamma(n = n, shape = 0.5*(kp + 1), rate = 0.5*(kp*sig2 + (y - rep(theta, nj))^2))
            # actualizar theta
            # sumas para la vtheta
            acumu <- cumsum(1/varsig2)
            sum1 <- acumu[lims] - c(0,acumu[lims[-length(lims)]])
            #vtheta
            vtheta <- 1/(1/tau2 + sum1)
            # sumas para mediatheta
            acumu <- cumsum(y/varsig2)
            sum2 <- acumu[lims] - c(0,acumu[lims[-length(lims)]])
            # los thetas
            theta  <- rnorm(n = m, mean = vtheta*(mu/tau2 + sum2), sd = sqrt(vtheta))
            # actualizar mu
            vmu <- 1/(1/g20 + m/tau2)
            mu  <- rnorm(n = 1, mean = vmu*(mu0/g20 + m*mean(theta)/tau2), sd = sqrt(vmu))
            # actualizar tau2
            tau2 <- 1/rgamma(n = 1, shape = 0.5*(eta0+m), rate = 0.5*(eta0*t20 + (m-1)*var(theta) + m*(mean(theta) - mu)^2))
            # actualizar sigma^2
            sig2 <- rgamma(n = 1, shape = 0.5*((kp*n)+ al0), rate = 0.5*((kp*sum(1/varsig2))+ be0))
            # almacenar
            THETA[b,] <- c(theta, sig2, mu, tau2)
            # log-verosimilitud
            LL[b] <- sum(dt.scaled(x = y, df = kp,mean = rep(theta, nj), sd = sqrt(sig2),log = T))
      }
      # fin de la cadena
      # salida
      THETA <- as.data.frame(THETA[-c(1:1000),])
      LL    <- as.data.frame(LL[-c(1:1000),])
      colnames(THETA) <- c(paste0("theta", 1:m), "sig2", "mu", "tau2")
      colnames(LL) <- c("ll")
      return(list(THETA = THETA, LL = LL))
}
