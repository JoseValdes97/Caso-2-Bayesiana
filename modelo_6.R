# Sexto modelo
#y_j <- aggregate(1/y ~ Data$Dominio, data = Data, sum)
MCMC_6 <- function(y, B, nj, yb){
      # tamaños
      n <- sum(nj)
      m <- length(nj)
      # hiperparametros
      mu0  <- 13.495 
      g20  <- 11.382
      eta0 <- 1  
      t20  <- 1.182
      al   <- 1  
      ab   <- 1
      bb   <- 1.182 
      k    <- 3
      # valores iniciales
      theta <- yb
      set.seed(2023)
      beta  <- rgamma(1,shape = ab*0.5,rate = bb*0.5)
      sig2  <- rgamma(n=m,shape= al/2,rate = beta/2) #sigma_j^2
      mu    <- mean(theta)
      tau2  <- var(theta)
      # almacenamiento
      THETA <- matrix(data = NA, nrow = B, ncol = 2*m+3)
      LL    <- matrix(data = NA, nrow = B, ncol = 1)
      # auxiliares sumas paralos parametros
      lims <- cumsum(nj)
      set.seed(2023)
      for (b in 1:B) {
            # Generar varsigma
            varsig2   <- 1/rgamma(n = n, shape = 0.5*(k+1), rate = 0.5*(k*rep(sig2,nj) + (y - rep(theta, nj))^2))
            # sumas para la vtheta
            acumu <- cumsum(1/varsig2)
            sum1 <- acumu[lims] - c(0,acumu[lims[-length(lims)]])
            #vtheta
            vtheta <- 1/(1/tau2 + sum1)
            # sumas para mediatheta
            acumu <- cumsum(y/varsig2)
            sum2 <- acumu[lims] - c(0,acumu[lims[-length(lims)]])
            # actualizar theta
            theta     <- rnorm(n = m, mean = vtheta*(mu/tau2 + sum2), sd = sqrt(vtheta))
            # actualizar mu
            vmu       <- 1/(1/g20 + m/tau2)
            mu        <- rnorm(n = 1, mean = vmu*((mu0/g20)+(m*mean(theta)/tau2)), sd=sqrt(vmu))
            # actualizar tau2
            tau2      <- 1/rgamma(n = 1, shape = 0.5*(eta0 + m), rate = 0.5*(eta0*t20+(m-1)*var(theta)+m*(mean(theta)-mu)^2))
            # actualizar sigma_j^2
            sig2      <- rgamma(n = m, shape = 0.5*(k*nj+al), rate = (0.5*(k*sum1 + beta)))
            # crear beta y actualizar
            beta      <- rgamma(n = 1, shape=0.5*(ab+al*m),rate = 0.5*(sum(sig2)+bb))
            # almacenar
            THETA[b,] <- c(theta, sig2, mu, tau2, beta)
            # log-verosimilitud
            LL[b]     <- sum(dt.scaled(x = y, df = k, mean = rep(theta, nj), sd = sqrt(rep(sig2, nj)),log = T))
      }
      
      # fin de la cadena
      # salida
      THETA <- as.data.frame(THETA[-c(1:1000),])
      LL    <- as.data.frame(LL[-c(1:1000)])
      colnames(THETA) <- c(paste0("theta", 1:m), paste0("sig2", 1:m),"mu", "tau2","beta")
      colnames(LL) <- c("ll")
      return(list(THETA = THETA, LL = LL))
}
