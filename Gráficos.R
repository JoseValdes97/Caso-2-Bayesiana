# grafico  1
# Gráfico
#par(mfrow = c(1,1))
#modelos normal
yrange <- range(cadena1_M1$LL$ll, cadena1_M2$LL$ll,cadena1_M3$LL$ll,cadena1_M4$LL$ll, cadena1_M5$LL$ll,cadena1_M6$LL$ll)
plot(cadena1_M1$LL$ll, type = "p", pch = 20, cex = 0.5, col = "darkblue", ylim = yrange, xlab = "Iteracion",ylab = "Log-verosimilitud", main = "Modelos", xlim= c(0,11000))
lines(cadena1_M2$LL$ll, type = "p",pch = 20, cex = 0.5,col = "cyan")
lines(cadena1_M3$LL$ll, type = "p", pch = 20, cex = 0.5,col = "brown")
#modelos t
#yrange <- range(cadena1_M4$LL$ll, cadena1_M5$LL$ll,cadena1_M6$LL$ll)
#plot(cadena1_M4$LL$ll, type = "p", pch = 20, cex = 0.6, col = "lightblue", ylim = yrange,  xlab = "Iteración",ylab = "Log-verosimilitud", main = "Modelos t")
lines(cadena1_M4$LL$ll, type = "p", pch = 20, cex = 0.5,col = "gold")
lines(cadena1_M5$LL$ll, type = "p", pch = 20, cex = 0.5,col = "magenta")
lines(cadena1_M6$LL$ll, type = "p", pch = 20, cex = 0.5,col = "darkgreen")
# cuadro legendario
legend(locator(1),legend = c("M1","M2","M3","M4","M5","M6"),
       bty = "n",col= c("darkblue","cyan","brown","gold","magenta","darkgreen"),
       lty=c(0,0,0,0,0,0),lwd=2, pch =  17, pt.cex = 1, cex= 0.9)
#Gráfico 2
par(mfrow = c(1,2))
boxplot(Bog$Ingtot, main = "Bogotá", ylab= "Ingresos")
plot(x = NA, y = NA, ylab = "Densidad", xlab = "Ingresos totales", cex.axis = 0.7,
      xlim = range(Bog$Ingtot),ylim = c(0,0.5), main= "Histograma de Bogotá")
hist(Bog$Ingtot, freq = F, add = T, col = "mistyrose", border = "mistyrose")
lines(density(Bog$Ingtot), col = "red")
#Gráifco 3
#ranking
THETA    <- cadena1_M6$THETA
Dominios <- Estadisticos$Dominio
#Estimaciones puntuales
that <- colMeans(THETA[,1:m])
# Int. credibilidad al 95%
IC95 <- apply(X = THETA[,1:m], MARGIN = 2, FUN = function(x) quantile(x, c(0.025,0.975)))
# ordernar el ranking (menor a mayor)
orden    <- order(that)
Dominios <- Dominios[orden]
that     <- that[orden]
IC95     <- IC95[,orden]
# Colores
colo <- rep(2,m)
colo[which(IC95[1,]>13.830)] <- 1
colo[which(IC95[2,]<13.830)] <- 3
colo <- c("green3","black","red3")[colo]
#Gráfico
par(mfrow = c(1,1), mar = c(4,7,1.5,1), mgp = c(2.5,0.75,0))
plot(NA, NA, xlab = "Ingreso total", ylab = "", main = "Ranking Bayesisano: Modelo 6", xlim = c(12.8,14.3), ylim = c(1,m), cex.axis = 0.75, yaxt = "n")
axis(side = 2, at = 1:m, labels = Dominios, las = 2, cex.axis = 0.6)
abline(v = 13.830,  col = "gray", lwd = 3)
abline(h = 1:m, col = "lightgray", lwd = 1)
# intervalo y estimación puntual
for (j in 1:m) {
      segments(x0 = IC95[1,j], y0 = j, x1 = IC95[2,j], y1 = j, col = colo[j])
      lines(x = that[j], y = j, type = "p", pch = 16, cex = 0.8, col = colo[j])
}

