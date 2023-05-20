#install.packages("dendextend")
#install.packages("data.table")
library(dendextend)
library(data.table)
library(factoextra)
library(cluster)
library(grid)
library(datawizard)

theta  <- cadena1_M6$THETA[c(1:25)]
sigma2 <- sqrt(cadena1_M6$THETA[c(26:50)])
theta_est  <- colMeans(theta)
sigma2_est <- colMeans(sigma2)
arreglo <- data.frame(theta_est=normalize(theta_est),
                      sigma2_est=normalize(sigma2_est))
rownames(arreglo) <- unique(Data$Dominio)
arreglo

# Paso 1: Calcular la matriz de distancias entre los dominios
dist_matrix <- get_dist(arreglo, method = "euclidian", stand = F)
                        
# Paso 2: Realizar la agrupación jerárquica
set.seed(2013)
hclust_result <- hclust(dist_matrix, method = "complete")
# Paso 3: Obtener la segmentación en cuatro grupos
groups <- cutree(hclust_result, k = 4)

# Paso 4: Mostrar los resultados
results <- data.frame(Grupo = as.factor(groups), arreglo)
results 

# Paso 5: Ordenar los grupos y los colores según el orden de las ramas
order <- order.dendrogram(as.dendrogram(hclust_result))
groups <- groups[order]
colores <- c("red", "blue", "green", "darkorange")
colores_grupo <- colores[groups]

# Paso 6: Graficar el dendrograma con los nombres de las hojas
plot(hclust_result, main = "Dendrograma de Segmentación", cex= 0.45,
     xlab= "Distancias entre individuos", ylab="Pesos",check = T,ann = T, hang = -1)
rect.hclust(hclust_result, k = 4, border = unique(colores_grupo))
# Paso 7: Ajustar la leyenda
legend("topright", legend = c("Grupo 1", "Grupo 2", "Grupo 3", "Grupo 4"),
       col = colores, lty = 1, lwd = 2)

# fviz_cluster(list(data = arreglo, cluster = groups))
library(ggrepel)
clases_aj <- cutree(hclust_result, k = 4)
results$Grupo <- clases_aj
results$Etiqueta <- rownames(results)
ggplot() + geom_text_repel(aes(x = theta_est, y = sigma2_est, color = Grupo, label = Etiqueta), 
                           data = results, size = 2) +
      geom_point(aes(x = theta_est, y = sigma2_est, color = Grupo), data = results, size = 0.5) +
      scale_colour_gradientn(colours= colores) +
      ggtitle('Agrupamiento Jerárquico en 4 grupos') + 
      xlab(expression(hat(theta))) + ylab(expression(hat(sigma))) 
