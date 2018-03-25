

# Load packages -----------------------------------------------------------

library(propagate)
library(psych)
library(elasticnet)
library(rospca)
library(PMA)

# Créer différents scenarios -----------------------------------------------

# Structure en blocks
vectors_blocks <- matrix(c(0.096, 0.082, 0.08, 0.594, 0.584, 0.533, -0.537, -0.565, -0.608, 0.085, 0.096, 0.074, 0.759, -0.599, -0.119, -0.074, -0.114, 0.18, -0.12, 0.231, -0.119, -0.308, -0.418, 0.805, 0.335, 0.511, -0.771, 0.069, 0.052, -0.157, -0.021, -0.013, 0.016, 0.731, -0.678, -0.069), nrow = 6, ncol = 6)
valeurs_blocks <- c(1.8367, 1.64, 0.751, 0.659, 0.607, 0.506)

# Structure intermediaire
vectors_intermediaires <- matrix(c(0.224, 0.253, 0.227, 0.553, 0.521, 0.507, -0.509, -0.519, -0.553, 0.249, 0.254, 0.199, 0.604, -0.361, -0.246, -0.249, -0.258, 0.561, 0.297, -0.644, 0.377, -0.052, 0.451, -0.384, -0.327, -0.341, 0.608, 0.262, -0.509, 0.281, 0.361, -0.064, -0.267, 0.706, -0.367, -0.402), nrow = 6, ncol = 6)
valeurs_intermediaires <- c(1.795, 1.674, 0.796, 0.618, 0.608, 0.510)

# Structure platte
vectors_plattes <- matrix(c(-0.455, -0.439, -0.415, 0.434, 0.301, 0.385, 0.336, 0.370, 0.422, 0.458, 0.435, 0.416, -0.087, -0.212, 0.378, 0.040, -0.697, 0.563, 0.741, -0.630, -0.110, -0.136, 0.114, 0.104, -0.328, -0.445, 0.697, -0.167, 0.356, -0.234, 0.125, -0.175, 0.099, 0.744, -0.306, -0.545), nrow = 6, ncol = 6)
valeurs_plattes <- c(1.841, 1.709, 0.801, 0.649, 0.520, 0.480)


# Simuler les données -----------------------------------------------------


cov_function <- function(vecteurs_propres, valeurs_propres, cov = TRUE){
  
  n <- nrow(vecteurs_propres)
  S = vecteurs_propres
  M = matrix(rep(0, n^2), ncol = n, nrow = n)
  diag(M) <- valeurs_propres
  
  cor = round(S %*% M %*% solve(S), 6)
  
  for (j in 1:ncol(cor)) {
    for (i in 1:j) {
      if (cor[i, j] != diag(cor)[j]){
        cor[i, j] = cor[j, i]
      }
    }
  }
  
  if (cov){
    cov = round(propagate::cor2cov(C = cor, var = valeurs_propres), 6)
  } else {
    round(cor, 6)
  }
}

cov_blocks <- cov_function(vecteurs_propres = vectors_blocks, valeurs_propres = valeurs_blocks)
cov_intermediaire <- cov_function(vecteurs_propres = vectors_intermediaires, valeurs_propres = valeurs_intermediaires)
cov_platte <- cov_function(vecteurs_propres = vectors_plattes, valeurs_propres = valeurs_plattes)

data_blocks <- rmvnorm(10000, mean = rep(0, 6), sigma = cov_blocks)
data_intermediaire <- rmvnorm(10000, mean = rep(0, 6), sigma = cov_intermediaire)
data_platte <- rmvnorm(10000, mean = rep(0, 6), sigma = cov_platte)


# Appliquer les différentes méthodes --------------------------------------

acp_blocks <- prcomp(data_blocks, scale. = TRUE)
rpca_blocks <- principal(cor(data_blocks), nfactors = 4, rotate = "varimax")
spca_blocks <- spca(cor(data_intermediaire), K = 4, para = rep(3, 4), sparse = "varnum", type = "Gram")

acp_intermediaire <- prcomp(cor(data_intermediaire), center = TRUE, scale. = TRUE)
rpca_intermediaire <- principal(cor(data_intermediaire), nfactors = 4, rotate = "varimax")
spca_intermediaire <- spca(cor(data_intermediaire), K = 4, para = rep(6, 4), sparse = "varnum")

acp_platte <- prcomp(cor(data_platte), center = TRUE, scale. = TRUE)
rpca_platte <- principal(cor(data_platte), nfactors = 4, rotate = "varimax")
spca_platte <- spca(cor(data_platte), K = 4, para = rep(6, 4), sparse = "varnum")

# Analyser les performances -----------------------------------------------

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

theta <- acos( sum(a*b) / ( sqrt(sum(a * a)) * sqrt(sum(b * b)) ) )

sapply(1:4, function(i) {
  angle(rpca_blocks$loadings[,i], vectors_blocks[,i]) 
})

sapply(1:4, function(x) {
  angle(-acp_blocks$rotation[,x], vectors_blocks[,x]) 
})

sapply(1:4, function(x) {
  angle(spca_blocks$loadings[,x], vectors_blocks[,x]) 
})


sapply(1:4, function(x) {
  angle(rpca_intermediaire$loadings[,x], vectors_intermediaires[,x]) 
})

sapply(1:4, function(x) {
  angle(acp_intermediaire$rotation[,x], vectors_intermediaires[,x]) 
})

sapply(1:4, function(x) {
  angle(spca_intermediaire$loadings[,x], vectors_intermediaires[,x]) 
})

sapply(1:4, function(x) {
  angle(rpca_platte$loadings[,x], vectors_platte[,x]) 
})

sapply(1:4, function(x) {
  angle(acp_platte$rotation[,x], vectors_platte[,x]) 
})

sapply(1:4, function(x) {
  angle(spca_platte$loadings[,x], vectors_platte[,x]) 
})




