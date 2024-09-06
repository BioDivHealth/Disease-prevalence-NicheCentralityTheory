library(adehabitatHS)# Run the ENFA analysis

#









#
#
#
#
#
data(lynxjura)

map <- lynxjura$map

## We keep only "wild" indices.
locs <- lynxjura$locs
locs <- locs[slot(locs, "data")[,2]!="D",]

hist(map, type = "l")
## The variable artif is far from symetric

## We perform a square root transformation
## of this variable
## We therefore normalize the variable 'artif'
slot(map,"data")[,4] <- sqrt(slot(map,"data")[,4])
hist(map, type = "l")

## We prepare the data for the ENFA
tab <- slot(map, "data")
pr <- slot(count.points(locs, map), "data")[,1]

## We then perform the PCA before the ENFA
pc <- dudi.pca(tab, scannf = FALSE)

## The object 'pc' contains the transformed table (i.e.
## centered so that all columns have a mean of 0
## and scaled so that all columns have a variance of 1
## 'pc' also contains the weights of the habitat variables,
## and the weights of the pixels in the analysis

# The FULL ENFA FUNCTION to addapt to make it work with spat raster data


#######################################################################
#######################################################################
#######                                                          ######
#######                                                          ######
#######                  Generalized ENFA                        ######
#######                                                          ######
#######                                                          ######
#######################################################################
#######################################################################


enfa <- function(dudi, pr, scannf = TRUE, nf = 1)
{
  ## Verifications
  if (!inherits(dudi, "dudi"))
    stop("object of class dudi expected")
  call <- match.call()
  if (any(is.na(dudi$tab)))
    stop("na entries in table")
  if (!is.vector(pr))
    stop("pr should be a vector")
  
  ## Bases of the function
  prb <- pr
  pr <- pr/sum(pr)
  row.w <- dudi$lw/sum(dudi$lw)
  col.w <- dudi$cw
  Z <- as.matrix(dudi$tab)
  n <- nrow(Z)
  f1 <- function(v) sum(v * row.w)
  center <- apply(Z, 2, f1)
  Z <- sweep(Z, 2, center)
  
  
  ## multiply with the square root of the column weights
  Ze <- sweep(Z, 2, sqrt(col.w), "*")
  
  ## Inertia matrices S and G
  DpZ <- apply(Ze, 2, function(x) x*pr)
  
  ## Marginality computation
  mar <- apply(Z,2,function(x) sum(x*pr))
  me <- mar*sqrt(col.w)
  Se <- crossprod(Ze, DpZ)
  Ge <- crossprod(Ze, apply(Ze,2,function(x) x*row.w))
  
  ## Computation of S^(-1/2)
  eS <- eigen(Se)
  kee <- (eS$values > 1e-9)     ## keep only the positive values
  S12 <- eS$vectors[,kee] %*% diag(eS$values[kee]^(-0.5)) %*% t(eS$vectors[,kee])
  
  ## Passage to the third problem
  W <- S12 %*% Ge %*% S12
  x <- S12%*%me
  b <- x / sqrt(sum(x^2))
  
  ## Eigenstructure of H
  H <- (diag(ncol(Ze)) - b%*%t(b)) %*% W %*% (diag(ncol(Ze)) - b%*%t(b))
  s <- eigen(H)$values[-ncol(Z)]
  
  ## Number of eigenvalues
  if (scannf) {
    barplot(s)
    cat("Select the number of specialization axes: ")
    nf <- as.integer(readLines(n = 1))
  }
  if (nf <= 0 | nf > (ncol(Ze) - 1))
    nf <- 1
  
  ## coordinates of the columns on the specialization axes
  co <- matrix(nrow = ncol(Z), ncol = nf + 1)
  tt <- data.frame((S12 %*% eigen(H)$vectors)[, 1:nf])
  ww <- apply(tt, 2, function(x) x/sqrt(col.w))
  norw <- sqrt(diag(t(as.matrix(tt))%*%as.matrix(tt)))
  co[, 2:(nf + 1)] <- sweep(ww, 2, norw, "/")
  
  ## coordinates of the columns on the marginality axis
  m <- me/sqrt(col.w)
  co[, 1] <- m/sqrt(sum(m^2))
  
  ## marginality
  m <- sum(m^2)
  
  ## Coordinates of the rows on these axes
  li <- Z %*% apply(co, 2, function(x) x*col.w)
  
  ## Output
  co <- as.data.frame(co)
  li <- as.data.frame(li)
  names(co) <- c("Mar", paste("Spe", (1:nf), sep = ""))
  row.names(co) <- dimnames(dudi$tab)[[2]]
  names(li) <- c("Mar", paste("Spe", (1:nf), sep = ""))
  enfa <- list(call = call, tab = data.frame(Z), pr = prb, cw = col.w,
               nf = nf, m = m, s = s, lw = row.w, li = li,
               co = co, mar = mar)
  class(enfa) <- "enfa"
  return(invisible(enfa))
}








