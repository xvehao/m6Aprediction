mylasso <- function(X, y, lambda=0.1, tol=1e-6, iter=100){
  n = nrow(X)
  p = ncol(X)
  w = solve(crossprod(X) + diag(lambda, ncol(X))) %*% crossprod(X,y)
  w0 = mean(y) - colMeans(X)%*%w
  tol_curr = Inf
  X2 = X^2
  X2_cm = colMeans(X2)
  i=1
  while(tol < tol_curr && i < iter){
    w_old = w
    for (j in 1:p){
      if (w[j] == 0){
        next
      }
      y_pred = X%*%w
      rho = (1/n)*t(X[,j])%*%(y - matrix(w0, length(y),1) - y_pred + w[j]*X[,j])
      w[j] = soft_thresh(rho, lambda)/X2_cm[j]
    }
    w0 = mean(y) - colMeans(X)%*%w
    tol_curr = crossprod(w - w_old)
    i = i + 1
  }
  rbind(w0,w)
}

soft_thresh <- function(a, lambda) {
  if (a > lambda){
    out = a - lambda
  }
  else if (a < -lambda){
    out = a + lambda
  }
  else{
    out = 0
  }
}

myGIREN <- function(X, y, A, lambda=0.1, gamma = 0.03, alpha=1, tol=1e-6, iter=100){
  n = nrow(X)
  p = ncol(X)
  w = solve(crossprod(X) + diag(lambda, ncol(X))) %*% crossprod(X,y)
  w0 = mean(y) - colMeans(X)%*%w
  tol_curr = Inf
  X2 = X^2
  X2_cm = colMeans(X2)
  A_cm = colMeans(A)
  i = 1
  while(tol < tol_curr && i < iter && (sum(abs(w))+abs(w0))<1e6){
    w_old = w
    for (j in 1:p){
      if (w[j] == 0){
        next
      }
      y_pred = X%*%w
      rho = (1/n)*t(X[,j])%*%(y - matrix(w0, length(y),1) - y_pred + w[j]*X[,j]) + 2*gamma*(crossprod(A[,j],w)-A[j,j]*w[j])
      w[j] = soft_thresh(rho, lambda*alpha)/(X2_cm[j]+lambda*(1-alpha)-2*gamma*A[j,j])
    }
    w0 = mean(y) - colMeans(X)%*%w
    tol_curr = crossprod(w - w_old)
    i = i + 1
  }
  rbind(w0,w)
}

