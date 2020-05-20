#Original IRLS method
# M --- trajectory matrix
# k --- rank of signal
# maxITER --- number of iterations N_IRLS
# maxiter --- number of iterations N_alpha

IRLS_orig<-function(M, k, alpha=4.685, eps=1e-5, maxITER=20, maxiter=10){
  m<-nrow(M)
  n<-ncol(M)
  initial<-svd(M) #initialization U and V using svd
  U<-initial$u[1:nrow(initial$u),1:k]
  Lambda<-initial$d[1:k]
  V<-initial$v[1:nrow(initial$v),1:k]
  ITER<-0
  iter<-0
  
  repeat {
    R<-M-U%*%t(V) #residuals matrix
    r<-as.vector(t(R))
    sigma<-1.4826*median(abs(r-median(abs(r))))
    W<-weights(R,sigma,alpha,m,n)
    
    repeat{
      # updating U using QR-decomposition
      for (i in (1:m)){
        Wi<-diag(W[i,1:ncol(W)])
        mi<-M[i,1:ncol(M)]
        QR <- qr(t(V)%*%Wi%*%V)
        Q <- qr.Q( QR )
        R <- qr.R( QR )
        beta<-backsolve(R,t(Q)%*%t(V)%*%Wi%*%mi)
        U[i,1:ncol(U)]<-beta
      }
      U<-U[1:nrow(U),1:k]
      
      # updating V using QR-decomposition
      for (j in (1:n)){
        Wj<-diag(W[1:nrow(W),j])
        mj<-M[1:nrow(M),j]
        QR <- qr(t(U)%*%Wj%*%U)
        Q <- qr.Q( QR )
        R <- qr.R( QR )
        beta<-backsolve(R,t(Q)%*%t(U)%*%Wj%*%mj)
        V[j,1:ncol(V)]<-beta
      }
      V<-V[1:nrow(V),1:k]
      
      iter<-iter+1
      if ( ((frobenius.norm(W^{1/2}*(M-U%*%t(V))))^2<eps) | (iter > maxiter) ) {break}
    }
    
    ITER<-ITER+1
    if ( ((frobenius.norm(W^{1/2}*(M-U%*%t(V))))^2<eps) | (ITER > maxITER) ) { break}
  }
  M_est<-U%*%t(V)
  return(M_est)
}