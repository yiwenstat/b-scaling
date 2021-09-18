##############################
## input: coef
##        ith data point
##        knots interval
## output: evaluated value
##############################
NT <- function(coef, t, t.knots.seq) {
  seq.pit = NULL
  for( ite.interval in 2:length(t.knots.seq) ){
    if( (t < t.knots.seq[ite.interval] && t>= t.knots.seq[ite.interval-1]) ) {
      seq.pit <- c(seq.pit, c(1, t, t^2, t^3))
    } else{
      seq.pit <- c(seq.pit, c(0,0,0,0))
    }
  }
  return(seq.pit%*%coef)
}

##############################
## input: number of knots k0
##        n by K data matrix
## output: basis matrix: n by K*l
##############################
generate.basis <- function(data, k0, newdata=NULL){
  n = dim(data)[1]
  K = dim(data)[2]

  ## cubic spline
  m = 3
  num.constr = (k0-2)*m
  num.para = (m+1)*(k0-1)
  # knots.seq = apply(data, 2, function(t){seq(min(t), max(t), length=k0)})
  knots.seq = apply(data, 2, function(t){quantile(t, prob=seq(0, 1, length=k0))})
  ## generate basis
  q1=function(t){ return( c(1, t, t^2, t^3, -1, -t, -t^2, -t^3) ) }
  q2=function(t){ return( c(0, 1, 2*t, 3*t^2, 0, -1, -2*t, -3*t^2) ) }
  q3=function(t){ return( c(0, 0, 2, 6*t, 0, 0, -2, -6*t) ) }

  ## format U matrix
  U = matrix( 0, nrow = K*num.constr, ncol = K*num.para )

  for(ite.K in 1:K){

    con.knots <- knots.seq[2:(k0-1),ite.K]

    id.row = c( ((ite.K-1)*num.constr+1) : (ite.K*num.constr) )
    id.col = c( ((ite.K-1)*num.para+1) : (ite.K*num.para) )

    tempU = NULL
    for(ite.knot in 1:length(con.knots)){
      id.start = (ite.knot-1)*(m+1)+1
      seqcon1 <- seqcon2 <- seqcon3 <- rep(0, num.para)
      seqcon1[id.start:(id.start+(m+1)*2-1)] = q1( con.knots[ite.knot] )
      seqcon2[id.start:(id.start+(m+1)*2-1)] = q2( con.knots[ite.knot] )
      seqcon3[id.start:(id.start+(m+1)*2-1)] = q3( con.knots[ite.knot] )

      # print(id.start)
      tempU <- rbind(tempU, seqcon1, seqcon2, seqcon3)
    }

    U[id.row, id.col] = tempU
  }
  ## coeffcient of spline basis
  eigU <- eigen(t(U)%*%U)
  c.basis <- eigU$vec[, (K*num.constr+1):(K*num.para) ]

  ## evaluate the basis at those data points
  if(is.null(newdata)) newdata <- data
  p.basis <- NULL
  for(ite.K in 1:K){
    id.col = c( ((ite.K-1)*num.para+1) : (ite.K*num.para) )
    solution.id = apply( c.basis[id.col,], 2, function(t){ sum(t!=0) > 0 } )

    temp.p.basis = NULL
    for(ite.n in 1:dim(newdata)[1]){
      temp.p.basis <- rbind(temp.p.basis, NT(c.basis[id.col, solution.id], newdata[ite.n, ite.K], knots.seq[,ite.K]))
    }
    p.basis <- cbind(p.basis, temp.p.basis )
  }
  return(p.basis)
}

bscaling = function(data, num.knots=3){
    ## sample size
    n = dim(data)[1]
    ## measurements
    K = dim(data)[2]
    ## the number of knots
    k0 <- num.knots
    m <- 3
    l <- num.basis <- k0+m-1

    ## initialize basis matrix
    basis <- generate.basis(data, k0)

    ## Sigma matrix
    sigma <- cov(basis)/(K^2)

    ## eigenvalue
    # eigcovbasis <- eigen(sigma)
    # eigval = eigcovbasis$val
    # id = K*l

    svdbasis = svd( scale(basis/K,scale=FALSE) )
    eigvec = svdbasis$v
    eigval = (svdbasis$d)^2/n
    # ## numerical adjust eigval
    # if( sum(eigval<=.Machine$double.eps)>=1 ){
    #     rep.id <- (eigval<=.Machine$double.eps)
    #     eigval[rep.id] <- eigval[ K*l-sum(rep.id) ]
    #     eigval[rep.id] <- .Machine$double.eps
    # }

    id = K*l

    ## sigma -1/2
    sigma_inv <- eigvec[,1:id] %*% diag( (eigval[1:id])^{-1/2} ) %*% t(eigvec[,1:id])
    if(!isSymmetric(sigma_inv)){
        sigma_inv = (sigma_inv + t(sigma_inv))/2
    }

    ## Delta matrix
    Delta = matrix(0, nrow = K*num.basis, ncol = K*num.basis)
    for( i in 1:n){
        G = matrix(0, nrow = K*num.basis, ncol = K)
        for(j in 1:K){
          G[ ((j-1)*num.basis+1):(j*num.basis), j] = basis[i, (num.basis*(j-1)+1):(num.basis*j) ]
        }
        one = matrix(1, nrow = K, 1)
        Q = diag(1, K) - one %*% t(one)/K
        Delta <- (G %*% Q %*% t(G)) + Delta
    }
    Delta = Delta/n
    if(!isSymmetric(Delta)){
        Delta = (Delta + t(Delta))/2
    }
    ## obj.mat = sigma_inv %*% Delta %*% sigma_inv
    ## if(!isSymmetric(obj.mat)){
    ##    obj.mat = (obj.mat + t(obj.mat))/2
    ## }
    ##Delta=-Delta
    ############################
    # mat.b = eigen(obj.mat)
    # b = mat.b$vec[,id]
    ############################
    svdDelta = svd(Delta)
    Delta_half = svdDelta$u %*% diag(sqrt(svdDelta$d)) %*% t(svdDelta$v)
    svd_objmat_half = svd(sigma_inv %*% Delta_half)
    # svd_objmat_half$val
    b = svd_objmat_half$u[,id]

    a = sigma_inv %*% b
    est = (basis %*% a)/K

    ## b-variance
    NAmat=t(apply(basis, 1, function(t){t*a}))
    mybreaks=ceiling(seq(1,length(a), length=K+1))
    fwmat=NULL
    for(i in 1:(length(mybreaks)-1)){
    	fwmat=cbind(fwmat, rowSums(NAmat[,(mybreaks[i]):(mybreaks[i]+length(a)/K-1)]))
    }
    bvar=apply(fwmat,1,function(t){var(t)*(length(t)-1)/length(t) })
    
    return( list(Delta, sigma, basis, b, a, est, bvar, svd_objmat_half$d) )
}

my.normalize <- function(mydata){
  mydata=as.matrix(mydata)
  for(j in 1:(dim(mydata)[2])) {
    mydata[,j] <- mydata[,j]-min(mydata[,j])
    mydata[,j] <- mydata[,j]/max(mydata[,j])
  }
  return(mydata)
}
