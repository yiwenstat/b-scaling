# B-scaling method

## Example

### Generate example data
```
set.seed(12345)
n=1000
K=10

input.x=sort(runif(n))
F_logit <- function(x){
   1/( 1+exp(20*(x-0.5)) )
  }

nu <- 2; unif.par <- 3; H <- 5
Zmat <- matrix(runif(H*K, -sqrt(unif.par), sqrt(unif.par)), K, H)
zeta <- sapply(1:H, function(h)(-1)^(h+1)*h^(-nu/2))
s <- runif(K,-10,10)

input.fx = NULL
for (i in 1:K){
    val <- 0
    temp.input.x <- input.x+rnorm(n,0,0.1)
    for (h in 1:H) val <- val + zeta[h] * Zmat[i,h] * F_logit(temp.input.x)
    input.fx = cbind(input.fx, val*s[i])
}

x <- input.x
wmat <- input.fx
```

### Calculate B-mean
```
source("bscaling.R")
k0=3
bsobj=bscaling(wmat, num.knots=k0)
bmean=bsobj[[6]]
cor(bmean, x) 
```

