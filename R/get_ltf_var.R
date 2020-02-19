get_ltf_var_inner <- function(vart,icc,rho.c,rho.s,J,K,L1,L2,G1,G2){
  varct <- vart*icc*(1-rho.c)
  vars <- vart*(1-icc)*rho.s
  varst <- vart*(1-icc)*(1-rho.s)
  lambda1 <- L1/K
  lambda2 <- L2/K
  gamma1 <- G1/K
  gamma2 <- G2/K
  eta <- (1-lambda1)/(1-lambda1+gamma1) + (1-lambda2)/(1-lambda2+gamma2) - 1
  rho.s.star <- rho.s - 0.25*(1/(1-lambda1+gamma1) + 1/(1-lambda2+gamma2) - 2*(eta*vars + varst)/(vars+varst) )
  4*(varct/J + (1-rho.s.star)*(vars+varst)/(J*K))
}

get_ltf_var <- function(vart,icc,rho.c,rho.s,J,K,L1,L2,G1,G2){
  ltfnr <- get_ltf_var_inner(vart,icc,rho.c,rho.s,J,K,L1,L2,G1=0,G2=0)
  cohort <- get_ltf_var_inner(vart,icc,rho.c,rho.s,J,K,L1=0,L2=0,G1=0,G2=0)
  xsec <- get_ltf_var_inner(vart,icc,rho.c,rho.s,J,K,L1=K,L2=K,G1=K,G2=K)
  mix <- get_ltf_var_inner(vart,icc,rho.c,rho.s,J,K,L1,L2,G1=L1,G2=L2)
  Kreduced <- K - max(c(L1,L2))
  reduced <- get_ltf_var_inner(vart,icc,rho.c,rho.s,J,K=Kreduced,L1=0,L2=0,G1=0,G2=0)
  varlist <- list(ltfnr=ltfnr,
                  mix=mix,
                  reduced=reduced,
                  xsec=xsec,
                  cohort=cohort)
  as.data.frame(t(unlist(varlist)))
}

get_var <- function(vart,icc,rho.c,rho.s,J,K,L1,L2,G1,G2){
  tempdf <- expand.grid(vart=vart,icc=icc,rho.c=rho.c,rho.s=rho.s,J=J,K=K,L1=L1,L2=L2,G1=G1,G2=G2)
  vardf <- pmap_dfr(tempdf, get_ltf_var)
  cbind(tempdf, vardf) %>%
    mutate(K1 = K - L1 + G1, K2 = K - L2 + G2)
}


get_did_var <- function(varc,varct,vars,varst,J,K,ltf.0,ltf.1,gtf.0,gtf.1){
  lambda1 <- ltf.0/K
  lambda2 <- ltf.1/K
  gamma1 <- gtf.0/K
  gamma2 <- gtf.1/K
  eta <- (1-lambda1)/(1-lambda1+gamma1) + (1-lambda2)/(1-lambda2+gamma2) - 1
  rho.s.star <- rho.s - 0.25*(1/(1-lambda1+gamma1) + 1/(1-lambda2+gamma2) - 2*(eta*vars + varst)/(vars+varst) )
  4*(varct/J + (1-rho.s.star)*(vars+varst)/(J*K))
}

