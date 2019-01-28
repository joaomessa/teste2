# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

hello <- function() {
  print("Hello, world!")
}

somalength<-function(x){
  func<-do.call(sumC,list(x))
  return(func+length(x))
}

maisvalor<-function(x,y,z){
  valor=do.call(signC,list(x,y))
  return(valor+z)
}

outerquad<-function(x){
  yu=do.call(rcppeigen_outerproduct,list(x))
  return(yu^2)
}

##' Fast version of Matrix :: .bdiag() -- for the case of *many*  (k x k) matrices:
##' @param lmat list(<mat1>, <mat2>, ....., <mat_N>)  where each mat_j is a  k x k 'matrix'
##' @return a sparse (N*k x N*k) matrix of class  \code{"\linkS4class{dgCMatrix}"}.
bdiag_m <- function(lmat) {
  ## Copyright (C) 2016 Martin Maechler, ETH Zurich
  if(!length(lmat)) return(new("dgCMatrix"))
  stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
            (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
            all(vapply(lmat, dim, integer(2)) == k)) # all of them
  N <- length(lmat)
  if(N * k > .Machine$integer.max)
    stop("resulting matrix too large; would be  M x M, with M=", N*k)
  M <- as.integer(N * k)
  ## result: an   M x M  matrix
  new("dgCMatrix", Dim = c(M,M),
      ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
      i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
      p = k * 0L:M,
      x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

bdiag_john <- function(lmat) {
k = length(lmat[[1]])
N = length(lmat)
M <- as.integer(N * k)
teste33=new("dgCMatrix", Dim = c(M,N),
            ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
            i = 0L:(M-1L),
            p = as.integer(k * 0L:N),
            x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
}

PWIGLSADCV <-function(vresp, listx, listz, IDCluster, wj_rep, wi_j  ){
  start<-Sys.time()

  y <- as.matrix(cbind(get(vresp)))
  nsubjc_t <- NROW(y)

  cons <- matrix(1,nsubjc_t,1)

  cluster_var<-as.matrix(cbind(IDCluster))
  cluster<-as.matrix(unique(cluster_var) )

  varlistz <-eval(parse(text=paste("cbind(",listz,")")))
  z <- as.matrix(cbind(cons,varlistz))

  varlist <-eval(parse(text=paste("cbind(",listx,")")))
  x <- as.matrix(cbind(cons,varlist))

  ##############

  mc <- cbind(y,x,z,IDCluster,wj_rep,wi_j)
  mc <- mc[order(IDCluster),]
  col_mc <- dim(mc)[2]

  #################

  nvar <- ncol(x)
  ncluster <- nrow(cluster)
  q <- ncol(z)
  s <-  ((q*(q+1))/2)+1

  y <- mc[,1]
  x <- mc[,2:(1+nvar)]
  z <- mc[,(nvar+2):(nvar+1+q)]
  wj_rep <- mc[,(col_mc-1)]
  wi_j <- mc[,col_mc]

  y<-as.matrix(y)
  x<-as.matrix(x)
  z<-as.matrix(z)
  wj_rep<-matrix(wj_rep)
  wi_j<-matrix(wi_j)


  #((I(s)[vec(makesymmetric(invvech(1::s-1))), ]))'
  h_matrix<- t(diag(s)[ks::invvech(1:(s-1)),])
  name1 <- matrix(cbind("intercept",t(unlist(noquote(strsplit(listx, ","))))))
  #name3 <- matrix(cbind("sigma2_u0","sigma_u10","sigma2_u1" ,"sigma2_e"))
  #name3 <- matrix( paste("teta_", 0:(s-1), sep=""))
  name3 <- matrix(paste ("sigma_", rep((0:(q-1)), ((q):1)), ((sequence(q:1))-1)+ (rep((0:(q-1)), ((q):1))), sep=""))
  name3<-rbind(name3, "sigma2_e")
  lambidaj <- matrix(0,nsubjc_t,1)

  info_j <-as.matrix(cbind( (tapply(cluster_var, cluster_var, function(x) NROW(x))), cumsum(tapply(cluster_var, cluster_var, function(x) NROW(x))) ))
  info_j <- cbind(info_j[,2]-info_j[,1] +1, info_j[,2])

  cluster_wgt<- matrix(0,ncluster,1)

  #-------------- Calculating the Scaled Weights --------------
  for (i in seq(along=cluster)){

    a<-info_j[i,1]
    b<-info_j[i,2]

    nsubjc<- cluster_var[ (a:b),]

    np <- NROW(nsubjc)

    k<-mean(wi_j[(a:b),])

    parte <- matrix(k,np,1)

    lambidaj[a:b,]<-parte

    wj_rep_j=wj_rep[(a:b),]

    cluster_wgt[i,]<-unique(wj_rep_j)
  }

  wi_j_star <- wi_j / lambidaj

  wj_star = cluster_wgt / mean(cluster_wgt)

  wjesc_r = wj_rep / mean(cluster_wgt)



  # Voltar com a quest?o do rescalonamento do peso
  # Retirei a parte re escalona os pesos
  # (wi_j_star <- (matrix(1,nrow(x),1)))
  # (wj_star <- (matrix(1,ncluster,1)))
  # (wjesc_r <- (matrix(1,nrow(x),1)))




  #Calcula Beta zero para cada J
  mat_t1j = matrix(0,ncluster,nvar^2)
  mat_t3j = matrix(0,ncluster,nvar)


  for (i in seq(along=cluster)){
    nsubjc=cluster_var[(info_j[i,1]:info_j[i,2]),]
    yj=y[(info_j[i,1]:info_j[i,2]),]
    xj=x[(info_j[i,1]:info_j[i,2]),]
    np = NROW(nsubjc)

    #peso
    wi_j_starj=wi_j_star[(info_j[i,1]:info_j[i,2]),]
    diag_ = diag(wi_j_starj)
    wj = wj_star[i,]

    mat_t1j[i,] = t(as.vector(crossprod(xj)))
    mat_t3j[i,] = t(as.vector(crossprod(yj,xj)))
  }


  somat1 <- matrix(colSums(mat_t1j),nvar)
  somat3 <- colSums(mat_t3j)
  beta0 = solve(somat1) %*% t(t(somat3))

  # ---------------- Parte aleat?ria
  vec_t6 = matrix(0,ncluster,1)
  vec_aux = matrix(0,ncluster,1)
  #teta0 =  matrix(1:s,s,1) #Aux
  teta0 =  matrix(rep(c(0.5,0),s),s,1)

  for (i in seq(along=cluster)){
    nsubjc=cluster_var[(info_j[i,1]:info_j[i,2]),]

    yj=y[(info_j[i,1]:info_j[i,2]),]
    xj=x[(info_j[i,1]:info_j[i,2]),]
    np = NROW(nsubjc)

    #peso
    wi_j_starj=wi_j_star[(info_j[i,1]:info_j[i,2]),]
    wj = wj_star[i,]


    #res?duos MQO
    resid = yj - xj %*% beta0
    ebarj = colSums(resid)/np
    vec_t6[i,] = colSums((resid - ebarj) ^ 2)
    vec_aux[i,]= np - 1
  }

  teta0[s,]=sum(vec_t6)/sum(vec_aux)


  #--------------------------------------------------------------------------
  #                     IGLS - ITERATIVO
  #--------------------------------------------------------------------------

  itera = 1

  H = h_matrix

  ## Truques
  beta_ant = beta0
  beta = beta_ant*2
  teta_ant=teta0
  teta=teta_ant * 2

  ss=s*s
  sncluster = s*ncluster
  ssncluster = s*sncluster
  sq=s*q
  qncluster = q*ncluster
  xxncluster = nvar^2*ncluster

  T5_res=list()
  Q1=list()
  T1=list()
  T2_nw=list()
  T2_ww=list()
  T3=NULL
  T4=NULL
  T5=list()
  T5_TETA=list()
  Nj=NULL
  Yj=list()
  Xj=list()
  T5_S_trans=list()
  NP=NULL
  for (i in seq(along=cluster)){
    #para cada j calculamos tudo isso
    yj=y[(info_j[i,1]:info_j[i,2]),]
    xj=x[(info_j[i,1]:info_j[i,2]),]
    zj=z[(info_j[i,1]:info_j[i,2]),]

    #peso
    v= (wi_j_star[(info_j[i,1]:info_j[i,2]),])
    diag_ = diag(v)
    wj = wj_star[i,]

    np = NROW(yj)

    if(np != 1){
    q1j <- crossprod(xj,diag_)
    t5j_S <- crossprod(zj,diag_)
    Xj[[i]] <- xj
    } else {
    q1j <- tcrossprod(xj,diag_)
    t5j_S <- tcrossprod(zj,diag_)
    Xj[[i]] <- t(xj)
    }
    t1j <- q1j %*% xj
    t2j <- q1j %*% zj
    t3j <- q1j %*% yj
    t4j <- t5j_S %*% yj
    t5j <- t5j_S %*% zj

    T5_res[[i]] <- t(zj) #Parte dos resíduos
    Q1[[i]] <- q1j       #Parte das variâncias
    T1[[i]] <- wj*t1j
    T2_nw[[i]] <- t(t2j)
    T2_ww[[i]] <- wj*t2j
    T3 <- c(T3,wj*t3j)
    T4 <- c(T4,t4j)
    T5_S_trans[[i]] <- t5j_S
    T5[[i]] <- t5j
    T5_TETA=c(T5_TETA,rep(list(t5j),s))
    Nj=c(Nj,rep(np,ss))
    Yj[[i]] <- yj
    NP=c(NP,rep(1 - 1/np, np)) #Parte dos resíduos

  }
  T5_res <- bdiag(T5_res)
  Q1 <- bdiag(Q1)
  T1 <- bdiag_m(T1)           #Seria melhor utilizar bdiag?
  T2_nw <- bdiag_m(T2_nw)
  T2_ww <- bdiag_m(T2_ww)
  T5 <- bdiag_m(T5)
  T5_TETA = bdiag_m(T5_TETA)
  Nj=Matrix::Diagonal(x = Nj)
  Yj <- Matrix::bdiag(Yj)
  Xj <- Matrix::bdiag(Xj)
  T5_S_trans <-Matrix::bdiag(T5_S_trans)


  T5_TETA_2=vector_eigen(T5_TETA, sncluster*q^2)
  T5_TETA_2=split(T5_TETA_2,rep(1:sncluster,each=q^2))
  T5_TETA_2=lapply(T5_TETA_2, function(i) matrix(i,q))
  T5_TETA_2=rep(T5_TETA_2,each=s)
  T5_TETA_2=bdiag_m(T5_TETA_2)

  #q = ncol(z)
  #s = ((q*(q+1))/2)+1
  DELTA=rep(c(rep(0,sq-q),rep(1,q)),ncluster)
  DELTA = Matrix::Diagonal(x = DELTA)
  DELTA = Matrix::drop0(DELTA, tol = 0, is.Csparse = NA)

  DELTA_Rl=rep(c(rep(0,s-1),1),sncluster)
  DELTA_Rl = Matrix::Diagonal(x = DELTA_Rl)
  DELTA_Rl = Matrix::drop0(DELTA_Rl, tol = 0, is.Csparse = NA)
  DELTA_Rl=as(DELTA_Rl,"dgCMatrix")

  DELTA_Rk=rep(c(rep(0,ss-s),rep(1,s)),ncluster)
  DELTA_Rk = Matrix::Diagonal(x = DELTA_Rk)
  DELTA_Rk = Matrix::drop0(DELTA_Rk, tol = 0, is.Csparse = NA)
  DELTA_Rk=as(DELTA_Rk,"dgCMatrix")

  DELTA_Rk_2 = rep(c(rep(0,s-1),1),ncluster)
  DELTA_Rk_2 = Matrix::Diagonal(x = DELTA_Rk_2)
  DELTA_Rk_2 = Matrix::drop0(DELTA_Rk_2, tol = 0, is.Csparse = NA)
  DELTA_Rk_2 = as(DELTA_Rk_2,"dgCMatrix")

  H_kj=list()
  for (k in 1:s) {
    H_kj[[k]] <- Matrix::Matrix(H[k,],q,sparse = T)
  }
  H_KJ = Matrix::bdiag(H_kj)
  H_kj=rep(list(H_KJ),ncluster)
  H_kj=Matrix::bdiag(H_kj)
  H_kj <- as(H_kj, "dgCMatrix")
  H_lj=rep(list(H_KJ),sncluster)
  H_lj=Matrix::bdiag(H_lj)
  H_lj <- as(H_lj, "dgCMatrix")

  wj_star_1=as.numeric(wj_star)
  wj_star_1_quad = wj_star_1^2
  wj_star_2=rep(wj_star_1,each=ss)
  wj_star_2=Matrix::Diagonal(x = wj_star_2)
  wj_star_3=rep(wj_star_1,each=s)
  wj_star_3=Matrix::Diagonal(x = wj_star_3)
  wj_star_4=rep(wj_star_1_quad,each=nvar)
  wj_star_4=Matrix::Diagonal(x = wj_star_4)
  wj_star_5=rep(wj_star_1_quad,each=s)
  wj_star_5=Matrix::Diagonal(x = wj_star_5)

  DIAG = Matrix::Diagonal(x=wi_j_star)
  R_1 = DELTA_Rk%*%DELTA_Rl%*%Nj

  while (itera<= 200  & (any(abs((teta-teta_ant))> 0.000001) | any(abs((beta-beta_ant))> 0.000001) )){

    if (itera != 1) {
      sol_invvech = teta[s,]*(solve(ks::invvech(teta[1:(s-1),])))
    } else {
      sol_invvech = teta0[s,]*(solve(ks::invvech(teta0[1:(s-1),])))
    }
    SOL_BETA=rep(list(sol_invvech),ncluster)
    SOL_BETA=bdiag_m(SOL_BETA)

    #Matrix::Diagonal(dim(T5)[1])
    AJ=Matrix::solve(T5+SOL_BETA, sparse = T) #Eigen

    T0 = T2_ww %*% AJ

    matp = T1 - T0 %*% T2_nw

    matq = T3 - T0 %*% T4

    matp = vector_eigen(matp, ncluster*nvar^2)
    s_matp = matrix(colSums(matrix(matp,ncluster,nvar^2,byrow = T)),nvar)

    s_matq = colSums(matrix(matq,ncluster,nvar,byrow = T))

    if (itera != 1){
      beta_ant = beta
    }

    beta = solve(s_matp) %*% (s_matq)

    #--------- looping for teta=inv(R)x S ------------------

    SOL_TETA=rep(list(sol_invvech), sncluster)
    SOL_TETA=bdiag_m(SOL_TETA)

    AJ_TETA=vector_eigen(AJ, q*qncluster)
    AJ_TETA=split(AJ_TETA,rep(1:ncluster,each=q^2))
    AJ_TETA=lapply(AJ_TETA, function(i) matrix(i,q))
    AJ_TETA=rep(AJ_TETA,each=s)
    AJ_TETA=bdiag_m(AJ_TETA)

    DELTA_AJ = DELTA%*%AJ_TETA
    B_1 = AJ_TETA%*%SOL_TETA%*%H_kj - DELTA_AJ
    C_1 = - DELTA_AJ + B_1 - B_1%*%T5_TETA%*%AJ_TETA #simetrica (supomos)

    #Evitando a criação de 'C_2'
    R_2 = T5_TETA%*%C_1

    R_2 = vector_eigen(R_2, sncluster*q^2)
    R_2 = split(R_2,rep(1:sncluster,each=q^2))
    R_2 = lapply(R_2, function(i) matrix(i,q))
    R_2 = rep(R_2,each=s)
    R_2 = bdiag_m(R_2)

    R_3=T5_TETA_2%*%H_lj #Não consegui melhorar (expandindo)
    R_4=R_2%*%R_3

    R_21=diag_eigen(R_2)
    R_21=matrix(R_21, nrow=q)
    R_21=Matrix::Diagonal(x=colSums(R_21))

    R_31=diag_eigen(R_3)
    R_32=Matrix::Matrix(R_31, nrow=q, sparse = T) #matriz esparsa?
    R_33=Matrix::Diagonal(x=Matrix::colSums(R_32))

    R_41=diag_eigen(R_4)
    R_42=Matrix::Matrix(R_41, nrow=q, sparse = T)
    R_43=Matrix::Diagonal(x=Matrix::colSums(R_42))

    R_ = wj_star_2%*%(R_1 + DELTA_Rl%*%R_21 + DELTA_Rk%*%R_33 + R_43) #Posso multiplicar sem precisar de Matriz
    R_ = diag_eigen(R_)
    R = matrix(R_, ncluster, ss, byrow=T)

    ###################

    BETAJ=rep(list(beta),ncluster)
    BETAJ=bdiag_john(BETAJ)

    E_IJ = Yj-Xj%*%BETAJ

    S_1 = Matrix::t(E_IJ)%*%DIAG%*%E_IJ
    S_1 = diag_eigen(S_1)
    S_1 = rep(S_1,each = s)
    S_1 = Matrix::Diagonal(x = S_1)

    S_2_1 = T5_S_trans%*%E_IJ
    S_2 = vector_eigen(S_2_1, qncluster)
    S_2 = split(S_2,rep(1:ncluster,each = q))
    S_2 = rep(S_2,each=s)
    S_2 = bdiag_john(S_2)

    S_ = wj_star_3%*%(DELTA_Rk_2%*%S_1 + Matrix::t(S_2)%*%C_1%*%S_2)
    S_ = diag_eigen(S_)
    S = matrix(S_, ncluster, s, byrow=T)

    #########################
    matr=colSums(R)
    mats=colSums(S)

    r_mat = matrix(matr,s)
    s_mat = matrix(mats,s)

    #-----teta ----

    if (itera != 1) {
      teta_ant = teta
    }

    teta = solve(r_mat) %*% s_mat


    #------End of iterative process----------
    itera = itera + 1
  }

  # Number of Iterations
  n_it = itera - 1

  #-------------------------------------------------------------------------


  #-------------------------------------------------------------------------
  # Vari?ncias
  #-------------------------------------------------------------------------

  inv_teta=ks::invvech(teta[1:(s-1),])
  sol_invvech = teta[s,]*(solve(inv_teta))

  SOL_BETA=rep(list(sol_invvech),ncluster)
  SOL_BETA=bdiag_m(SOL_BETA)

  #Matrix::Diagonal(dim(T5)[1])
  AJ=Matrix::solve(T5+SOL_BETA, sparse = T) #Eigen

  CJ = Q1%*%E_IJ - Matrix::t(T2_nw)%*%AJ%*%S_2_1
  mat_c = wj_star_4%*%CJ%*%Matrix::t(CJ)              #valores levemente diferentes!
  mat_c = matrix(vector_eigen(mat_c, xxncluster), ncluster, nvar^2, byrow = T)

  R_ = split(R_, rep(1:ncluster,each=ss))
  R_ = lapply(R_, function(i) matrix(i,s))
  R_ = bdiag_m(R_)

  S_ = split(S_, rep(1:ncluster,each=s))
  S_ = bdiag_john(S_)

  TETA = rep(list(teta),ncluster)
  TETA = bdiag_john(TETA)

  DKJ = S_ - R_%*%TETA
  DKJ = DKJ%*%Matrix::t(DKJ)

  mat_d = wj_star_5%*%DKJ
  mat_d = matrix(vector_eigen(mat_d, ssncluster), ncluster, ss, byrow = T)

  #----------Variances-----------------
  var_beta = solve(s_matp)%*%((ncluster /( ncluster-1))*matrix(colSums(mat_c),nvar))%*%solve(s_matp)
  dp_beta = sqrt(diag(var_beta))
  dp_beta

  #var_teta = 2*solve(r_mat)
  var_teta = solve(r_mat)%*%(ncluster/(ncluster-1)* matrix(colSums(mat_d),s))%*% solve(r_mat)

  dp_teta = sqrt(diag(var_teta))


  z_star = beta/dp_beta
  z_star2= teta/dp_teta

  z_star_l= beta-abs(qnorm(0.025))*dp_beta
  z_star_u= beta+abs(qnorm(0.025))*dp_beta
  pz_star= 2*(1-pnorm(abs(z_star)))

  z_star2_l= teta-abs(qnorm(0.05))*dp_teta
  z_star2_u= teta+abs(qnorm(0.05))*dp_teta
  pz_star2= 2*(1-pnorm(abs(z_star2)))

  beta_mod <<- beta
  teta_mod <<- teta


  #-------------------------------------------------------------------------


  #/*-------------------------------------------------------------------------
  #  Residuals
  #-------------------------------------------------------------------------*/
  T5_res_trans=Matrix::t(T5_res)

  SIGMAU = rep(list(inv_teta),ncluster)
  SIGMAU = bdiag_m(SIGMAU)
  RHJ=SIGMAU%*%T5_res
  VJ = T5_res_trans%*%RHJ + Matrix::Diagonal(x = rep(teta[s,],nsubjc_t))

  #Matrix::Diagonal(dim(VJ)[1])
  AUX=RHJ %*% Matrix::solve(VJ, sparse = T)

  MM = AUX%*%E_IJ
  MM = as(MM,"dgCMatrix")
  U_1 = vector_eigen(MM, qncluster)
  u = matrix(U_1, ncluster, byrow = T) #nrow(inv_teta) = q ?

  dp_u = SIGMAU - AUX%*%Matrix::t(RHJ)
  dp_u = sqrt(diag_eigen(dp_u))       #nivel 2
  dp_u = matrix(dp_u, ncluster, byrow = T)

  U_2 = split(U_1, rep(1:ncluster,each=q))
  U_2 = bdiag_john(U_2)
  v = E_IJ - T5_res_trans%*%U_2    #nivel 1
  v = vector_eigen(v,nsubjc_t)

  var_v = rep(teta[s,], nsubjc_t)*NP

  dp_v <-  sqrt(var_v)
  u_pad = u /  dp_u
  v_pad = v / sqrt(var_v)
  #-------------------------------------------------------------------
  # N?vel 2: u dp_u u_pad
  # N?vel 1: v dp_v v_pad
  #-------------------------------------------------------------------

  residuosN2 <- data.frame(u, u_pad, dp_u)
  namesun2 <- c("u0", "u1", "upad0", "upad1", "dpu0", "dpu1")
  names(residuosN2) <- namesun2
  residuosN2$li <- residuosN2$u0-1.96*residuosN2$dpu0
  residuosN2$ls <- residuosN2$u0+1.96*residuosN2$dpu0
  residuosN2$idCluster <- cluster


  residuosN1 <- data.frame(v, v_pad, sqrt(var_v))
  namesun1 <- c("e0", "epad0", "dpe0")
  names(residuosN1) <- namesun1

  residuosN2 <<- residuosN2
  residuosN1 <<- residuosN1

  #finish= st_global("c(current_time)")
  finish<-Sys.time()

  #printing results

  testeprint2<-function(){
    cat("------------------------------------------------------------------------\n")
    cat("*-*-*-*-*-*-*-*-*-*-*-*-*-*-*PWIGLS*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-* \n")
    cat("------------------------------------------------------------------------\n")
    cat("Vari?vel Resposta =", vresp, "\n" )
    cat("Come?ou a rodar: ", sprintf("%19s\n",start))
    cat("N?mero de Itera??es    = " , n_it, "\n")
    cat("N?mero de unidades n?vel 1 = ", nsubjc_t, "\n")
    cat("N?mero de unidades n?vel 2 = ", ncluster, "\n")
    cat("------------------------------------------------------------------------\n")
    cat(" Ef.Fixos    Coef.    S.E.      z   P>|z|   [95%Interval.Conf]  Val.Ini.\n")
    cat("------------------------------------------------------------------------\n")
    for (mi in 1:nvar ) {
      cat(sprintf("%9s %8.3f %6.3f  %6.2f  %6.4f  %8.3f   %8.3f  %8.3f\n",name1[mi], beta[mi,], dp_beta[mi],z_star[mi,], pz_star[mi,] , z_star_l[mi,],z_star_u[mi,],beta0[mi,] ))
    }
    cat("------------------------------------------------------------------------\n")
    cat(("Comp.Var.    Coef.    S.E.     z   P>|z|   [95%Interval.Conf]  Val.Ini.\n"))
    cat("------------------------------------------------------------------------\n")
    for (mi in 1:s) {
      cat(sprintf("%9s %8.3f %6.3f  %6.2f  %6.4f  %8.3f   %8.3f  %8.3f\n",name3[mi,], teta[mi,], dp_teta[mi],z_star2[mi,], pz_star2[mi,] , z_star2_l[mi,],z_star2_u[mi,],teta0[mi,] ))
    }
    cat("------------------------------------------------------------------------\n")
    cat("Nota: Erros Padr?o Robustos \n")
    cat("Terminou de rodar em: ", sprintf("%19s\n",finish))
    cat("?Alinne Veiga - Janeiro 2018\n")
    cat("------------------------------------------------------------------------\n")

  }


  testeprint2()

}
