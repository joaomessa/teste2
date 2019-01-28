
#include <RcppEigen.h>
#include <iostream.h>
#include <Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SparseMatrixXd;

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
Eigen::MatrixXd invvech_eigen(const Eigen::VectorXd & x){
  int k=0;
  int leng = x.size();
  int tam=(sqrt(8*leng+1)-1)/2;
  Eigen::MatrixXd mat(tam, tam);
  for(int i=0; i<tam; ++i){
    for(int j=i; j<tam; ++j){
      mat(i,j)=x[k];
      k=k+1;
    }
  }
  return mat.selfadjointView<Eigen::Upper>();
}

// [[Rcpp::export]]
Eigen::MatrixXd rcppeigen_hello_world() {
    Eigen::MatrixXd m1 = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd m2 = Eigen::MatrixXd::Random(3, 3);

    return m1 + 3 * (m1 + m2);
}


// another simple example: outer product of a vector,
// returning a matrix
//

// and we can use Rcpp::List to return both at the same time
//
// [[Rcpp::export]]
Rcpp::List rcppeigen_bothproducts(const Eigen::VectorXd & x) {
    Eigen::MatrixXd op = x * x.transpose();
    double          ip = x.transpose() * x;
    return Rcpp::List::create(Rcpp::Named("outer")=op,
                              Rcpp::Named("inner")=ip);
}

//Multiplicando matrix 3x3 com vetor [1 2 3]
// [[Rcpp::export]]
Eigen::VectorXd functest1(const Eigen::MatrixXd & x)
{
  Eigen::VectorXd v(3);
  v << 1, 2, 3;
  Eigen::VectorXd u = x * v;
  return u;
}

//Multiplicando matrix nxn com vetor nx1
// [[Rcpp::export]]
Eigen::VectorXd functest2(const Eigen::MatrixXd & x, const Eigen::VectorXd & i)
{
  Eigen::VectorXd u = x * i;
  return u;
}

//Multiplicando matrix nxn com vetor nx1
// [[Rcpp::export]]
Rcpp::List autos(const Eigen::MatrixXd & x){
  Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(x);
  return Rcpp::List::create(Rcpp::Named("valores")=eigensolver.eigenvalues(),
                            Rcpp::Named("vetores")=eigensolver.eigenvectors());
}

// [[Rcpp::export]]
Rcpp::List pw(int nvar,
              int itera,
              int ncluster,
              const Eigen::MatrixXd & info_j,
              const Eigen::MatrixXd & x,
              const Eigen::VectorXd & y,
              const Eigen::MatrixXd & z,
              const Eigen::VectorXd & wi_j_star,
              const Eigen::VectorXd & wj_star,
              const Eigen::VectorXd & teta,
              int s,
              const Eigen::VectorXd & teta0){
  int nvar_quad = std::pow(nvar,2);
  Eigen::MatrixXd matp(ncluster,nvar_quad); //mudar o 4
  Eigen::MatrixXd matq(ncluster,nvar); //mudar o 2

  for(int i = 0; i < ncluster; ++i){
    int a = info_j(i,0);
    int b = info_j(i,1);
    Eigen::VectorXd yj = y.segment(a-1,b-a+1); //Outro formato de matriz 'info_j' era melhor
    Eigen::MatrixXd xj = x.block(a-1,0,b-a+1,x.cols());
    Eigen::MatrixXd zj = z.block(a-1,0,b-a+1,x.cols());

    Eigen::VectorXd v = wi_j_star.segment(a-1,b-a+1);
    Eigen::MatrixXd diag = v.asDiagonal();
    double wj = wj_star(i);


    Eigen::MatrixXd t1j = xj.transpose()*diag*xj;
    Eigen::MatrixXd t2j = xj.transpose()*diag*zj;
    Eigen::VectorXd t3j = xj.transpose()*diag*yj;
    Eigen::VectorXd t4j = zj.transpose()*diag*yj;
    Eigen::MatrixXd t5j = zj.transpose()*diag*zj;

    Eigen::MatrixXd aj;
    if(itera!=1){
      aj = (t5j + teta(s-1)*invvech_eigen(teta.head(s-1)).inverse()).inverse();
    }
    else {
      aj = (t5j + teta0(s-1)*invvech_eigen(teta0.head(s-1)).inverse()).inverse();
    }

    Eigen::MatrixXd t0 = t2j*aj;
    Eigen::MatrixXd t1 = wj*(t1j - t0*t2j.transpose());
    Eigen::VectorXd t11(Eigen::Map<Eigen::VectorXd>(t1.data(), t1.cols()*t1.rows()));
    Eigen::VectorXd t2 = wj*(t3j - t0*t4j);
    matp.row(i) = t11;
    matq.row(i) = t2;
   }

  return Rcpp::List::create(Rcpp::Named("matp")=matp,
                            Rcpp::Named("matq")=matq);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> teste_sparse(){
  SparseMatrixXd matA(10, 10);
  matA.coeffRef(0, 0) = 1;
  matA.coeffRef(1, 1) = 1;
  matA.makeCompressed();
  return matA;
}

// [[Rcpp::export]]
Eigen::VectorXd vector_eigen(const Eigen::SparseMatrix<double> & A,
                             int a){
  Eigen::VectorXd v(a);
  int i = 0;
  for (int k=0; k<A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
    {
      v[i] = it.value();
      i = i+1;
    }
  return v;
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> diagonalizar(const Eigen::SparseMatrix<double> & AJ_TETA,
                                         const Eigen::SparseMatrix<double> & SOL_TETA,
                                         const Eigen::SparseMatrix<double> & H_kj,
                                         const Eigen::VectorXd & delta){
  return AJ_TETA*SOL_TETA*H_kj;
}

// [[Rcpp::export]]
Eigen::VectorXd diag_eigen(const Eigen::SparseMatrix<double> & A){
  Eigen::VectorXd v = A.diagonal();
    return v;
}

// [[Rcpp::export]]
Eigen::VectorXd vecto_eigen(const Eigen::SparseMatrix<double> & A,
                             int a){
  Eigen::VectorXd v(a);
  int i = 0;
  for (int k=0; k<A.outerSize(); ++k)
    for (Eigen::SparseMatrix<double>::InnerIterator it(A,k); it; ++it)
    {
      v[i] = it.value();
      i = i+1;
    }
    return v;
}
