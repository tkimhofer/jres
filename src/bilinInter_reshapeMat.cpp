#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace arma;

// [[Rcpp::export]]
arma::mat biliInter(arma::vec x, arma::vec y, arma::mat z, arma::vec xout, arma::vec yout) {
  
  arma::mat zout;
  arma::interp2(x, y, z, xout, yout, zout);
  
  return zout;
  
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// /*** R
// biliInter(t1.temp, t2.temp, df)
// */



// [[Rcpp::export]]
arma::mat reshape_submatToMat(int si_f2, int si_f1, int xdim_f2, int xdim_f1, arma::rowvec spec) {
  
  
  mat spec_m(1, spec.size());
  spec_m.row(0) = spec;
  
  // calculate number of submatrices
  int nsub_f1 = si_f1 / xdim_f1;
  int nsub_f2 = si_f2 / xdim_f2;
  
  // calc nb of entries in submatrix
  int nsub_n = xdim_f2 * xdim_f1;
  
  // create output matrix
  mat out(si_f1, si_f2);
  //  mat smat(xdim_f1, xdim_f2);
  
  int count = 1;
  Rcpp::IntegerVector idx_m_row;
  Rcpp::IntegerVector idx_m_col;
  int idx_s_start;
  int idx_s_end;
  Rcpp::IntegerVector idx_sub;
  
  for ( int j=0; j < nsub_f1; j++ ) {
    idx_m_row = Rcpp::IntegerVector::create( (((j) * xdim_f1 )), ((j+1)*xdim_f1)-1 );
    
    for ( int i=0; i < nsub_f2; i++ ) {
      idx_m_col = Rcpp::IntegerVector::create( (((i) * xdim_f2 )), ((i+1)*xdim_f2)-1 );
      
      
      // Rcpp::Rcout << idx_m_col << "!!!!!!!!!!!\n---";
      
      idx_s_start = (count - 1) * nsub_n;
      idx_s_end = ( count * nsub_n ) - 1;
      
      // Rcpp::Rcout << spec_m.size() << "\n---";
      // Rcpp::Rcout << idx_s_start << "\n";
      // Rcpp::Rcout << idx_s_end << "\n";
      
      mat inter = spec_m.submat(0,  idx_s_start, 0, idx_s_end);
      // Rcpp::Rcout << "fine\n";
      // Rcpp::Rcout << inter.size() << "\n";
      
      inter.reshape(xdim_f2, xdim_f1);
      inter=inter.t();
      
      // Rcpp::Rcout << inter.n_rows() << "\n";
      // Rcpp::Rcout << inter.n_cols() << "\n";
      // 
      // Rcpp::Rcout << "---\n";
      // Rcpp::Rcout << idx_m_row << "\n";
      // Rcpp::Rcout << idx_m_col << "\n";
      
      // return inter;
      out.submat(idx_m_row(0),  idx_m_col(0), idx_m_row(1), idx_m_col(1)) = inter;
      
      count +=1;
    }
    
  }
  return out;
}



/*** R
reshape_submatToMat(si_f2 = 20, si_f1 = 8, xdim_f2 = 6, xdim_f1 = 4, spec = seq(160, by=-2,  length.out = 160 ))
*/