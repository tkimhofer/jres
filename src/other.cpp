#include <Rcpp.h>
using namespace Rcpp;

//' @title List all combinations of two vectors
//' @param x Vector A
//' @param y Vector B
//' @return all combinations, iterating first on B, then on A
// [[Rcpp::export]]
NumericMatrix expandGrid_rcpp(NumericVector x, NumericVector y) {
  
  int le_x = x.size();
  int le_y = y.size();
  int count = 0;
  NumericMatrix out(le_x*le_y, 2);
  
  
  for(int j=0; j < le_y; ++j){
    for(int i = 0; i < le_x; ++i) {
      
      out(count, 0) = x(i);
      out(count, 1) = y(j);
      
      count += 1;
    }
  }
  
  return out;
}


//' @title Symmetrise signal using f2 as mirror axis and using min values
//' @param A Matrix A 
//' @return summetrised A
// [[Rcpp::export]]
NumericMatrix matsymminf2f1(NumericMatrix A, int idx_cent_row, int idx_cent_col) {
  
  NumericMatrix out(A.rows(), A.cols());
  
  int indx_row = std::min( (idx_cent_row), (A.rows()-idx_cent_row) );
  int indx_col = std::min( (idx_cent_col), (A.cols()-idx_cent_col) );
  
// Rcout << indx ;
  
  for(int i = 0; i < indx_row; ++i) {
    for(int j = 0; j < indx_col; ++j) {
      out(idx_cent_row - i, idx_cent_col - j) =  out(idx_cent_row + i, idx_cent_col + j) = std::min( A(idx_cent_row - i, idx_cent_col - j), A(idx_cent_row + i, idx_cent_col + j) ); 
      out(idx_cent_row - i, idx_cent_col + j) =  out(idx_cent_row + i, idx_cent_col - j) = std::min( A(idx_cent_row - i, idx_cent_col + j), A(idx_cent_row + i, idx_cent_col - j) ); 
    }
  }
  
  return out;
}

// [[Rcpp::export]]
NumericMatrix matsymminf2(NumericMatrix A, int idx_cent) {
  
  NumericMatrix out(A.rows(), A.cols());
  
  int indx = std::min( (idx_cent), (A.cols()-idx_cent) );
  Rcout << indx ;
  
  for(int i = 0; i < A.rows(); ++i) {
    for(int j = 0; j < indx; ++j) {
      out(i, indx - j) =  out(i, indx + j) = std::min( A(i, indx - j), A(i, indx + j) ); 
    }
  }
  
  return out;
}



// test=matrix(c(1,2,3,1,2,3,1,2,3), nrow=3, byrow=T)
// matsym_minf2(test, 1) 
//   




//' @title Vector boundary finder when adding number to index 
//' @param i index of vector
//' @param add_int index count added / subtracted from i
//' @param mat_end length of the vector 
//' @param type addition or subtraction (specified by + or - symbol)
//' @return all combinations, iterating first on B, then on A
// [[Rcpp::export]]
int getBoundary(int i, int add_int, int mat_end, char type='-') {
  // indexing starts at zero, so mat_end is acutally ncol-1 in index terms
  
  
  if(type == '-'){
    int low_mat = i - add_int;
    IntegerVector bound_low = IntegerVector::create(low_mat,  0);
    int out = max(bound_low);
    return out;
  }
  
  if(type == '+'){
    int high_mat = i + add_int;
    IntegerVector bound_high = IntegerVector::create(high_mat,  mat_end);
    int out = min(bound_high);
    return out;
  }
  
}


//' @title Summation of two matrices of the same dimensions
//' @param A Matrix A 
//' @param A Matrix B 
//' @return A + B
// [[Rcpp::export]]
NumericMatrix matsum(NumericMatrix A,NumericMatrix B) {
  
  NumericMatrix out(A.rows(), A.cols());
  for(int i = 0; i < A.rows(); ++i) {
    for(int j = 0; j < A.cols(); ++j) {
      out(i,j) = A(i,j) + B(i,j); 
    }
  }
  
  
  return out;
}

//' @title Power transform each matrix element
//' @param Matrix A
//' @param exp power
//' @return A^exp
// [[Rcpp::export]]
NumericMatrix matpow(NumericMatrix A, const double exp) {
  
  NumericMatrix out(A.rows(), A.cols());
  for(int i = 0; i < A.rows(); ++i) {
    for(int j = 0; j < A.cols(); ++j) {
      out(i,j) = pow(A(i,j), exp); 
    }
  }
  return out;
}


//' @title Summation of all matrix elements to a double
//' @param Matrix A
//' @return sum(A)
// [[Rcpp::export]]
double matsumup(NumericMatrix A) {
  
  double out = 0;
  for(int i = 0; i < A.rows(); ++i) {
    for(int j = 0; j < A.cols(); ++j) {
      out += A(i,j); 
    }
  }
  return out;
}

//' @title Min-max normalisation for matrix
//' @param Matrix A
//' @return min-max-scaled A (values range from 0-1)
// [[Rcpp::export]]
NumericMatrix minmaxMat(NumericMatrix A) {
  double mins = min(A);
  double mdif = max(A) - mins;
  
  NumericMatrix out(A.rows(), A.cols());
  
  for(int i = 0; i < A.rows(); ++i) {
    for(int j = 0; j < A.cols(); ++j) {
      out(i,j) = ( A(i,j) - mins ) / mdif; 
    }
  }
  return out;
  
}

//' @title Find indices of points on axis that are closest after adding/subtracting a constant from a center point
//' @param ppm Axis/ scale vector
//' @param i index for ppm, determining point where constand should be added/subtracted
//' @param add Constand to be added/subtracted (on ppm scale)
//' @return indieces of ppm (lower, upper)
// [[Rcpp::export]]
IntegerVector getBBdim_rcpp(NumericVector ppm, int i, double add) {
  
  
  int idx_ppm_min = which_min( abs(( ppm[i] - add ) - ppm ) );
  int idx_ppm_max = which_min( abs(( ppm[i] + add ) - ppm ) );
  
  IntegerVector idx_ppm = { idx_ppm_min , idx_ppm_max };
  
  return(idx_ppm);
  
}


//' @title Calculate Laplacian of Gaussians to determine peak size / bounding box dimension
//' @param sub Matrix with signal
//' @param f1hz Row scale of sub (Jres f1, expresed in Hz)
//' @param f2hz Col scale of sub (Jres f2, expressed in Hz)
//' @param cent_f1hz Location of peak maximum in f1 dimension
//' @param cent_f2hz Location of peak maximum in f2 dimension
//' @param sf Spectrometer frequency
//' @param npix SD parameter values used to create Gaussians - should span expected Jres signal fwhm range 
//' @return List, 1st element: Best matching Gaussian SD (peak size), 2nd element: Evaluation criterium for all tested SD
// [[Rcpp::export]]
List lapOfG(NumericMatrix sub, NumericVector f1hz, NumericVector f2hz, double cent_f1hz, double cent_f2hz, double sf, NumericVector npix) {
  
  //npix seq(0.01, 3, by=0.05)
  
  sub = sub + abs(min(sub));
  NumericMatrix te=expandGrid_rcpp(f1hz, f2hz);
  
  // smoothing vector
  NumericVector fconv = {0.3333333, 0.7453560, 0.3333333,0.7453560, 1, 0.7453560, 0.3333333, 0.7453560, 0.3333333};
  fconv.attr("dim") = Dimension(3, 3);
  NumericMatrix fm = as<NumericMatrix>(fconv);
  
  
  // apply slight smoothing
  // grab the R function F1
  Function convolve3cpp( "convolve3cpp" );
  NumericMatrix sub_smooth = as<NumericMatrix>(convolve3cpp(sub, fm));
  //sub_smooth = minmaxMat(sub_smooth);
  // perform laplacian of gaussians
  NumericVector lapl = {0.5, 1, 0.5, 1, -6, 1, 0.5, 1, 0.5};
  lapl.attr("dim") = Dimension(3, 3);
  NumericMatrix laplm = as<NumericMatrix>(lapl);
  NumericVector out;
  
  List retout;
  int count = 0;
  //List siml;
  double sdiff;
  
  for(int i = 0; i < npix.size(); ++i) {
    
    NumericVector sig = {npix(i), npix(i)*1.3};
    NumericVector mus = {cent_f1hz, cent_f2hz};
    NumericVector ta = te( _ , 0 );
    
    Function d2gauss_cpp( "d2gauss_cpp" );
    NumericVector gay1 = as<NumericVector>(d2gauss_cpp(te(_,0),te(_,1), mus, sig));
    
    //NumericVector gay1 = d2gauss_cpp(x = te(_,0), y = te(_,1), mu = mus, sigma = sig);
    gay1.attr("dim") = Dimension(f1hz.size(), f2hz.size());
    
    NumericMatrix gaula2 = as<NumericMatrix>(convolve3cpp(gay1, laplm)) * pow(npix(i), 2);
    //siml.push_back(gaula2);
    
    NumericMatrix inter = gaula2 * (-1/min(gaula2));
    NumericMatrix ssmomm = minmaxMat(sub_smooth);
    NumericMatrix sim =  matsum(ssmomm, inter); // export this
    
    NumericMatrix inter1 = matpow(sim, 2);
    
    double simadd = matsumup(inter1);
    
    out.push_back(simadd);
    
    
    // stop fitting if sum of LoG is increasing again (by > 1% in three steps)
    if( i > 1)
      
      sdiff = ( out[i] - out[i-1] ) / out[i-1] ;
      
      if( sdiff > 0.01)
        ++count;
      
    if( count == 3 )
      break;
    
    
    
    
    //out.push_back(List::create(gaula2, inter, sub_smooth, ssmomm, sim, inter1, simadd));
  }
  
  int idx_min = which_min(out);
  retout = List::create(npix(idx_min), out, npix);
  return retout;
  
}


//' @title Perform peak picking and get bounding box dimensions
//' @param jr Jres matrix, f1 in rows and f2 in cols
//' @param f1hz F1 scale of jr (Hz)
//' @param f2hz F2 scale of jr (ppm)
//' @param noise Intensity threshold for noise (no peaks are detected below this value)
//' @param boundary Initial bounding box estimate for determining peak size (one side estimate: x +/- boundary), this should be large enought to capture large signals (expressed in Hz)
//' @param sf Spectrometer frequency
//' @return List of dataframes summarising the detected peaks/features
// [[Rcpp::export]]
List pickPeaks_rcpp(NumericMatrix jr, NumericVector f1hz, NumericVector f2ppm, double noise, double boundary, double sf) {
  
  NumericVector f1ppm = f1hz / sf;
  NumericVector f2hz = f2ppm * sf;
  
  NumericVector dfb = diff(f1hz); // 0.3 Hz in one increment
  double doub = boundary / median(dfb); // if boundary is 10, then 34 indices increments amounts to about 10 Hz
  double add_row = round( doub ) +1 ; // if boundary is 10, then 34 indices increments amounts to about 10 Hz
  
  NumericVector df2b = diff(f2hz);
  double doub1 = boundary / median(df2b);
  int add_col = round( doub1 ) +1; // if boundary is 10, then 24 indices increments amounts to about 10 Hz
  
  int n = f1ppm.size();
  int m = f2ppm.size();
  
  // generate sd for LoG -> upper limit sd is determined by boundary (max sd = boundary / 10) -> related to gaussian
  // boundary is added and subtracted from mean, that means total interval captured in sub is 2 * boundary
  NumericVector npix;
  
  double val = 0.01;
  double stepsize = 0.05;
  double upperlim = (boundary*2/10); // restrict size of gaussian to boundary size (otherwise sideband effects after convolution)
  npix.push_back(val);
  
  while (val < upperlim) {
    val += stepsize;
    npix.push_back(val);
  }
  // for(int i = 1; i < upperlim; ++i) {
  //   val += stepsize;
  //   npix.push_back(val);
  // }
  
  String sym;
  List res;
  NumericMatrix sub_sym;
  // define npix for LoG (sigma)
  
  for(int i = 1; i < (n-1); ++i) {
    for(int j = 1; j < (m-1); ++j) {
      
      if(jr(i,j) > jr(i-1,j) & jr(i,j) > jr(i+1,j) & jr(i,j) > jr(i,j-1) & jr(i,j) > jr(i,j+1) & jr(i,j) > noise) {
        
        IntegerVector out = {i, j};
        sym = "no";
        // define center in hz and ppm
        double cent_f1hz_sub = f1hz(i);
        double cent_f2hz_sub = f2hz(j);
        
        
        // create sub matrix that can be input for LoG
        // get lower and upper boundary (in case add_col/add_row exeed matrix dimennsions)
        int f1l_idx = getBoundary(i, add_row, jr.nrow()-1, '-');
        int f1h_idx = getBoundary(i, add_row, jr.nrow()-1, '+');
        
        int f2l_idx = getBoundary(j, add_col, jr.ncol()-1, '-');
        int f2h_idx = getBoundary(j, add_col, jr.ncol()-1, '+');
        
        // IntegerVector rowidx = {f1l_idx, f1h_idx};
        // IntegerVector colidx = {f2l_idx, f2h_idx};
        
        // define submatrix that can be inputed into LoG
        NumericMatrix sub = jr( Range(f1l_idx, f1h_idx), Range(f2l_idx, f2h_idx)); // submatrix
        
        // update f1 hz and f2 ppm range for sub
        NumericVector f1hz_upd = f1hz[Range(f1l_idx, f1h_idx)];
        NumericVector f2hz_upd = f2hz[Range(f2l_idx, f2h_idx)];
        
        // define indices for center position in sub
        //IntegerVector out_upd = {(out(0) - (f1l_idx)), (out(1) - (f2l_idx))}; // indices of max in submatrix
        
        // get center index of sub and symmetrise
        int cent_f1_idx = i - f1l_idx ;
        int cent_f2_idx = j - f2l_idx ;
        
        // define center in hz and ppm
        double cent_f1hz_upd = f1hz_upd(cent_f1_idx);
        double cent_f2hz_upd = f2hz_upd(cent_f2_idx);
       
        // include function here to perform LoG for peak matching
        // List logres = lapOfG(sub, f1hz_upd, f2hz_upd, cent_f1hz_sub, cent_f2hz_sub, sf, npix);
        
        // re-do lapOfG after f2 axis symmetrisation (peak overlap results in too small bb) 
        
        
       
        
        List logres = lapOfG(sub, f1hz_upd, f2hz_upd, cent_f1hz_sub, cent_f2hz_sub, sf, npix);
        double gsd_ori = logres[0];
        double add_gsd = gsd_ori;
        
       
        
        // 
        // NumericMatrix sub_symf2 = matsymminf2(sub, cent_f2_idx);
        // List logres_symf2 = lapOfG(sub_symf2, f1hz_upd, f2hz_upd, cent_f1hz_sub, cent_f2hz_sub, sf, npix);
        // double gsd_sumf2 = logres_symf2[0];
       
        if (gsd_ori < 0.25 ) {
          sub_sym = matsymminf2f1(sub, cent_f1_idx, cent_f2_idx);
          List logres_sym = lapOfG(sub_sym, f1hz_upd, f2hz_upd, cent_f1hz_sub, cent_f2hz_sub, sf, npix);
          add_gsd = logres_sym[0];
          sym = "yes";
        }
          
        
        // define bounding box with log sd value
        
        // characterise signal
        //double gsd_ori = logres[0]; 
        // Rcout << gsd_ori << "\n";
        
        // // peaks unequally extend in f1 and f2, therefore applying adjustment factors for f1 and f2
        double add1_f1 = gsd_ori * 1;
        double add1_f2 = gsd_ori * 1;
        
        doub = add1_f1 / median(dfb); // if boundary is 10, then 34 indices increments amounts to about 10 Hz
        add_row = round( doub ) +1; 
        // 
        doub1 = add1_f2 / median(df2b);
        add_col = round( doub1 ) +1;
        
        
        IntegerVector d_f1 ;
        IntegerVector d_f2 ;
        //
        d_f1[0] = getBoundary(i, add_row, jr.nrow(), '-');
        d_f1[1] = getBoundary(i, add_row, jr.nrow(), '+');

        d_f2[0] = getBoundary(j, add_col, jr.ncol(), '-');
        d_f2[1] = getBoundary(j, add_col, jr.ncol(), '+');

        // // 
        // //Integer d_f2 = {(which_min( abs(( f2hz[j] - add1_f2 ) - f2hz ))),  (which_max( abs(( f2hz[i] - add1_f2 ) - f1hz )))};
        // 
        // // Rcout << d_f1 << "\n";
        // // Rcout << d_f2 << "\n";
        NumericMatrix feat = jr( Range(d_f1[0], d_f1[1]) , Range(d_f2[0], d_f2[1]) );
        
        
        
        DataFrame outs = DataFrame::create( _["P.peakOri"] = 1, 
                                            _["cent.f1"] =cent_f1hz_sub,  
                                            _["cent.f2"] = cent_f2hz_sub / sf, 
                                            _["cent.f1_sub"] =cent_f1hz_upd,  
                                            _["cent.f2_sub"] = cent_f2hz_upd / sf, 
                                            _["LoG.hz_ori"] = gsd_ori,
                                            // _["LoG.hz_symf1f2"] = gsd_sumf2f1,
                                            // _["LoG.hz_sumf2"] = gsd_sumf2,
                                            //_["npix_max"] = max( npix ),
                                            // _["addrow"] = add_row,
                                            // _["addrow"] = add_col,
                                           // _["bb.width"] = add1,
                                           _["bb.width.f1"] = add_gsd,
                                           _["bb.width.f2"] = add_gsd,
                                            _["f1.idx"] = i,
                                            _["f2.idx"] = j,
                                            _["f1.idx_sub"] = cent_f1_idx,
                                           _["f2.idx_sub"] = cent_f2_idx,
                                            // _["f1.add.low"] = d_f1[0],
                                            // _["f1.add.high"] = d_f1[1],
                                            // _["f2.add.low"] = d_f2[0] ,
                                            // _["f2.add.high"] = d_f2[1],
                                             _["Int"] = jr(i,j),
                                             _["sym"] = sym
        );
        
        // 
        List sblist =  List::create(_["info"]=outs,
                                    _["sub"]=feat
                                    );
        res.push_back(sblist);
        
      }
      
    }
    
  }
  
  return(res);
  
}




















































