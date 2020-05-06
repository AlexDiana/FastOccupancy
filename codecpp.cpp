#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

const double TRUNC = .64;
const double TRUNC_RECIP = 1.0 / .64;
const double log2pi = std::log(2.0 * M_PI);

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

double aterm(int n, double x, double t) {
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }    
  return (double)exp(f);
}

double exprnd(double mu) {
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

double truncgamma() {
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

double randinvg(double mu) {
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

double tinvgauss(double z, double t) {
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }    
  return X;
}

double samplepg(double z) {
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;
  
  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);
  
  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q); 
  
  double u, X;
  
  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1) 
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss(z, t);
    }
    
    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;
    
    while(1) 
    {
      Sn = Sn + asgn * aterm(i, X, t);
      
      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      
      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

// [[Rcpp::export]]
double rpg(int n, double z){
  
  double x = 0;
  for(int i = 0; i < n; i++){
    x += samplepg(z);
  }
  
  return(x);
}

// [[Rcpp::export]]
double ldmvnorm_cpp(arma::vec data, arma::vec m, arma::mat Sigma){
  int xdim = data.size();
  arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(Sigma))));
  double rootisum = arma::sum(log(rooti.diag()));
  
  double constants = -(xdim/2) * log2pi;
  arma::vec z = rooti * ( data - m) ;  
  
  return constants - 0.5 * arma::sum(z%z) + rootisum;   
}

arma::vec mvrnormArma(arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mu + arma::chol(sigma) * Y;
}

arma::vec mvrnormArmaQuick(arma::vec mu, arma::mat cholsigma) {
  int ncols = cholsigma.n_cols;
  arma::vec Y = arma::randn(ncols);
  return mu + cholsigma * Y;
}

// [[Rcpp::export]]
arma::mat diagMatrixProd(arma::mat& X, arma::vec& D){
  // this is slow
  arma::mat result(X.n_rows, D.size());
  for(int j = 0; (unsigned)j < result.n_cols; j++){
    result.col(j) = X.col(j) * D(j);
  }
  
  // RMatrix<double>::Row rowi = A.row(i);
  
  // out(i,1) = std::inner_product(rowi.begin(), rowi.end(), x.begin(), 0.0);
  
  return(result);
}

// [[Rcpp::export]]
arma::vec sample_Omega_cpp(arma::mat& X, arma::vec& beta, arma::vec& n){
  
  int nsize = n.size();
  arma::vec Omega_vec(nsize);
  
  arma::vec Xbeta = X * beta;
  
  for(int i = 0; i < nsize; i++){
    
    Omega_vec[i] = rpg(n[i], Xbeta[i]);
    
  }
  
  return(Omega_vec);
}

// [[Rcpp::export]]
arma::vec sample_Omega_cpp_noXb(arma::vec& Xbeta, arma::vec& n){
  
  int nsize = n.size();
  arma::vec Omega_vec(nsize);
  
  for(int i = 0; i < nsize; i++){
    
    Omega_vec[i] = rpg(n[i], Xbeta[i]);
    
  }
  
  return(Omega_vec);
}

double log_L_gamma_cpp(arma::vec gamma, arma::mat X, arma::vec indexes_covariates,
                   arma::vec b, arma::mat B, arma::vec Omega, arma::vec k){


  IntegerVector index_present(indexes_covariates.size());
  int l = 0;
  for(int i = 0; (unsigned)i < indexes_covariates.size(); i++){
    if(gamma[indexes_covariates[i]-1] == 1){
      index_present[l] = i;
      l += 1;
    }
  }
  
  arma::mat X_gamma(X.n_rows, l);
  arma::vec b_gamma(l);
  arma::mat B_gamma(l, l);

  for(int i = 0; i < l; i++){
    X_gamma.col(i) = X.col(index_present[i]);
    b_gamma[i] = b[index_present[i]];
    for(int j = 0; j < l; j++){
      B_gamma(i,j) = B(index_present[i], index_present[j]);
    }
  }
  
  arma::mat tX = arma::trans(X_gamma);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat cholXgOmX = arma::chol(tXOmega * X_gamma + arma::inv(B_gamma));
  
  double firstTerm = (.5) * log(det(arma::inv(B_gamma))) - log(det(cholXgOmX));
  
  arma::mat tXKbplusBb = arma::trans(X_gamma) * k + arma::inv(B_gamma) * b_gamma;
  
  arma::vec v = solve(arma::trimatl(arma::trans(cholXgOmX)),tXKbplusBb);
  arma::mat vtv = arma::trans(v) * v;
  
  arma::mat secondTerm = - .5 * ( (arma::trans(b_gamma) * arma::inv(B_gamma) * b_gamma) - vtv);
      
  // double firstTerm = (.5) * log(det(arma::inv(B_gamma))) - 
  //   (.5) * log(det(arma::trans(X_gamma) * Omega * X_gamma + arma::inv(B_gamma)));
  
  // arma::mat tXKbplusBb = arma::trans(X_gamma) * k + arma::inv(B_gamma) * b_gamma;
  
  // arma::mat secondTerm = - .5 * ( (arma::trans(b_gamma) * arma::inv(B_gamma) * b_gamma) -
  //     arma::trans(tXKbplusBb) * arma::inv(arma::trans(X_gamma) * Omega * X_gamma + 
  //     arma::inv(B_gamma)) * tXKbplusBb);
  
  return(firstTerm + secondTerm(0,0));
}

int sample_int(IntegerVector samples) {
  // Rcpp::IntegerVector pool = Rcpp::seq(1, 10);
  std::random_shuffle(samples.begin(), samples.end());
  return samples[0];
} 

arma::vec sample_gamma_cpp(arma::vec gamma, arma::mat X, arma::vec Omega, 
                       arma::vec b, arma::mat B, int ncov, arma::vec k, 
                       arma::vec indexes_covariates, int fixedIndexes, 
                       double d_bar){

  arma::vec gamma_star = gamma;
  
  double h_ratio = 1;

  if(R::runif(0, 1) < .33333){ // add

    if((sum(gamma) - fixedIndexes) != ncov){

      // Find zero covariates
      int numOfZeroCov = ncov - (sum(gamma) - fixedIndexes);
      IntegerVector zeroCov(numOfZeroCov);
      int i = 0;
      for(int l = fixedIndexes; (unsigned)l < gamma.size(); l++){
        if(gamma[l] == 0) {
          zeroCov[i] = l;
          i += 1;
        }
      }
      
      int covariate_to_update = sample_int(zeroCov);
      
      gamma_star[covariate_to_update] = 1;

      h_ratio = (ncov -(sum(gamma) - fixedIndexes)) / (ncov - (sum(gamma) - fixedIndexes) - 1 + ((ncov - d_bar) / (d_bar)) );

    }

  } else if(R::runif(0, 1) < .5){ // delete

    if((sum(gamma) - fixedIndexes) != 0){
      
      // Find non zero covariates
      int numOfNonZeroCov = sum(gamma) - fixedIndexes;
      IntegerVector nonZeroCov(numOfNonZeroCov);
      int i = 0;
      for(int l = fixedIndexes; (unsigned)l < gamma.size(); l++){
        if(gamma[l] == 1) {
          nonZeroCov[i] = l;
          i += 1;
        }
      }
      
      int covariate_to_update = sample_int(nonZeroCov);
      
      gamma_star[covariate_to_update] = 0;

      h_ratio = (ncov - (sum(gamma) - fixedIndexes) + ((ncov - d_bar) / (d_bar)) ) / (ncov - (sum(gamma) - fixedIndexes) + 1);

    }

  } else { // swap

    if(((sum(gamma) - fixedIndexes) != 0) & ((sum(gamma) - fixedIndexes) != ncov)){

      // Find zero covariates
      int numOfZeroCov = ncov - (sum(gamma) - fixedIndexes);
      IntegerVector zeroCov(numOfZeroCov);
      int i = 0;
      for(int l = fixedIndexes; (unsigned)l < gamma.size(); l++){
        if(gamma[l] == 0) {
          zeroCov[i] = l;
          i += 1;
        }
      }
      int covariates2_to_swap = sample_int(zeroCov);
      
      // Find non zero covariates
      int numOfNonZeroCov = sum(gamma) - fixedIndexes;
      IntegerVector nonZeroCov(numOfNonZeroCov);
      i = 0;
      for(int l = fixedIndexes; (unsigned)l < gamma.size(); l++){
        if(gamma[l] == 1) {
          nonZeroCov[i] = l;
          i += 1;
        }
      }
      
      int covariates1_to_swap = sample_int(nonZeroCov);

      gamma_star[covariates1_to_swap] = 0;
      gamma_star[covariates2_to_swap] = 1;

      h_ratio = 1;

    }

  }
  
  double L_gamma_star = log_L_gamma_cpp(gamma_star, X, indexes_covariates, b, B, Omega, k);
  
  double L_gamma = log_L_gamma_cpp(gamma, X, indexes_covariates, b, B, Omega, k);
  
  h_ratio = h_ratio * exp(L_gamma_star - L_gamma);
  
  if(R::runif(0, 1) < h_ratio){
    gamma = gamma_star;
  }

  return(gamma);
}

// [[Rcpp::export]]
arma::vec sample_z_cpp(arma::vec psi, arma::vec p, arma::mat k_s){
  
  arma::vec z = arma::zeros(k_s.n_rows);
  
  // this loops over p
  int l = 0;
  
  for(int i = 0; (unsigned)i < z.size(); i++){
    
    if(k_s(i, 0) == 1){
      
      z[i] = 1;
      l += k_s(i,1);
      
    } else {
      
      double prod1mp = 1; //prod(1 - p[i + seq_len(k_s[i,4])])
      for(int k = 0; k < k_s(i,1); k++){
        prod1mp *= (1 - p(l + k));
      }
      
      double p_zsequal1 = (psi[i] * prod1mp) / (psi[i] * prod1mp + (1 - psi[i])) ;
      z[i] = R::rbinom(1, p_zsequal1);
      l += k_s(i,1);
      
    }
     
  }
  
  return(z);
}

// [[Rcpp::export]]
arma::mat vec_subset_mat(const arma::mat& x, const arma::uvec& idx) {
  return x.cols(idx);
}

// [[Rcpp::export]]
arma::vec subset_vec(const arma::vec& x, const arma::uvec& idx) {
  return x.elem(idx);
}

// [[Rcpp::export]]
arma::vec sample_beta_cpp(arma::mat& X, arma::mat& B, arma::vec& b, 
                          arma::vec& Omega, arma::vec& k){
  
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat XtOmegaX = tXOmega * X;
  
  // arma::mat cov_matrix = arma::inv(tXOmega * X + arma::inv(B));
  // arma::vec result = mvrnormArma(cov_matrix * (arma::trans(X) * k + arma::inv(B) * b), cov_matrix);
  
  
  arma::mat L = arma::trans(arma::chol(XtOmegaX + arma::inv(B)));
  arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
  
  arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
  
  return(result);
}

// [[Rcpp::export]]
List sample_gamma_beta_cpp(arma::vec gamma, arma::vec beta, arma::mat X, 
                           arma::vec b, arma::mat B, int D, arma::vec n, arma::vec k, 
                           arma::vec indexes_covariates, int fixedIndexes, double d_bar){
  
  // resize beta and X
  IntegerVector index_present(indexes_covariates.size());
  int l = 0;
  for(int i = 0; (unsigned)i < indexes_covariates.size(); i++){
    if(gamma[indexes_covariates[i]-1] == 1){
      index_present[l] = i;
      l += 1;
    }
  }
  
  arma::mat X_gamma(X.n_rows, l);
  arma::vec beta_gamma(l);
  
  for(int i = 0; i < l; i++){
    X_gamma.col(i) = X.col(index_present[i]);
    beta_gamma[i] = beta[index_present[i]];
  }
  
  // sample Omega
  arma::vec Omega = sample_Omega_cpp(X_gamma, beta_gamma, n);
  
  // sample gamma
  gamma = sample_gamma_cpp(gamma, X, Omega, b, B, D, k, indexes_covariates, fixedIndexes, d_bar);
  
  // resize X, b and B
  l = 0;
  for(int i = 0; (unsigned)i < indexes_covariates.size(); i++){
    if(gamma[indexes_covariates[i]-1] == 1){
      index_present[l] = i;
      l += 1;
    }
  }
  
  arma::mat X_gamma2(X.n_rows, l);
  arma::vec b_gamma(l);
  arma::mat B_gamma(l, l);
  
  for(int i = 0; i < l; i++){
    X_gamma2.col(i) = X.col(index_present[i]);
    b_gamma[i] = b[index_present[i]];
    for(int j = 0; j < l; j++){
      B_gamma(i,j) = B(index_present[i], index_present[j]);
    }
  }
  
  // sample beta
  arma::vec beta_new = sample_beta_cpp(X_gamma2, B_gamma, b_gamma, Omega, k);
  
  for(int i = 0; i < l; i++){
    beta[index_present[i]] = beta_new[i];
  }  
  
  return List::create(_["gamma"] = gamma,
                      _["beta"] = beta,
                      _["Omega"] = Omega);
}

// [[Rcpp::export]]
arma::mat matrixProductXtOmegaX(arma::mat& X, int Y, int ncov_psi, arma::vec Omega,
                                arma::vec X_y_index){
  
  arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
  
  XtOmegaX2(0, 0) = sum(Omega);
  
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    XtOmegaX2(X_y_index[i], X_y_index[i]) += Omega[i];
    
    for(int j = 0; j < ncov_psi; j++){
      
      XtOmegaX2(Y + 1 + j, X_y_index[i]) +=  X(i, Y + 1 + j) * Omega[i];
      
    }
  }
  
  for (int i = 1; i <= Y; i++) {
    
    XtOmegaX2(0, i) = XtOmegaX2(i, i);
    XtOmegaX2(i, 0) = XtOmegaX2(i, i);
    
    for(int j = 0; j < ncov_psi; j++){
      
      XtOmegaX2(i, Y + 1 + j) = XtOmegaX2(Y + 1 + j, i);
      
    }
    
  }
  
  for(int i = 0; i < ncov_psi; i++){
    
    double result = arma::as_scalar(Omega.t() * X.col(Y + 1 + i));
    
    XtOmegaX2(Y + 1 + i, 0) = result;
    XtOmegaX2(0, Y + 1 + i) = XtOmegaX2(Y + 1 + i, 0);
  }
  
  for(int i = 0; i < ncov_psi; i++){
    for (int j = 0; j <= i; j++) {
      arma::vec firstProduct = Omega % X.col(Y + 1 + i);
      arma::vec secondProduct = firstProduct % X.col(Y + 1 + j);
      XtOmegaX2(Y + 1 + i, Y + 1 + j) = sum(secondProduct);
    }
  }
  
  for (int i = 0; i < (ncov_psi- 1); i++) {
    for (int j = i; j < ncov_psi; j++) {
      XtOmegaX2(Y + 1 + i, Y + 1 + j) = XtOmegaX2(Y + 1+ j, Y + 1 + i);
    }
  }
  
  return(XtOmegaX2);
}

// [[Rcpp::export]]
arma::vec sample_beta_cpp_fast(arma::mat& X, arma::mat& B, arma::vec& b, 
                               arma::vec& k, arma::mat XtOmegaX){
  
  arma::mat tX = arma::trans(X);
  
  arma::mat L = arma::trans(arma::chol(XtOmegaX + arma::inv(B)));
  arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + arma::inv(B) * b);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
  
  arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
  
  return(result);
}

// [[Rcpp::export]]
List sampler_beta(arma::vec beta,
                  arma::vec a_s,
                  arma::mat X, 
                  arma::vec b, 
                  arma::mat B, 
                  arma::vec n, 
                  arma::vec k, 
                  int Y, 
                  int ncov_psi,
                  arma::vec X_y_index){
  
  arma::vec Xbeta = X * beta + a_s;
  arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
  
  // wrong
  arma::mat XtOmegaX = matrixProductXtOmegaX(X, Y, ncov_psi, Omega,
                                              X_y_index);
  
  // arma::mat tX = arma::trans(X);
  // arma::mat tXOmega = diagMatrixProd(tX, Omega);
  // arma::mat XtOmegaX = tXOmega * X;

  arma::vec knew = k - Omega % a_s;

  beta = sample_beta_cpp_fast(X, B, b, knew, XtOmegaX);
  
  return(List::create(_["beta"] = beta,
                      _["Omega"] = Omega,
                      _["XtOmegaX"] = XtOmegaX));
}

// [[Rcpp::export]]
List sampler_beta_integrated(arma::vec beta,
                             arma::vec a_s,
                             arma::vec mu_i,
                             arma::vec sigmasq_i,
                             arma::mat X, 
                             arma::vec b, 
                             arma::mat B, 
                             arma::vec n, 
                             arma::vec k, 
                             int Y, 
                             int ncov_psi,
                             arma::vec X_y_index){
  
  arma::vec Xbeta = X * beta + a_s;
  arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
  
  arma::mat XtOmegaX = matrixProductXtOmegaX(X, Y, ncov_psi, Omega,
                                                  X_y_index);
  
  // create new variables for integrated posterior
  arma::vec ktilde = arma::zeros(k.size());
  for(int i = 0; i < ktilde.size(); i++){
    
    ktilde[i] = k[i] - (Omega[i] / (Omega[i] + (1 / sigmasq_i[i]))) * (k[i] + mu_i[i] / sigmasq_i[i]);
    
  }
  
  arma::vec Omegatilde = arma::zeros(Omega.size());
  for(int i = 0; i < Omegatilde.size(); i++){
    
    Omegatilde[i] = Omega[i] * (1 - 1 / (Omega[i] + (1 / sigmasq_i[i])));
    
  } 
    
  arma::mat XtOmegaTildeX = matrixProductXtOmegaX(X, Y, ncov_psi, Omegatilde,
                                             X_y_index);
  
  beta = sample_beta_cpp_fast(X, B, b, ktilde, XtOmegaTildeX);
  
  return(List::create(_["beta"] = beta,
                      _["Omega"] = Omega,
                      _["XtOmegaX"] = XtOmegaX));
}

// [[Rcpp::export]]
List sample_beta_omega_cpp(arma::vec beta,
                           arma::mat X, 
                           arma::vec b, 
                           arma::mat B, 
                           arma::vec n, 
                           arma::vec k){
  
  // sample Omega
  arma::vec Omega = sample_Omega_cpp(X, beta, n);
  
  arma::mat tX = arma::trans(X);
  arma::mat tXOmega = diagMatrixProd(tX, Omega);
  arma::mat XtOmegaX = tXOmega * X;
  
  beta = sample_beta_cpp_fast(X, B, b, k, XtOmegaX);
  
  return(List::create(_["beta"] = beta,
                      _["Omega"] = Omega));
}

// SAMPLER L (SCALE PARAMETER OF GAUSSIAN PROCESS)

// [[Rcpp::export]]

double k_cpp(double x1, double x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-(x1-x2)*(x1-x2)/(2*pow(l,2)));
  // return 1;
}

// [[Rcpp::export]]
arma::mat K(arma::vec x1, arma::vec x2, double a, double l){
  arma::mat res(x1.size(), x2.size());
  
  for(int i = 0; (unsigned)i < x1.size(); i++){
    for(int j = 0; (unsigned)j < x2.size(); j++){
      res(i,j) = k_cpp(x1[i],x2[j], a, l);
    }  
  }
  
  return res;
}

// [[Rcpp::export]]
double loglikelihood_l_gp_cpp(double l, 
                              double a, 
                              int Y, 
                              arma::mat XtOmegaX, 
                              arma::mat Xtz){
  
  arma::vec years(Y);
  std::iota(years.begin(), years.end(), 1);
  arma::mat K_l = K(years, years, a, l);
  
  arma::mat invKl = arma::inv(K_l);
  
  arma::mat XtOmegaXpsolveK_l = XtOmegaX + invKl;
  arma::mat identityMatrix = arma::eye(invKl.n_rows, invKl.n_rows);
  arma::mat invKlXtOmegaX = invKl * XtOmegaX;
  arma::mat IpXtOmegaXpsolveK_l = identityMatrix + invKlXtOmegaX;
  
  double logdet1det2;
  double sign;
  
  log_det(logdet1det2, sign, IpXtOmegaXpsolveK_l);
  
  arma::vec term = arma::trans(Xtz) * arma::inv(XtOmegaXpsolveK_l) * Xtz;
  
  return (- .5 * logdet1det2 + .5 * term[0]);
}

// [[Rcpp::export]]
double sampler_l(double l, double a, arma::vec beta, arma::mat X,
                 double a_l, double b_l, double sd_l,
                 arma::vec a_s, int Y,
                 arma::mat XtOmegaX, 
                 arma::vec k, arma::vec Omega,
                 arma::vec X_y_index){
  
  arma::vec idxes = arma::zeros(X.n_cols - Y);
  idxes[0] = 0;
  for(int i = 0; (unsigned)i < (idxes.size() - 1); i++) idxes[i + 1] = 1 + Y + i;
  arma::uvec uidxes = arma::conv_to<arma::uvec>::from(idxes);
  arma::mat X_no_y = vec_subset_mat(X, uidxes);
  arma::vec beta_no_y = subset_vec(beta, uidxes);
  
  arma::vec c_i = X_no_y * beta_no_y + a_s;
  
  arma::vec y = k - Omega % c_i;
  
  arma::vec Xtz = arma::zeros(Y);
  for(int i = 0; (unsigned)i < X_y_index.size(); i++){
    Xtz[X_y_index[i] - 1] += y[i];
  }
  
  double l_star = R::rnorm(l, sd_l);
  
  arma::mat XtOmegaX_subset = arma::zeros(Y, Y);
  for(int i = 0; i < Y; i++){
    XtOmegaX_subset(i,i) = XtOmegaX(i + 1, i + 1);
  }
  
  if(l_star > 0){
    
    double loglikelihood_star = loglikelihood_l_gp_cpp(l_star, a, Y, XtOmegaX_subset, Xtz);
    double loglikelihood_current = loglikelihood_l_gp_cpp(l, a, Y, XtOmegaX_subset, Xtz);
    
    double logprior_star = R::dgamma(l_star, a_l, 1 / b_l, 1);
    double logprior_current = R::dgamma(l, a_l, 1 / b_l, 1);
    
    double logposterior_star = loglikelihood_star + logprior_star;
    double logposterior_current = loglikelihood_current + logprior_current;
    
    if(R::runif(0, 1) < exp(logposterior_star - logposterior_current)){
      
      l = l_star;
      
    }
    
  }
  
  return l;
}

// SAMPLE AS (CLUSTERING SITE VARIABLES)

// [[Rcpp::export]]
arma::vec sample_as_cpp(arma::vec k_s, arma::vec sites, 
                        arma::vec beta_psi, arma::mat X_psi,
                        arma::vec z, arma::vec Omega,
                        double sigma_a){
  
  arma::vec a_s = arma::zeros(k_s.size());
  
  int index_site = 0;
  
  for(int i = 0; (unsigned)i < sites.size(); i++){
    
    int site = sites[i];
    
    int l = 1;
    
    // find rows associated with current site
    if((unsigned)i != (sites.size() - 1)){
      
      while(k_s[index_site + l] == site){
        l += 1;
      }
      
    } else {
      
      l = k_s.size() - index_site;
      
    }
    
    IntegerVector indexes_site(l);
    for(int j = 0; j < l; j++){
      indexes_site[j] = index_site + j;
    }
    index_site += l;  
    
    arma::mat data_psi_site(l, X_psi.n_cols); 
    arma::vec z_site = arma::zeros(l);
    arma::vec w_site = arma::zeros(l);
    for(int j = 0; j < l; j++){
      data_psi_site.row(j) = X_psi.row(indexes_site[j]);
      z_site[j] = z[indexes_site[j]];
      w_site[j] = Omega[indexes_site[j]];
    }  
    
    // sample x
    arma::vec k_site = (z_site - .5);
    arma::vec n_site = arma::ones(l);   
    arma::vec c_site = data_psi_site * beta_psi;
    
    double a = 1 / (sum(w_site) + 1 / (sigma_a*sigma_a));
    double b = sum(k_site - w_site % c_site);
    
    double x = R::rnorm(b * a, sqrt(a));
    
    for(int j = 0; j < l; j++){
      a_s[indexes_site[j]] = x;
    }
    
  }
  
  return(a_s);
}

// [[Rcpp::export]]
double loglikelihood_binom(double z, double a, double c){
  return(z * (a + c) - log(1 + exp(a + c)));
}

// [[Rcpp::export]]
NumericVector GenerateWeightsDP(arma::vec& data, IntegerVector& clusters,
                                IntegerVector& nec, arma::vec& a_s_clusters,
                                arma::vec& c_site, arma::vec& w_site, arma::vec& k_site,
                                double alpha, double sigma_as){
  int k = clusters.length();
  
  NumericVector weights(k + 1);
  
  for(int j = 0; j < clusters.length(); j++){

    int n_ij = nec[clusters[j] - 1];

    double a_s = a_s_clusters[clusters[j] - 1];

    double loglikelihood = 0;
    for(int i = 0; (unsigned)i < data.size(); i++){
      loglikelihood += loglikelihood_binom(data[i], c_site[i], a_s);
    }

    weights[j] = n_ij * exp(loglikelihood);
  }
  
  double term1 = 1 / (sigma_as * sqrt(2 * M_PI));
  arma::vec c_site_squared = c_site % c_site;
  double kc = arma::as_scalar(k_site.t() * c_site);
  double wcc = arma::as_scalar(w_site.t() * c_site_squared);
  double term2 = exp(kc - .5 * wcc);
  double term3 = 1 / sqrt((1 / (sigma_as * sigma_as) + sum(w_site)));
  double kpwc = sum(k_site + w_site % c_site);
  double term4 = exp(.5 * kpwc / (1 / (sigma_as * sigma_as)) + sum(w_site));
  double loglikelihood = term1 * term2 * term3 * term4;

  weights(k) = alpha * loglikelihood;
  
  weights = weights / sum(weights);
  
  return weights;
}

// [[Rcpp::export]]
NumericVector GenerateWeightsDP_neal(arma::vec& data, IntegerVector& clusters,
                                     IntegerVector& nec, arma::vec& a_s_clusters,
                                     arma::vec& c_site, double alpha, double sigma_as, int M){
  int k = clusters.length();
  
  NumericVector weights(k + 1);
  
  for(int j = 0; j < clusters.length(); j++){

    int n_ij = nec[clusters[j] - 1];

    double a_s = a_s_clusters[clusters[j] - 1];

    double loglikelihood = 0;
    for(int i = 0; (unsigned)i < data.size(); i++){
      loglikelihood += loglikelihood_binom(data[i], c_site[i], a_s);
    }

    weights[j] = n_ij * exp(loglikelihood);
  }
  
  NumericVector weights_new(M);
  for(int j = 0; j < M; j++){
    
    double a_s = R::rnorm(0, sigma_as);
    
    double loglikelihood = 0;
    for(int i = 0; (unsigned)i < data.size(); i++){
      loglikelihood += loglikelihood_binom(data[i], c_site[i], a_s);
    }
    
    weights_new[j] = (alpha / M) * exp(loglikelihood);
  }
  
  weights[k] = sum(weights_new);

  weights = weights / sum(weights);
  
  return weights;
}

// [[Rcpp::export]]
List sample_as_cpp_clustering_old(IntegerVector c_s,
                              arma::vec k_s, arma::vec sites,
                              arma::vec beta_psi, arma::mat X_psi,
                              arma::vec z, arma::vec Omega,
                              double sigma_as, double alpha_dp){
  
  // create vector of cluster sizes
  IntegerVector n_clusters(sites.size());
  for(int i = 0; i < c_s.size(); i++){
    n_clusters[c_s[i] - 1] += 1;
  }  
  
  // create vector containing all the clusters labels
  IntegerVector clusters1(sites.size());
  int l = 0;
  for(int i = 0; i < n_clusters.size(); i++){
    if(n_clusters[i] > 0){
      clusters1[l] = i + 1;
      l += 1;
    }
  }
  
  // resize it to the right size
  IntegerVector clusters(l);
  for(int i = 0; i < l; i++){
    clusters[i] = clusters1[i];
  }
  
  // sample cluster locations
  arma::vec a_s_clusters = arma::zeros(sites.size()); // value of random effects for each cluster
  Rcout << clusters << std::endl;
  for(int j = 0; j < clusters.length(); j++){
    
    // find rows associated with cluster
    IntegerVector indexes_cluster(k_s.n_rows);
    int l = 0;
    for(int k = 0; (unsigned)k < k_s.n_rows; k++){
      if(c_s[k_s[k] - 1] == clusters[j]){
        indexes_cluster[l] = k;
        l += 1;
      }
    }
    
    arma::mat data_psi_cluster(l, X_psi.n_cols);
    arma::vec z_cluster = arma::zeros(l);
    arma::vec w_cluster = arma::zeros(l);Rcout << l << "- ne" << std::endl;
    for(int k = 0; k < l; k++){
      Rcout << k << std::endl;
      data_psi_cluster.row(k) = X_psi.row(indexes_cluster[k]);
      z_cluster[k] = z[indexes_cluster[k]];
      w_cluster[k] = Omega[indexes_cluster[k]];
    }
    
    // sample cluster location
    arma::vec k_cluster = (z_cluster - .5);
    arma::vec c_cluster = data_psi_cluster * beta_psi;
    
    double a = 1 / (sum(w_cluster) + 1 / (sigma_as*sigma_as));
    double b = sum(k_cluster - w_cluster % c_cluster);
    
    a_s_clusters[clusters[j] - 1] = R::rnorm(b * a, sqrt(a));
    
  }
  
  // sample cluster allocations
  int index_site = 0; // keep track of where we are in the array k_s
  arma::vec a_s = arma::zeros(k_s.size()); // value of random effect for each site
  
  for(int i = 0; (unsigned)i < sites.size(); i++){
    
    int site = sites[i];
    
    // find rows associated with current site
    int l = 1;
    
    if((unsigned)i != (sites.size() - 1)){
      
      while(k_s[index_site + l] == site){
        l += 1;
      }
      
    } else {
      
      l = k_s.size() - index_site;
      
    }
    
    // create indexes of observation in the cluster
    IntegerVector indexes_site(l);
    for(int j = 0; j < l; j++){
      indexes_site[j] = index_site + j;
    }
    index_site += l;
    
    arma::mat data_psi_site(l, X_psi.n_cols);
    arma::vec z_site = arma::zeros(l);
    arma::vec w_site = arma::zeros(l);
    for(int j = 0; j < l; j++){
      data_psi_site.row(j) = X_psi.row(indexes_site[j]);
      z_site[j] = z[indexes_site[j]];
      w_site[j] = Omega[indexes_site[j]];
    }
    
    arma::vec k_site = (z_site - .5);
    arma::vec n_site = arma::ones(l);
    arma::vec c_site = data_psi_site * beta_psi;
    
    n_clusters[c_s[i] - 1] -= 1;
    Rcout << n_clusters << std::endl;
    // generate weights
    NumericVector weights = GenerateWeightsDP(z_site, clusters, n_clusters, a_s_clusters,
                                              c_site, w_site, k_site, alpha_dp, sigma_as);
    
    // Find Index New Cluster
    int k;
    for(k = 0; k < c_s.size(); k++){
      if(n_clusters[k] == 0) break;
    }
    
    // create vector of clusters to sample from
    IntegerVector clusters_new(clusters.size() + 1);
    for(k = 0; k < clusters.size(); k++){
      clusters_new[k] = clusters[k];
    }
    clusters_new[clusters.size()] = k + 1;
    
    c_s[i] = as<int>(RcppArmadillo::sample(clusters_new, 1, 1, weights));
    
    // if it is a new cluster
    if (c_s[i] == (k+1)) {

      double a = 1 / (sum(w_site) + 1 / (sigma_as*sigma_as));
      double b = sum(k_site - w_site % c_site);

      a_s_clusters[c_s[i] - 1] = R::rnorm(b * a, sqrt(a));

      clusters = clusters_new;
    }
    
    n_clusters[c_s[i] - 1] += 1;
    
    // update the random effect variable for all the elements in the cluster
    for(int j = 0; j < l; j++){
      a_s[indexes_site[j]] = a_s_clusters[c_s[i] - 1];
    }
    
  }
  
  return(List::create(_["a_s"] = a_s,
                      _["c_s"] = c_s));
}

// [[Rcpp::export]]
arma::vec sample_cluster_locations(IntegerVector c_s,
                              arma::vec k_s, arma::vec sites,
                              arma::vec beta_psi, arma::mat X_psi,
                              arma::vec z, arma::vec Omega,
                              double sigma_as){
  
  // sample cluster locations
  arma::vec a_s_clusters = arma::zeros(sites.size()); // value of random effects for each cluster
  
  for(int j = 0; j < sites.size(); j++){
    
    // find rows associated with cluster
    IntegerVector indexes_cluster(k_s.n_rows);
    int l = 0;
    for(int k = 0; (unsigned)k < k_s.n_rows; k++){
      if(c_s[k_s[k] - 1] == (j + 1)){
        indexes_cluster[l] = k;
        l += 1;
      }
    }
    
    arma::mat data_psi_cluster(l, X_psi.n_cols);
    arma::vec z_cluster = arma::zeros(l);
    arma::vec w_cluster = arma::zeros(l);
    for(int k = 0; k < l; k++){
      data_psi_cluster.row(k) = X_psi.row(indexes_cluster[k]);
      z_cluster[k] = z[indexes_cluster[k]];
      w_cluster[k] = Omega[indexes_cluster[k]];
    }
    
    // sample cluster location
    arma::vec k_cluster = (z_cluster - .5);
    arma::vec c_cluster = data_psi_cluster * beta_psi;
    
    double a = 1 / (sum(w_cluster) + 1 / (sigma_as*sigma_as));
    double b = sum(k_cluster - w_cluster % c_cluster);
    
    a_s_clusters[j] = R::rnorm(b * a, sqrt(a));
    
  }
  
  return(a_s_clusters);
}

// [[Rcpp::export]]
List sample_as_cpp_clustering(IntegerVector c_s,
                              arma::vec k_s, arma::vec sites,
                              arma::vec beta_psi, arma::mat X_psi,
                              arma::vec z, arma::vec Omega,
                              double sigma_as, double alpha_dp, int M_neal){
  
  // create vector of cluster sizes
  IntegerVector n_clusters(sites.size());
  for(int i = 0; i < c_s.size(); i++){
    n_clusters[c_s[i] - 1] += 1;
  }  
  
  // create vector containing all the clusters labels
  IntegerVector clusters1(sites.size());
  int l = 0;
  for(int i = 0; i < n_clusters.size(); i++){
    if(n_clusters[i] > 0){
      clusters1[l] = i + 1;
      l += 1;
    }
  }
  
  // resize it to the right size
  IntegerVector clusters(l);
  for(int i = 0; i < l; i++){
    clusters[i] = clusters1[i];
  }
  
  // sample cluster locations
  arma::vec a_s_clusters = arma::zeros(sites.size()); // value of random effects for each cluster
  
  for(int j = 0; j < clusters.length(); j++){
    
    // find rows associated with cluster
    IntegerVector indexes_cluster(k_s.n_rows);
    int l = 0;
    for(int k = 0; (unsigned)k < k_s.n_rows; k++){
      if(c_s[k_s[k] - 1] == clusters[j]){
        indexes_cluster[l] = k;
        l += 1;
      }
    }
    
    arma::mat data_psi_cluster(l, X_psi.n_cols);
    arma::vec z_cluster = arma::zeros(l);
    arma::vec w_cluster = arma::zeros(l);
    for(int k = 0; k < l; k++){
      data_psi_cluster.row(k) = X_psi.row(indexes_cluster[k]);
      z_cluster[k] = z[indexes_cluster[k]];
      w_cluster[k] = Omega[indexes_cluster[k]];
    }
    
    // sample cluster location
    arma::vec k_cluster = (z_cluster - .5);
    arma::vec c_cluster = data_psi_cluster * beta_psi;
    
    double a = 1 / (sum(w_cluster) + 1 / (sigma_as*sigma_as));
    double b = sum(k_cluster - w_cluster % c_cluster);
    
    a_s_clusters[clusters[j] - 1] = R::rnorm(b * a, sqrt(a));
    
  }
  
  // sample cluster allocations
  int index_site = 0; // keep track of where we are in the array k_s
  arma::vec a_s = arma::zeros(k_s.size()); // value of random effect for each site
  
  for(int i = 0; (unsigned)i < sites.size(); i++){
    
    int site = sites[i];
    
    // find rows associated with current site
    int l = 1;
    
    if((unsigned)i != (sites.size() - 1)){
      
      while(k_s[index_site + l] == site){
        l += 1;
      }
      
    } else {
      
      l = k_s.size() - index_site;
      
    }
    
    // create indexes of observation in the cluster
    IntegerVector indexes_site(l);
    for(int j = 0; j < l; j++){
      indexes_site[j] = index_site + j;
    }
    index_site += l;
    
    arma::mat data_psi_site(l, X_psi.n_cols);
    arma::vec z_site = arma::zeros(l);
    arma::vec w_site = arma::zeros(l);
    for(int j = 0; j < l; j++){
      data_psi_site.row(j) = X_psi.row(indexes_site[j]);
      z_site[j] = z[indexes_site[j]];
      w_site[j] = Omega[indexes_site[j]];
    }
    
    arma::vec k_site = (z_site - .5);
    arma::vec c_site = data_psi_site * beta_psi;
    
    n_clusters[c_s[i] - 1] -= 1;
    
    // generate weights
    NumericVector weights = GenerateWeightsDP_neal(z_site, clusters, n_clusters, a_s_clusters,
                                                   c_site, alpha_dp, sigma_as, M_neal);
    
    // Find Index New Cluster
    int k;
    for(k = 0; k < c_s.size(); k++){
      if(n_clusters[k] == 0) break;
    }
    
    // create vector of clusters to sample from
    IntegerVector clusters_new(clusters.size() + 1);
    for(k = 0; k < clusters.size(); k++){
      clusters_new[k] = clusters[k];
    }
    clusters_new[clusters.size()] = k + 1;
    
    c_s[i] = as<int>(RcppArmadillo::sample(clusters_new, 1, 1, weights));
    
    // if it is a new cluster
    if (c_s[i] == (k+1)) {

      double a = 1 / (sum(w_site) + 1 / (sigma_as*sigma_as));
      double b = sum(k_site - w_site % c_site);

      a_s_clusters[c_s[i] - 1] = R::rnorm(b * a, sqrt(a));

      clusters = clusters_new;
    }
    
    n_clusters[c_s[i] - 1] += 1;
    
    // update the random effect variable for all the elements in the cluster
    for(int j = 0; j < l; j++){
      a_s[indexes_site[j]] = a_s_clusters[c_s[i] - 1];
    }
    
  }
  
  return(List::create(_["a_s"] = a_s,
                      _["c_s"] = c_s));
}

// [[Rcpp::export]]
List sample_as_cpp_clustering_integrated(IntegerVector c_s,
                                         arma::vec k_s, arma::vec sites,
                                         arma::vec beta_psi, arma::mat X_psi,
                                         arma::vec z, arma::vec Omega,
                                         double sigma_as, double alpha_dp, int M_neal){
  
  // create vector of cluster sizes
  IntegerVector n_clusters(sites.size());
  for(int i = 0; i < c_s.size(); i++){
    n_clusters[c_s[i] - 1] += 1;
  }  
  
  // create vector containing all the clusters labels
  IntegerVector clusters1(sites.size());
  int l = 0;
  for(int i = 0; i < n_clusters.size(); i++){
    if(n_clusters[i] > 0){
      clusters1[l] = i + 1;
      l += 1;
    }
  }
  
  // resize it to the right size
  IntegerVector clusters(l);
  for(int i = 0; i < l; i++){
    clusters[i] = clusters1[i];
  }
  
  // sample cluster locations
  arma::vec a_s_clusters = arma::zeros(sites.size()); // value of random effects for each cluster
  
  for(int j = 0; j < clusters.length(); j++){
    
    // find rows associated with cluster
    IntegerVector indexes_cluster(k_s.n_rows);
    int l = 0;
    for(int k = 0; (unsigned)k < k_s.n_rows; k++){
      if(c_s[k_s[k] - 1] == clusters[j]){
        indexes_cluster[l] = k;
        l += 1;
      }
    }
    
    arma::mat data_psi_cluster(l, X_psi.n_cols);
    arma::vec z_cluster = arma::zeros(l);
    arma::vec w_cluster = arma::zeros(l);
    for(int k = 0; k < l; k++){
      data_psi_cluster.row(k) = X_psi.row(indexes_cluster[k]);
      z_cluster[k] = z[indexes_cluster[k]];
      w_cluster[k] = Omega[indexes_cluster[k]];
    }
    
    // sample cluster location
    arma::vec k_cluster = (z_cluster - .5);
    arma::vec c_cluster = data_psi_cluster * beta_psi;
    
    double a = 1 / (sum(w_cluster) + 1 / (sigma_as*sigma_as));
    double b = sum(k_cluster - w_cluster % c_cluster);
    
    a_s_clusters[clusters[j] - 1] = R::rnorm(b * a, sqrt(a));
    
  }
  
  // sample cluster allocations
  int index_site = 0; // keep track of where we are in the array k_s
  arma::vec a_s = arma::zeros(k_s.size()); // value of random effect for each site
  
  for(int i = 0; (unsigned)i < sites.size(); i++){
    
    int site = sites[i];
    
    // find rows associated with current site
    int l = 1;
    
    if((unsigned)i != (sites.size() - 1)){
      
      while(k_s[index_site + l] == site){
        l += 1;
      }
      
    } else {
      
      l = k_s.size() - index_site;
      
    }
    
    // create indexes of observation in the cluster
    IntegerVector indexes_site(l);
    for(int j = 0; j < l; j++){
      indexes_site[j] = index_site + j;
    }
    index_site += l;
    
    arma::mat data_psi_site(l, X_psi.n_cols);
    arma::vec z_site = arma::zeros(l);
    arma::vec w_site = arma::zeros(l);
    for(int j = 0; j < l; j++){
      data_psi_site.row(j) = X_psi.row(indexes_site[j]);
      z_site[j] = z[indexes_site[j]];
      w_site[j] = Omega[indexes_site[j]];
    }
    
    arma::vec k_site = (z_site - .5);
    arma::vec c_site = data_psi_site * beta_psi;
    
    n_clusters[c_s[i] - 1] -= 1;
    
    // generate weights
    NumericVector weights = GenerateWeightsDP_neal(z_site, clusters, n_clusters, a_s_clusters,
                                                   c_site, alpha_dp, sigma_as, M_neal);
    
    // Find Index New Cluster
    int k;
    for(k = 0; k < c_s.size(); k++){
      if(n_clusters[k] == 0) break;
    }
    
    // create vector of clusters to sample from
    IntegerVector clusters_new(clusters.size() + 1);
    for(k = 0; k < clusters.size(); k++){
      clusters_new[k] = clusters[k];
    }
    clusters_new[clusters.size()] = k + 1;
    
    c_s[i] = as<int>(RcppArmadillo::sample(clusters_new, 1, 1, weights));
    
    // if it is a new cluster
    if (c_s[i] == (k+1)) {

      double a = 1 / (sum(w_site) + 1 / (sigma_as*sigma_as));
      double b = sum(k_site - w_site % c_site);

      a_s_clusters[c_s[i] - 1] = R::rnorm(b * a, sqrt(a));

      clusters = clusters_new;
    }
    
    n_clusters[c_s[i] - 1] += 1;
    
    // update the random effect variable for all the elements in the cluster
    for(int j = 0; j < l; j++){
      a_s[indexes_site[j]] = a_s_clusters[c_s[i] - 1];
    }
    
  }
  
  // resample cluster parameters
  
  // create vector containing all the clusters labels
  IntegerVector clusters3(sites.size());
  l = 0;
  for(int i = 0; i < n_clusters.size(); i++){
    if(n_clusters[i] > 0){
      clusters3[l] = i + 1;
      l += 1;
    }
  }
  
  // resize it to the right size
  IntegerVector clusters4(l);
  for(int i = 0; i < l; i++){
    clusters4[i] = clusters3[i];
  }
  
  // sample cluster locations
  arma::vec mu_clusters = arma::zeros(sites.size()); // value of random effects for each cluster
  arma::vec sigma_clusters = arma::zeros(sites.size()); // value of random effects for each cluster
  
  for(int j = 0; j < clusters4.length(); j++){
    
    // find rows associated with cluster
    IntegerVector indexes_cluster(k_s.n_rows);
    int l = 0;
    for(int k = 0; (unsigned)k < k_s.n_rows; k++){
      if(c_s[k_s[k] - 1] == clusters4[j]){
        indexes_cluster[l] = k;
        l += 1;
      }
    }
    
    arma::mat data_psi_cluster(l, X_psi.n_cols);
    arma::vec z_cluster = arma::zeros(l);
    arma::vec w_cluster = arma::zeros(l);
    for(int k = 0; k < l; k++){
      data_psi_cluster.row(k) = X_psi.row(indexes_cluster[k]);
      z_cluster[k] = z[indexes_cluster[k]];
      w_cluster[k] = Omega[indexes_cluster[k]];
    }
    
    // sample summary statistics
    arma::vec k_cluster = (z_cluster - .5);
    arma::vec c_cluster = data_psi_cluster * beta_psi;
    
    double a = 1 / (sum(w_cluster) + 1 / (sigma_as*sigma_as));
    double b = sum(k_cluster - w_cluster % c_cluster);
    
    mu_clusters[clusters4[j] - 1] = b * a;
    sigma_clusters[clusters4[j] - 1] = sqrt(a);
    
  }
  
  return(List::create(_["a_s"] = a_s,
                      _["c_s"] = c_s,
                      _["mu_clusters"] = mu_clusters,
                      _["sigma_clusters"] = sigma_clusters));
}

// [[Rcpp::export]]
List sample_as_cpp_clustering_try(IntegerVector c_s,
                              arma::vec k_s, arma::vec sites,
                              arma::vec beta_psi, arma::mat X_psi,
                              arma::vec z, arma::vec Omega,
                              double sigma_as, double alpha_dp, int M_neal){
  
  // create vector of cluster sizes
  IntegerVector n_clusters(sites.size());
  for(int i = 0; i < c_s.size(); i++){
    n_clusters[c_s[i] - 1] += 1;
  }  
  
  // create vector containing all the clusters labels
  IntegerVector clusters1(sites.size());
  int l = 0;
  for(int i = 0; i < n_clusters.size(); i++){
    if(n_clusters[i] > 0){
      clusters1[l] = i + 1;
      l += 1;
    }
  }
  
  // resize it to the right size
  IntegerVector clusters(l);
  for(int i = 0; i < l; i++){
    clusters[i] = clusters1[i];
  }
  
  // sample cluster locations
  arma::vec a_s_clusters = arma::zeros(sites.size()); // value of random effects for each cluster
  
  for(int j = 0; j < clusters.length(); j++){
    
    // find rows associated with cluster
    IntegerVector indexes_cluster(k_s.n_rows);
    int l = 0;
    for(int k = 0; (unsigned)k < k_s.n_rows; k++){
      if(c_s[k_s[k] - 1] == clusters[j]){
        indexes_cluster[l] = k;
        l += 1;
      }
    }
    
    arma::mat data_psi_cluster(l, X_psi.n_cols);
    arma::vec z_cluster = arma::zeros(l);
    arma::vec w_cluster = arma::zeros(l);
    for(int k = 0; k < l; k++){
      data_psi_cluster.row(k) = X_psi.row(indexes_cluster[k]);
      z_cluster[k] = z[indexes_cluster[k]];
      w_cluster[k] = Omega[indexes_cluster[k]];
    }
    
    // sample cluster location
    arma::vec k_cluster = (z_cluster - .5);
    arma::vec c_cluster = data_psi_cluster * beta_psi;
    
    double a = 1 / (sum(w_cluster) + 1 / (sigma_as*sigma_as));
    double b = sum(k_cluster - w_cluster % c_cluster);
    
    a_s_clusters[clusters[j] - 1] = R::rnorm(b * a, sqrt(a));
    
  }
  
  // sample cluster allocations
  int index_site = 0; // keep track of where we are in the array k_s
  arma::vec a_s = arma::zeros(k_s.size()); // value of random effect for each site
  
  for(int i = 0; (unsigned)i < sites.size(); i++){
    
    int site = sites[i];
    
    // find rows associated with current site
    int l = 1;
    
    if((unsigned)i != (sites.size() - 1)){
      
      while(k_s[index_site + l] == site){
        l += 1;
      }
      
    } else {
      
      l = k_s.size() - index_site;
      
    }
    
    // create indexes of observation in the cluster
    IntegerVector indexes_site(l);
    for(int j = 0; j < l; j++){
      indexes_site[j] = index_site + j;
    }
    index_site += l;
    
    arma::mat data_psi_site(l, X_psi.n_cols);
    arma::vec z_site = arma::zeros(l);
    arma::vec w_site = arma::zeros(l);
    for(int j = 0; j < l; j++){
      data_psi_site.row(j) = X_psi.row(indexes_site[j]);
      z_site[j] = z[indexes_site[j]];
      w_site[j] = Omega[indexes_site[j]];
    }
    
    arma::vec k_site = (z_site - .5);
    arma::vec c_site = data_psi_site * beta_psi;
    
    n_clusters[c_s[i] - 1] -= 1;

    // // generate weights
    // NumericVector weights = GenerateWeightsDP_neal(z_site, clusters, n_clusters, a_s_clusters,
    //                                                c_site, alpha_dp, sigma_as, M_neal);
    // generate weights
    NumericVector weights = 0;//GenerateWeightsDP_neal(z_site, clusters, n_clusters, a_s_clusters,
                                //                   c_site, alpha_dp, sigma_as, M_neal);

    // Find Index New Cluster
    int k;
    for(k = 0; k < c_s.size(); k++){
      if(n_clusters[k] == 0) break;
    }

    // create vector of clusters to sample from
    IntegerVector clusters_new(clusters.size() + 1);
    for(k = 0; k < clusters.size(); k++){
      clusters_new[k] = clusters[k];
    }
    clusters_new[clusters.size()] = k + 1;

    c_s[i] = as<int>(RcppArmadillo::sample(clusters_new, 1, 1, weights));

    // if it is a new cluster
    if (c_s[i] == (k+1)) {

      double a = 1 / (sum(w_site) + 1 / (sigma_as*sigma_as));
      double b = sum(k_site - w_site % c_site);

      a_s_clusters[c_s[i] - 1] = R::rnorm(b * a, sqrt(a));

      clusters = clusters_new;
    }

    n_clusters[c_s[i] - 1] += 1;
    
    // update the random effect variable for all the elements in the cluster
    for(int j = 0; j < l; j++){
      a_s[indexes_site[j]] = a_s_clusters[c_s[i] - 1];
    }
    
  }
  
  return(List::create(_["a_s"] = a_s,
                      _["c_s"] = c_s));
}

// // [[Rcpp::export]]
// double sampleSingleObs(IntegerVector& n_clusters,
//                        IntegerVector& clusters,
//                        arma::vec& a_s,
//                        arma::vec& a_s_clusters, //
//                        IntegerVector& c_s,
//                        arma::vec& k_s, arma::vec& sites,
//                        arma::vec& beta_psi, arma::mat& X_psi,
//                        arma::vec& z, arma::vec& Omega,
//                        double sigma_as, double alpha_dp){
//   
//   int index_site = 0;
//   
//   int i = 0;
//   
//   int site = sites[i];
//   
//   int l = 1;
//   
//   // find rows associated with current site
//   if((unsigned)i != (sites.size() - 1)){
//     
//     while(k_s[index_site + l] == site){
//       l += 1;
//     }
//     
//   } else {
//     
//     l = k_s.size() - index_site;
//     
//   }
//   
//   IntegerVector indexes_site(l);
//   for(int j = 0; j < l; j++){
//     indexes_site[j] = index_site + j;
//   }
//   index_site += l;
//   
//   arma::mat data_psi_site(l, X_psi.n_cols);
//   arma::vec z_site = arma::zeros(l);
//   arma::vec w_site = arma::zeros(l);
//   for(int j = 0; j < l; j++){
//     data_psi_site.row(j) = X_psi.row(indexes_site[j]);
//     z_site[j] = z[indexes_site[j]];
//     w_site[j] = Omega[indexes_site[j]];
//   }
//   
//   arma::vec k_site = (z_site - .5);
//   arma::vec n_site = arma::ones(l);
//   arma::vec c_site = data_psi_site * beta_psi;
//   
//   n_clusters[c_s[i] - 1] -= 1;
//   
//   // generate weights
//   NumericVector weights = GenerateWeightsDP(z_site, clusters, n_clusters, a_s_clusters,
//                                             c_site, w_site, k_site, alpha_dp, sigma_as);
//   
//   // Find Index New Cluster
//   int k;
//   for(k = 0; k < c_s.size(); k++){
//     if(n_clusters[k] == 0) break;
//   }
//   IntegerVector clusters_new(clusters.size() + 1);
//   for(k = 0; k < clusters.size(); k++){
//     clusters_new[k] = clusters[k];
//   }
//   clusters_new[clusters.size()] = k + 1;
//   
//   c_s[i] = as<int>(RcppArmadillo::sample(clusters_new, 1, 1, weights));
//   
//   if (c_s[i] == (k+1)) {
//     
//     double a = 1 / (sum(w_site) + 1 / (sigma_as*sigma_as));
//     double b = sum(k_site - w_site % c_site);
//     
//     a_s_clusters[c_s[i] - 1] = R::rnorm(b * a, sqrt(a));
//     
//     clusters = clusters_new;
//   }
//   
//   n_clusters[c_s[i] - 1] += 1;
//   
//   for(int j = 0; j < l; j++){
//     a_s[indexes_site[j]] = a_s_clusters[c_s[i] - 1];
//   }
//   
//   return(0);
// }
// 
// // [[Rcpp::export]]
// double sampleSingleObs2(IntegerVector& n_clusters,
//                        IntegerVector& clusters,
//                        arma::vec& a_s,
//                        arma::vec& a_s_clusters, //
//                        IntegerVector& c_s,
//                        arma::vec& k_s, arma::vec& sites,
//                        arma::vec& beta_psi, arma::mat& X_psi,
//                        arma::vec& z, arma::vec& Omega,
//                        double sigma_as, double alpha_dp){
//   
//   int index_site = 0;
//   
//   int i = 0;
//   
//   int site = sites[i];
//   
//   int l = 1;
//   
//   // find rows associated with current site
//   if((unsigned)i != (sites.size() - 1)){
// 
//     while(k_s[index_site + l] == site){
//       l += 1;
//     }
// 
//   } else {
// 
//     l = k_s.size() - index_site;
// 
//   }
// 
//   IntegerVector indexes_site(l);
//   for(int j = 0; j < l; j++){
//     indexes_site[j] = index_site + j;
//   }
//   index_site += l;
// 
//   arma::mat data_psi_site(l, X_psi.n_cols);
//   arma::vec z_site = arma::zeros(l);
//   arma::vec w_site = arma::zeros(l);
//   for(int j = 0; j < l; j++){
//     data_psi_site.row(j) = X_psi.row(indexes_site[j]);
//     z_site[j] = z[indexes_site[j]];
//     w_site[j] = Omega[indexes_site[j]];
//   }
// 
//   arma::vec k_site = (z_site - .5);
//   arma::vec n_site = arma::ones(l);
//   arma::vec c_site = data_psi_site * beta_psi;
//   
//   n_clusters[c_s[i] - 1] -= 1;
//   
//   // generate weights
//   NumericVector weights = 0;
// 
//   // Find Index New Cluster
//   int k;
//   for(k = 0; k < c_s.size(); k++){
//     if(n_clusters[k] == 0) break;
//   }
//   IntegerVector clusters_new(clusters.size() + 1);
//   for(k = 0; k < clusters.size(); k++){
//     clusters_new[k] = clusters[k];
//   }
//   clusters_new[clusters.size()] = k + 1;
// 
//   c_s[i] = as<int>(RcppArmadillo::sample(clusters_new, 1, 1, weights));
// 
//   if (c_s[i] == (k+1)) {
// 
//     double a = 1 / (sum(w_site) + 1 / (sigma_as*sigma_as));
//     double b = sum(k_site - w_site % c_site);
// 
//     a_s_clusters[c_s[i] - 1] = R::rnorm(b * a, sqrt(a));
// 
//     clusters = clusters_new;
//   }
// 
//   n_clusters[c_s[i] - 1] += 1;
// 
//   for(int j = 0; j < l; j++){
//     a_s[indexes_site[j]] = a_s_clusters[c_s[i] - 1];
//   }
//   
//   return(0);
// }