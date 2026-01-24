// STRIDER: Fused Gromov-Wasserstein Optimal Transport Implementation
// Author: Zaoqu Liu (liuzaoqu@163.com)
// 
// Implementation of Fused Gromov-Wasserstein (FGW) distance for spatial
// transcriptomics sample integration, based on:
// 
// Vayer, T., et al. (2020). Fused Gromov-Wasserstein distance for structured 
// objects. Algorithms.
//
// Peyré, G., et al. (2016). Gromov-Wasserstein averaging of kernel and 
// distance matrices. ICML.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// ============================================================================
// Sinkhorn Algorithm for Entropic Optimal Transport
// ============================================================================

//' Sinkhorn-Knopp Algorithm
//'
//' @description Compute optimal transport plan using the Sinkhorn-Knopp 
//'   algorithm with entropic regularization.
//'
//' @param M Cost matrix (n x m)
//' @param p Source distribution (n), must sum to 1
//' @param q Target distribution (m), must sum to 1
//' @param reg Entropic regularization parameter (default: 0.1)
//' @param max_iter Maximum iterations (default: 1000)
//' @param tol Convergence tolerance (default: 1e-9)
//'
//' @return Transport plan matrix (n x m)
//' @keywords internal
// [[Rcpp::export]]
arma::mat sinkhorn_cpp(const arma::mat& M, 
                       const arma::vec& p, 
                       const arma::vec& q,
                       double reg = 0.1,
                       int max_iter = 1000,
                       double tol = 1e-9) {
    
    int n = M.n_rows;
    int m = M.n_cols;
    
    // Gibbs kernel: K = exp(-M / reg)
    arma::mat K = exp(-M / reg);
    
    // Initialize scaling vectors
    arma::vec u = ones<vec>(n);
    arma::vec v = ones<vec>(m);
    
    arma::vec Kv, Ktu;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // Save previous u for convergence check
        arma::vec u_prev = u;
        
        // Update u: u = p / (K * v)
        Kv = K * v;
        Kv = clamp(Kv, 1e-300, datum::inf);
        u = p / Kv;
        
        // Update v: v = q / (K^T * u)
        Ktu = K.t() * u;
        Ktu = clamp(Ktu, 1e-300, datum::inf);
        v = q / Ktu;
        
        // Check convergence every 10 iterations
        if ((iter + 1) % 10 == 0) {
            double err = norm(u - u_prev, 1) / (norm(u, 1) + 1e-10);
            if (err < tol) {
                break;
            }
        }
    }
    
    // Compute transport plan: T = diag(u) * K * diag(v)
    arma::mat T = diagmat(u) * K * diagmat(v);
    
    return T;
}


// ============================================================================
// Gromov-Wasserstein Distance Functions
// ============================================================================

//' Compute GW Gradient (Tensor Product)
//'
//' @description Compute the gradient of GW loss w.r.t. transport plan.
//'   For squared loss: L(C1, C2, T) = sum_{ijkl} (C1_ik - C2_jl)^2 * T_ij * T_kl
//'   Gradient: dL/dT_ij = -4 * sum_kl C1_ik * T_kl * C2_jl = -4 * (C1 * T * C2)_ij
//'
//' @param C1 Source distance matrix (n x n)
//' @param C2 Target distance matrix (m x m)
//' @param T Current transport plan (n x m)
//'
//' @return Gradient matrix (n x m)
//' @keywords internal
// [[Rcpp::export]]
arma::mat gwgrad_cpp(const arma::mat& C1, 
                     const arma::mat& C2, 
                     const arma::mat& T) {
    // Gradient for squared loss: -4 * C1 * T * C2
    return -4.0 * C1 * T * C2;
}

//' Initialize Constant Term for GW
//'
//' @description Initialize the constant term in GW loss for squared loss.
//'   constC = f1(C1) * p * 1_m^T + 1_n * q^T * f2(C2)^T
//'   where f1(a) = f2(a) = a^2 for squared loss
//'
//' @param C1 Source distance matrix (n x n)
//' @param C2 Target distance matrix (m x m)
//' @param p Source distribution (n)
//' @param q Target distribution (m)
//'
//' @return Constant cost matrix (n x m)
//' @keywords internal
// [[Rcpp::export]]
arma::mat init_gw_cost_cpp(const arma::mat& C1, 
                           const arma::mat& C2,
                           const arma::vec& p,
                           const arma::vec& q) {
    
    int n = C1.n_rows;
    int m = C2.n_rows;
    
    // For squared loss: f(a) = a^2
    // f1(C1) * p gives sum over k of C1_ik^2 * p_k for each i
    arma::vec f1_C1_p = sum(square(C1), 1) % p;  // Element-wise: (sum_k C1_ik^2) * p_i
    
    // Actually, the constant term is:
    // sum_k C1_ik^2 * p_k (weighted by p over columns of C1^2)
    arma::vec row_term = square(C1) * p;  // n x 1
    arma::vec col_term = square(C2) * q;  // m x 1
    
    // constC_ij = sum_k C1_ik^2 * p_k + sum_l C2_jl^2 * q_l
    arma::mat constC = row_term * ones<rowvec>(m) + ones<vec>(n) * col_term.t();
    
    return constC;
}

//' Gromov-Wasserstein Distance
//'
//' @description Compute Gromov-Wasserstein distance between two metric spaces
//'   using conditional gradient descent with Sinkhorn projections.
//'
//' @param C1 Source distance matrix (n x n)
//' @param C2 Target distance matrix (m x m)
//' @param p Source distribution (n)
//' @param q Target distribution (m)
//' @param reg Entropic regularization (default: 0.1)
//' @param max_iter Maximum iterations (default: 100)
//' @param tol Convergence tolerance (default: 1e-7)
//'
//' @return List with transport plan T and GW distance
//' @keywords internal
// [[Rcpp::export]]
List gw_cpp(const arma::mat& C1, 
            const arma::mat& C2,
            const arma::vec& p,
            const arma::vec& q,
            double reg = 0.1,
            int max_iter = 100,
            double tol = 1e-7) {
    
    // Initialize transport plan with outer product of marginals
    arma::mat T = p * q.t();
    
    // Constant term in GW loss
    arma::mat constC = init_gw_cost_cpp(C1, C2, p, q);
    
    double prev_loss = datum::inf;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // Compute gradient
        arma::mat tens = gwgrad_cpp(C1, C2, T);
        
        // Total cost matrix
        arma::mat M = constC + tens;
        
        // Sinkhorn projection
        T = sinkhorn_cpp(M, p, q, reg, 100, 1e-8);
        
        // Compute GW loss
        double loss = accu(M % T);  // Element-wise product and sum
        
        // Check convergence
        if (std::abs(prev_loss - loss) < tol) {
            break;
        }
        prev_loss = loss;
    }
    
    // Final loss
    arma::mat tens = gwgrad_cpp(C1, C2, T);
    arma::mat M = constC + tens;
    double gw_dist = accu(M % T);
    
    return List::create(
        Named("T") = T,
        Named("gw_dist") = gw_dist
    );
}


// ============================================================================
// Fused Gromov-Wasserstein Distance
// ============================================================================

//' Fused Gromov-Wasserstein Distance
//'
//' @description Compute Fused Gromov-Wasserstein (FGW) distance, which combines
//'   feature-based Wasserstein distance and structure-based Gromov-Wasserstein
//'   distance.
//'
//' @param M Feature distance/cost matrix (n x m)
//' @param C1 Source structural distance matrix (n x n)
//' @param C2 Target structural distance matrix (m x m)
//' @param p Source distribution (n), must sum to 1
//' @param q Target distribution (m), must sum to 1
//' @param alpha Trade-off parameter between 0 and 1 (default: 0.5)
//'   \itemize{
//'     \item alpha = 0: pure Wasserstein (only feature cost)
//'     \item alpha = 1: pure GW (only structural cost)
//'   }
//' @param reg Entropic regularization (default: 0.1)
//' @param max_iter Maximum iterations (default: 100)
//' @param tol Convergence tolerance (default: 1e-7)
//'
//' @return List with:
//'   \itemize{
//'     \item T: Optimal transport plan (n x m)
//'     \item fgw_dist: FGW distance value
//'   }
//'
//' @details
//' The FGW distance is defined as:
//' FGW(p,q) = min_T (1-alpha) * <M, T> + alpha * GW(C1, C2, T)
//' subject to T * 1 = p, T^T * 1 = q
//'
//' @references
//' Vayer, T., et al. (2020). Fused Gromov-Wasserstein distance for structured 
//' objects. Algorithms.
//'
//' @export
// [[Rcpp::export]]
List fused_gromov_wasserstein_cpp(const arma::mat& M,
                                   const arma::mat& C1,
                                   const arma::mat& C2,
                                   const arma::vec& p,
                                   const arma::vec& q,
                                   double alpha = 0.5,
                                   double reg = 0.1,
                                   int max_iter = 100,
                                   double tol = 1e-7) {
    
    // Initialize transport plan
    arma::mat T = p * q.t();
    
    // Constant term for GW
    arma::mat constC = init_gw_cost_cpp(C1, C2, p, q);
    
    double prev_loss = datum::inf;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // GW gradient
        arma::mat tens = gwgrad_cpp(C1, C2, T);
        
        // Combined cost: (1-alpha)*M + alpha*(GW cost)
        arma::mat cost = (1.0 - alpha) * M + alpha * (constC + tens);
        
        // Sinkhorn projection
        T = sinkhorn_cpp(cost, p, q, reg, 100, 1e-8);
        
        // Compute loss
        double fgw_loss = accu(cost % T);
        
        // Check convergence
        if (std::abs(prev_loss - fgw_loss) < tol) {
            break;
        }
        prev_loss = fgw_loss;
    }
    
    // Final loss computation
    arma::mat tens = gwgrad_cpp(C1, C2, T);
    arma::mat cost = (1.0 - alpha) * M + alpha * (constC + tens);
    double fgw_dist = accu(cost % T);
    
    return List::create(
        Named("T") = T,
        Named("fgw_dist") = fgw_dist
    );
}


// ============================================================================
// KL Divergence Computation
// ============================================================================

//' Compute KL Divergence Matrix
//'
//' @description Compute pairwise KL divergence between probability distributions
//'   (rows of P and Q).
//'
//' @param P First matrix (n x k), rows are distributions
//' @param Q Second matrix (m x k), rows are distributions
//'
//' @return KL divergence matrix (n x m) where entry at row i, col j is KL(P_i || Q_j)
//'
//' @details
//' KL(p||q) = sum_k p_k * log(p_k / q_k)
//'
//' @export
// [[Rcpp::export]]
arma::mat kl_divergence_matrix_cpp(const arma::mat& P, const arma::mat& Q) {
    
    int n = P.n_rows;
    int m = Q.n_rows;
    int k = P.n_cols;
    
    // Normalize rows to sum to 1
    arma::mat P_norm = P;
    arma::mat Q_norm = Q;
    
    for (int i = 0; i < n; i++) {
        double sum_p = accu(P.row(i));
        if (sum_p > 0) {
            P_norm.row(i) = P.row(i) / sum_p;
        }
    }
    
    for (int j = 0; j < m; j++) {
        double sum_q = accu(Q.row(j));
        if (sum_q > 0) {
            Q_norm.row(j) = Q.row(j) / sum_q;
        }
    }
    
    // Add small constant for numerical stability
    double eps = 1e-300;
    P_norm = clamp(P_norm, eps, 1.0);
    Q_norm = clamp(Q_norm, eps, 1.0);
    
    // Compute KL divergence matrix
    arma::mat KL = zeros<mat>(n, m);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double kl = 0.0;
            for (int d = 0; d < k; d++) {
                kl += P_norm(i, d) * (log(P_norm(i, d)) - log(Q_norm(j, d)));
            }
            KL(i, j) = kl;
        }
    }
    
    return KL;
}


// ============================================================================
// Distance Matrix Computation
// ============================================================================

//' Compute Euclidean Distance Matrix
//'
//' @description Compute pairwise Euclidean distances between points.
//'
//' @param X First point set (n x d)
//' @param Y_nullable Second point set (m x d), optional. If NULL, computes
//'   pairwise distances within X.
//'
//' @return Distance matrix (n x m or n x n)
//' @export
// [[Rcpp::export]]
arma::mat euclidean_distance_cpp(const arma::mat& X, 
                                  Rcpp::Nullable<Rcpp::NumericMatrix> Y_nullable = R_NilValue) {
    
    // Squared norms of X
    arma::vec XX = sum(square(X), 1);
    
    if (Y_nullable.isNull()) {
        // Self-distance: D[i,j] = ||X[i] - X[j]||
        // = sqrt(||X[i]||^2 + ||X[j]||^2 - 2*X[i]'*X[j])
        int n = X.n_rows;
        arma::mat D2 = XX * ones<rowvec>(n) + ones<vec>(n) * XX.t() - 2.0 * X * X.t();
        D2 = clamp(D2, 0.0, datum::inf);  // Handle numerical errors
        return sqrt(D2);
    } else {
        Rcpp::NumericMatrix Y_rcpp(Y_nullable);
        arma::mat Y = as<arma::mat>(Y_rcpp);
        arma::vec YY = sum(square(Y), 1);
        
        int n = X.n_rows;
        int m = Y.n_rows;
        arma::mat D2 = XX * ones<rowvec>(m) + ones<vec>(n) * YY.t() - 2.0 * X * Y.t();
        D2 = clamp(D2, 0.0, datum::inf);
        return sqrt(D2);
    }
}


// ============================================================================
// Procrustes Analysis
// ============================================================================

//' Compute Optimal Rotation via SVD
//'
//' @description Find the optimal rotation matrix R to align X to Y 
//'   using Procrustes analysis (SVD of cross-covariance).
//'
//' @param X Source points (n x d)
//' @param Y Target points (n x d)
//' @param weights Point weights (n)
//'
//' @return Rotation matrix R (d x d) such that X * R aligns with Y
//'
//' @details
//' Minimizes weighted sum of squared alignment errors.
//' Solution uses SVD: R = U * V' where H = Y' * W * X = U * S * V'
//'
//' @export
// [[Rcpp::export]]
arma::mat compute_rotation_svd_cpp(const arma::mat& X, 
                                    const arma::mat& Y,
                                    const arma::vec& weights) {
    
    // Weighted cross-covariance matrix
    arma::mat W = diagmat(weights);
    arma::mat H = Y.t() * W * X;
    
    // SVD: H = U * S * V^T
    arma::mat U, V;
    arma::vec s;
    svd(U, s, V, H);
    
    // Optimal rotation: R = U * V^T
    arma::mat R = U * V.t();
    
    // Handle reflection (ensure det(R) = 1)
    if (det(R) < 0) {
        V.col(V.n_cols - 1) *= -1;
        R = U * V.t();
    }
    
    return R;
}


//' Center Points by Weighted Mean
//'
//' @description Center point coordinates by subtracting the weighted centroid.
//'
//' @param X Points (n x d)
//' @param weights Point weights (n), should sum to 1
//'
//' @return Centered points (n x d)
//' @export
// [[Rcpp::export]]
arma::mat center_points_cpp(const arma::mat& X, const arma::vec& weights) {
    
    // Weighted centroid: sum_i w_i * X[i]
    arma::rowvec centroid = weights.t() * X;
    centroid /= accu(weights);
    
    // Center by subtracting centroid
    arma::mat X_centered = X.each_row() - centroid;
    
    return X_centered;
}
