// STRIDER: Fused Gromov-Wasserstein Optimal Transport Implementation
// Author: Zaoqu Liu (liuzaoqu@163.com)
// 
// This file implements the Fused Gromov-Wasserstein (FGW) distance computation
// for spatial transcriptomics sample integration.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;
using namespace arma;

// ============================================================================
// Sinkhorn Algorithm for Optimal Transport
// ============================================================================

//' Sinkhorn-Knopp Algorithm for Optimal Transport
//'
//' @description Compute optimal transport plan using Sinkhorn-Knopp algorithm.
//'
//' @param M Cost matrix (n x m)
//' @param p Source distribution (n)
//' @param q Target distribution (m)
//' @param reg Regularization parameter (default: 0.1)
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
    
    // Gibbs kernel
    arma::mat K = exp(-M / reg);
    
    // Initialize dual variables
    arma::vec u = ones<vec>(n);
    arma::vec v = ones<vec>(m);
    
    arma::vec Kv, Ktu;
    double err;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // Update u
        Kv = K * v;
        Kv = arma::clamp(Kv, 1e-300, datum::inf);  // Numerical stability
        u = p / Kv;
        
        // Update v
        Ktu = K.t() * u;
        Ktu = arma::clamp(Ktu, 1e-300, datum::inf);
        v = q / Ktu;
        
        // Check convergence every 10 iterations
        if ((iter + 1) % 10 == 0) {
            arma::mat T = diagmat(u) * K * diagmat(v);
            arma::vec row_sums = sum(T, 1);
            err = norm(row_sums - p, 1);
            
            if (err < tol) {
                break;
            }
        }
    }
    
    // Compute transport plan
    arma::mat T = diagmat(u) * K * diagmat(v);
    
    return T;
}


// ============================================================================
// Gromov-Wasserstein Distance
// ============================================================================

//' Compute Tensor-Matrix Product for GW
//'
//' @description Efficiently compute the tensor-matrix product for GW update.
//' @keywords internal
// [[Rcpp::export]]
arma::mat gwgrad_cpp(const arma::mat& C1, 
                     const arma::mat& C2, 
                     const arma::mat& T) {
    
    // Compute gradient: -4 * C1 * T * C2
    // This is the key operation in GW optimization
    
    arma::mat grad = C1 * T * C2;
    grad *= -4.0;
    
    return grad;
}

//' Initialize Cost Matrix for GW
//'
//' @description Initialize the cost matrix using squared loss.
//' @keywords internal
// [[Rcpp::export]]
arma::mat init_gw_cost_cpp(const arma::mat& C1, 
                           const arma::mat& C2,
                           const arma::vec& p,
                           const arma::vec& q) {
    
    int n = C1.n_rows;
    int m = C2.n_rows;
    
    // f1(a) = a^2, f2(b) = b^2
    // h1(a) = a, h2(b) = 2*b
    
    arma::vec f1_p = sum(square(C1), 1);  // C1^2 * p
    arma::vec f2_q = sum(square(C2), 1);  // C2^2 * q
    
    // Constant term: f1(C1) * p * 1_m^T + 1_n * q^T * f2(C2)^T
    arma::mat constC = f1_p * ones<rowvec>(m) + ones<vec>(n) * f2_q.t();
    
    return constC;
}

//' Gromov-Wasserstein Distance
//'
//' @description Compute Gromov-Wasserstein distance between two metric spaces.
//'
//' @param C1 Distance matrix for source (n x n)
//' @param C2 Distance matrix for target (m x m)
//' @param p Source distribution (n)
//' @param q Target distribution (m)
//' @param loss_fun Loss function: "square" (default) or "kl"
//' @param reg Sinkhorn regularization (default: 0.1)
//' @param max_iter Maximum iterations (default: 100)
//' @param tol Convergence tolerance (default: 1e-7)
//'
//' @return List with transport plan and distance
//' @keywords internal
// [[Rcpp::export]]
List gw_cpp(const arma::mat& C1, 
            const arma::mat& C2,
            const arma::vec& p,
            const arma::vec& q,
            double reg = 0.1,
            int max_iter = 100,
            double tol = 1e-7) {
    
    int n = C1.n_rows;
    int m = C2.n_rows;
    
    // Initialize with product of marginals
    arma::mat T = p * q.t();
    
    // Constant cost matrix
    arma::mat constC = init_gw_cost_cpp(C1, C2, p, q);
    
    double prev_loss = datum::inf;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // Compute gradient
        arma::mat tens = gwgrad_cpp(C1, C2, T);
        
        // Cost matrix for this iteration
        arma::mat M = constC + tens;
        
        // Sinkhorn update
        T = sinkhorn_cpp(M, p, q, reg, 100, 1e-7);
        
        // Compute loss
        double loss = accu(M % T);
        
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
//' @description Compute Fused Gromov-Wasserstein distance combining feature
//'   similarity and structural similarity.
//'
//' @param M Feature distance matrix (n x m)
//' @param C1 Spatial distance matrix for source (n x n)
//' @param C2 Spatial distance matrix for target (m x m)
//' @param p Source distribution (n)
//' @param q Target distribution (m)
//' @param alpha Trade-off parameter (default: 0.5)
//'   - 0: pure OT (only feature similarity)
//'   - 1: pure GW (only structural similarity)
//' @param reg Sinkhorn regularization (default: 0.1)
//' @param max_iter Maximum iterations (default: 100)
//' @param tol Convergence tolerance (default: 1e-7)
//'
//' @return List with transport plan and distance
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
    
    int n = C1.n_rows;
    int m = C2.n_rows;
    
    // Initialize with product of marginals
    arma::mat T = p * q.t();
    
    // Constant GW cost matrix
    arma::mat constC = init_gw_cost_cpp(C1, C2, p, q);
    
    double prev_loss = datum::inf;
    
    for (int iter = 0; iter < max_iter; iter++) {
        // Compute GW gradient
        arma::mat tens = gwgrad_cpp(C1, C2, T);
        
        // Combined cost: (1-alpha)*M + alpha*(GW cost)
        arma::mat cost = (1.0 - alpha) * M + alpha * (constC + tens);
        
        // Sinkhorn update
        T = sinkhorn_cpp(cost, p, q, reg, 100, 1e-7);
        
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
// KL Divergence
// ============================================================================

//' Compute KL Divergence Matrix
//'
//' @description Compute pairwise KL divergence between rows of two matrices.
//'
//' @param P Matrix 1 (n x k)
//' @param Q Matrix 2 (m x k)
//'
//' @return KL divergence matrix (n x m)
//' @export
// [[Rcpp::export]]
arma::mat kl_divergence_matrix_cpp(const arma::mat& P, const arma::mat& Q) {
    
    int n = P.n_rows;
    int m = Q.n_rows;
    
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
    P_norm = clamp(P_norm, 1e-300, 1.0);
    Q_norm = clamp(Q_norm, 1e-300, 1.0);
    
    // Log matrices
    arma::mat logP = log(P_norm);
    arma::mat logQ = log(Q_norm);
    
    // KL(P||Q) = sum_k P_k * (log(P_k) - log(Q_k))
    arma::mat KL = zeros<mat>(n, m);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            double kl = 0.0;
            for (int k = 0; k < (int)P_norm.n_cols; k++) {
                kl += P_norm(i, k) * (logP(i, k) - logQ(j, k));
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
//' @description Fast computation of pairwise Euclidean distance matrix.
//'
//' @param X Matrix (n x d)
//' @param Y Matrix (m x d), optional. If NULL, compute self-distance.
//'
//' @return Distance matrix (n x m or n x n)
//' @export
// [[Rcpp::export]]
arma::mat euclidean_distance_cpp(const arma::mat& X, 
                                  Rcpp::Nullable<Rcpp::NumericMatrix> Y_nullable = R_NilValue) {
    
    arma::mat XX = sum(square(X), 1);
    
    if (Y_nullable.isNull()) {
        // Self-distance
        arma::mat D2 = XX * ones<rowvec>(X.n_rows) + ones<vec>(X.n_rows) * XX.t() - 2 * X * X.t();
        D2 = clamp(D2, 0.0, datum::inf);  // Ensure non-negative
        return sqrt(D2);
    } else {
        Rcpp::NumericMatrix Y_rcpp(Y_nullable);
        arma::mat Y = as<arma::mat>(Y_rcpp);
        arma::mat YY = sum(square(Y), 1);
        
        arma::mat D2 = XX * ones<rowvec>(Y.n_rows) + ones<vec>(X.n_rows) * YY.t() - 2 * X * Y.t();
        D2 = clamp(D2, 0.0, datum::inf);
        return sqrt(D2);
    }
}


// ============================================================================
// SVD-based Rotation
// ============================================================================

//' Compute Optimal Rotation Matrix via SVD
//'
//' @description Compute the optimal rotation matrix to align two point sets
//'   using singular value decomposition (Procrustes analysis).
//'
//' @param X Source points (n x d)
//' @param Y Target points (n x d)
//' @param weights Point weights (n)
//'
//' @return Rotation matrix (d x d)
//' @export
// [[Rcpp::export]]
arma::mat compute_rotation_svd_cpp(const arma::mat& X, 
                                    const arma::mat& Y,
                                    const arma::vec& weights) {
    
    // Weighted cross-covariance matrix
    arma::mat W = diagmat(weights);
    arma::mat H = Y.t() * W * X;
    
    // SVD
    arma::mat U, V;
    arma::vec s;
    svd(U, s, V, H);
    
    // Rotation matrix
    arma::mat R = U * V.t();
    
    // Handle reflection (ensure proper rotation)
    if (det(R) < 0) {
        V.col(V.n_cols - 1) *= -1;
        R = U * V.t();
    }
    
    return R;
}


// ============================================================================
// Center Points
// ============================================================================

//' Center Points by Weighted Mean
//'
//' @description Center point coordinates by subtracting the weighted mean.
//'
//' @param X Points (n x d)
//' @param weights Point weights (n)
//'
//' @return Centered points (n x d)
//' @export
// [[Rcpp::export]]
arma::mat center_points_cpp(const arma::mat& X, const arma::vec& weights) {
    
    // Weighted mean
    arma::rowvec mean = (weights.t() * X);
    mean /= accu(weights);
    
    // Center
    arma::mat X_centered = X.each_row() - mean;
    
    return X_centered;
}
