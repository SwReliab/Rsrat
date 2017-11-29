#include <Rcpp.h>
#include <R_ext/Lapack.h>

using namespace Rcpp;

//' dggglm
//'
//' Solve the least square problem
//'
//' Solve the following least square problem:
//' \deqn{\min_x || B^{-1} (d - A x) ||_2}
//' \describe{
//'   \item{WORK}{(workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK)).
//'               On exit, if INFO = 0, WORK(1) returns the optimal LWORK.}
//'   \item{LWORK}{(input) INTEGER.
//'               The dimension of the array WORK. LWORK >= max(1,N+M+P).
//'               For optimum performance, LWORK >= M+min(N,P)+max(N,P)*NB,
//'               where NB is an upper bound for the optimal blocksizes for
//'               DGEQRF, SGERQF, DORMQR and SORMRQ.
//'               If LWORK = -1, then a workspace query is assumed; the routine
//'               only calculates the optimal size of the WORK array, returns
//'               this value as the first entry of the WORK array, and no error
//'               message related to LWORK is issued by XERBLA.}
//'   \item{INFO}{(output) INTEGER.
//'     \describe{
//'       \item{= 0}{successful exit.}
//'       \item{< 0}{if INFO = -i, the i-th argument had an illegal value.}
//'       \item{= 1}{the upper triangular factor R associated with A in the
//'                  generalized QR factorization of the pair (A, B) is
//'                  singular, so that rank(A) < M; the least squares
//'                  solution could not be computed.}
//'       \item{= 2}{the bottom (N-M) by (N-M) part of the upper trapezoidal
//'                  factor T associated with B in the generalized QR
//'                  factorization of the pair (A, B) is singular, so that
//'                  rank( A B ) < N; the least squares solution could not
//'                  be computed.}
//'       }}
//' }
//'
//' @param A A matrix A (dimension n,m)
//' @param B A matrix B (dimension n,p)
//' @param d A vector d. dimension n
//' @param NB An integer indicating an upper bound for the optimal blocksizes. The default is 64.
//' @return A list with components;
//' \item{x}{A vector m. dimension m}
//' \item{y}{A vector p. dimension p}
//' \item{info}{A status. See details.}
//'
//' @export
// [[Rcpp::export]]

List dggglm(NumericMatrix A, NumericMatrix B, NumericVector d, const int NB = 64) {
  const int n = A.nrow();
  const int m = A.ncol();
  const int p = B.ncol();
  const int lwork = m + ((n<p)?n:p) + ((n>p)?n:p) * NB;
  int i, info;

  NumericMatrix tmpA = clone(A);
  NumericMatrix tmpB = clone(B);
  NumericVector tmpd = clone(d);

  NumericVector x(m);
  NumericVector y(p);
  NumericVector work(lwork);

  dggglm_(&n, &m, &p, &tmpA[0], &n, &tmpB[0], &n, &tmpd[0],
    &x[0], &y[0], &work[0], &lwork, &info);

  return List::create(
    Named("x") = x,
    Named("y") = y,
    Named("info") = info
  );
}
