// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// em_exp_emstep
List em_exp_emstep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_exp_emstep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_exp_emstep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_gamma_emstep
List em_gamma_emstep(NumericVector params, List data, int divide, double eps);
RcppExport SEXP _Rsrat_em_gamma_emstep(SEXP paramsSEXP, SEXP dataSEXP, SEXP divideSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type divide(divideSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(em_gamma_emstep(params, data, divide, eps));
    return rcpp_result_gen;
END_RCPP
}
// em_llogis_emstep
List em_llogis_emstep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_llogis_emstep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_llogis_emstep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_llogis_estep
List em_llogis_estep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_llogis_estep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_llogis_estep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_llogis_pllf
double em_llogis_pllf(NumericVector params, List data, double w1);
RcppExport SEXP _Rsrat_em_llogis_pllf(SEXP paramsSEXP, SEXP dataSEXP, SEXP w1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type w1(w1SEXP);
    rcpp_result_gen = Rcpp::wrap(em_llogis_pllf(params, data, w1));
    return rcpp_result_gen;
END_RCPP
}
// em_lnorm_emstep
List em_lnorm_emstep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_lnorm_emstep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_lnorm_emstep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_lxvmax_estep
List em_lxvmax_estep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_lxvmax_estep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_lxvmax_estep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_lxvmax_pllf
double em_lxvmax_pllf(NumericVector params, List data, double w1);
RcppExport SEXP _Rsrat_em_lxvmax_pllf(SEXP paramsSEXP, SEXP dataSEXP, SEXP w1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type w1(w1SEXP);
    rcpp_result_gen = Rcpp::wrap(em_lxvmax_pllf(params, data, w1));
    return rcpp_result_gen;
END_RCPP
}
// em_lxvmin_estep
List em_lxvmin_estep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_lxvmin_estep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_lxvmin_estep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_lxvmin_pllf
double em_lxvmin_pllf(NumericVector params, List data, double w1);
RcppExport SEXP _Rsrat_em_lxvmin_pllf(SEXP paramsSEXP, SEXP dataSEXP, SEXP w1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type w1(w1SEXP);
    rcpp_result_gen = Rcpp::wrap(em_lxvmin_pllf(params, data, w1));
    return rcpp_result_gen;
END_RCPP
}
// em_pareto_emstep
List em_pareto_emstep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_pareto_emstep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_pareto_emstep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_tlogis_emstep_mo
List em_tlogis_emstep_mo(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_tlogis_emstep_mo(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_tlogis_emstep_mo(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_tlogis_estep
List em_tlogis_estep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_tlogis_estep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_tlogis_estep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_tlogis_pllf
double em_tlogis_pllf(NumericVector params, List data, double w0, double w1);
RcppExport SEXP _Rsrat_em_tlogis_pllf(SEXP paramsSEXP, SEXP dataSEXP, SEXP w0SEXP, SEXP w1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< double >::type w1(w1SEXP);
    rcpp_result_gen = Rcpp::wrap(em_tlogis_pllf(params, data, w0, w1));
    return rcpp_result_gen;
END_RCPP
}
// em_tnorm_emstep
List em_tnorm_emstep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_tnorm_emstep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_tnorm_emstep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_txvmax_emstep_mo
List em_txvmax_emstep_mo(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_txvmax_emstep_mo(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_txvmax_emstep_mo(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_txvmax_estep
List em_txvmax_estep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_txvmax_estep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_txvmax_estep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_txvmax_pllf
double em_txvmax_pllf(NumericVector params, List data, double w0, double w1);
RcppExport SEXP _Rsrat_em_txvmax_pllf(SEXP paramsSEXP, SEXP dataSEXP, SEXP w0SEXP, SEXP w1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< double >::type w1(w1SEXP);
    rcpp_result_gen = Rcpp::wrap(em_txvmax_pllf(params, data, w0, w1));
    return rcpp_result_gen;
END_RCPP
}
// em_txvmin_estep
List em_txvmin_estep(NumericVector params, List data);
RcppExport SEXP _Rsrat_em_txvmin_estep(SEXP paramsSEXP, SEXP dataSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    rcpp_result_gen = Rcpp::wrap(em_txvmin_estep(params, data));
    return rcpp_result_gen;
END_RCPP
}
// em_txvmin_pllf
double em_txvmin_pllf(NumericVector params, List data, double w0, double w1);
RcppExport SEXP _Rsrat_em_txvmin_pllf(SEXP paramsSEXP, SEXP dataSEXP, SEXP w0SEXP, SEXP w1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type w0(w0SEXP);
    Rcpp::traits::input_parameter< double >::type w1(w1SEXP);
    rcpp_result_gen = Rcpp::wrap(em_txvmin_pllf(params, data, w0, w1));
    return rcpp_result_gen;
END_RCPP
}
// dgumbel
NumericVector dgumbel(NumericVector x, double loc, double scale, bool log, bool min);
RcppExport SEXP _Rsrat_dgumbel(SEXP xSEXP, SEXP locSEXP, SEXP scaleSEXP, SEXP logSEXP, SEXP minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type loc(locSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    Rcpp::traits::input_parameter< bool >::type min(minSEXP);
    rcpp_result_gen = Rcpp::wrap(dgumbel(x, loc, scale, log, min));
    return rcpp_result_gen;
END_RCPP
}
// pgumbel
NumericVector pgumbel(NumericVector q, double loc, double scale, bool lower, bool log, bool min);
RcppExport SEXP _Rsrat_pgumbel(SEXP qSEXP, SEXP locSEXP, SEXP scaleSEXP, SEXP lowerSEXP, SEXP logSEXP, SEXP minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type q(qSEXP);
    Rcpp::traits::input_parameter< double >::type loc(locSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    Rcpp::traits::input_parameter< bool >::type min(minSEXP);
    rcpp_result_gen = Rcpp::wrap(pgumbel(q, loc, scale, lower, log, min));
    return rcpp_result_gen;
END_RCPP
}
// qgumbel
NumericVector qgumbel(NumericVector p, double loc, double scale, bool lower, bool log, bool min);
RcppExport SEXP _Rsrat_qgumbel(SEXP pSEXP, SEXP locSEXP, SEXP scaleSEXP, SEXP lowerSEXP, SEXP logSEXP, SEXP minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type loc(locSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type lower(lowerSEXP);
    Rcpp::traits::input_parameter< bool >::type log(logSEXP);
    Rcpp::traits::input_parameter< bool >::type min(minSEXP);
    rcpp_result_gen = Rcpp::wrap(qgumbel(p, loc, scale, lower, log, min));
    return rcpp_result_gen;
END_RCPP
}
// rgumbel
NumericVector rgumbel(int n, double loc, double scale, bool min);
RcppExport SEXP _Rsrat_rgumbel(SEXP nSEXP, SEXP locSEXP, SEXP scaleSEXP, SEXP minSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type loc(locSEXP);
    Rcpp::traits::input_parameter< double >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< bool >::type min(minSEXP);
    rcpp_result_gen = Rcpp::wrap(rgumbel(n, loc, scale, min));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_Rsrat_em_exp_emstep", (DL_FUNC) &_Rsrat_em_exp_emstep, 2},
    {"_Rsrat_em_gamma_emstep", (DL_FUNC) &_Rsrat_em_gamma_emstep, 4},
    {"_Rsrat_em_llogis_emstep", (DL_FUNC) &_Rsrat_em_llogis_emstep, 2},
    {"_Rsrat_em_llogis_estep", (DL_FUNC) &_Rsrat_em_llogis_estep, 2},
    {"_Rsrat_em_llogis_pllf", (DL_FUNC) &_Rsrat_em_llogis_pllf, 3},
    {"_Rsrat_em_lnorm_emstep", (DL_FUNC) &_Rsrat_em_lnorm_emstep, 2},
    {"_Rsrat_em_lxvmax_estep", (DL_FUNC) &_Rsrat_em_lxvmax_estep, 2},
    {"_Rsrat_em_lxvmax_pllf", (DL_FUNC) &_Rsrat_em_lxvmax_pllf, 3},
    {"_Rsrat_em_lxvmin_estep", (DL_FUNC) &_Rsrat_em_lxvmin_estep, 2},
    {"_Rsrat_em_lxvmin_pllf", (DL_FUNC) &_Rsrat_em_lxvmin_pllf, 3},
    {"_Rsrat_em_pareto_emstep", (DL_FUNC) &_Rsrat_em_pareto_emstep, 2},
    {"_Rsrat_em_tlogis_emstep_mo", (DL_FUNC) &_Rsrat_em_tlogis_emstep_mo, 2},
    {"_Rsrat_em_tlogis_estep", (DL_FUNC) &_Rsrat_em_tlogis_estep, 2},
    {"_Rsrat_em_tlogis_pllf", (DL_FUNC) &_Rsrat_em_tlogis_pllf, 4},
    {"_Rsrat_em_tnorm_emstep", (DL_FUNC) &_Rsrat_em_tnorm_emstep, 2},
    {"_Rsrat_em_txvmax_emstep_mo", (DL_FUNC) &_Rsrat_em_txvmax_emstep_mo, 2},
    {"_Rsrat_em_txvmax_estep", (DL_FUNC) &_Rsrat_em_txvmax_estep, 2},
    {"_Rsrat_em_txvmax_pllf", (DL_FUNC) &_Rsrat_em_txvmax_pllf, 4},
    {"_Rsrat_em_txvmin_estep", (DL_FUNC) &_Rsrat_em_txvmin_estep, 2},
    {"_Rsrat_em_txvmin_pllf", (DL_FUNC) &_Rsrat_em_txvmin_pllf, 4},
    {"_Rsrat_dgumbel", (DL_FUNC) &_Rsrat_dgumbel, 5},
    {"_Rsrat_pgumbel", (DL_FUNC) &_Rsrat_pgumbel, 6},
    {"_Rsrat_qgumbel", (DL_FUNC) &_Rsrat_qgumbel, 6},
    {"_Rsrat_rgumbel", (DL_FUNC) &_Rsrat_rgumbel, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_Rsrat(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
