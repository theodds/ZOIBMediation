// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// rcpp_hello
List rcpp_hello();
RcppExport SEXP _ZOIBMediation_rcpp_hello() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello());
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_stan_fit4bayes_zoib_mod();

static const R_CallMethodDef CallEntries[] = {
    {"_ZOIBMediation_rcpp_hello", (DL_FUNC) &_ZOIBMediation_rcpp_hello, 0},
    {"_rcpp_module_boot_stan_fit4bayes_zoib_mod", (DL_FUNC) &_rcpp_module_boot_stan_fit4bayes_zoib_mod, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_ZOIBMediation(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
