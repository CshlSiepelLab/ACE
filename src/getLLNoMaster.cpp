#include <Rcpp.h>
#include "LLfunctions.h"
using namespace Rcpp;

// Helper FUNCTIONS
// get_pois_muVec
// get_nb_cpp (if var_model == 2)

// MAIN FUNCTION
//
//' Function getLLNoMaster
//' Version of getLL.cpp modified to return likelihood without master library
//' information (the infection distribution).
//' Call get_nb and get_poiss within cpp code, returning a likelihood value to optimize.
//' @param gene_essentiality NumericVector (for flexibility) by-guide essentiality to test.
//' @param guide_efficiency NumericVector by-guide efficiency.
//' @param sample_effects NumericVector by-sample effects on gene_essentiality.
//' @param init_counts NumericMatrix
//' @param dep_counts NumericMatrix
//' @param var_model int specifying which mean-var model to use for count data.
//' @param init_scaling NumericVector of scaling factors to normalize total
//'         read depth of initial sample.
//' @param dep_scaling NumericVector of scaling factors to normalize total
//'         read depth and LFC ratio shifts due to guide drop-out.
//' @param nsg_vals Vector of percent of cells infected (unobserved, integrating over.)
//' @param var_params List of mean~variance model parameters (fit in ModelObj).
//' @param step_size Intervals of infected cells for discrete integral evaluation.
//' @export
// [[Rcpp::export]]
double getLLNoMaster(NumericVector gene_essentiality, NumericVector guide_efficiency,
             NumericVector sample_effects,
             NumericMatrix init_counts, NumericMatrix dep_counts, int var_model,
             NumericVector init_scaling, NumericVector dep_scaling,
             NumericVector nsg_vals, NumericVector var_params,
             int step_size) {

  // init.
  NumericVector ll_per_datum = no_init(guide_efficiency.length() * sample_effects.length()),
  infected_guides = no_init(nsg_vals.length()),
  init_seq_mean = no_init(nsg_vals.length()),
  dep_seq_mean = no_init(nsg_vals.length()),
  inf_ll = no_init(nsg_vals.length()),
  init_ll = no_init(nsg_vals.length()),
  dep_ll = no_init(nsg_vals.length()),
  ll_per_nsg = no_init(nsg_vals.length()),
  init_size = no_init(nsg_vals.length()),
  dep_size = no_init(nsg_vals.length());
  double guide_effect, init_seq_count, dep_seq_count, q_max;
  int i = 0;

  // extract & match relative vectors.
  for (int guide = 0; guide < guide_efficiency.length(); guide++) {
    for (int sample = 0; sample < sample_effects.length(); sample++) {
      guide_effect = 1 - guide_efficiency[guide] + guide_efficiency[guide] *
        gene_essentiality[guide] * sample_effects[sample];
      //NumericVector::iterator guide_max = std::max_element(infected_guides.begin(), infected_guides.end());
      init_seq_mean = nsg_vals * init_scaling[sample];
      dep_seq_mean = nsg_vals * guide_effect * dep_scaling[sample];
      init_seq_count = init_counts(guide, sample);
      dep_seq_count = dep_counts(guide, sample);

      //var==mean (poisson)
      if (var_model == 1) {
        init_ll = get_pois_muVec(init_seq_count, init_seq_mean);
        dep_ll = get_pois_muVec(dep_seq_count, dep_seq_mean);
        // use size as 1/k = 1/exp(par1 + par2/mu)
      } else if (var_model == 2) {
        init_size = 1/exp(var_params[0] + var_params[1]*log(1/init_seq_mean));
        dep_size = 1/exp(var_params[0] + var_params[1]*log(1/dep_seq_mean));
        init_ll = get_nb_cpp(init_seq_mean, init_seq_count, init_size);
        dep_ll = get_nb_cpp(dep_seq_mean, dep_seq_count, dep_size);
        // use ribodiff variance model (separately iteratively optimized shape parameter 1/k)
      } else if (var_model == 3) {
        init_size = var_params[0];
        dep_size = var_params[1];
        init_ll = get_nb_cpp(init_seq_mean, init_seq_count, init_size);
        dep_ll = get_nb_cpp(dep_seq_mean, dep_seq_count, dep_size);
      } else {
        Rcout << "Shape parameter for this model not implemented." << std::endl;
      }

      // debugging.
      //        Rcout << nsg_vals.length() << " is nsg_vals length; all vectors should match." << std::endl;
      //        Rcout << "input lengths to muVec were:" << init_seq_mean.length() << std::endl;
      //        Rcout << "initial LL " << init_ll.length() << std::endl;
      //        Rcout << "dep LL " << dep_ll.length() << std::endl;
      //        Rcout << "max LL mean val, count, nsgval fraction for init, & dep LL:" << std::endl;
      //        NumericVector::iterator initLlMax = std::max_element(init_ll.begin(), init_ll.end());
      //        NumericVector::iterator depLlMax = std::max_element(dep_ll.begin(), dep_ll.end());
      //        int depMax = which_max(dep_ll);
      //        int initMax = which_max(init_ll);
      //
      //        Rcout << init_seq_mean[initMax] << ", " << init_seq_count << ", " << nsg_vals[initMax] << std::endl;
      //        Rcout << dep_seq_mean[depMax] <<", " << dep_seq_count << ", " << nsg_vals[depMax] << std::endl;
      //        for (int i = 0; i < 5; i++) {
        //          Rcout << " " << init_ll[i] << " " << dep_ll[i] << std::endl;
        //          Rcout << " " << init_ll[endI] << " " << dep_ll[endI] << std::endl;
        //          Rcout << "====" << std::endl;
        //        }
      if (any(is_na(init_ll))) {
        Rcout << "inital ll is na" << std::endl;
        Rcout << guide << std::endl << sample << std::endl << guide_effect << std::endl;
      }
      if (any(is_na(dep_ll))) {
        Rcout << "depleted ll contains na" << std::endl;
        Rcout << guide << std::endl << sample << std::endl << guide_effect << std::endl;
      }

      ll_per_nsg = init_ll + dep_ll + log(step_size);

      // sum over all possible number of infected cells (nsg_vals)
      q_max = max(ll_per_nsg);
      ll_per_datum[i++] = q_max + log(sum(exp(ll_per_nsg - q_max)));
    }
  }
  // Rcout << "ll per sample per guide for this gene is: " << ll_per_datum << std::endl;
  return sum(ll_per_datum);
}


////// END OF SCRIPT
