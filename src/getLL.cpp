#include <Rcpp.h>
#include "LLfunctions.h"
using namespace Rcpp;

// Helper FUNCTIONS
// get_nb_cpp (if var_model ==2 for nb)
// get_pois_countVec
// get_pois_muVec

// MAIN FUNCTION
//' Function getLL
//' For experiments with master library, initial infection, AND depleted samples
//' sequenced.
//' Call get_nb and get_poiss within cpp code, returning a likelihood value to optimize.
//' @param gene_essentiality NumericVector (for flexibility) by-guide essentiality to test.
//' @param guide_efficiency NumericVector by-guide efficiency.
//' @param sample_effects NumericVector by-sample effects on gene_essentiality.
//' @param init_counts NumericMatrix
//' @param dep_counts NumericMatrix
//' @param var_model int specifying which count mean~variation model to use.
//' @param master_freq NumericMatrix; log fraction of cells infected by each guide;
//'         rows correspond to guides, columns to different masterlibraries.
//' @param masterlib_key NumericVector of columns of master_freq to use per sample.
//' @param cells_infected NumericVector number of cells infected per sample in experiment.
//' @param init_scaling NumericVector of scaling factors to normalize total
//'         read depth of initial sample.
//' @param dep_scaling NumericVector of scaling factors to normalize total
//'         read depth and LFC ratio shifts due to guide drop-out.
//' @param nsg_vals Vector of percent of cells infected (unobserved, integrating over.)
//' @param var_params List of mean~variance model parameters (fit in ModelObj).
//' @export
// [[Rcpp::export]]
double getLL(NumericVector gene_essentiality, NumericVector guide_efficiency,
             NumericVector sample_effects,
             NumericMatrix init_counts, NumericMatrix dep_counts, int var_model,
             NumericMatrix master_freq, NumericVector masterlib_key,
             NumericVector cells_infected,
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
  double guide_effect, infected_prior, init_seq_count, dep_seq_count, q_max;
  int i = 0;

  // eval L(f) = Poiss(nsg*τs; mg*τs)NB(xsg; nsg*λs)NB(ysg; nsg*λs*[(fg*eg)+(1-eg)])
  // extract & match relative vectors.
  for (int guide = 0; guide < guide_efficiency.length(); guide++) {
    for (int sample = 0; sample < sample_effects.length(); sample++) {
      guide_effect = 1 - guide_efficiency[guide] + guide_efficiency[guide] *
        gene_essentiality[guide] * sample_effects[sample];
      infected_guides = nsg_vals;
      // master_freq is log.
      int masterlib_idx = masterlib_key[sample];
      infected_prior = exp(master_freq(guide, masterlib_idx) + log(cells_infected[sample]));
      init_seq_mean = nsg_vals * init_scaling[sample];
      dep_seq_mean = nsg_vals * guide_effect * dep_scaling[sample];
      init_seq_count = init_counts(guide, sample);
      dep_seq_count = dep_counts(guide, sample);

      // get infection likelihood (poisson in all models; beta & gamma dist suggested).
      inf_ll = get_pois_countVec(infected_guides, infected_prior);

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
      //NumericVector::iterator infLlMax = std::max_element(inf_ll.begin(), inf_ll.end());
      //NumericVector::iterator initLlMax = std::max_element(init_ll.begin(), init_ll.end());
      //NumericVector::iterator depLlMax = std::max_element(dep_ll.begin(), dep_ll.end());
      //int depMax = which_max(dep_ll);
      //int initMax = which_max(init_ll);
      //int infMax = which_max(inf_ll);

      //Rcout << infected_prior << ", " << infected_guides[infMax] << ", " << nsg_vals[infMax] << std::endl;
      //Rcout << init_seq_mean[initMax] << ", " << init_seq_count << ", " << nsg_vals[initMax] << std::endl;
      //Rcout << dep_seq_mean[depMax] <<", " << dep_seq_count << ", " << nsg_vals[depMax] << std::endl;

      if (any(is_na(inf_ll))) {
        Rcout << "infected likelihood is na." << std::endl;
        Rcout << guide << std::endl << sample << std::endl << guide_effect << std::endl;
        NumericVector::iterator infLlMax = std::max_element(inf_ll.begin(), inf_ll.end());
        Rcout << "opt infected_guides, prior, ll are:" << infected_guides[infLlMax - inf_ll.begin()] << std::endl;
        Rcout << infected_prior << std::endl;
        Rcout << inf_ll[infLlMax - inf_ll.begin()] << std::endl;
      }
      if (any(is_na(init_ll))) {
        Rcout << "inital ll is na" << std::endl;
        Rcout << guide << std::endl << sample << std::endl << guide_effect << std::endl;
      }
      if (any(is_na(dep_ll))) {
        Rcout << "depleted ll contains na" << std::endl;
        Rcout << guide << std::endl << sample << std::endl << guide_effect << std::endl;
      }

      ll_per_nsg = inf_ll + init_ll + dep_ll + log(step_size);

      // sum over all possible number of infected cells (nsg_vals)
      q_max = max(ll_per_nsg);
      // Rcout << "q_max is " << q_max << std::endl;
      ll_per_datum[i++] = q_max + log(sum(exp(ll_per_nsg - q_max)));
    }
  }
  // Rcout << "ll per sample per guide for this gene is: " << ll_per_datum << std::endl;
  return sum(ll_per_datum);
}


////// END OF SCRIPT
