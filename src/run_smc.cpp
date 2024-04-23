#include <RcppArmadillo.h>
#include <vector>

#include "parallel_utils.h"
#include "smc_classes.h"
#include "particles.h"

using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List  run_smc(
  Rcpp::List data,
  Rcpp::List new_data,
  Rcpp::List model_options,
  Rcpp::List smc_options,
  Rcpp::List compute_options,
  Rcpp::List priors,
  Rcpp::Nullable<arma::mat> pfun_values,
  Rcpp::Nullable<arma::mat> pfun_estimate) {

  SMCData dat{data};
  SMCParameters pars{model_options, compute_options, smc_options};
  Priors pris{priors};
  SMCAugmentation aug{compute_options, smc_options};

  vec alpha_init = randg(pars.n_particles, distr_param(pris.gamma, 1 / pris.lambda));
  Rcpp::List initial_values;
  initial_values["alpha_init"] = alpha_init;
  mat rho_init(dat.n_items, pars.n_particles);
  for(size_t i{}; i < pars.n_particles; i++) {
    rho_init.col(i) = conv_to<vec>::from(randperm(dat.n_items) + 1);
  }
  initial_values["rho_init"] = rho_init;
  initial_values["aug_init"] = R_NilValue;

  std::vector<StaticParticle> particle_vector =
    initialize_particles(initial_values, pars.n_particles, pars.n_particle_filters, dat);

  auto pfun = choose_partition_function(
    dat.n_items, pars.metric, pfun_values, pfun_estimate);
  auto distfun = choose_distance_function(pars.metric);
  auto resampler = choose_resampler(pars.resampling_method);

  size_t T = new_data.size();
  vec ess(T);
  mat alpha_trace(particle_vector.size(), T);
  for(size_t t{}; t < T; t++) {
    Rcpp::Rcout << t << std::endl;
    dat.update(new_data[t]);
    par_for_each(particle_vector.begin(), particle_vector.end(),
                 [&dat](StaticParticle& p){ p.prepare_particle_filter(dat); });
    aug.run_particle_filter(particle_vector, dat, pfun, distfun, resampler, t);

    vec resampling_probs = normalize_probs(particle_vector);
    ess(t) = 1 / accu(pow(resampling_probs, 2));

    if(ess(t) < pars.resampling_threshold) {
      resample(particle_vector, resampling_probs, resampler);
      vec current_alphas(particle_vector.size());
      std::transform(particle_vector.cbegin(), particle_vector.cend(),
                     current_alphas.begin(),
                     [](const StaticParticle& p) { return p.alpha; });
      double alpha_sd = stddev(current_alphas);
      Rcpp::Rcout << "alpha_sd = " << alpha_sd << std::endl;
      uvec initial_alphas = find_unique(current_alphas);
      unsigned int initial_unique_alphas = initial_alphas.size();
      Rcpp::Rcout << "initial unique = " << initial_unique_alphas << std::endl;

      for(size_t i{}; i < pars.mcmc_steps; i++) {
        double acceptance_rate{};
        for(auto& p : particle_vector) {
         RankProposal rho_proposal = pars.rho_proposal_function->propose(p.rho);
         double alpha_proposal = R::rlnorm(std::log(p.alpha), alpha_sd);

         Rcpp::List initial_values = Rcpp::List::create(
           Rcpp::Named("alpha_init") = vec{alpha_proposal, p.alpha},
           Rcpp::Named("rho_init") = join_rows(rho_proposal.rankings, p.rho),
           Rcpp::Named("aug_init") = R_NilValue
         );

         SMCData proposal_dat{data};
         std::vector<StaticParticle> proposal_particle =
           initialize_particles(initial_values, 2, pars.n_particle_filters, proposal_dat);

         for(size_t s{}; s < t; s++) {
           proposal_dat.update(new_data[s]);
           std::for_each(proposal_particle.begin(), proposal_particle.end(),
                         [&proposal_dat](StaticParticle& p){ p.prepare_particle_filter(proposal_dat); });
           aug.run_particle_filter(proposal_particle, proposal_dat, pfun, distfun, resampler, s);
         }

         double log_ratio = proposal_particle[0].marginal_log_likelihood - proposal_particle[1].marginal_log_likelihood +
           std::log(alpha_proposal) - std::log(p.alpha) - std::log(rho_proposal.prob_forward) +
           std::log(rho_proposal.prob_backward) +
           pris.lambda * (p.alpha - alpha_proposal) +
           pris.gamma * (std::log(alpha_proposal) - std::log(p.alpha));

         bool accept = std::log(R::runif(0, 1)) < log_ratio;

         if(accept) {
           p.alpha = alpha_proposal;
           p.rho = rho_proposal.rankings;
         }
         acceptance_rate += accept;
       }
        Rcpp::Rcout << "acceptance rate " << acceptance_rate / particle_vector.size() << std::endl;

        std::transform(particle_vector.cbegin(), particle_vector.cend(),
                       current_alphas.begin(),
                       [](const StaticParticle& p) { return p.alpha; });
        uvec unique_alphas = find_unique(current_alphas);
        unsigned int current_unique_alphas = unique_alphas.size();
        Rcpp::Rcout << "unique particles " << current_unique_alphas << std::endl;
        if(initial_unique_alphas < (2 * current_unique_alphas - particle_vector.size())) {
          Rcpp::Rcout << "stopping the rejuvenation" << std::endl;
          break;
        }
     }
    }

    std::transform(
      particle_vector.cbegin(), particle_vector.cend(), alpha_trace.col(t).begin(),
      [](const StaticParticle& p) { return p.alpha; }
    );
  }

  Rcpp::List particle_history = Rcpp::List::create(
    Rcpp::Named("rho_samples") = wrapup_rho(particle_vector),
    Rcpp::Named("alpha_samples") = wrapup_alpha(particle_vector),
    Rcpp::Named("augmented_rankings") = wrapup_augmented_data(particle_vector),
    Rcpp::Named("ess") = ess,
    Rcpp::Named("alpha_trace") = alpha_trace
  );
  return particle_history;
}
