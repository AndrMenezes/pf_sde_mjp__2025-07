#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

void get_ws_lml(std::vector<double> &lw, std::vector<double> &w, double &s,
                int n) {
  double lw_max = *std::max_element(lw.begin(), lw.end());
  s = 0.0;
  for (int j=0; j < n; j++) {
    w[j] = std::exp(lw[j] - lw_max);
    s += w[j];
  }
  // normalise the weights
  for (int j=0; j < n; j++) w[j] /= s;
  // log-likelihood
  s = lw_max + std::log(s) - std::log(n);
}

int sample_discrete(const std::vector<double> &probs, const int &k) {
  double u = R::unif_rand();
  double cdf = 0.0;
  for(int j=0; j < k; j++) {
    cdf += probs[j];
    if (u < cdf) return(j);
  }
  return(k - 1);
}

std::vector<int> resample_indices(std::vector<int>& ancestors,
                                  const std::vector<double>& probs,
                                  int n_particles) {
  std::vector<int> new_ancestors(n_particles);
  for (int j = 0; j < n_particles; ++j)
    new_ancestors[j] = ancestors[sample_discrete(probs, n_particles)];
  return new_ancestors;
}


class BDI {
public:
  // (de)constructor
  BDI(std::vector<double> y, double x0, double dt, double shape,
      int n_particles, int target_success, int max_trials);
  ~BDI();
  // fields
  std::vector<double> y;
  int n_obs, n_particles, n_dt, max_trials, target_success;

  // Model parameters
  double mu, lambda, gamma, dt, shape;

  // log-likelihood of y_{1:t}
  double lml;

  // Simulate from the initial prior (x_0)
  // double SimX0(); // prior time t = 0
  // initial state (fixed for now)
  double x0;

  // Simulate from the transition x_t | x_{t-1}
  double SimXt(double x);
  void RunBootstrapFilter(std::vector<double> parameters);
  void RunFrankenFilter(std::vector<double> parameters);
  void RunParticleMCMC_mu(double lambda_fix, double gamma_fix, double sd_prior,
                          double sd_proposal, int ndpost);
  void RunParticleMCMCFranken_mu(double lambda_fix, double gamma_fix, double sd_prior,
                          double sd_proposal, int ndpost);
  void RunParticleMCMC_mu_gamma(double lambda_fix, double sd_prior,
                                double sd_proposal, int ndpost);
  // int AncestorSample(int m, int k, std::vector<double> w);

  // void RunParticleMCMC_mu_gamma(std::vector<double> parameters);

  std::vector<double> hazard(double x);

  // index offset to keep the bridge trajectory
  int idx_offset;

  // Fields for the BF
  arma::mat x_bridge;
  std::vector<double> x_particles, x_weights, log_uweights, weights;
  std::vector<int> x_ancestors, ancestors, aux_ancestors;
  double pct_m_max_reached;
  std::vector<int> number_trials;

  // Fields for the particleMCMC
  std::vector<double> draws_mu, draws_gamma, draws_lml;
  int accept_rate;
};

BDI::BDI(std::vector<double> y, double x0, double dt, double shape,
         int n_particles, int target_success, int max_trials)
  : y(y), x0(x0), n_particles(n_particles), dt(dt), shape(shape),
    target_success(target_success), max_trials(max_trials) {
  n_obs = y.size();
  n_dt = static_cast<int>(1.0 / dt);
  idx_offset = (n_dt-1);
}

// Destructor
BDI::~BDI() {}

// Hazard function
std::vector<double> BDI::hazard(double x) {
  return {x * mu, x * lambda, gamma};
}

// Construct the bridge, simulating from Poisson-leap
double BDI::SimXt(double x) { //, int t_cur, int j_cur

  for (int i = 0; i < n_dt; i++) {
    std::vector<double> h = hazard(x);
    std::vector<double> r(3);
    r[0] = R::rpois(h[0] * dt);
    r[1] = R::rpois(h[1] * dt);
    r[2] = R::rpois(h[2] * dt);
    x += r[0] - r[1] + r[2];
    if (x <= 0) { x = 0.001;} //; //std::cout << x << "\n";
    // Save the trajectory
    // int row_index = (t_cur)*n_dt + i - idx_offset;
    // x_bridge(row_index, j_cur) = x;
  }
  return x;
}

void BDI::RunBootstrapFilter(std::vector<double> parameters) {

  // Set parameters
  mu = parameters[0];
  lambda = parameters[1];
  gamma = parameters[2];

  // Aux to keep the information for a given time t
  log_uweights.resize(n_particles);
  weights.resize(n_particles);
  ancestors.resize(n_particles);
  aux_ancestors.resize(n_particles);

  // To keep the particles, ancestors and the weights,
  // Access: x[j * m + i], ith observation jth particle
  x_particles.resize(n_obs * n_particles);
  // x_bridge = arma::zeros<arma::mat>(n_obs * n_dt - idx_offset, n_particles);
  // x_ancestors.resize(n_obs * n_particles);
  // x_weights.resize(n_obs * n_particles);

  // Initialise the ancestors: this tells us which particle is "alive", hence
  // we don't need to copy or change the x_particles
  for (int j=0; j < n_particles; j++) {
    aux_ancestors[j] = j;
    ancestors[j] = j;
  }

  // Set the log-likelihood at 0
  lml = 0.0;
  double l_tmp = 0.0;

  // Set the initial values at time t=0
  for (int j=0; j < n_particles; j++) {
    x_particles[j] = x0;
    log_uweights[j] = R::dgamma(y[0], shape, x0 / shape, 1);
    // weights[j] = (double)1.0/n_particles;
  }

  // Normalise the weights and get the log-likelihood contribution into l_tmp
  get_ws_lml(log_uweights, weights, l_tmp, n_particles);
  // std::cout << l_tmp << "\n";
  lml += l_tmp;
  // Resample the ancestors
  ancestors = resample_indices(aux_ancestors, weights, n_particles);

  // Save the weights and the ancestors
  // for (int j=0; j < n_particles; j++) {
  //   x_weights[j] = (double) 1.0 / n_particles;
  //   x_ancestors[j] = j;//ancestors[j];
  // }
  int idx;

  for (int t=1; t < n_obs; t++) {
    // std::cout<< t << "\n";
    // Propagate
    for (int j=0; j < n_particles; j++) {
      idx = t*n_particles + j;
      // Sample X_t | x^{a_{t-1}}_{t-1}
      x_particles[idx] = SimXt(x_particles[(t-1)*n_particles + ancestors[j]]);
      // Compute particle weights
      log_uweights[j] = R::dgamma(y[t], shape, x_particles[idx] / shape, 1);
    }
    // Normalise the weights and get the log-likelihood contribution into l_tmp
    l_tmp = 0.0;
    get_ws_lml(log_uweights, weights, l_tmp, n_particles);
    // Increase the likelihood
    lml += l_tmp;
    // Resample the ancestors
    ancestors = resample_indices(aux_ancestors, weights, n_particles);

    // Save the weights and the ancestors
    // for (int j=0; j < n_particles; j++) {
    //   x_weights[t*n_particles + j] = weights[j];
    //   x_ancestors[t*n_particles + j] = ancestors[j];
    // }
  }
  // std::cout << lml << "\n";
}

void BDI::RunFrankenFilter(std::vector<double> parameters) {

  // Set parameters
  mu = parameters[0];
  lambda = parameters[1];
  gamma = parameters[2];

  // int max_particles = max_trials - 1;
  // Create container to keep maximum number particles for a given time t
  std::vector<double> x_particles_cur(max_trials);
  std::vector<double> x_particles_prev(max_trials, x0);

  // Ancestors
  // std::vector<int> ancestors(max_trials);

  // Initialise number of trials per observation time
  number_trials.resize(n_obs);

  // Weights
  log_uweights.resize(max_trials);
  weights.resize(max_trials);
  //std::vector<double> ln_w(max_trials, -1000);
  // std::vector<double> w(max_trials);
  for (int j=0; j<max_trials; j++) {
    weights[j] = (double) 1.0/max_trials;
    log_uweights[j] = -5000.0;
  }
  // Auxiliary variables
  int m = 0, anc = 0;
  double k = 0.0, l_tmp = 0.0, ln_sup_g;

  pct_m_max_reached = 0.0;

  // Set log-like to 0
  lml = 0.0;

  // Loop over observations
  for (int t=0; t < n_obs; t++) {
    m=0;
    k=0.0;
    ln_sup_g = R::dgamma(y[t], shape, y[t]/shape, 1);
    while (m < max_trials & k < target_success) {
      // Sample ancestors
      anc = sample_discrete(weights, max_trials);
      // Simulate transition: x_t | x_{t-1}
      x_particles_cur[m] = SimXt(x_particles_prev[anc]);
      // Compute log-weights
      log_uweights[m] = R::dgamma(y[t], shape, x_particles_cur[m]/shape, 1);
      // Compute number of success (k)
      k += exp(log_uweights[m] - ln_sup_g);
      // double inc = exp(log_uweights[m] - ln_sup_g);
      // std::cout << "t: " << t << " m " << m << " " << inc << "\n";// " x[m] " << x_particles_cur[m] << " ln w: " << ln_w[m] << " ln g(x*) " << ln_sup_g << "\n";
      // Increment number of trials
      m++;
    }

    // use only m_t - 1 particles, i.e., set to -1000.
    //if (k == target_success) ln_w[m-1] = -1000;

    // Normalise the weights and compute the log-likelihood contribution
    std::fill(weights.begin(), weights.end(), 0.0);
    l_tmp = 0.0;
    std::vector<double> ln_w(log_uweights.begin() + 1, log_uweights.begin() + m-1);
    get_ws_lml(ln_w, weights, l_tmp, ln_w.size());
    lml += l_tmp;

    // This compute the "correct" log-likelihood when the target success isn't reached
    // if (k < target_success) {
    //   // Normalise the weights and compute the log-likelihood contribution
    //   get_ws_lml(ln_w, w, l_tmp, max_trials);
    //   lml += l_tmp;
    //   //anc_type = 0;
    // } else {
    //   // Normalise the weights and compute the log-likelihood contribution
    //   ln_w[m]=-1000;
    //   get_ws_lml(ln_w, w, l_tmp, max_trials);
    //   lml += l_tmp; // sum(log(w)) for j=1, m_t-1
    //   //anc_type = 1;
    // }

    // Update the particles, by copying _cur into _prev. I think this is not
    // efficient might need to change
    x_particles_prev = x_particles_cur;
    // Set the log-weights equal to -Inf for next time y[t]
    std::fill(log_uweights.begin(), log_uweights.end(), -5000.0);

    // save number of trials
    number_trials[t] = m;
    if (m == max_trials) pct_m_max_reached++;
    //m_prev = m;
  }
  pct_m_max_reached /= n_obs;
}


void BDI::RunParticleMCMC_mu(double lambda_fix, double gamma_fix,
                             double sd_prior, double sd_proposal, int ndpost) {

  // Resise the vector to keep the draws
  draws_lml.resize(ndpost);
  draws_mu.resize(ndpost);
  // Initialise \mu by drawing from the prior
  double log_mu_cur = R::norm_rand()*sd_prior;
  double mu_cur = exp(log_mu_cur);
  double log_mu_prop, mu_prop, lr;
  accept_rate = 0;
  // Initialise the parameter vector of the model
  std::vector<double> parms = {mu_cur, lambda_fix, gamma_fix};
  // Get estimate of log-likelihood
  RunBootstrapFilter(parms);
  double lml_cur = lml;
  // std::cout << lml << " " << parms[0] << " " << parms[1] << " " << parms[2];
  // Start MCMC
  double progress = 0.0;
  for (int k=0; k < ndpost; k++) {
    progress = (double) 100 * k / ndpost;
    Rprintf("%3.2f%% Sampling completed", progress);
    Rprintf("\r");
    // if (k % 100 == 0) std::cout << k << "\n";
    // RW proposal
    log_mu_prop = log_mu_cur + R::norm_rand()*sd_proposal;
    mu_prop = exp(log_mu_prop);
    // Compute estimate of the log-likelihood
    parms[0] = mu_prop;
    RunBootstrapFilter(parms);
    double lml_prop = lml;
    // MH ratio
    lr = lml_prop - lml_cur;
    // std::cout << k << " " << lr << "\n";

    // Add the log-prior
    lr += 0.5*std::pow(log_mu_cur/sd_prior,2.0) - 0.5*std::pow(log_mu_prop/sd_prior,2.0);
    // MH step
    if (log(R::unif_rand()) < lr) {
      mu_cur = mu_prop;
      log_mu_cur = log_mu_prop;
      lml_cur = lml_prop;
      accept_rate++;
    }
    draws_mu[k] = mu_cur;
    draws_lml[k] = lml_cur;
  }
}

void BDI::RunParticleMCMCFranken_mu(double lambda_fix, double gamma_fix,
                                    double sd_prior, double sd_proposal, int ndpost) {

  // Resise the vector to keep the draws
  draws_lml.resize(ndpost);
  draws_mu.resize(ndpost);
  // Initialise \mu by drawing from the prior
  double log_mu_cur = R::norm_rand()*sd_prior;
  double mu_cur = exp(log_mu_cur);
  double log_mu_prop, mu_prop, lr;
  accept_rate = 0;
  // Initialise the parameter vector of the model
  std::vector<double> parms = {mu_cur, lambda_fix, gamma_fix};
  // Get estimate of log-likelihood
  RunFrankenFilter(parms);
  double lml_cur = lml;
  // std::cout << lml << " " << parms[0] << " " << parms[1] << " " << parms[2];
  // Start MCMC
  double progress = 0.0;
  for (int k=0; k < ndpost; k++) {
    progress = (double) 100 * k / ndpost;
    Rprintf("%3.2f%% Sampling completed", progress);
    Rprintf("\r");
    // if (k % 100 == 0) std::cout << k << "\n";
    // RW proposal
    log_mu_prop = log_mu_cur + R::norm_rand()*sd_proposal;
    mu_prop = exp(log_mu_prop);
    // Compute estimate of the log-likelihood
    parms[0] = mu_prop;
    RunFrankenFilter(parms);
    double lml_prop = lml;
    // MH ratio
    lr = lml_prop - lml_cur;
    // std::cout << k << " " << lr << "\n";

    // Add the log-prior
    lr += 0.5*std::pow(log_mu_cur/sd_prior,2.0) - 0.5*std::pow(log_mu_prop/sd_prior,2.0);
    // MH step
    if (log(R::unif_rand()) < lr) {
      mu_cur = mu_prop;
      log_mu_cur = log_mu_prop;
      lml_cur = lml_prop;
      accept_rate++;
    }
    draws_mu[k] = mu_cur;
    draws_lml[k] = lml_cur;
  }
}

// void BDI::RunParticleMCMC_mu_gamma(double lambda_fix, double sd_prior,
//                                    double sd_proposal, int ndpost) {
//
//   // Resise the vector to keep the draws
//   draws_mu.resize(n_particles);
//   draws_gamma.resize(n_particles);
//   // Initialise \mu by drawing from the prior
//   double log_mu_cur = R::norm_rand()*sd_prior;
//   double mu_cur = exp(log_mu_cur);
//   double log_mu_prop, mu_prop, lml_prop, lr;
//   accept_rate = 0;
//   // Initialise the parameter vector of the model
//   std::vector<double> parms = {mu_cur, lambda_fix, gamma_fix};
//   // Get estimate of log-likelihood
//   RunBootstrapFilter(parms);
//   double lml_cur = lml;
//   // Start MCMC
//   for (int k=0; k < ndpost; k++) {
//     // RW proposal
//     log_mu_prop = log_mu_prop + R::norm_rand()*sd_prior;
//     mu_prop = exp(log_mu_prop);
//     // Compute log-likelihood
//     parms[0] = mu_prop;
//     RunBootstrapFilter(parms);
//     lml_prop = lml;
//     // MH ratio
//     lr = lml_prop - lml_cur;
//     // Add the log-prior
//     lr += 0.5*std::pow(log_mu_cur/sd_prior,2) - 0.5*std::pow(log_mu_prop/sd_prior,2);
//     // MH step
//     if (log(R::unif_rand()) < lr) {
//       mu_cur = mu_prop;
//       log_mu_cur = log_mu_prop;
//       lml_cur = lml_prop;
//       accept_rate++;
//     }
//     draws_mu[k] = mu_cur;
//   }
// }

// Expose above class in R
RCPP_MODULE(BDI) {

  Rcpp::class_<BDI>("BDI")

  .constructor<std::vector<double>, double, double, double, int, int, int>()

  .method("RunBootstrapFilter", &BDI::RunBootstrapFilter)
  .method("RunFrankenFilter", &BDI::RunFrankenFilter)
  .method("RunParticleMCMC_mu", &BDI::RunParticleMCMC_mu)
  .method("RunParticleMCMCFranken_mu", &BDI::RunParticleMCMCFranken_mu)

   // parameters fields
  .field("y", &BDI::y)

   // BF fields
  .field("lml", &BDI::lml)
  .field("particles", &BDI::x_particles)
  .field("bridge", &BDI::x_bridge)
  .field("ancestors", &BDI::x_ancestors)
  .field("weights", &BDI::x_weights)
  .field("weights_curr", &BDI::weights)
  .field("number_trials", &BDI::number_trials)
  .field("pct_m_max_reached", &BDI::pct_m_max_reached)

   // particleMCMC fields
  .field("draws_mu", &BDI::draws_mu)


;
}











