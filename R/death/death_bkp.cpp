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


class Death {
public:
  // (de)constructor
  Death(std::vector<double> y, double x0, double dt, int n_particles,
        int target_success, double max_trials);
  ~Death();
  // fields
  std::vector<double> y;
  int n_obs, n_particles, target_success, n_dt;

  // Model parameters
  double theta, dt;
  double max_trials;

  // log-likelihood of y_{1:t}
  double lml;

  // Simulate from the initial prior (x_0)
  // double SimX0(); // prior time t = 0
  // initial state (fixed for now)
  double x0;

  // Simulate from the transition x_t | x_{t-1}
  double SimXt(double x);
  void RunBootstrapFilter(double par);
  void RunAliveFilter(double par);
  void RunFrankenFilter(double par);
  void RunParticleMCMC(double sd_prior, double sd_proposal, int ndpost);

  double hazard(double x);

  // index offset to keep the bridge trajectory
  int idx_offset;

  // Fields for the BF
  std::vector<double> x_particles, x_weights, log_uweights, weights;
  std::vector<int> x_ancestors, ancestors, aux_ancestors;

  // Fields for the particleMCMC
  std::vector<double> draws, draws_lml;
  int accept_rate;
};

Death::Death(std::vector<double> y, double x0, double dt, int n_particles,
             int target_success, double max_trials)
  : y(y), x0(x0), dt(dt), n_particles(n_particles), target_success(target_success), max_trials(max_trials) {
  n_obs = y.size();
  n_dt = static_cast<int>(1.0 / dt);
  idx_offset = (n_dt-1);
}

// Destructor
Death::~Death() {}

// double Death::SimX0() { return 0.0}

double Death::hazard(double x) {
  return theta * x;
}

// Construct the bridge, simulating from Poisson-leap
double Death::SimXt(double x) {
  for (int i = 0; i < n_dt; i++) x -= R::rpois(x * theta * dt);
  if (x < 0) x = 0.0;
  return x;
}

void Death::RunBootstrapFilter(double par) {

  // Set parameters
  theta = par;

  // Aux to keep the information for a given t
  log_uweights.resize(n_particles);
  weights.resize(n_particles);
  ancestors.resize(n_particles);
  aux_ancestors.resize(n_particles);

  // To keep the particles, ancestors and the weights
  x_particles.resize(n_obs * n_particles);
  x_ancestors.resize(n_obs * n_particles);
  x_weights.resize(n_obs * n_particles);
  // Access: x[t * m + j], t-th observation j-th particle

  // Initialise the ancestors: this tells us which particle is "alive", hence
  // we don't need to copy or change the x_particles
  for (int j=0; j < n_particles; j++) aux_ancestors[j] = j;

  // Initialise the log-likelihood
  lml = -(n_obs-1) * log(n_particles);

  // Set the initial values at time t=0
  for (int j=0; j < n_particles; j++) {
    //x_bridge(0, j) = x0;
    x_particles[j] = x0;
    // lml += log(1.0 / n_particles);
    ancestors[j] = j;
    weights[j] = 1.0/n_particles;
    x_weights[j] = 1.0/n_particles;
  }

  double total_weights=0.0;
  for (int t=1; t < n_obs; t++) {
    // std::cout<< t << "\n";
    total_weights = 0.0;
    // Propagate
    for (int j=0; j < n_particles; j++) {
      int idx = t*n_particles + j;
      //std::cout << idx << "\n";
      // std::cout << x_particles[(t-1)*n_particles + ancestors[j]] << "\n";
      x_particles[idx] = SimXt(x_particles[(t-1)*n_particles + ancestors[j]]);
      // Compute particle weights (likelihood is an indicator function)
      weights[j] = 1.0*(y[t] == x_particles[idx]);
      total_weights += weights[j];
    }
    // Normalise the weights
    if (total_weights == 0.0) {
      // extreme case, all observations are 0, leading to issues.
      // std::cout << "total weights is zero\n";
      lml = -1000;
      return;
      // total_weights = 1.0;
      // for (int j=0; j < n_particles; j++) weights[j] = 1.0 / n_particles;
    } else for (int j=0; j < n_particles; j++) weights[j] /= total_weights;
    lml += std::log(total_weights);

    // Resample the ancestors
    ancestors = resample_indices(aux_ancestors, weights, n_particles);

    // Save the weights and the ancestors
    for (int j=0; j < n_particles; j++) {
      x_weights[t*n_particles + j] = weights[j];
      // x_ancestors[t*n_particles + j] = ancestors[j];
    }
  }
  // std::cout << lml << "\n";

}

void Death::RunAliveFilter(double par) {

  // Set parameters
  theta = par;

  // the actual number of particles to keep is target_success - 1
  int n_particles_keep = target_success - 1;

  // Aux to keep the information for a given t
  log_uweights.resize(n_particles_keep);
  weights.resize(n_particles_keep);
  ancestors.resize(n_particles_keep);
  aux_ancestors.resize(n_particles_keep);

  // To keep the particles, ancestors and the weights
  x_particles.resize(n_obs * n_particles_keep);
  x_ancestors.resize(n_obs * n_particles_keep);
  x_weights.resize(n_obs * n_particles_keep);
  // Access: x[j * m + i], ith observation jth particle

  // Set the log-likelihood at 0
  lml = 0.0;

  // Set the initial values at time t=0
  for (int j=0; j < n_particles_keep; j++) {
    // x_bridge(0, j) = x0;
    x_particles[j] = x0;
    // lml += log(1.0 / n_particles_keep);
    // ancestors[j] = j;
  }
  // Auxiliary
  int m, k, s, anc;
  double x_prop;
  double log_sf = std::log((double)target_success - 1.0);
  std::vector<double> probs(n_particles_keep, (double)1.0/n_particles_keep);

  // Start
  for (int t=1; t < n_obs; t++) {
    m = 0;
    k = 0;
    while (m < max_trials && k < target_success) {
      // std::cout << m << "\n";
      m++;
      anc = sample_discrete(probs, n_particles_keep);
      // std::cout << "ancestors " << anc << "\n";
      // x_prop = SimXt(x_particles[(t-1)*n_particles_keep + anc]);
      x_prop = SimXt(y[t-1]);
      // s = 1*(y[t] == x_prop);
      if (y[t] == x_prop) {
        x_particles[t * n_particles_keep + k] = x_prop;
        k++;
      }
      // x_particles[t * n_particles + k] = x_prop;
      // k += s;
    }
    if (m == max_trials) {
      // filter collapsed, stop loop
      lml = -1000;//-INFINITY;
      // std::cout << "Filter collapsed: max_trials reached at time " << t << " with total number of success: " << k << "\n";
      return;
      // throw std::runtime_error(
      //     "Filter collapsed: max_trials reached at time " +
      //       std::to_string(t) +
      //       " with total number of success: " +
      //       std::to_string(k) + "\n"
      //     );
    }
    else lml += log_sf - std::log((double) m - 1.0);

    // Save the weights and the ancestors
    for (int j=0; j < n_particles_keep; j++) {
      x_weights[t*n_particles_keep + j] = weights[j];
      // x_ancestors[t*n_particles_keep + j] = ancestors[j];
    }
  }
}

void Death::RunFrankenFilter(double par) {
  // Set parameters
  theta = par;

  // the maximum number of particles we can use
  int n_particles_max = max_trials - 1;

  // To keep the particles
  x_particles.resize(n_obs * n_particles_max);
  //x_weights.resize(n_obs * n_particles_max);
  weights.resize(n_particles_max);

  // Set the log-likelihood at 0
  lml = 0.0;

  // Initialise particles at intial condition
  for (int j=0; j < n_particles_max; j++) {
    // x_bridge(0, j) = x0;
    x_particles[j] = x0;
    weights[j] = 1.0 / n_particles_max;
    // lml += log(1.0 / n_particles_max);
    // ancestors[j] = j;
  }
  // Auxiliary
  int m, anc;
  double log_sf = std::log((double)target_success - 1.0);
  int total_success;

  for (int t=1; t < n_obs; t++) {
    m = 0;
    total_success = 0;
    while (m < max_trials && total_success < target_success) {
      m++;
      int idx = t * n_particles_max + m;
      // x_particles[idx] = SimXt(x_particles[(t-1)*n_particles + anc]);
      x_particles[idx] = SimXt(y[t-1]);
      // Compute particle weights (likelihood is an indicator function)
      weights[m] = 1.0*(y[t] == x_particles[idx]);
      total_success += weights[m];
    }
    if (total_success < target_success) lml += std::log((double) total_success) - std::log((double) m);
    else lml += std::log((double) total_success) - std::log((double) m - 1.0);

    // Normalise weights
    // for (int j=0; j < n_particles_max; j++) weights[j] /= total_success;

  }

}


// Expose above class in R
RCPP_MODULE(Death) {

  Rcpp::class_<Death>("Death")

  .constructor<std::vector<double>, double, double, int, int, double>()

  .method("RunBootstrapFilter", &Death::RunBootstrapFilter)
  .method("RunAliveFilter", &Death::RunAliveFilter)

   // parameters fields
  .field("y", &Death::y)

   // BF fields
  .field("lml", &Death::lml)
  .field("particles", &Death::x_particles)
  // .field("ancestors", &Death::x_ancestors)
  .field("weights", &Death::x_weights)

   // particleMCMC fields
  .field("draws", &Death::draws)


;
}











