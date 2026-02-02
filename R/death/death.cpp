#include <Rcpp.h>

// Pure death process with exact observations y_t = x_t

class Death {
public:
  // (de)constructor
  Death(std::vector<double> y, double x0, double dt, int n_particles,
        int target_success, int max_trials);
  ~Death();
  // fields
  std::vector<double> y;
  int n_obs, n_particles, target_success, max_trials, n_dt;

  // Model parameters
  double theta, dt;

  // log-likelihood of y_{1:T} = (y_1, ..., y_T)
  double lml;

  // initial state (fixed)
  double x0;

  // Simulate from the transition x_t | x_{t-1}
  double SimXt(double x);

  // Different particle filter
  double RunBootstrapFilter(double par);
  double RunAliveFilter(double par);
  double RunFrankenFilter(double par);
  // void RunParticleMCMC(double sd_prior, double sd_proposal, int ndpost);

  double x_prev, sum_ws;

  // Fields for the BF
  std::vector<double> x_particles, x_weights, log_uweights, weights;
  std::vector<int> x_ancestors, ancestors, aux_ancestors;

  // Fields for the particleMCMC
  std::vector<double> draws_theta, draws_lml;
  int accept_rate;
};

Death::Death(std::vector<double> y, double x0, double dt, int n_particles,
             int target_success, int max_trials)
  : y(y), x0(x0), dt(dt), n_particles(n_particles), target_success(target_success), max_trials(max_trials) {
  n_obs = y.size();
  n_dt = static_cast<int>(1.0 / dt);
}

// Destructor
Death::~Death() {}

// Construct the path using Poisson(tau)-leap
double Death::SimXt(double x) {
  for (int i = 0; i < n_dt; i++) x -= R::rpois(x * theta * dt);
  // if (x < 0) x = 0.0;
  return x;
}

double Death::RunBootstrapFilter(double par) {

  // Set parameters
  theta = par;

  // Initialise the log-likelihood
  lml = -(n_obs-1) * log(n_particles);

  // Set the initial values at time t=0
  sum_ws = 0.0;

  for (int t=1; t < n_obs; t++) {
    // std::cout<< t << "\n";
    sum_ws = 0.0;
    // Propagate
    for (int j=0; j < n_particles; j++) {
      // Simulate using the exact observations
      x_prev = SimXt(y[t-1]);
      // Likelihood is an indicator function
      sum_ws += 1.0*(y[t] == x_prev);
    }
    // Extreme case, all observations are 0, return -Inf.
    if (sum_ws == 0.0) return -1000;
    lml += std::log(sum_ws);
  }
  return lml;
}

double Death::RunAliveFilter(double par) {

  // Set parameters
  theta = par;

  // Set the log-likelihood at 0
  lml = 0.0;

  // Auxiliary
  int m, k;
  double log_sf = std::log((double)target_success - 1.0);

  // Start
  for (int t=1; t < n_obs; t++) {
    m = 0;
    k = 0;
    while (m < max_trials && k < target_success) {
      // std::cout << m << "\n";
      m++;
      x_prev = SimXt(y[t-1]);
      k += 1*(y[t] == x_prev);
    }
    // filter collapsed, stop loop
    if (m == max_trials) return -1000;
    else lml += log_sf - std::log((double) m - 1.0);
  }
  return lml;
}

double Death::RunFrankenFilter(double par) {
  // Set parameters
  theta = par;

  // Set the log-likelihood at 0
  lml = 0.0;

  // Auxiliary
  int m, k;

  for (int t=1; t < n_obs; t++) {
    m = 0;
    k = 0;
    while (m < max_trials && k < target_success) {
      m++;
      x_prev = SimXt(y[t-1]);
      k += 1*(y[t] == x_prev);
    }
    if (k < target_success) lml += std::log((double) k) - std::log((double) m);
    else lml += std::log((double) k) - std::log((double) m - 1.0);

  }
  return lml;
}


// Expose above class in R
RCPP_MODULE(Death) {

  Rcpp::class_<Death>("Death")

  .constructor<std::vector<double>, double, double, int, int, double>()

  .method("RunBootstrapFilter", &Death::RunBootstrapFilter)
  .method("RunAliveFilter", &Death::RunAliveFilter)
  .method("RunFrankenFilter", &Death::RunFrankenFilter)

   // parameters fields
  .field("y", &Death::y)

   // BF fields
  .field("lml", &Death::lml)
  // .field("particles", &Death::x_particles)
  // .field("ancestors", &Death::x_ancestors)
  // .field("weights", &Death::x_weights)

   // particleMCMC fields
  // .field("draws", &Death::draws)


;
}











