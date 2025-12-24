#include <Rcpp.h>

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


class StateSpaceAR1 {
public:
  // (de)constructor
  StateSpaceAR1(std::vector<double> y, std::vector<double> parameters,
                std::vector<double> initial_ss, int n_particles);
  ~StateSpaceAR1();
  // fields
  std::vector<double> y;
  int n_obs, n_particles;
  // Model parameters
  double phi, sigma, tau, sigma2, tau2;

  // log-likelihood of y_{1:t}
  double lml;

  // Parameters (mean and variance) of the current state
  double mean_state, var_state, sd_initial, mean_initial;

  // Functions to run the KF and BF
  void RunKalmanFilter();
  // Vector to keep the mean and variance of the states
  std::vector<double> mean_ss;
  std::vector<double> var_ss;
  std::vector<double> var_pred;


  // void BootstrapFilter();

  // Simulate from the initial prior (x_0), the transition x_t | x_{t-1}, and the
  // observation y_{t} | x_{t}
  double SimX0(); // prior time t = 0
  double SimXt(double x); // Transition x_t | x_{t-1}
  void RunBootstrapFilter();
  // State filtering
  std::vector<double> x_particles, x_weights, log_uweights, weights;
  std::vector<int> x_ancestors, ancestors, aux_ancestors;


};

StateSpaceAR1::StateSpaceAR1(std::vector<double> y, std::vector<double> parameters,
                             std::vector<double> initial_ss,
                             int n_particles) : y(y), n_particles(n_particles) {
  n_obs = y.size();
  phi = parameters[0];
  sigma = parameters[1];
  sigma2 = sigma*sigma;
  tau = parameters[2];
  tau2 = tau*tau;

  mean_initial = initial_ss[0];
  sd_initial = initial_ss[1];

  mean_state = mean_initial;
  var_state = sd_initial*sd_initial;

}

StateSpaceAR1::~StateSpaceAR1() {}

void StateSpaceAR1::RunKalmanFilter() {

  // Initialise the objects to keep the mean and variance of state
  mean_ss.resize(n_obs);
  var_ss.resize(n_obs);
  var_pred.resize(n_obs);

  // Run the Kalman filter
  for (int t=0; t < n_obs; t++) {

    std::cout << t << "\n";

    // Evolve to the prior x_{t} | y_{1:{t-1}}
    double at = phi * mean_state;
    double rt = std::pow(phi, 2.0) * var_state + tau2;
    // Predictive y_{t} | y_{1:{t-1}} : same mean, but variance increases
    double qt = rt + sigma2;
    // Posterior x_{t} | y_{1:t} : Kalman update
    double kg = rt / qt;
    mean_state = at + kg * (y[t] - at);
    var_state = (1 - kg) * rt;//rt - std::pow(rt, 2.0) / qt;
    // Save the state and observational mean and variances
    mean_ss[t] = mean_state;
    var_ss[t] = var_state;
    var_pred[t] = qt;
  }
}

double StateSpaceAR1::SimX0() {
  return R::norm_rand()*sd_initial + mean_initial;
}

double StateSpaceAR1::SimXt(double x) {
  return phi * x + R::norm_rand()*tau ;
}

void StateSpaceAR1::RunBootstrapFilter() {

  // Aux to keep the information for a given t
  log_uweights.resize(n_particles);
  weights.resize(n_particles);
  ancestors.resize(n_particles);
  aux_ancestors.resize(n_particles);

  // To keep the particles, ancestors and the weights
  x_particles.resize(n_obs * n_particles);
  x_ancestors.resize(n_obs * n_particles);
  x_weights.resize(n_obs * n_particles);
  // Access: x[j * m + i], ith observation jth particle

  // Initialise the ancestors: this tells us which particle is "alive", hence
  // we don't need to copy or change the x_particles
  for (int j=0; j < n_particles; j++) aux_ancestors[j] = j;

  double l_tmp = 0.0;
  lml = 0.0;

  // At time t=0 simulate from the prior x_0 and compute its weight using the likelihood
  for (int j=0; j < n_particles; j++) {
    x_particles[j] = SimX0();
    log_uweights[j] = R::dnorm4(y[0], x_particles[j], sigma, 1);
  }

  // Compute weights and the log-likelihood
  get_ws_lml(log_uweights, weights, l_tmp, n_particles);
  lml += l_tmp;
  // Resample the ancestors
  ancestors = resample_indices(aux_ancestors, weights, n_particles);

  // Save the weights and the ancestors
  for (int j=0; j < n_particles; j++) {
    x_weights[j] = weights[j];
    x_ancestors[j] = ancestors[j];
  }

  for (int t=1; t < n_obs; t++) {
    // std::cout<< t << "\n";
    // Propagate
    for (int j=0; j < n_particles; j++) {
      int idx = t*n_particles + j;
      //std::cout << idx << "\n";
      x_particles[idx] = SimXt(x_particles[(t-1)*n_particles + j]);
      // Get the particle-specific weight
      log_uweights[j] = R::dnorm4(y[t], x_particles[idx], sigma, 1);
    }
    // Resample
    get_ws_lml(log_uweights, weights, l_tmp, n_particles);
    // Increase the likelihood
    lml += l_tmp;
    // Resample the ancestors
    ancestors = resample_indices(aux_ancestors, weights, n_particles);

    // Save the weights and the ancestors
    for (int j=0; j < n_particles; j++) {
      x_weights[t*n_particles + j] = weights[j];
      x_ancestors[t*n_particles + j] = ancestors[j];
    }
  }
}


// Expose above class in R
RCPP_MODULE(StateSpaceAR1) {

  Rcpp::class_<StateSpaceAR1>("StateSpaceAR1")

  .constructor<std::vector<double> , std::vector<double> , std::vector<double>, int>()

  .method("RunKalmanFilter", &StateSpaceAR1::RunKalmanFilter)
  .method("RunBootstrapFilter", &StateSpaceAR1::RunBootstrapFilter)

   // parameters fields
  .field("phi", &StateSpaceAR1::phi)
  .field("sigma2", &StateSpaceAR1::sigma2)
  .field("tau2", &StateSpaceAR1::tau2)
  .field("y", &StateSpaceAR1::y)

   // KF fields
  .field("mean_ss", &StateSpaceAR1::mean_ss)
  .field("var_ss", &StateSpaceAR1::var_ss)
  .field("var_pred", &StateSpaceAR1::var_pred)
  .field("mean_state", &StateSpaceAR1::mean_state)
  .field("var_state", &StateSpaceAR1::var_state)

   // BF fields
  .field("lml", &StateSpaceAR1::lml)
  .field("particles", &StateSpaceAR1::x_particles)
  .field("ancestors", &StateSpaceAR1::x_ancestors)
  .field("weights", &StateSpaceAR1::x_weights)
  .field("weights_curr", &StateSpaceAR1::weights)

;
}











