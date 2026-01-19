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
  BDI(std::vector<double> y, double x0, int n_particles, double dt);
  ~BDI();
  // fields
  std::vector<double> y;
  int n_obs, n_particles, n_dt;

  // Model parameters
  double mu, lambda, gamma, dt;

  // log-likelihood of y_{1:t}
  double lml;

  // Simulate from the initial prior (x_0)
  // double SimX0(); // prior time t = 0
  // initial state (fixed for now)
  double x0;

  // Simulate from the transition x_t | x_{t-1}
  double SimXt(double x, int t_cur, int j_cur);
  void RunBootstrapFilter(std::vector<double> parameters);

  std::vector<double> hazard(double x);


  // Fields for the BF
  arma::mat x_path;
  std::vector<double> x_particles, x_weights, log_uweights, weights;
  std::vector<int> x_ancestors, ancestors, aux_ancestors;

};

BDI::BDI(std::vector<double> y, double x0, int n_particles, double dt)
  : y(y), x0(x0), n_particles(n_particles), dt(dt) {
  n_obs = y.size();
  n_dt = static_cast<int>(1.0 / dt);

}

// Destructor
BDI::~BDI() {}

// double BDI::SimX0() { return 0.0}

std::vector<double> BDI::hazard(double x) {
  return {x * mu, x * lambda, gamma};
}

// Simulate the path trajectory
double BDI::SimXt(double x, int t_cur, int j_cur) {

  for (int i = 0; i < n_dt; i++) {
    std::vector<double> h = hazard(x);
    std::vector<double> r(3);
    r[0] = R::rpois(h[0] * dt);
    r[1] = R::rpois(h[1] * dt);
    r[2] = R::rpois(h[2] * dt);
    x += r[0] - r[1] + r[2];
    // Save the path
    int row_index = (t_cur-1)*n_dt + i;
    x_path(row_index, j_cur) = x;
  }
  return x;
}

void BDI::RunBootstrapFilter(std::vector<double> parameters) {

  // Set parameters
  mu = parameters[0];
  lambda = parameters[1];
  gamma = parameters[2];

  // Aux to keep the information for a given t
  log_uweights.resize(n_particles);
  weights.resize(n_particles);
  ancestors.resize(n_particles);
  aux_ancestors.resize(n_particles);

  // To keep the particles, ancestors and the weights
  x_path = arma::zeros<arma::mat>(n_obs * n_dt, n_particles);
  x_particles.resize((n_obs+1) * n_particles);
  x_ancestors.resize((n_obs+1) * n_particles);
  x_weights.resize((n_obs+1) * n_particles);
  // Access: x[j * m + i], ith observation jth particle

  // Initialise the ancestors: this tells us which particle is "alive", hence
  // we don't need to copy or change the x_particles
  for (int j=0; j < n_particles; j++) aux_ancestors[j] = j;

  double l_tmp = 0.0;
  lml = 0.0;

  // At time t=0 simulate from the prior x_0 and compute its weight using the likelihood
  for (int j=0; j < n_particles; j++) {
    x_particles[j] = x0;
    // x_particles[j] = SimXt(x0, 0, j);
    log_uweights[j] = R::dpois(y[0], x_particles[j], 1);
  }
  // Normalise the weights and get the log-likelihood contribution into l_tmp
  get_ws_lml(log_uweights, weights, l_tmp, n_particles);
  lml += l_tmp;
  // Resample the ancestors
  ancestors = resample_indices(aux_ancestors, weights, n_particles);

  // Save the weights and the ancestors
  for (int j=0; j < n_particles; j++) {
    x_weights[j] = (double) 1.0 / n_particles;
    x_ancestors[j] = ancestors[j];
  }

  for (int t=1; t <= n_obs; t++) {
    // std::cout<< t << "\n";
    // Propagate
    for (int j=0; j < n_particles; j++) {
      int idx = t*n_particles + j;
      //std::cout << idx << "\n";
      // std::cout << x_particles[(t-1)*n_particles + ancestors[j]] << "\n";
      x_particles[idx] = SimXt(x_particles[(t-1)*n_particles + ancestors[j]], t, j);
      // Compute particle weights
      log_uweights[j] = R::dpois(y[t], x_particles[idx], 1);
    }
    // Normalise the weights and get the log-likelihood contribution into l_tmp
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
RCPP_MODULE(BDI) {

  Rcpp::class_<BDI>("BDI")

  .constructor<std::vector<double>, double, int, double>()
//std::vector<double> , double , int ,  dt


  .method("RunBootstrapFilter", &BDI::RunBootstrapFilter)

   // parameters fields
  .field("y", &BDI::y)

   // BF fields
  .field("lml", &BDI::lml)
  .field("particles", &BDI::x_particles)
  .field("path", &BDI::x_path)
  .field("ancestors", &BDI::x_ancestors)
  .field("weights", &BDI::x_weights)
  .field("weights_curr", &BDI::weights)

;
}











