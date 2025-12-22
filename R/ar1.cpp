#include <Rcpp.h>

double log_sum_exp(std::vector<double> &x) {
  double x_max = *std::max_element(x.begin(), x.end());
  double sum = 0.0;
  for (double xi : x) {
    sum += std::exp(xi - x_max);
  }
  return x_max + std::log(sum);
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
  std::vector<double> x_filter, log_uweights;


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
    mean_state = at + rt / (rt + sigma2) * (y[t] - at);
    var_state = rt - std::pow(rt, 2.0) / qt;
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

  x_filter.resize(n_obs * n_particles);
  log_uweights.resize(n_particles);
  // Access: x[i * m + j], ith observation jth particle

  lml = 0.0;

  // STOP HERE: NEED TO CHECK HOW TO NORMALISE AND COMPUTE THE LIKELIHOOD

  // At time t=0 simulate from the prior x_0 and compute its weight using the likelihood
  for (int j=0; j < n_particles; j++) {
    x_filter[j] = SimX0();
    log_uweights[j] = R::dnorm4(y[0], x_filter[j], sigma, 1);
  }

  // Compute the log-likelihood
  lml += log_sum_exp(log_uweights);
  // Resample

  for (int t=1; t < n_obs; t++) {

    for (int j=0; j < n_particles; j++) {

    }

  }



}


// Expose above class in R
RCPP_MODULE(StateSpaceAR1) {

  Rcpp::class_<StateSpaceAR1>("StateSpaceAR1")

  .constructor<std::vector<double> , std::vector<double> , std::vector<double>>()

  .method("RunKalmanFilter", &StateSpaceAR1::RunKalmanFilter)

  .field("mean_ss", &StateSpaceAR1::mean_ss)
  .field("var_ss", &StateSpaceAR1::var_ss)
  .field("var_pred", &StateSpaceAR1::var_pred)
  .field("mean_state", &StateSpaceAR1::mean_state)
  .field("var_state", &StateSpaceAR1::var_state)

  .field("phi", &StateSpaceAR1::phi)
  .field("sigma2", &StateSpaceAR1::sigma2)
  .field("tau2", &StateSpaceAR1::tau2)
  .field("y", &StateSpaceAR1::y)

;
}











