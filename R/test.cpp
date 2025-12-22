#include <Rcpp.h>
class World {
public:
  World();
  void set(std::string msg);
  std::string greet();
private:
  std::string msg;
};

World::World(): msg("hello") {};

void World::set(std::string msg) { this->msg = msg;};
std::string World::greet( ) { return msg;};

RCPP_MODULE(yada){

  using namespace Rcpp;
class_<World>("World")
  .constructor()
  .method("greet", &World::greet)
  .method("set", &World::set)
;
}
