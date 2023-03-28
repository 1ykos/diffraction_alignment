#include "wmath.hpp"
#include "wmath_color.hpp"
#include "wmath_optimisation.hpp"
#include "patchmap.hpp"

using dlib::abs;
using dlib::cholesky_decomposition;
using dlib::diag;
using dlib::diagm;
using dlib::dot;
using dlib::eigenvalue_decomposition;
using dlib::identity_matrix;
using dlib::inv;
using dlib::length;
using dlib::length_squared;
using dlib::make_symmetric;
using dlib::matrix;
using dlib::matrix_exp;
using dlib::matrix_op;
using dlib::normalize;
using dlib::op_make_symmetric;
using dlib::round;
using dlib::set_colm;
using dlib::tmp;
using dlib::trace;
using dlib::zeros_matrix;
using std::abs;
using std::array;
using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::fill;
using std::fixed;
using std::floor;
using std::get;
using std::getline;
using std::ifstream;
using std::isnan;
using std::istream;
using std::lower_bound;
using std::map;
using std::max_element;
using std::numeric_limits;
using std::ofstream;
using std::pow;
using std::round;
using std::setprecision;
using std::setw;
using std::sort;
using std::copy;
using std::stod;
using std::streamsize;
using std::nth_element;
using std::string;
using std::stringstream;
using std::swap;
using std::tuple;
using std::unordered_map;
using std::vector;
using wmath::long_mul;
using wmath::clz;
using wmath::digits;
using wmath::hsv_rgb;

using whash::patchmap;

const double pi = 3.14159265;

int main(int argc,char **argv) {
  vector<vector<matrix<double,2,1>>> features(2);
  for (size_t i=0;i!=2;++i) {
    ifstream file(argv[i+1]);
    for (string line;getline(file,line);) {
      stringstream ss(line);
      double x0,x1,v;
      ss >> x0 >> x1 >> v;
      matrix<double,2,1> x{x0,x1};
      features[i].push_back(x);
    }
  }
  const size_t n = max(features[0].size(),features[1].size());
  matrix<int64_t> cost(n,n);
  for (size_t j=0;j!=n;++j) {
    for (size_t i=0;i!=n;++i) {
      if ((i<features[0].size())&&(j<features[1].size())) {
        cost(i,j) = -16*sqrt(length(features[0][i]-features[1][j]));
        // yes really square root, to bias the assignment towards closer
        // features
      } else {
        cost(i,j) = -2147483647ll;
      }
    }
  }
  const auto assignment = max_cost_assignment(cost);
  // TODO given assignments find the best shift
  return 0;
}
