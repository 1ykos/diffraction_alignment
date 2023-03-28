#include "wmath.hpp"

#include <dlib/matrix.h>

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
using std::stod;
using std::streamsize;
using std::string;
using std::stringstream;
using std::swap;
using std::tuple;
using std::unordered_map;
using std::vector;
using wmath::long_mul;
using wmath::clz;
using wmath::digits;

int main() {
  vector<array<double,3>> data;
  for (double x0,x1,v;cin>>x0>>x1>>v;data.push_back({x0,x1,v}));
  // A * x = y
  matrix<double,0,6> A(data.size(),6);
  matrix<double,0,1> y(data.size(),1);
  size_t n=0;
  for (auto [x0,x1,v] : data) {
    A(n,0) = 0.5*x0*x0;
    A(n,1) = 1.0*x0*x1;
    A(n,2) = 0.5*x1*x1;
    A(n,3) = 1.0*x0;
    A(n,4) = 1.0*x1;
    A(n,5) = 1.0;
    y(n) = v;
    ++n;
  }
  const matrix<double,6,1> c = inv(trans(A)*A)*trans(A)*y;
  const matrix<double,2,2> D2 = {c(0),c(1),c(1),c(2)};
  const matrix<double,2,1> D1 = {c(3),c(4)};
  const matrix<double,2,1> m = inv(D2)*D1;
  cout << trans(-m);
}
