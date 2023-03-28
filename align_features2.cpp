#include "wmath.hpp"
#include "wmath_color.hpp"
#include "wmath_optimisation.hpp"
#include "patchmap.hpp"

#include <dlib/optimization/max_cost_assignment.h>

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
using dlib::max_cost_assignment;
using dlib::assignment_cost;
using std::abs;
using std::max;
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

template< class ForwardIt, class UnaryPredicate >
constexpr ForwardIt move_if(
    ForwardIt first,
    ForwardIt last,
    UnaryPredicate p ) {
  for (auto it=first;it!=last;) {
    if (p(*it)) {
      --last;
      swap(*it,*last);
    } else {
      ++it;
    }
  }
  return last;
}

struct filter_with_pivot{
  const matrix<double,2,1>& x; // the piovt
  const matrix<double,2,1>& d; // the derivative
  const bool operator()(
    const tuple<
      double,
      matrix<double,2,1>,
      matrix<double,2,1>,
      double
    > & v
      ) const {
    const matrix<double,2,1>& y = get<1>(v);
    return (trans(d)*(y-x))>0;
  }
};

tuple<double,matrix<double,2,1>> fast_centerpoint(
    const vector<tuple<double,matrix<double,2,1>>>& points
    ) {
  if (points.size()==0) {
    return {0.0,zeros_matrix<double>(2,1)};
  }
  if (points.size()==1) {
    return *points.begin();
  }
  vector<
    tuple<
      double,
      matrix<double,2,1>,
      matrix<double,2,1>,
      double
    >
  > data;
  for (auto it=points.begin();it!=points.end();++it) {
    const double             w = get<0>(*it);
    const matrix<double,2,1> x = get<1>(*it);
    const matrix<double,2,1> d{0.0,0.0};
    data.emplace_back(w,x,d,0.0);
  }
  auto pivot = data.begin();
  auto pivot_end = data.begin()+1;
  auto end = data.end();
  while (true) {
    get<2>(*pivot) = zeros_matrix<double>(2,1);
    for (auto it=data.begin()+1;it!=data.end();++it) {
      get<2>(*pivot) += get<0>(*it)
                       *(get<1>(*pivot)-get<1>(*it))
                       /length(get<1>(*pivot)-get<1>(*it));
      get<3>(*pivot) +=get<0>(*it)*length(get<1>(*pivot)-get<1>(*it));
    }
    // filter the previous pivot elements and throw out the ones that are
    // not suitable pivots any more because they are clearly not median elements
    pivot_end = move_if(
        data.begin()+1,
        pivot_end,
        filter_with_pivot{get<1>(*data.begin()),get<2>(*data.begin())}
        );
    if (pivot_end==end) break;
    // filter all the other elements with the new pivot
    end = move_if(
        pivot_end,
        end,
        filter_with_pivot{get<1>(*data.begin()),get<2>(*data.begin())}
        );
    if (pivot_end==end) break;
    // move the pivot element to the front
    swap(*pivot_end,*data.begin());
    ++pivot_end;
    pivot = data.begin();
  }
  auto median = *min_element(
      data.begin(),
      pivot_end,[](auto & l,auto & r){
        return get<3>(l)<get<3>(r);
      });
  return {get<0>(median),get<1>(median)};
}

matrix<double,2,1> geometric_median(
    const vector<tuple<double,matrix<double,2,1>>>& points
    ) {
  if (points.size()==0) {
    return zeros_matrix<double>(2,1);
  }
  if (points.size()==1) {
    return get<1>(*points.begin());
  }
  matrix<double,2,1> x = get<1>(fast_centerpoint(points));
  // cout << trans(x) ;
  for (size_t i=0;i!=1024;++i) {
    matrix<double,2,1> d1 = zeros_matrix<double>(2,1);
    matrix<double,2,2> d2 = zeros_matrix<double>(2,2);
    double sum = 0;
    double w = 0;
    double sumw = 0;
    for (auto it=points.begin();it!=points.end();++it) {
      const matrix<double,2,1> y = get<1>(*it);
      const double l2 = length_squared(x-y);
      const double l  = sqrt(l2);
      sum += get<0>(*it)*l;
      sumw+= get<0>(*it);
      if (l<1e-16) {
        w+= get<0>(*it);
        continue;
      }
      d1  += get<0>(*it)*(x-y)/l;
      d2  += get<0>(*it)*(identity_matrix<double>(2)-(x-y)*trans(x-y)/l2)/l;
    }
    if (length(d1)<w) break;
    matrix<double,2,1> s = inv(d2)*(d1-w*d1/length(d1));
    if (length(s)<1e-16) break;
    //s = inv(d2)*(d1-w*s/length(s));
    //if (length(s)<1e-16) break;
    matrix<double,2,1> t = x-s;
    //cout << length(t-x) << endl;
    //cout << trans(t-x);
    bool breakthrough = false;
    while (true) {
      double sum2 = 0;
      for (auto it=points.begin();it!=points.end();++it) {
        const matrix<double,2,1> y = get<1>(*it);
        sum2 += get<0>(*it)*length(t-y);
      }
      if (sum2<sum) {
        x = t;
        break;
      } else {
        t = x-(d1-d1*w/sumw)/sumw;
        d1*=exp(-1);
        //cout << trans(x-t);
      }
      if (length(x-t)<1e-16) {
        breakthrough = true;
        break;
      }
    }
    if (breakthrough) break;
  }
  return x;
}

int main(int argc,char **argv) {
  vector<vector<tuple<double,double,double>>> features(2);
  for (size_t i=0;i!=2;++i) {
    ifstream file(argv[i+1]);
    for (string line;getline(file,line);) {
      stringstream ss(line);
      double x0,x1,v;
      ss >> x0 >> x1 >> v;
      v = -v/(2*5.3);
      features[i].emplace_back(x0,x1,v);
    }
  }
  // cerr << features[0].size() << " " << features[1].size() << endl;
  const size_t n = 2*max(features[0].size(),features[1].size());
  matrix<int64_t> cost(n,n);
  for (size_t j=0;j!=n;++j) {
    for (size_t i=0;i!=n;++i) {
      if ((i<features[0].size())&&(j<features[1].size())) {
        const matrix<double,3,1> x0{
          get<0>(features[0][i]),
          get<1>(features[0][i]),
          get<2>(features[0][i])
        };
        const matrix<double,3,1> x1{
          get<0>(features[1][j]),
          get<1>(features[1][j]),
          get<2>(features[1][j])
        };
        const double dist = length(x0-x1);
        cost(j,i) = floor(-16*dist);
      } else {
        cost(j,i) = -16*8ll;
      }
    }
  }
  const auto assignment = max_cost_assignment(cost);
  /*
  for (size_t j=0;j!=n;++j) {
    const size_t i = assignment[j];
    if ((i<features[0].size())&&(j<features[1].size())) {
      cout << get<0>(features[0][i]) << " "
           << get<1>(features[0][i]) << " "
           << get<2>(features[0][i]) << endl;
      cout << get<0>(features[1][j]) << " "
           << get<1>(features[1][j]) << " "
           << get<2>(features[1][j]) << endl;
      cout << endl;
    }
  }
  */
  vector<tuple<double,matrix<double,2,1>>> points;
  for (size_t j=0;j!=n;++j) {
    const size_t i = assignment[j];
    if ((i<features[0].size())&&(j<features[1].size())) {
      const matrix<double,2,1> x0{
        get<0>(features[0][i]),
        get<1>(features[0][i])
      };
      const matrix<double,2,1> x1{
        get<0>(features[1][j]),
        get<1>(features[1][j])
      };
      const matrix<double,2,1> d = x1-x0;
      const double w = sqrt(get<2>(features[1][j])*get<2>(features[0][i]));
      points.emplace_back(w,d);
    }
  }
  cout << trans(geometric_median(points));
  return 0;
}
