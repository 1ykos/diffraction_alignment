#include "wmath.hpp"
#include "wmath_color.hpp"
#include "wmath_optimisation.hpp"
#include "patchmap.hpp"

// phi=1.61803398874989;print sqrt(phi)*exp(-I*phi*2*pi)
// awk 'BEGIN{a=-0.886157625224292;b=0.811792974752093;x=1.0;y=1.0;for(i=0;i!=256;++i) {print x,y;_x=x*a-y*b;_y=x*b+y*a;x=_x;y=_y}}'

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
using wmath::destructive_weighted_median_variance;

using whash::patchmap;

const double pi = 3.14159265;

int main(int argc,char **argv) {
  const size_t nfs = 1032;
  const size_t nss = 1032;
  double width = 5.0;
  int32_t* data = new int32_t[nfs*nss];
  ifstream(argv[1]).read(reinterpret_cast<char*>(data),nfs*nss*sizeof(int32_t));
  uint8_t* mask = new uint8_t[nfs*nss]();
  for (size_t i=0;i!=nfs*nss;++i) {
    if (data[i]<0) {
      data[i] = 0;
      mask[i] = 1;
    } else {
      data[i] = sqrt(data[i]);
    }
  }
  patchmap<size_t,vector<tuple<double,double>>> tmp;
  for (double mx,my;cin>>mx>>my;) {
    for (auto it=tmp.begin();it!=tmp.end();++it) {
      it->second.clear();
    }
    patchmap<size_t,double> v;
    for (size_t ss=0;ss!=nss;++ss) {
      for (size_t fs=0;fs!=nfs;++fs) {
        const size_t i = ss*nfs+fs;
        if (mask[i]) continue;
        const double d = sqrt(pow(fs-mx,2)+pow(ss-my,2))/width;
        const size_t r = round(d);
        tmp[r].emplace_back(data[i],1.0);
      }
    }
    for (auto it=tmp.begin();it!=tmp.end();++it) {
      v[it->first]=get<0>(destructive_weighted_median_variance(
            it->second.begin(),
            it->second.end()
            ));
      it->second.clear();
    }
    for (size_t i=0;i!=2;++i) {
      for (size_t ss=0;ss!=nss;++ss) {
        for (size_t fs=0;fs!=nfs;++fs) {
          const size_t i = ss*nfs+fs;
          if (mask[i]) continue;
          const double d = sqrt(pow(fs-mx,2)+pow(ss-my,2))/width;
          const size_t f = floor(d);
          const size_t c =  ceil(d);
          if (f==c) {
            tmp[f].emplace_back(data[i]-v[f],1.0);
          } else {
            const double e = v[c]*(d-f)+v[f]*(c-d);
            //cerr << data[i] << " " << e << " " << (d-f) << " " << v[f] << " " << v[c] << endl;
            const double n2 = pow(d-f,2)+pow(c-d,2);
            tmp[c].emplace_back((v[f]-data[i])*(c-d)/n2+data[i],(d-f));
            tmp[f].emplace_back((v[c]-data[i])*(d-f)/n2+data[i],(c-d));
          }
        }
      }
      for (auto it=tmp.begin();it!=tmp.end();++it) {
        v[it->first]+=get<0>(destructive_weighted_median_variance(
              it->second.begin(),
              it->second.end()
              ));
        it->second.clear();
      }
    }
    size_t n = 0;
    double score = 0;
    for (size_t ss=0;ss!=nss;++ss) {
      for (size_t fs=0;fs!=nfs;++fs) {
        const size_t i = ss*nfs+fs;
        if (mask[i]) continue;
        const double d = sqrt(pow(fs-mx,2)+pow(ss-my,2))/width;
        const size_t f = floor(d);
        const size_t c =  ceil(d);
        const double e = v[c]*(d-f)+v[f]*(c-d);
        score += abs(data[i]-e);
        ++n;
      }
    }
    cout << mx << " " << my << " " << score/n << endl;
  }
}
