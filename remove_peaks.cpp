#include <execution>
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
using std::array;
using std::cerr;
using std::cin;
using std::copy;
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
using std::max;
using std::max_element;
using std::nth_element;
using std::numeric_limits;
using std::ofstream;
using std::pow;
using std::random_device;
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
using std::uniform_int_distribution;
using std::unordered_map;
using std::vector;
using wmath::long_mul;
using wmath::clz;
using wmath::digits;
using wmath::hsv_rgb;
using wmath::destructive_weighted_median_variance;

using whash::patchmap;

const double pi = 3.14159265;

void remove_peaks(
  const size_t nfs,
  const size_t nss,
  const uint8_t* mask,
        int32_t* data,
  const double sigma = 3.0
    )
{
  random_device gen{};
  unordered_map<size_t,vector<tuple<float,float>>> accumulate_neighbours;
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      const size_t minfs = fs>2*sigma?fs-2*sigma:0;
      const size_t minss = ss>2*sigma?ss-2*sigma:0;
      const size_t maxfs = fs+2*sigma<nfs?fs+2*sigma+1:nfs; // one more than max
      const size_t maxss = ss+2*sigma<nfs?ss+2*sigma+1:nss; // one more than max
      cerr << minfs << " " << minss << " " << maxfs << " " << maxss << endl;
      for (size_t iss=minss;iss!=maxss;++iss) {
        for (size_t ifs=minfs;ifs!=maxfs;++ifs) {
          const size_t i=nfs*iss+ifs;
          if (mask[i]) continue;
          const size_t j=nfs*ss+fs;
          if (mask[j]) continue;
          const double r2 = pow(1.0*ss-1.0*iss,2)+pow(1.0*fs-1.0*ifs,2);
          if (r2>pow(2*sigma,2)) continue;
          const double e = exp(-0.5*(pow(1.0*ss-1.0*iss,2)
                                    +pow(1.0*fs-1.0*ifs,2))
                                    /pow(sigma,2));
          accumulate_neighbours[j].emplace_back(sqrt(data[i]),e);
        }
      }
    }
  }
  for (auto it =accumulate_neighbours.begin();
            it!=accumulate_neighbours.end();
          ++it) {
    if (mask[it->first]) continue;
    const double median = get<0>(destructive_weighted_median_variance(
        it->second.begin(),
        it->second.end()));
    for (auto it0=it->second.begin();it0!=it->second.end();++it0) {
      get<0>(*it0)=abs(get<0>(*it0)-median);
      //cerr << "# " << get<0>(*it0) << " " << get<1>(*it0) << endl;
    }
    const double mad = get<0>(destructive_weighted_median_variance(
        it->second.begin(),
        it->second.end()));
    //cerr << sqrt(data[it->first]) << " " << median-mad << " " << median+mad << " "
    //     << it->second.size() << endl; 
    if (abs(sqrt(data[it->first])-median)>mad) {
      uniform_int_distribution<int32_t> dist(-mad,mad);
      data[it->first] = pow(median+(mad+dist(gen))/2,2);
    }
  }
}

int main(int argc,char **argv) {
  const size_t nfs = 1032;
  const size_t nss = 1032;
  ifstream input_data(argv[1]);
  int32_t* data = new int32_t[nfs*nss];
  input_data.read(reinterpret_cast<char*>(data),nfs*nss*sizeof(int32_t));
  uint8_t* mask = new uint8_t[nfs*nss];
  for (size_t i=0;i!=nfs*nss;++i) mask[i]|=(data[i]<0);
  remove_peaks(nfs,nss,mask,data,10.0);
  cout.write(reinterpret_cast<char*>(data),nfs*nss*sizeof(int32_t));
}
