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
using std::numeric_limits;
using std::ofstream;
using std::pow;
using std::round;
using std::setprecision;
using std::setw;
using std::sort;
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
using wmath::weighted_quartile_value;
using wmath::lower_weighted_quartile_element;
using wmath::primitive_weighted_quartile_value;

using whash::patchmap;

const double pi = 3.14159265;

double radial_symmetry(
  const size_t nfs,
  const size_t nss,
  const uint8_t* mask,
  const int32_t* data,
  const double width,
  const double mx,
  const double my
    )
{
  size_t n = 0;
  n = max(n,size_t(ceil(sqrt(pow(  0-mx,2)+pow(  0-my,2))/width)));
  n = max(n,size_t(ceil(sqrt(pow(nfs-mx,2)+pow(  0-my,2))/width)));
  n = max(n,size_t(ceil(sqrt(pow(  0-mx,2)+pow(nss-my,2))/width)));
  n = max(n,size_t(ceil(sqrt(pow(nfs-mx,2)+pow(nss-my,2))/width)));
  ++n;
  vector<vector<tuple<double,double>>> tmp(n);
  vector<double> v(n,0.0);
  for (size_t i=0;i!=n;++i) tmp[i].reserve(ceil(pow(width,2)*pi*(i+1)));
  for (size_t i=0;i!=16;++i) {
    for (size_t ss=0;ss!=nss;++ss) {
      for (size_t fs=0;fs!=nfs;++fs) {
        const size_t j = ss*nfs+fs;
        if (mask[j]) continue;
        const double d = sqrt(pow(fs-mx,2)+pow(ss-my,2))/width;
        const size_t f = floor(d);
        const size_t c =  ceil(d);
        //if (c>32) continue;
        if (f==c) {
          tmp[f].emplace_back(data[j]-v[f],1.0);
        } else {
          if (i==0) { // move more in first step
            tmp[f].emplace_back(data[j],(c-d));
            tmp[c].emplace_back(data[j],(d-f));
          } else {
            // oh yes this is where the magic happens.
            // This is the smallest change to v[c] and v[f] to make the linear
            // interpolation between the two coincide with data[j]
            const double n2 = pow(d-f,2)+pow(c-d,2);
            tmp[f].emplace_back((v[c]-data[j])*(d-f)/n2+data[j],(c-d));
            tmp[c].emplace_back((v[f]-data[j])*(c-d)/n2+data[j],(d-f));
          }
        }
      }
    }
    for (size_t i=0;i!=n;++i) {
      double sumw = 0;
      for (auto it=tmp[i].begin();it!=tmp[i].end();++it) {
        sumw+=get<1>(*it);
      }
      /*
      v[i]=primitive_weighted_quartile_value(
          tmp[i].begin(),
          sumw*0.5,
          tmp[i].end()
          );
          */
      v[i]=weighted_quartile_value(
            tmp[i].begin(),
            sumw*0.5,
            tmp[i].end());
      tmp[i].clear();
    }
  }
  vector<double> mad;
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      const size_t i = ss*nfs+fs;
      if (mask[i]) continue;
      const double d = sqrt(pow(fs-mx,2)+pow(ss-my,2))/width;
      const size_t f = floor(d);
      const size_t c =  ceil(d);
      double e = 0.0;
      if (f==c) {
        e = v[f];
      } else {
        e = v[c]*(d-f)+v[f]*(c-d);
      }
      mad.push_back(abs(data[i]-e));
    }
  }
  double score = 0;
  constexpr bool skip_peaks = false; // not worth it I think
  if constexpr (skip_peaks) {
    const auto cutoff = mad.begin()+15*(mad.end()-mad.begin())/16;
    nth_element(mad.begin(),cutoff,mad.end());
    n=0;
    for (auto it=mad.begin();it!=cutoff;++it) {
      score+=*it;
      ++n;
    }
    for (auto it=cutoff;it!=mad.end();++it) {
      score+=*cutoff;
      ++n;
    }
  } else {
    for (auto it=mad.begin();it!=mad.end();++it) {
      score+=*it;
      ++n;
    }
  }
  return score/n;
}

double maximize_radial_symmetry(
  const size_t nfs,
  const size_t nss,
  const uint8_t* mask,
  const int32_t* data,
  const double width,
  size_t& mx,
  size_t& my
    )
{
  constexpr size_t r = 2; // search radius
  patchmap<tuple<size_t,size_t>,double> memory;
  double best_score = numeric_limits<double>::infinity();
  for (size_t i=0;i!=256;++i) {
    bool stop = true;
    const size_t _mx = mx;
    const size_t _my = my;
    for (size_t tmy = _my-r;tmy<=_my+r;++tmy) { // lazily search a squre
      for (size_t tmx = _mx-r;tmx<=_mx+r;++tmx) { // to find a circle within
        // yes I know there is a better way, but it's not worth the arithmetic
        const size_t r2 = pow(tmy>_my?tmy-_my:_my-tmy,2)
                         +pow(tmx>_mx?tmx-_mx:_mx-tmx,2);
        if (r2>r*r) continue; // now we are within the radius we want to search
        // employ mnemonic to reduce computations
        if (memory.count({tmx,tmy})==0) {
          memory[{tmx,tmy}]=radial_symmetry(nfs,nss,mask,data,width,tmx,tmy);
          //cout << tmx << " " << tmy << " " << memory[{tmx,tmy}] << endl;
        }
        // keep track of the best position and score so far
        if (memory[{tmx,tmy}]<best_score) {
          mx=tmx;
          my=tmy;
          best_score = memory[{tmx,tmy}];
          stop = false;
          //cerr << "new best score is: " << best_score << endl;
        }
      }
    }
    // if we did not find a more radially symmetric position within the search
    // radius than the one we already had, we stop
    // This search procedure is most similar to numerical gradient descent with
    // a step size of one pixel
    if (stop) break;
  }
  return best_score;
}

int main(int argc,char **argv) {
  const size_t nfs = 1032;
  const size_t nss = 1032;
  double width = 4.0;
  ifstream input_data(argv[1]);
  int32_t* data = new int32_t[nfs*nss];
  input_data.read(reinterpret_cast<char*>(data),nfs*nss*sizeof(int32_t));
  uint8_t* mask = new uint8_t[nfs*nss];
  for (size_t i=0;i!=nfs*nss;++i) {
    if (data[i]<0) { // we are using negative values to encode the pixel mask
      data[i] = 0;
      mask[i]|= 1;
    } else {
      // this makes the noise in the image more uniform, but it is not necessary
      data[i] = sqrt(data[i]);
    }
  }
  size_t mx,my;
  cin >> mx >> my;
  double score = maximize_radial_symmetry(nfs,nss,mask,data,width,mx,my);
  cout << mx << " "
       << my << " "
       << score << endl;
}
