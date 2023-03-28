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

// mask the edges of the panels
bool is_valid(const int64_t i) {
  if (i<   0)     return false;
  if (i>1023)     return false;
  if (i%256==  0) return false;
  if (i%256==255) return false;
  return true;
}

// convolve with the following kernel
// 1  2  1
// 2  4  2
// 1  2  1
void smooth_image(
  const size_t nfs,
  const size_t nss,
  const uint8_t* mask,
  const double* a,
        double* b
)
{
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      if (is_valid(fs)&&is_valid(ss)) {
        double sum = 4*a[ss*nfs+fs];
        double n = 4;
        if (is_valid(ss-1)) {
          n+=2;
          sum+=2*a[(ss-1)*nfs+fs];
          if (is_valid(fs-1)) {
            n+=1;
            sum+=a[(ss-1)*nfs+(fs-1)];
          }
          if (is_valid(fs+1)) {
            n+=1;
            sum+=a[(ss-1)*nfs+(fs+1)];
          }
        }
        if (is_valid(ss+1)) {
          n+=2;
          sum+=2*a[(ss+1)*nfs+fs];
          if (is_valid(fs-1)) {
            n+=1;
            sum+=a[(ss+1)*nfs+(fs-1)];
          }
          if (is_valid(fs+1)) {
            n+=1;
            sum+=a[(ss+1)*nfs+(fs+1)];
          }
        }
        b[nfs*ss+fs]=sum/n;
      } else {
        b[nfs*ss+fs] = 0;
      }
    }
  }
}

// compute the discrete Laplacian of the image
// i.e. add up the discrete second derivatives
// essentially convolve the image with the following kernel:
// 0  1  0
// 1 -4  1
// 0  1  0
void extract_features(
    const size_t nfs,
    const size_t nss,
    const uint8_t* mask,
    const double* data,
          double* features
) {
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      if (!is_valid(fs)) continue;
      if (!is_valid(ss)) continue;
      const size_t i = ss*nfs+fs;
      if (mask[i])       continue;
      features[i] = 0;
      const size_t j0 = (ss-1)*nfs+fs;
      const size_t j1 = (ss+1)*nfs+fs;
      if (is_valid(ss-1) && is_valid(ss+1)) {
        if ((!mask[j0])&&(!mask[j1])) {
          features[i]+= data[j0] + data[j1] - 2*data[i];
        }
      }
      const size_t k0 = ss*nfs+fs-1;
      const size_t k1 = ss*nfs+fs+1;
      if (is_valid(fs-1) && is_valid(fs+1)) {
        if ((!mask[k0])&&(!mask[k1])) {
          features[i]+= data[k0] + data[k1] - 2*data[i];
        }
      }
      continue; // this skipped part may improve the result or may not
                // but it complicates the algorithm so it is skipped
      const size_t l0 = (ss-1)*nfs+fs-1;
      const size_t l1 = (ss+1)*nfs+fs+1;
      const size_t m0 = (ss-1)*nfs+fs+1;
      const size_t m1 = (ss+1)*nfs+fs-1;
      if (is_valid(fs-1) && is_valid(fs+1)
       && is_valid(ss-1) && is_valid(ss+1)) {
        if ((!mask[l0])&&(!mask[l1])) {
          features[i]+= 0.5*data[l0] + 0.5*data[l1] - data[i];
        }
        if ((!mask[m0])&&(!mask[m1])) {
          features[i]+= 0.5*data[m0] + 0.5*data[m1] - data[i];
        }
      }
    }
  }
}

int main(int argc,char **argv) {
  const size_t nfs = 1024;
  const size_t nss = 1024;
  uint8_t* mask = new uint8_t[nfs*nss];
  cin.read(reinterpret_cast<char*>(mask),nfs*nss);
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      if (is_valid(fs)&&is_valid(ss)) continue;
      mask[ss*nfs+fs]|=1u;
    }
  }
  double * data = new double[nfs*nss];
  cin.read(reinterpret_cast<char*>(data),nfs*nss*sizeof(double));
  for (size_t i=0;i!=nfs*nss;++i) data[i] = data[i]>0?sqrt(data[i]):0;
  double * tmp = new double[nfs*nss]();
  //smooth_image(nfs,nss,mask,data,tmp);
  double * features = new double[nfs*nss];
  extract_features(nfs,nss,mask,data,features);
  //copy(features,features+nfs*nss,tmp);
  size_t j=0;
  for (size_t i=0;i!=nfs*nss;++i) if (!mask[i]) {
    tmp[j++]=features[i];
    //cerr << tmp[j-1] << endl;
  }
  // calculate median of curvature (expected 0 )
  nth_element(
      tmp,
      tmp+j/2,
      tmp+j);
  double median = *(tmp+j/2); 
  // robust estimate of background standard deviation at 83% percentile
  nth_element(
      tmp,
      tmp+j*41/256,
      tmp+j);
  double mad = *(tmp+j*41/256)-median;
  cerr << median << " "  << mad << endl;
  patchmap<tuple<uint32_t,uint32_t>,uint64_t>         pixel_color;
  patchmap<uint64_t,vector<tuple<uint32_t,uint32_t>>> color_pixel;
  size_t next_color = 0;
  for (uint32_t ss=0;ss!=nss;++ss) {
    for (uint32_t fs=0;fs!=nfs;++fs) {
      const size_t i = ss*nfs+fs;
      if (2.5*mad-(features[i]-median)<=0) continue;
      uint64_t color = next_color++;
      pixel_color[{fs,ss}] = color;
      color_pixel[color].push_back({fs,ss});
      const vector<tuple<uint32_t,uint32_t>> search({
          {fs-1,ss-1},
          {fs+0,ss-1},
          {fs-1,ss+0}
          });
      for (auto s : search) {
        if (pixel_color.count(s)==0) continue;
        const size_t old_color = pixel_color.at({fs,ss});
        const size_t new_color = pixel_color.at(s);
        if (old_color==new_color) continue;
        for (auto p : color_pixel[old_color]) {
          pixel_color[p]=new_color;
          color_pixel[new_color].push_back(p);
        }
        color_pixel.erase(old_color);
      }
    }
  }
  patchmap<uint64_t,tuple<double,double,double>> color_coordinates;
  for (auto entry : color_pixel) {
    const auto color = entry.first;
    if (color_pixel[color].size()<2) continue;
    for (auto fsss : color_pixel[color]) {
      const double dfs = get<0>(fsss) + 0.5;
      const double dss = get<1>(fsss) + 0.5;
      const double w   = features[get<1>(fsss)*nfs+get<0>(fsss)];
      get<0>(color_coordinates[color])+=w*dfs;
      get<1>(color_coordinates[color])+=w*dss;
      get<2>(color_coordinates[color])+=w;
    }
  }
  for (auto entry : color_coordinates) {
    const double i = 1.0/get<2>(entry.second);
    get<0>(entry.second)*=i;
    get<1>(entry.second)*=i;
    cout << get<0>(entry.second) << " "
         << get<1>(entry.second) << " "
         << get<2>(entry.second) << endl;  
  }
  return 0;
  // this is how I would export images if I want to view them quickly
  uint8_t* img = new uint8_t[3*nfs*nss]();
  for (uint32_t ss=0;ss!=nss;++ss) {
    for (uint32_t fs=0;fs!=nfs;++fs) {
      if (pixel_color.count({fs,ss})==0) continue;
      const uint64_t color = pixel_color.at({fs,ss});
      if (color_pixel[color].size()<2) continue;
      uint64_t h = whash::hash<uint64_t>{}(color);
      uint64_t s = ~uint64_t(0);
      uint64_t v = ~uint64_t(0);
      wmath::hsv_rgb(h,s,v);
      img[3*(ss*nfs+fs)+0] = wmath::apply_srgb_gamma(h);
      img[3*(ss*nfs+fs)+1] = wmath::apply_srgb_gamma(s);
      img[3*(ss*nfs+fs)+2] = wmath::apply_srgb_gamma(v);
    }
  }
  cout << "P6 " << nfs << " " << nss << " 255" << endl;
  cout.write(reinterpret_cast<char*>(img),3*nfs*nss);
  return 0;
}
