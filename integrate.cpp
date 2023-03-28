#include "geometry.hpp"
#include "wmath.hpp"
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
using std::stod;
using std::streamsize;
using std::string;
using std::stringstream;
using std::swap;
using std::to_string;
using std::tuple;
using std::unordered_map;
using std::vector;
using wmath::long_mul;
using wmath::clz;
using wmath::digits;
using wmath::clip;
using wmath::destructive_median;

using whash::patchmap;

typedef tuple<int32_t,int32_t,int32_t> IDX;

const double pi = 3.14159265;

const double gauss_cdf(
    const double x,
    const double m,
    const double v
    ) {
  return 0.5*(1+erf((x-m)/sqrt(2*v)));
}

const double gauss_integral(
    const double m,
    const double v,
    const double a,
    const double b
    ) {
  return gauss_cdf(b,m,v)-gauss_cdf(a,m,v);
}

inline matrix<double,3,3> const rotation_matrix(
    const double a,
    const matrix<double,3,1>& u
    ) {
  const double q0 = u(0)*sin(a/2);
  const double q1 = u(1)*sin(a/2);
  const double q2 = u(2)*sin(a/2);
  const double q3 = cos(a/2);
  return matrix<double,3,3>
  {
    q0*q0+q1*q1-q2*q2-q3*q3,         2*(q1*q2-q0*q3),         2*(q1*q3+q0*q2),
            2*(q1*q2+q0*q3), q0*q0-q1*q1+q2*q2-q3*q3,         2*(q2*q3-q0*q1),
            2*(q1*q3-q0*q2),         2*(q2*q3+q0*q1), q0*q0-q1*q1-q2*q2+q3*q3
  };
}

struct integration_entry{
  double foreground_sum = 0;
  double foreground_sumsquared = 0;
  size_t foreground_num = 0;
  double background_sum = 0;
  double background_sumsquared = 0;
  size_t background_num = 0;
};

int main(int argc,char **argv) {
  const size_t nfs=1024,nss=1024;
  const double angle_per_image = 0.1*pi/180;
  double angle = -32.0*pi/180; // start at 32 deg
  matrix<double,3,1> axis{0.0,1.0,0.0};
  ifstream axis_file("axis");
  if (axis_file) {
    axis_file >> axis(0) >> axis(1) >> axis(2);
  }
  matrix<double,3,1> direct_beam{0.0,0.0,284.5};
  ifstream beam_file("beam");
  if (beam_file) {
    beam_file >> direct_beam(0) >> direct_beam(1) >> direct_beam(2);
  }
  const double max_trusted_value = 1e5;
  geometry::geometry geom;
  {
    ifstream file("geom");
    read_geometry(file,geom);
  }
  uint32_t* pixel_to_panel = new uint32_t[nfs*nss];
  for (size_t i=0;i!=nfs*nss;++i) pixel_to_panel[i] = ~uint32_t(0);
  vector<tuple<size_t,size_t,size_t,size_t>> reorder_panels{
    { 256, 512, 256, 512},
    { 512, 768, 256, 512},
    { 768,1024, 256, 512},
    { 256, 512, 512, 768},
    { 512, 768, 512, 768},
    { 768,1024, 512, 768}
  };
  for (size_t p=0;p!=reorder_panels.size();++p) {
    const size_t min_fs = get<0>(reorder_panels[p]);
    const size_t max_fs = get<1>(reorder_panels[p]);
    const size_t min_ss = get<2>(reorder_panels[p]);
    const size_t max_ss = get<3>(reorder_panels[p]);
    for (size_t ss=min_ss+1;ss!=max_ss-1;++ss) {
      for (size_t fs=min_fs+1;fs!=max_fs-1;++fs) {
        pixel_to_panel[ss*nfs+fs]=p;
      }
    }
  }
  size_t nx=256,ny=256,nz=256;
  /*matrix<double,3,3> U{
     5.13181, 1.22186,-0.07982,
    -2.15376, 9.00584, 0.07675,
    -1.84814,-0.14887,20.73700
  };*/
  /*
  matrix<double,3,3> U{
       1.514237, 5.079067,-0.012401,
      -8.910385, 2.655986,-0.201858,
      -1.022791,-1.843335,20.893917
  };*/
  matrix<double,3,3> U{
       1.532545, 5.115983, 0.020794,
      -8.895800, 2.665134,-0.075360,
      -0.808783,-2.087564,21.644617
  };
  ifstream cell_file("cell");
  if (cell_file) {
    cell_file >> U(0,0) >> U(0,1) >> U(0,2);
    cell_file >> U(1,0) >> U(1,1) >> U(1,2);
    cell_file >> U(2,0) >> U(2,1) >> U(2,2);
  }
  matrix<double,3,3> R = inv(U);
  //cerr << U << endl;
  //cerr << R << endl;
  uint8_t* img  = new uint8_t[3*nfs*nss];
  uint32_t* buffer = new uint32_t[nfs*nss];
  uint8_t* mask = new uint8_t[nfs*nss]; // 0 for masked
  cin.read(reinterpret_cast<char*>(mask),nfs*nss);
  double total_integrated_intensity = 0;
  double total_background_intensity = 0;
  size_t num_integrated_pixels = 0;
  size_t num_background_pixels = 0;
  double score = 0;
  double norm = 0;
  constexpr bool write_images = true;
  constexpr bool integrate = true;
  for (size_t i=0;
       cin.read(reinterpret_cast<char*>(buffer),nfs*nss*sizeof(uint32_t));
       ++i
       ) {
    const matrix<double,3,3> rot   = rotation_matrix(-angle,axis);
    const matrix<double,3,3> rot_U = rot*U;
    const matrix<double,3,3> R_rot = inv(rot_U);
    //cerr << U_rot << endl;
    patchmap<IDX,integration_entry> hkl_intensities;
    for (size_t ss=256;ss!=768;++ss) {
      for (size_t fs=0;fs!=nfs;++fs) {
        const size_t j = ss*nfs+fs;
        if constexpr (write_images) img[3*j+0] = img[3*j+1] = img[3*j+1] = 0;
        if (mask[j]) continue;
        const size_t panel_number = pixel_to_panel[j];
        if (panel_number==~uint32_t(0)) continue;
        if (buffer[j]>max_trusted_value) continue;
        // gray for everything else
        if constexpr (write_images) {
          img[3*j+0] = clip(sqrt(16*buffer[j]),0,255);
          img[3*j+1] = clip(sqrt(16*buffer[j]),0,255);
          img[3*j+2] = clip(sqrt(16*buffer[j]),0,255);
        }
        //cerr << panel_number << endl;
        const size_t min_fs = get<0>(reorder_panels[panel_number]);
        const size_t max_fs = get<1>(reorder_panels[panel_number]);
        const size_t min_ss = get<2>(reorder_panels[panel_number]);
        const size_t max_ss = get<3>(reorder_panels[panel_number]);
        const size_t pfs    = fs-min_fs;
        const size_t pss    = ss-min_ss;
        //cout << pfs << " "  << pss << endl;
        if (!geom.panels[panel_number].isvalid(pfs,pss)) continue;
        const matrix<double,2,1> fsss{double(pfs+0.5),double(pss+0.5)};
        matrix<double,3,1> x = geom.panels[panel_number](fsss);
        const matrix<double,3,1> m = normalize(x)*length(direct_beam)-direct_beam;
        const matrix<double,3,1> rot_U_m = rot_U*m;
        matrix<double,3,1> dhkl=round(rot_U_m);
        const double w = exp(-8*length(dhkl-round(dhkl))); 
        norm += w;
        score += w*buffer[j];
        IDX hkl;
        get<0>(hkl) = dhkl(0);
        get<1>(hkl) = dhkl(1);
        get<2>(hkl) = dhkl(2);
        //cerr << get<0>(hkl) << " " << get<1>(hkl) << " " << get<2>(hkl) << endl;
        if ((get<0>(hkl)==0)&&(get<1>(hkl)==0)&&(get<2>(hkl)==0)) continue;
        if (((get<0>(hkl)+get<1>(hkl))%2)!=0) continue;
        const matrix<double,3,1> y = R_rot*dhkl;
        const matrix<double,3,3> S {
          0.001 ,0.0    ,0.0    ,
          0.0   ,0.001  ,0.0    ,
          0.0   ,0.0    ,0.00003
        };
        //cerr << trans(m) << trans(y);
        if (trans(m-y)*inv(S)*(m-y)<2) {
          auto & e = hkl_intensities[hkl];
          /*
          cerr << get<0>(hkl) << " "
               << get<1>(hkl) << " "
               << get<2>(hkl) << " "
               << buffer[j] << endl;*/
          //cerr << acos(trans(x)*y) << endl;
          if (trans(m-y)*inv(S)*(m-y)<1) {
            if constexpr (write_images) {
              img[3*j+0] = 0;
              img[3*j+1] = clip(sqrt(16*buffer[j]),0,255);
              img[3*j+2] = clip(sqrt(16*buffer[j]),0,255);
            }
            if constexpr (integrate) {
              e.foreground_sum+=buffer[j];
              e.foreground_sumsquared+=pow(buffer[j],2);
              e.foreground_num+=1;
            }
            num_integrated_pixels += 1;
            total_integrated_intensity += buffer[j];
            //cerr << i << " "
            //     << _x << " "
            //     << _y << " "
            //     << buffer[j] << endl;
          } else {
            if constexpr (write_images) {
              img[3*j+0] = clip(sqrt(16*buffer[j]),0,255);
              img[3*j+1] = 0;
              img[3*j+2] = 0;
            }
            if constexpr (integrate) {
              e.background_sum+=buffer[j];
              e.background_sumsquared+=pow(buffer[j],2);
              e.background_num+=1;
            }
            num_background_pixels += 1;
            total_background_intensity += buffer[j];
          }
        }
      }
    }
    if constexpr (write_images) {
      string name = to_string(i);
      name = string(3-name.length(),'0')+name; // prefix 0
      name = "tmp"+name+".pgm";
      ofstream img_file(name);
      img_file << "P6 " << nfs << " " << nss << " 255" << endl;
      img_file.write(reinterpret_cast<const char*>(img),3*nfs*nss);
    }
    /*
    vector<double> variance_ratios;
    for (auto it=hkl_intensities.begin();it!=hkl_intensities.end();++it) {
      const double variance =
          (it->second).background_sumsquared
          -pow((it->second).background_sum,2)/(it->second).background_num;
      //cerr << (it->second).background_sum << " " << variance << endl;
      variance_ratios.push_back(variance/abs((it->second).background_sum));
    }
    */
    const double r = 2.8;
    /*
    const double r = destructive_median(
          variance_ratios.begin(),
          variance_ratios.end()
        );
        */
    //cerr << "estimated ADU per quantum:" << sqrt(r) << endl;
    if constexpr (integrate) {
      for (auto it=hkl_intensities.begin();it!=hkl_intensities.end();++it) {
        if ((it->second).background_num<9) continue;
        if ((it->second).background_num<1) continue;
        cout << get<0>(it->first) << " ";
        cout << get<1>(it->first) << " ";
        cout << get<2>(it->first) << " ";
        const double c = (1.0*(it->second).foreground_num)
                        /(it->second).background_num;
        const double background_estimate = (it->second).background_sum*c;
        const double intensity = (it->second).foreground_sum-background_estimate;
        const double variance  = r*(it->second).foreground_sum;
                                +c*((it->second).background_sumsquared
                                   /(it->second).background_num
                                   -pow((it->second).background_sum
                                       /(it->second).background_num,2));
        cout << i << " ";
        cout << intensity << " ";
        cout << variance << endl;
      }
    }
    angle+=angle_per_image;
    //cerr << angle/pi*180 << endl;
  }
  score/=norm;
  //double score = total_integrated_intensity/num_integrated_pixels;
  //score -= total_background_intensity*num_integrated_pixels
  //  /num_background_pixels;
  cerr << score << endl;
}
