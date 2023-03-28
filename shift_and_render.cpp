#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

using std::abs;
using std::cerr;
using std::cin;
using std::cout;
using std::distance;
using std::endl;
using std::floor;
using std::fstream;
using std::get;
using std::ifstream;
using std::ios;
using std::nth_element;
using std::numeric_limits;
using std::ofstream;
using std::setw;
using std::stoull;
using std::streamsize;
using std::string;
using std::swap;
using std::to_string;
using std::tuple;
using std::vector;
using std::clamp;

const size_t inline map_down(const size_t i) {
  size_t r = (i/258)*256;
  size_t m = i%258;
  if (m==  0||m==  1) return r;
  if (m==256||m==257) return r+255;
  return r+m-1;
}

int main(int argc, char**argv) {
  std::ios::sync_with_stdio(false);
  const size_t nfs=1024,nss=1024;
  uint8_t* mask = new uint8_t[nfs*nss];
  cin.read(reinterpret_cast<char*>(mask),nfs*nss);
  double* a = new double[nfs*nss];
  int32_t* b = new int32_t[1032*1032]();
  cin.read(reinterpret_cast<char*>(a),nfs*nss*sizeof(double));
  for (size_t pss=0;pss!=1032;++pss) {
    for (size_t pfs=0;pfs!=1032;++pfs) {
      const size_t fs = map_down(pfs);
      const size_t ss = map_down(pss);
      //cerr <<  fs << " " <<  ss << " "
      //     << pfs << " " << pss << endl;
      const double c  = (((fs%256==0)||(fs%256==255))?0.5:1.0)
                       *(((ss%256==0)||(ss%256==255))?0.5:1.0);
      const size_t i = ss*nfs+fs;
      const size_t j = pss*1032+pfs;
      if (mask[i]) b[j] = -3;
      else b[j] = c*a[i];
    }
  }
  for (size_t ss=514;ss!=518;++ss) {
    for (size_t fs=0;fs!=1032;++fs) {
      b[ss*1032+fs] = -3;
    }
  }
  int shift_fs = 0, shift_ss = 0, dist = 0;
  if (argc>3) {
    shift_fs = stoull(argv[1]);
    shift_ss = stoull(argv[2]);
    dist = stoull(argv[3]);
  }
  cerr << shift_fs << " " << shift_ss << " " << dist << endl;
  int32_t* c = new int32_t[1032*1032]();
  for (size_t ss=0;ss!=1032;++ss) {
    for (size_t fs=0;fs!=1032;++fs) {
      const int shift_ss2 = shift_ss + (ss<516?(dist/2):-(dist+1)/2);
      if (ss<516) {
        if (516-ss<shift_ss+dist/2) {
          c[ss*1032+fs] = -3;
          continue;
        }
      } else {
        if (ss-516<-shift_ss+(dist+1)/2) {
          c[ss*1032+fs] = -3;
          continue;
        }
      }
      if ((ss+shift_ss2)*nfs+(fs+shift_fs)<nfs*nss) {
        c[ss*1032+fs] = b[(ss+shift_ss2)*1032+(fs+shift_fs)];
      }
    }
  }
  if constexpr (false) {
    uint8_t* img = new uint8_t[1032*1032]();
    for (size_t i=0;i!=1032*1032;++i) img[i]=
      clamp(
          int64_t(sqrt(abs(c[i]))),
          int64_t(0),
          int64_t(255)
          );
    cout << "P5 " << 1032 << " " << 1032 << " 255 " << endl;
    cout.write(reinterpret_cast<char*>(img),1032*1032);
    return 0;
  } else {
    cout.write(reinterpret_cast<char*>(c),1032*1032*sizeof(int32_t));
  }
}
