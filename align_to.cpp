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

bool is_valid(const int64_t i) {
  if (i<   0)     return false;
  if (i>1023)     return false;
  if (i%256==  0) return false;
  if (i%256==255) return false;
  return true;
}

bool is_valid(const int64_t i, const int64_t shift) {
  if (!is_valid(i))             return false;
  if (!is_valid(i+shift))       return false;
  if ((i/256)!=((i+shift)/256)) return false;
  return true;
}

inline const double map(double x) {
  return sqrt(x>1e4?1e4:(x<0?0:x));
}

const double median(
  double* begin,
  double* end
) {
  nth_element(
      begin,
      begin+(end-begin-1)/2,
      end);
  nth_element(
      begin,
      begin+(end-begin-1)/2,
      end);
  const double c =
      0.5*(*(begin+(end-begin-1)/2))
     +0.5*(*(begin+(end-begin-0)/2));
  return c;
}

void shift(
  const size_t nfs,
  const size_t nss,
  const double* a,
        double* b,
  const uint8_t* mask,
  const int64_t shift_fs,
  const int64_t shift_ss
) {
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      const size_t ib = ss*nfs+fs;
      const size_t ia = (ss+shift_ss)*nfs+(fs+shift_fs);
      bool is_valid_a = is_valid(fs,shift_fs)&&is_valid(ss,shift_ss);
      if (is_valid_a) is_valid_a=is_valid_a&&(mask[ia]==0);
      if (is_valid_a) b[ib] = a[ia];
      else            b[ib] = a[ib];
    }
  }
}

void scale(
  const size_t nfs,
  const size_t nss,
  const double* a,
        double* b,
        double* buffer,
  const uint8_t* mask
) {
  double* d = buffer;
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      const size_t i = ss*nfs+fs;
      bool valid = is_valid(fs)&&is_valid(ss);
      if (valid) valid=valid&&(mask[i]==0);
      if (!valid) continue;
      if ((a[i]==0)&&(b[i]==0)) continue;
      (*(d++))=a[i]/b[i];
    }
  }
  const double c = median(buffer,d);
  for (size_t ss=0;ss!=nss;++ss)
    for (size_t fs=0;fs!=nfs;++fs)
      b[ss*nfs+fs]*=c;
}


double smad(
  const size_t nfs,
  const size_t nss,
  const double* a,
  const double* b,
        double* buffer0,
        double* buffer1,
  const uint8_t* mask,
  const int64_t shift_fs,
  const int64_t shift_ss
) {
  size_t min_ss = 257;
  size_t max_ss = 766;
  size_t min_fs = 257;
  size_t max_fs = 766;
  shift(nfs,nss,b,buffer0,mask,shift_fs,shift_ss);
  scale(nfs,nss,a,buffer0,buffer1,mask);
  double n = 0;
  double s = 0;
  for (size_t ss=min_ss;ss<=max_ss;++ss) {
    for (size_t fs=min_fs;fs<=max_fs;++fs) {
      const size_t i = ss*nfs+fs;
      bool valid = is_valid(fs)&&is_valid(ss);
      if (valid) valid=valid&&(mask[i]==0);
      if (!valid) continue;
      s+= abs(a[i]-buffer0[i]);
    }
  }
  return s; 
}

int main(int argc, char**argv) {
  std::ios::sync_with_stdio(false);
  const size_t nfs=1024,nss=1024;
  double* a = new double[nfs*nss];
  double* b = new double[nfs*nss];
  uint8_t* mask = new uint8_t[nfs*nss]();
  if (argc>1) {
    ifstream(argv[1]).read(reinterpret_cast<char*>(mask),nfs*nss);
  }
  double* buffer0 = new double[nfs*nss];
  double* buffer1 = new double[nfs*nss];
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      mask[ss*nfs+fs]|=(fs-540)*(fs-540)+(ss-512)*(ss-512)<50*50;
      //mask[ss*nfs+fs]|=(fs-512)*(fs-512)+(ss-512)*(ss-512)>260*260;
    }
  }
  cin.read(reinterpret_cast<char*>(a),nfs*nss*sizeof(double));
  cin.read(reinterpret_cast<char*>(b),nfs*nss*sizeof(double));
  for (size_t i=0;i!=nfs*nss;++i) {
    a[i]=sqrt(abs(a[i]));
    b[i]=sqrt(abs(b[i]));
  }
  double best=1e300;
  //cerr << best_corr << endl;
  int64_t best_shift_fs=0,best_shift_ss=0;
  int64_t r = 10;
  for (int64_t ss = -r; ss<= r; ++ss) {
    for (int64_t fs = 0; fs*fs+ss*ss<=r*r; fs=(fs<=0?-fs+1:-fs)) {
      const double score = smad(nfs,nss,a,b,buffer0,buffer1,mask,fs,ss); 
      cout << fs << " " << ss << " " << score << endl;
      if (score<best) {
        best = score;
        best_shift_fs=fs;
        best_shift_ss=ss;
      }
    }
  }
  //cout << best_shift_fs << " " << best_shift_ss << " " << best << endl;
}
