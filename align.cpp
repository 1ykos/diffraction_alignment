#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

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
using std::numeric_limits;
using std::ofstream;
using std::setw;
using std::stoull;
using std::streamsize;
using std::string;
using std::to_string;
using std::tuple;
using std::vector;
using std::swap;
using std::abs;

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

template<class F>
double correlate(
  const size_t nfs,
  const size_t nss,
  const double* a,
  const double* b,
  const uint8_t* mask,
  const int64_t shift_fs,
  const int64_t shift_ss,
  F f
) {
  double s_xy = 0;
  size_t n = 0;
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      if (is_valid(fs)&&is_valid(ss)&&is_valid(fs+shift_fs)&&is_valid(ss+shift_ss)) {
        if ((mask[ss*nfs+fs]==0)&&(mask[(ss+shift_ss)*nfs+(fs+shift_fs)]==0)) {
          s_xy+=f(a[ss*nfs+fs])*f(b[(ss+shift_ss)*nfs+(fs+shift_fs)]);
          ++n;
        }
      }
    }
  }
  return s_xy/n;
}

void smooth_image(
  const size_t nfs,
  const size_t nss,
//  const uint8_t* mask,
  const double* a,
        double* b
)
{
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      if (is_valid(fs)&&is_valid(ss)) {
        double sum = 4*a[ss*nfs+fs];
        double n = 4;
        if (is_valid(ss,-1)) {
          n+=2;
          sum+=2*a[(ss-1)*nfs+fs];
          if (is_valid(fs,-1)) {
            n+=1;
            sum+=a[(ss-1)*nfs+(fs-1)];
          }
          if (is_valid(fs, 1)) {
            n+=1;
            sum+=a[(ss-1)*nfs+(fs+1)];
          }
        }
        if (is_valid(ss, 1)) {
          n+=2;
          sum+=2*a[(ss+1)*nfs+fs];
          if (is_valid(fs,-1)) {
            n+=1;
            sum+=a[(ss+1)*nfs+(fs-1)];
          }
          if (is_valid(fs, 1)) {
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

int main(int argc, char**argv) {
  std::ios::sync_with_stdio(false);
  const size_t nfs=1024,nss=1024;
  double* sum = new double[nfs*nss];
  double* norm = new double[nfs*nss];
  double* smooth = new double[nfs*nss];
  uint8_t* mask = new uint8_t[nfs*nss]();
  vector<double*> data;
  while (cin.read(reinterpret_cast<char*>(sum),nfs*nss*sizeof(double))) {
    for (size_t ss=0;ss!=nss;++ss) {
      for (size_t fs=0;fs!=nfs;++fs) {
        mask[ss*nfs+fs]|=sum[ss*nfs+fs]<0;
        mask[ss*nfs+fs]|=sum[ss*nfs+fs]>1e6;
        mask[ss*nfs+fs]|=(fs-512)*(fs-512)+(ss-512)*(ss-512)<64*64;
        //sum[ss*nfs+fs]=sum[ss*nfs+fs]>1e5?1e5:sum[ss*nfs+fs];
        //sum[ss*nfs+fs]=sum[ss*nfs+fs]<0?-1:sum[ss*nfs+fs];
      }
    }
    data.push_back(sum);
    sum = new double[nfs*nss];
  }
  vector<tuple<int64_t,int64_t>> shifts(data.size());
  for (size_t j=0;j!=2;++j) {
    /*for (size_t j=0;j!=data.size();++j) {
      double s = 0;
      for (size_t i=0;i!=nfs*nss;++i) s+=data[j][i];
      cerr << s << " ";
    }
    cerr << endl;*/
    for (size_t i=0;i!=nfs*nss;++i) {
      sum[i]=0;
      norm[i]=0;
    }
    for (size_t i=0;i!=data.size();++i) {
      for (size_t ss=0;ss!=nss;++ss) {
        for (size_t fs=0;fs!=nfs;++fs) {
          int64_t shift_fs = get<0>(shifts[i]);
          int64_t shift_ss = get<1>(shifts[i]);
          if (is_valid(ss,shift_ss)&&is_valid(fs,shift_fs)) {
            //if (mask[(ss+shift_ss)*nfs+(fs+shift_fs)]==0) {
              sum[ss*nfs+fs]+=data[i][(ss+shift_ss)*nfs+(fs+shift_fs)];
              norm[ss*nfs+fs]+=1.0;
            //}
          }
        }
      }
    }
    for (size_t i=0;i!=nfs*nss;++i) if (norm[i]>0.5) sum[i]/=norm[i];
    //smooth_image(nfs,nss,sum,smooth);
    //smooth_image(nfs,nss,smooth,sum);
    //smooth_image(nfs,nss,sum,smooth);
    for (size_t i=0;i!=data.size();++i) {
      int64_t shift_fs = get<0>(shifts[i]);
      int64_t shift_ss = get<1>(shifts[i]);
      double corr=0,best_corr=0;
      //cerr << best_corr << endl;
      int64_t best_shift_fs=shift_fs,best_shift_ss=shift_ss;
      for (int64_t ss = -2; ss<= 2; ++ss) {
        for (int64_t fs = -2; fs<= 2; ++fs) {
          shift_fs = get<0>(shifts[i])+fs;
          shift_ss = get<1>(shifts[i])+ss;
          corr = correlate(nfs,nss,sum,data[i],mask,shift_fs,shift_ss,map); 
          //cerr << corr << " ";
          if (corr>best_corr) {
            best_corr = corr;
            best_shift_fs=shift_fs;
            best_shift_ss=shift_ss;
          }
        }
        //cerr << endl;
      }
      //cerr << best_shift_fs << " " << best_shift_ss << " " << best_corr << endl;
      get<0>(shifts[i])=best_shift_fs;
      get<1>(shifts[i])=best_shift_ss;
    }
    //cerr << endl;
  }
  cout.write(reinterpret_cast<const char*>(sum),nfs*nss*sizeof(double));
}
