#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <iomanip>
#include <iterator>
#include <cmath>

using std::cerr;
using std::cin;
using std::cout;
using std::distance;
using std::endl;
using std::fstream;
using std::ifstream;
using std::ios;
using std::numeric_limits;
using std::ofstream;
using std::setw;
using std::stoull;
using std::streamsize;
using std::string;
using std::to_string;
using std::floor;

int main(int argc, char**argv){
  std::ios::sync_with_stdio(false);
  size_t nfs=1024,nss=1024;
  double* data = new double[nfs*nss*sizeof(double)];
  double* plus = new double[nfs*nss*sizeof(double)];
  double n = 0;
  if (cin.read(reinterpret_cast<char*>(data),nfs*nss*sizeof(double))) {
    ++n;
    while (cin.read(reinterpret_cast<char*>(plus),nfs*nss*sizeof(double))) {
      ++n;
      for (size_t ss=0;ss!=nss;++ss) {
        for (size_t fs=0;fs!=nfs;++fs) {
          data[nfs*ss+fs]+=(plus[nfs*ss+fs]-data[nfs*ss+fs])/n;
        }
      }
    }
    cout.write(reinterpret_cast<const char*>(data),nfs*nss*sizeof(double));
  }
}
