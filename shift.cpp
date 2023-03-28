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

int main(int argc, char**argv) {
  std::ios::sync_with_stdio(false);
  const size_t nfs=1024,nss=1024;
  double* a = new double[nfs*nss];
  double* b = new double[nfs*nss]();
  size_t shift_fs = 0, shift_ss = 0;
  if (argc>2) {
    shift_fs = stoull(argv[1]);
    shift_ss = stoull(argv[2]);
  }
  cin.read(reinterpret_cast<char*>(a),nfs*nss*sizeof(double));
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      if ((ss+shift_ss)*nfs+(fs+shift_fs)<nfs*nss) {
        b[ss*nfs+fs] = a[(ss+shift_ss)*nfs+(fs+shift_fs)];
      }
    }
  }
  cout.write(reinterpret_cast<char*>(b),nfs*nss*sizeof(double));
}
