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
  cin.read(reinterpret_cast<char*>(data),nfs*nss*sizeof(double));
  for (size_t ss=0;ss!=nss;++ss) {
    for (size_t fs=0;fs!=nfs;++fs) {
      uint32_t mask = 1;
      mask&=((fs%256)!=  0);
      mask&=((fs%256)!=255);
      mask&=((ss%256)!=  0);
      mask&=((ss%256)!=255);
      data[ss*nfs+fs]*=mask;
      data[ss*nfs+fs]-=(mask==0);
    }
  }
  cout.write(reinterpret_cast<const char*>(data),nfs*nss*sizeof(double));
}
