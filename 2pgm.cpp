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
  double* data = new double[nfs*nss];
  cin.read(reinterpret_cast<char*>(data),nfs*nss*sizeof(double));
  uint8_t* img = new uint8_t[nfs*nss];
  for (size_t i=0;i!=nfs*nss;++i) {
    data[i]=sqrt(abs(data[i]));
    data[i]=data[i]>255?255:data[i];
    data[i]=data[i]<  0?  0:data[i];
    img[i] = data[i];
  }
  cout << "P5 " << nfs << " " << nss << " 255 " << endl;
  cout.write(reinterpret_cast<char*>(img),nfs*nss);
}
