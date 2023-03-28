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
using std::floor;
using std::fstream;
using std::ifstream;
using std::ios;
using std::numeric_limits;
using std::ofstream;
using std::round;
using std::setw;
using std::stoull;
using std::streamsize;
using std::string;
using std::to_string;

int main(int argc, char**argv){
  double* d = new double[1024];
  uint32_t* b = new uint32_t[1024];
  while (cin.read(reinterpret_cast<char*>(d),8*1024)) {
    for (size_t i=0;i!=1024;++i) {
      b[i]=d[i]<0?0:floor(d[i]);
    }
    cout.write(reinterpret_cast<const char*>(b),4*1024);
  }
}
