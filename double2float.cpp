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
  double* d = new double[1024];
  float* f = new float[1024];
  while (cin.read(reinterpret_cast<char*>(d),8*1024)) {
    for (size_t i=0;i!=1024;++i) f[i]=d[i];
    cout.write(reinterpret_cast<const char*>(f),4*1024);
  }
}
