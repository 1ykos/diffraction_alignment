#include "wmath.hpp"
#include <bitset>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <string>
#include <vector>

using std::bitset;
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
using std::streamsize;
using std::string;
using std::to_string;
using std::vector;

using wmath::bswap;

int main(int argc, char**argv){
  vector<int32_t> data;
  for (uint32_t x;
       cin.read(reinterpret_cast<char*>(&x),sizeof(int32_t));
       data.push_back(x));
  vector<double> doubles;
  doubles.reserve(data.size());
  for (auto v : data) doubles.push_back(v);
  cout.write(
      reinterpret_cast<const char*>(doubles.data()),
      doubles.size()*sizeof(doubles)
      );
}
