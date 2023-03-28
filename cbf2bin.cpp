#include "wmath.hpp"
#include <iostream>
#include <fstream>
#include <bitset>
#include <string>
#include <limits>
#include <iomanip>
#include <iterator>

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

using wmath::bswap;

void inline cbf_data_decode(
          int32_t& v,
    const uint8_t& b,
          size_t & s
    ){
  switch (s){
    case 0:
      if (b==0x80) s=1;
      else v = static_cast<int8_t>(b);
      return;
    case 1:
      v=0;
      v^=b;
      s=2;
      return;
    case 2:
      v^=(int32_t(b)<<8);
      if(v==0x8000){
        s=3;
        return;
      }
      v=*reinterpret_cast<int16_t*>(&v);
      s=0;
      return;
    case 3:
      v=0;
      v^=b;
      s=4;
      return;
    case 4:
      v^=(int32_t(b)<<8);
      s=5;
      return;
    case 5:
      v^=(int32_t(b)<<16);
      s=6;
      return;
    case 6:
      v^=(int32_t(b)<<24);
      s=0;
      return;
    default:
      return;
  }
}

int main(int argc, char**argv){
  std::ios::sync_with_stdio(false);
  int32_t v=0,t=0;
  uint8_t b;
  size_t s=0;
  while (cin) {
    cin.read (reinterpret_cast<char*>(&b),1);
    switch (s){
      case 0:
        if (b==0x0C) s=1;
        break;
      case 1:
        if (b==0x1A) s=2;
        break;
      case 2:
        if (b==0x04) s=3;
        break;
      case 3:
        if (b==0xD5) s=4;
        break;
    }
    if (s==4) break;
  }
  v = 0;
  t = 0;
  s=0;
  while (cin){
    cin.read (reinterpret_cast<char*>(&b),1);
    cbf_data_decode(v,b,s);
    if (s==0){
      t+=v;
      cout.write(reinterpret_cast<char*>(&t),4);
      v=0;
    }
  }
}
