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

const void byte_offset_encode(const int32_t& v, char*& p){
  if ((v<127)&&(v>-127)){
    *p=int8_t(v);
    ++p;
    return;
  }
  if ((v<32767)&&(v>-32767)){
    (*p)=char(0x80);
    ++p;
    *p=0xFF&(int16_t(v));
    ++p;
    *p=0xFF&(int16_t(v)>>8);
    ++p;
    return;
  }
  *p=char(0x80);
  ++p;
  *p=char(0x00);
  ++p;
  *p=char(0x80);
  ++p;
  *p=0xFF&(int32_t(v));
  ++p;
  *p=0xFF&(int32_t(v)>>8);
  ++p;
  *p=0xFF&(int32_t(v)>>16);
  ++p;
  *p=0xFF&(int32_t(v)>>24);
  ++p;
  return;
}

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
      v<<=8;
      v^=b;
      if(v&0x8000){
        s=3;
        return;
      }
      v=int16_t((int16_t(0ll)^v));
      return;
    case 3:
      v=0;
      v^=b;
      s=4;
      return;
    case 4:
      v<<=8;
      v^=b;
      s=5;
      return;
    case 5:
      v<<=8;
      v^=b;
      s=6;
      return;
    case 6:
      v<<=8;
      v^=b;
      s=0;
      return;
    default:
      return;
  }
}

const string crlf = "\r\n";

string const cbf_header(
    const size_t& fs,
    const size_t& ss,
    const size_t& bs
    ){
  return  "###CBF: VERSION 1.5, CBFlib v0.7.8 - SLS/DECTRIS PILATUS detectors"+crlf+crlf
         +"_array_data.array_id image_1"+crlf
         +"_array_data.binary_id 1"+crlf
         +"_array_data.data"+crlf+";"+crlf
         +"--CIF-BINARY-FORMAT-SECTION--"+crlf
         +"Content-Type: application/octet-stream"+crlf
         +"     conversions=\"x-CBF_BYTE_OFFSET\""+crlf
         +"Content-Transfer-Encoding: BINARY"+crlf
         +"X-Binary-Size: "+to_string(bs)+crlf
         +"X-Binary-ID: 1"+crlf
         +"X-Binary-Element-Type: \"signed 32-bit integer\""+crlf
         +"X-Binary-Element-Byte-Order: LITTLE_ENDIAN"+crlf
         +"X-Binary-Number-of-Elements: "+to_string(fs*ss)+crlf
         +"X-Binary-Size-Fastest-Dimension: "+to_string(fs)+crlf
         +"X-Binary-Size-Second-Dimension: "+to_string(ss)+crlf
         +crlf;
}

const string cbf_footer = "--CIF-BINARY-FORMAT-SECTION----"+crlf
                          +";"+crlf+crlf;


int main(int argc, char**argv){
  std::ios::sync_with_stdio(false);
  size_t fs=0,ss=0,nc;
  if(argc<2) {
    cerr << "please specify the dimensions "
         << "I sadly cannot know from a binary inupt" << endl;
  } else {
    fs = stoull(argv[1]);
    ss = stoull(argv[2]);
  }
  char* p0 = new char[fs*ss*16];
  char* p1 = p0;
  double d = 0;
  int32_t v = 0;
  int32_t t;
  cin.read(reinterpret_cast<char*>(&d),8);
  v = floor(d+0.5);
  if (v<0) v=0; // CAREFULL
  *p1=0xFF&(int32_t(v));
  ++p1;
  *p1=0xFF&(int32_t(v)>>8);
  ++p1;
  *p1=0xFF&(int32_t(v)>>16);
  ++p1;
  *p1=0xFF&(int32_t(v)>>24);
  ++p1;
  t=v;
  while (cin){
    d=0;
    cin.read(reinterpret_cast<char*>(&d),8);
    v = floor(d+0.5);  
    byte_offset_encode(int32_t(v)-t,p1);
    t=v;
    if (p1+5-p0>fs*ss*16) break;
  }
  cout << cbf_header(fs,ss,distance(p0,p1)); 
  cout.put(0x0C);
  cout.put(0x1A);
  cout.put(0x04);
  cout.put(0xD5);
  cout.write(p0,distance(p0,p1));
  cout << cbf_footer;
}
