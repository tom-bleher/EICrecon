// SPDX-License-Identifier: LGPL-3.0-or-later
#include "B0ACLGADCalibration.h"
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
using namespace eicrecon::b0::aclgad;
int main(int argc, char** argv) {
  try {
    if(argc!=2) throw std::runtime_error("expected example path");
    std::ifstream stream(argv[1]);
    if(!stream) throw std::runtime_error("missing example");
    const std::string example((std::istreambuf_iterator<char>(stream)),{});
    auto parse=[](const std::string& text, bool allow) { std::istringstream in(text); return readCalibration(in,allow); };
    if(parse(example,true).config.pitchX!=.5) throw std::runtime_error("wrong calibration");
    auto reject=[&](const std::string& text,bool allow) {
      bool failed=false; try { parse(text,allow); } catch(const std::exception&) { failed=true; }
      if(!failed) throw std::runtime_error("invalid calibration accepted");
    };
    reject(example,false); reject(example+"\npitchX=.5",true);
    reject(example+"\nunknown=42",true); reject("schemaVersion=1",true);
    for(const auto& replacement: {std::string("-1"), std::string("1.2"), std::string("18446744073709551616")}) {
      auto text=example; const auto p=text.find("maxWork=20000000");
      text.replace(p,16,"maxWork="+replacement); reject(text,true);
    }
    auto text=example; text.replace(text.find("noiseSigma=500"),14,"noiseSigma=nan"); reject(text,true);
    std::cout<<"9 calibration checks passed\n";
  } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
