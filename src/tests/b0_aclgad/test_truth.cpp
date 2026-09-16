// SPDX-License-Identifier: LGPL-3.0-or-later
#include "B0MeasurementTruth.h"
#include <cmath>
#include <iostream>
#include <vector>
using namespace eicrecon::b0;
struct Object {
  struct Id { int collectionID, index; } id;
  Id getObjectID() const { return id; }
  Object getRawHit() const { return *this; }
};
struct Measurement {
  std::vector<Object> hits; std::vector<float> weights;
  const auto& getHits() const { return hits; }
  const auto& getWeights() const { return weights; }
  auto hits_size() const { return hits.size(); }
};
int main() {
  try {
    const TruthId a{9,1}, b{9,2};
    RawTruthVotes truth{{{1,1},{{a,.2}}},{{1,2},{{b,1}}}};
    Measurement m{{{{1,1}},{{1,2}}},{.75,.25}};
    auto result=measurementTruth(m,truth);
    if(std::abs(result[a]-.15)>1e-9 || std::abs(result[b]-.25)>1e-9) throw std::runtime_error("noise/charge weights lost");
    truth[TruthId{1,1}]={{a,2},{b,2}}; m.weights={};
    result=measurementTruth(m,truth);
    if(std::abs(result[a]-.25)>1e-9 || std::abs(result[b]-.75)>1e-9) throw std::runtime_error("legacy normalization");
    m.hits.push_back({{1,3}}); result=measurementTruth(m,truth);
    if(std::abs(result[a]+result[b]-2./3)>1e-9) throw std::runtime_error("unassociated hit disappeared");
    m.hits.push_back(m.hits[0]); bool failed=false;
    try { measurementTruth(m,truth); } catch(const std::invalid_argument&) { failed=true; }
    if(!failed) throw std::runtime_error("duplicate hit accepted");
    std::cout<<"4 truth-mixture checks passed\n";
  } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
