// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration
#include "B0ACLGADResponse.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace eicrecon::b0::aclgad;
namespace {
int checks = 0;
void check(bool value, const char* message) {
  ++checks; if (!value) throw std::runtime_error(message);
}
void near(double a, double b, double eps = 1e-9) { check(std::abs(a-b) < eps, "numeric mismatch"); }
template <typename F> void throws(F f) {
  bool caught = false; try { f(); } catch (const std::exception&) { caught = true; }
  check(caught, "invalid input did not throw");
}
Config clean() {
  Config c; c.noiseSigma=0; c.channelGainSigma=0; c.sensorTimeSigma=0;
  c.channelTimeSigma=0; c.clockTimeSigma=0; c.threshold=1; c.adcStep=1;
  c.sigmaSharingX=0.001; c.sigmaSharingY=0.001; c.gain=1; c.pairEnergy=1;
  c.noiseOnlySensors=false; return c;
}
bool equal(const Result& a, const Result& b) {
  if (a.channels.size()!=b.channels.size() || a.clusters.size()!=b.clusters.size()) return false;
  for (std::size_t i=0;i<a.channels.size();++i) {
    const auto& x=a.channels[i]; const auto& y=b.channels[i];
    if (std::tie(x.sensor,x.ix,x.iy,x.adc,x.time,x.truth)!=std::tie(y.sensor,y.ix,y.iy,y.adc,y.time,y.truth)) return false;
  }
  for (std::size_t i=0;i<a.clusters.size();++i) {
    const auto& x=a.clusters[i]; const auto& y=b.clusters[i];
    if (std::tie(x.sensor,x.channels,x.x,x.y,x.xx,x.yy,x.xy,x.time,x.timeVariance)!=
        std::tie(y.sensor,y.channels,y.x,y.y,y.xx,y.yy,y.xy,y.time,y.timeVariance)) return false;
  }
  return true;
}
}
int main() {
  try {
    std::vector<Sensor> sensors{{100,8,8},{200,8,8}};
    auto c=clean(); c.threshold=1500;
    auto sum=simulate(sensors,{{1,100,.25,.25,800,20},{2,100,.25,.25,1200,20}},c,10);
    check(sum.channels.size()==1 && sum.clusters.size()==1,"threshold must follow summation");
    near(sum.channels[0].charge,2000); near(sum.channels[0].truth[0].second,.4);
    near(sum.channels[0].truth[1].second,.6); near(sum.clusters[0].x,.25);
    c=clean();
    auto split=simulate(sensors,{{1,100,0,0,10000,20}},c,10);
    check(split.channels.size()==4 && split.clusters.size()==1,"four-pad boundary sharing");
    for(const auto& ch:split.channels) near(ch.charge,2500);
    near(split.clusters[0].x,0); near(split.clusters[0].y,0);
    check(split.clusters[0].xx*split.clusters[0].yy>std::pow(split.clusters[0].xy,2),"SPD covariance");
    auto edge=simulate(sensors,{{1,100,8,.25,10000,20}},c,10);
    near(edge.sharedSignal/edge.inputSignal,.5); check(edge.clusters[0].edge,"edge flag");
    check(edge.clusters[0].xx>split.clusters[0].xx,"edge uncertainty");
    auto different=simulate(sensors,{{1,100,.25,.25,1000,20},{2,200,.25,.25,1000,20}},c,10);
    check(different.clusters.size()==2,"must not merge distinct sensors");
    auto gate=simulate(sensors,{{1,100,.25,.25,1000,100},{2,100,.25,.25,0,20}},c,10);
    check(gate.outsideGate==1 && gate.zeroEnergy==1 && gate.channels.empty(),"gate and zero deposit");
    throws([&]{simulate(sensors,{{1,300,0,0,1,20}},c,10);});
    throws([&]{simulate(sensors,{{1,100,9,0,1,20}},c,10);});
    throws([&]{simulate(sensors,{{1,100,0,0,-1,20}},c,10);});
    throws([&]{simulate(sensors,{{1,100,0,0,1,20},{1,200,0,0,1,20}},c,10);});
    throws([&]{simulate({{100,8,8},{100,8,8}},{},c,10);});
    auto bad=c; bad.offsetX=0; throws([&]{simulate(sensors,{},bad,10);});
    bad=c; bad.gain=std::numeric_limits<double>::infinity(); throws([&]{simulate(sensors,{},bad,10);});
    bad=c; bad.maxWork=1; throws([&]{simulate(sensors,{{1,100,0,0,1,20}},bad,10);});
    bad=c; bad.maxChannels=1; throws([&]{simulate(sensors,{{1,100,0,0,10000,20}},bad,10);});
    c=clean(); c.maxADC=100;
    auto sat=simulate(sensors,{{1,100,.25,.25,1000,20}},c,10);
    check(sat.channels.size()==1 && sat.clusters.empty() && sat.rejectedSaturatedClusters==1,"saturation not silently fitted");
    c.rejectSaturatedClusters=false;
    check(simulate(sensors,{{1,100,.25,.25,1000,20}},c,10).clusters[0].saturated,"explicit saturated diagnostic");
    c=clean(); c.sensorTimeSigma=.03; c.clockTimeSigma=.02;
    auto timed=simulate(sensors,{{1,100,0,0,10000,20}},c,10);
    check(timed.clusters[0].timeVariance>=.03*.03+.02*.02,"common time errors incorrectly averaged down");
    near(timed.channels[0].time,timed.channels[3].time);
    c=Config{};
    std::vector<Deposit> ds{{3,200,1,.3,18000,20},{1,100,.3,.2,17000,20},{2,100,.35,.3,20000,20.1}};
    auto before=simulate(sensors,ds,c,22); std::reverse(ds.begin(),ds.end()); std::reverse(sensors.begin(),sensors.end());
    auto after=simulate(sensors,ds,c,22); check(equal(before,after),"input ordering changed digitization");
    check(!equal(before,simulate(sensors,ds,c,23)),"different event seeds unchanged");
    for(const auto& ch:before.channels) {
      double total=0; for(const auto& w:ch.truth) { check(w.second>=0,"negative truth weight"); total+=w.second; }
      check(total<=1+1e-12,"noise attributed to signal");
    }
    c=clean(); c.noiseOnlySensors=true; c.noiseSigma=1000; c.threshold=2000;
    auto noise=simulate(sensors,{},c,31); check(!noise.channels.empty(),"noise-only pads not generated");
    for(const auto& ch:noise.channels) check(ch.truth.empty(),"noise hit has fabricated truth");
    c=clean(); c.channelGainSigma=.2;
    auto fixed=simulate(sensors,{{1,100,.25,.25,1000,20}},c,31);
    auto fixed2=simulate(sensors,{{1,100,.25,.25,1000,20}},c,32);
    near(fixed.channels[0].charge,fixed2.channels[0].charge);
    c=clean(); c.clusterTimeWindow=.2;
    auto times=simulate(sensors,{{1,100,.25,.25,1000,20},{2,100,.75,.25,1000,22}},c,31);
    check(times.clusters.size()==2,"time-separated neighbors merged");
    c=clean(); c.timeWalk=100;
    auto walk=simulate(sensors,{{1,100,.25,.25,1000,20}},c,31);
    near(walk.clusters[0].time,20);
    std::cout<<checks<<" checks passed\n";
  } catch(const std::exception& e) { std::cerr<<e.what()<<'\n'; return 1; }
}
