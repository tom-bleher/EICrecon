// SPDX-License-Identifier: LGPL-3.0-or-later
#include "B0StationMap.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

void check(bool value) { if (!value) throw std::runtime_error("station-map assertion failed"); }
int main() {
  using Map = eicrecon::b0::StationMap;
  std::vector<Map::Surface> surfaces{{91, 6000}, {23, 6007}, {5, 6276}, {1, 6283},
                                    {17, 6525}, {18, 6532}, {48, 6750}, {49, 6757}};
  Map map(surfaces);
  check(map.stations().size() == 4 && map.surfaceCount() == 8);
  check(map.stationForSurface(91) == 0 && map.stationForSurface(49) == 3);
  check(!map.stationForSurface(999));
  check(!map.stationForZ(6200));
  check(!map.stationForZ(std::numeric_limits<double>::quiet_NaN()));
  std::reverse(surfaces.begin(), surfaces.end());
  Map reordered(surfaces);
  for (const auto& s : surfaces) check(map.stationForSurface(s.id) == reordered.stationForSurface(s.id));
  surfaces.push_back(surfaces.front());
  check(Map(surfaces).surfaceCount() == 8);
  bool rejected = false;
  try { Map bad({{1, 1}, {1, 2}}); } catch (const std::invalid_argument&) { rejected = true; }
  check(rejected);
  rejected = false;
  try { Map bad({}, 0); } catch (const std::invalid_argument&) { rejected = true; }
  check(rejected);
  check(Map({}).stations().empty());
  check(Map({{9, 5902}, {8, 6172}, {7, 6442}, {6, 6712}}).stations().size() == 4);
  check(Map({{1, 0}, {2, 50}, {3, 100.01}}).stations().size() == 2);
  check(std::abs(eicrecon::b0::ionFrameZ(2, 3, 0) - 3) < 1e-12);
  std::cout << "station-map tests passed\n";
}
