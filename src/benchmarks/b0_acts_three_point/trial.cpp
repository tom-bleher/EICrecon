#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
#include <ActsPlugins/DD4hep/DD4hepFieldAdapter.hpp>
#include <DD4hep/Detector.h>
#include <DD4hep/Fields.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <string>
using V = Acts::Vector3;
int main(int argc, char** argv) try {
  if (argc < 4 || argc > 5) {
    std::cerr << "Usage: trial GEOMETRY.xml INPUT.txt OUTPUT.txt [STEP_MM]\n";
    return 3;
  }
  const double stepSize = argc > 4 ? std::stod(argv[4]) : 1.;
  if (!(stepSize > 0.) || !std::isfinite(stepSize))
    throw std::invalid_argument("Step must be finite and positive");
  // Known exact circles for both charges, including a rotated field.
  for (double sign : {-1., 1.})
    for (double angle : {0., 0.7}) {
      double R            = 20 / (sign * 0.000299792458 * 1.2);
      Eigen::Matrix3d rot = Eigen::AngleAxisd(angle, V::UnitZ()).toRotationMatrix();
      auto point          = [&](double s) -> V {
        return rot * V(R * (std::cos(s / R) - 1), 0, R * std::sin(s / R));
      };
      auto check = Acts::estimateTrackParamsFromSeed(point(0), point(500), point(1000),
                                                     rot * V(0, 1.2 * Acts::UnitConstants::T, 0));
      std::cerr << std::setprecision(16) << "uniform expected q/p=" << sign * .05
                << " actual=" << check[Acts::eFreeQOverP] << "\n";
      if (!check.allFinite() || std::abs(check[Acts::eFreeQOverP] - sign * .05) > 1e-9 ||
          (check.segment<3>(Acts::eFreeDir0) - rot * V(0, 0, 1)).norm() > 1e-9)
        return 2;
    }
  auto& det = dd4hep::Detector::getInstance();
  det.fromCompact(argv[1]);
  ActsPlugins::DD4hepFieldAdapter field(det.field());
  auto cache = field.makeCache(Acts::MagneticFieldContext{});
  auto bf    = [&](const V& p) -> V {
    V value = field.getField(p, cache).value();
    if (!value.allFinite())
      throw std::runtime_error("Nonfinite geometry field");
    return value;
  };
  std::ifstream in(argv[2]);
  if (!in)
    throw std::runtime_error("Cannot open input " + std::string(argv[2]));
  std::ofstream out(argv[3]);
  if (!out)
    throw std::runtime_error("Cannot open output " + std::string(argv[3]));
  out << std::setprecision(17);
  int id;
  V a, b, c, r, d, fieldpoint;
  double q;
  std::size_t rows = 0;
  while (in >> id) {
    if (id != static_cast<int>(rows))
      throw std::runtime_error("Nonsequential input row ID");
    for (V* v : {&a, &b, &c, &r, &d, &fieldpoint}) {
      for (int i = 0; i < 3; ++i)
        in >> (*v)[i];
      if (!in || !v->allFinite())
        throw std::runtime_error("Malformed/nonfinite vector in input");
    }
    in >> q;
    if (!in || !std::isfinite(q) || d.z() <= 0 || r.z() > a.z())
      throw std::runtime_error("Invalid initial state");
    V B      = bf(fieldpoint);
    auto est = Acts::estimateTrackParamsFromSeed(a, b, c, B);
    // RK4 in global z, 1 mm steps, no material. Current state at origin perigee -> first hit plane.
    auto deriv = [&](V x, V u) {
      if (!u.allFinite() || u.z() <= 1e-8)
        throw std::runtime_error("RK4 state is not forward-going");
      return std::pair<V, V>{V(u / u.z()), V(q * u.cross(bf(x)) / u.z())};
    };
    std::size_t steps = 0;
    while (r.z() < a.z() - 1e-8) {
      if (++steps > 1000000)
        throw std::runtime_error("RK4 step limit reached");
      double h      = std::min(stepSize, a.z() - r.z());
      auto [k1, l1] = deriv(r, d);
      auto [k2, l2] = deriv(r + h / 2 * k1, d + h / 2 * l1);
      auto [k3, l3] = deriv(r + h / 2 * k2, d + h / 2 * l2);
      auto [k4, l4] = deriv(r + h * k3, d + h * l3);
      r += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
      d += h / 6 * (l1 + 2 * l2 + 2 * l3 + l4);
      d.normalize();
    }
    out << id << ' ' << est[Acts::eFreeQOverP] << ' ' << est.segment<3>(Acts::eFreeDir0).transpose()
        << ' ' << d.transpose() << ' ' << (B / Acts::UnitConstants::T).transpose() << ' '
        << (r - a).norm() << '\n';
    if (!out)
      throw std::runtime_error("Output write failed");
    ++rows;
  }
  if (!in.eof())
    throw std::runtime_error("Malformed input row ID");
  if (rows == 0)
    throw std::runtime_error("Input contains no seeds");
  out.close();
  if (!out)
    throw std::runtime_error("Output close failed");
  std::cerr << "Completed " << rows << " seeds\n";
  return 0;
} catch (const std::exception& error) {
  std::cerr << "trial: " << error.what() << '\n';
  return 1;
}
