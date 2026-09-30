// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2024 Souvik Paul

#pragma once

namespace eicrecon {

struct SiliconChargeSharingConfig {
  // Parameters of Silicon signal generation
  // determines the meaning of sigma_sharingx and y.
  // rel means relative, so charge sharing range = sigma_sharingx * cell_width_x, etc for y
  // abs means absolute, charge sharing range = sigma_shargeingx directly
  enum class ESigmaMode { abs = 0, rel = 1 } sigma_mode = ESigmaMode::abs;
  float sigma_sharingx;
  float sigma_sharingy;
  float min_edep;
  std::string readout;
  // Pad response: gaussian spreads the deposit with the sigmas above; trapezoid is the AC-LGAD
  // response of TrapezoidPadResponse.h for metal electrodes of electrode_x by electrode_y, with
  // a fraction gap_loss of the signal lost at mid-gap; table reads a PadResponseTable file
  enum class EModel { gaussian = 0, trapezoid = 1, table = 2 } model = EModel::gaussian;
  float electrode_x                                                  = 0;
  float electrode_y                                                  = 0;
  float gap_loss                                                     = 0;
  std::string table_file;
};

inline std::istream& operator>>(std::istream& in, SiliconChargeSharingConfig::EModel& model) {
  std::string s;
  in >> s;
  if (s == "gaussian" or s == "0") {
    model = SiliconChargeSharingConfig::EModel::gaussian;
  } else if (s == "trapezoid" or s == "1") {
    model = SiliconChargeSharingConfig::EModel::trapezoid;
  } else if (s == "table" or s == "2") {
    model = SiliconChargeSharingConfig::EModel::table;
  } else {
    in.setstate(std::ios::failbit);
  }
  return in;
}
inline std::ostream& operator<<(std::ostream& out,
                                const SiliconChargeSharingConfig::EModel& model) {
  switch (model) {
  case SiliconChargeSharingConfig::EModel::gaussian:
    out << "gaussian";
    break;
  case SiliconChargeSharingConfig::EModel::trapezoid:
    out << "trapezoid";
    break;
  case SiliconChargeSharingConfig::EModel::table:
    out << "table";
    break;
  default:
    out.setstate(std::ios::failbit);
  }
  return out;
}

std::istream& operator>>(std::istream& in, SiliconChargeSharingConfig::ESigmaMode& sigmaMode) {
  std::string s;
  in >> s;
  // stringifying the enums causes them to be converted to integers before conversion to strings
  if (s == "abs" or s == "0") {
    sigmaMode = SiliconChargeSharingConfig::ESigmaMode::abs;
  } else if (s == "rel" or s == "1") {
    sigmaMode = SiliconChargeSharingConfig::ESigmaMode::rel;
  } else {
    in.setstate(std::ios::failbit); // Set the fail bit if the input is not valid
  }

  return in;
}
std::ostream& operator<<(std::ostream& out,
                         const SiliconChargeSharingConfig::ESigmaMode& sigmaMode) {
  switch (sigmaMode) {
  case SiliconChargeSharingConfig::ESigmaMode::abs:
    out << "abs";
    break;
  case SiliconChargeSharingConfig::ESigmaMode::rel:
    out << "rel";
    break;
  default:
    out.setstate(std::ios::failbit);
  }
  return out;
}

} // namespace eicrecon
