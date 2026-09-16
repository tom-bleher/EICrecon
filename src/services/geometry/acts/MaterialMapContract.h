// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#pragma once

#include <array>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace eicrecon::acts_material {

// Small dependency-free SHA-256 implementation used only to establish material-map
// identity at initialization. This is an integrity/provenance check, not a cryptographic
// authentication mechanism.
class Sha256 {
public:
  Sha256() = default;

  void update(const std::uint8_t* data, std::size_t size) {
    if (size > (std::numeric_limits<std::uint64_t>::max() - m_bits) / 8U) {
      throw std::overflow_error("SHA-256 input length overflow");
    }
    m_bits += static_cast<std::uint64_t>(size) * 8U;
    while (size != 0U) {
      const std::size_t take = std::min(size, m_block.size() - m_used);
      for (std::size_t i = 0; i < take; ++i) {
        m_block[m_used + i] = data[i];
      }
      m_used += take;
      data += take;
      size -= take;
      if (m_used == m_block.size()) {
        transform(m_block.data());
        m_used = 0;
      }
    }
  }

  std::array<std::uint8_t, 32> digest() const {
    Sha256 copy = *this;
    copy.m_block[copy.m_used++] = 0x80U;
    if (copy.m_used > 56U) {
      while (copy.m_used < copy.m_block.size()) {
        copy.m_block[copy.m_used++] = 0U;
      }
      copy.transform(copy.m_block.data());
      copy.m_used = 0;
    }
    while (copy.m_used < 56U) {
      copy.m_block[copy.m_used++] = 0U;
    }
    for (unsigned i = 0; i < 8U; ++i) {
      copy.m_block[63U - i] = static_cast<std::uint8_t>(copy.m_bits >> (8U * i));
    }
    copy.transform(copy.m_block.data());

    std::array<std::uint8_t, 32> result{};
    for (std::size_t i = 0; i < copy.m_state.size(); ++i) {
      result[4U * i]     = static_cast<std::uint8_t>(copy.m_state[i] >> 24U);
      result[4U * i + 1] = static_cast<std::uint8_t>(copy.m_state[i] >> 16U);
      result[4U * i + 2] = static_cast<std::uint8_t>(copy.m_state[i] >> 8U);
      result[4U * i + 3] = static_cast<std::uint8_t>(copy.m_state[i]);
    }
    return result;
  }

private:
  static constexpr std::array<std::uint32_t, 64> k{
      0x428a2f98U, 0x71374491U, 0xb5c0fbcfU, 0xe9b5dba5U, 0x3956c25bU, 0x59f111f1U,
      0x923f82a4U, 0xab1c5ed5U, 0xd807aa98U, 0x12835b01U, 0x243185beU, 0x550c7dc3U,
      0x72be5d74U, 0x80deb1feU, 0x9bdc06a7U, 0xc19bf174U, 0xe49b69c1U, 0xefbe4786U,
      0x0fc19dc6U, 0x240ca1ccU, 0x2de92c6fU, 0x4a7484aaU, 0x5cb0a9dcU, 0x76f988daU,
      0x983e5152U, 0xa831c66dU, 0xb00327c8U, 0xbf597fc7U, 0xc6e00bf3U, 0xd5a79147U,
      0x06ca6351U, 0x14292967U, 0x27b70a85U, 0x2e1b2138U, 0x4d2c6dfcU, 0x53380d13U,
      0x650a7354U, 0x766a0abbU, 0x81c2c92eU, 0x92722c85U, 0xa2bfe8a1U, 0xa81a664bU,
      0xc24b8b70U, 0xc76c51a3U, 0xd192e819U, 0xd6990624U, 0xf40e3585U, 0x106aa070U,
      0x19a4c116U, 0x1e376c08U, 0x2748774cU, 0x34b0bcb5U, 0x391c0cb3U, 0x4ed8aa4aU,
      0x5b9cca4fU, 0x682e6ff3U, 0x748f82eeU, 0x78a5636fU, 0x84c87814U, 0x8cc70208U,
      0x90befffaU, 0xa4506cebU, 0xbef9a3f7U, 0xc67178f2U};

  static std::uint32_t rotateRight(std::uint32_t value, unsigned shift) {
    return (value >> shift) | (value << (32U - shift));
  }

  void transform(const std::uint8_t* block) {
    std::array<std::uint32_t, 64> w{};
    for (std::size_t i = 0; i < 16U; ++i) {
      w[i] = (static_cast<std::uint32_t>(block[4U * i]) << 24U) |
             (static_cast<std::uint32_t>(block[4U * i + 1]) << 16U) |
             (static_cast<std::uint32_t>(block[4U * i + 2]) << 8U) |
             static_cast<std::uint32_t>(block[4U * i + 3]);
    }
    for (std::size_t i = 16; i < w.size(); ++i) {
      const std::uint32_t s0 = rotateRight(w[i - 15], 7U) ^ rotateRight(w[i - 15], 18U) ^
                               (w[i - 15] >> 3U);
      const std::uint32_t s1 = rotateRight(w[i - 2], 17U) ^ rotateRight(w[i - 2], 19U) ^
                               (w[i - 2] >> 10U);
      w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }

    std::uint32_t a = m_state[0];
    std::uint32_t b = m_state[1];
    std::uint32_t c = m_state[2];
    std::uint32_t d = m_state[3];
    std::uint32_t e = m_state[4];
    std::uint32_t f = m_state[5];
    std::uint32_t g = m_state[6];
    std::uint32_t h = m_state[7];

    for (std::size_t i = 0; i < w.size(); ++i) {
      const std::uint32_t big1 = rotateRight(e, 6U) ^ rotateRight(e, 11U) ^ rotateRight(e, 25U);
      const std::uint32_t choose = (e & f) ^ ((~e) & g);
      const std::uint32_t t1 = h + big1 + choose + k[i] + w[i];
      const std::uint32_t big0 = rotateRight(a, 2U) ^ rotateRight(a, 13U) ^ rotateRight(a, 22U);
      const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
      const std::uint32_t t2 = big0 + majority;
      h = g;
      g = f;
      f = e;
      e = d + t1;
      d = c;
      c = b;
      b = a;
      a = t1 + t2;
    }

    m_state[0] += a;
    m_state[1] += b;
    m_state[2] += c;
    m_state[3] += d;
    m_state[4] += e;
    m_state[5] += f;
    m_state[6] += g;
    m_state[7] += h;
  }

  std::array<std::uint32_t, 8> m_state{0x6a09e667U, 0xbb67ae85U, 0x3c6ef372U,
                                        0xa54ff53aU, 0x510e527fU, 0x9b05688cU,
                                        0x1f83d9abU, 0x5be0cd19U};
  std::array<std::uint8_t, 64> m_block{};
  std::size_t m_used{0};
  std::uint64_t m_bits{0};
};

inline std::string hexDigest(const std::array<std::uint8_t, 32>& digest) {
  std::ostringstream stream;
  stream << std::hex << std::setfill('0');
  for (const auto byte : digest) {
    stream << std::setw(2) << static_cast<unsigned>(byte);
  }
  return stream.str();
}

inline std::string sha256File(const std::filesystem::path& path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) {
    throw std::runtime_error("Cannot open material map for SHA-256: " + path.string());
  }
  Sha256 hasher;
  std::array<char, 65536> buffer{};
  while (input) {
    input.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
    const auto count = input.gcount();
    if (count > 0) {
      hasher.update(reinterpret_cast<const std::uint8_t*>(buffer.data()),
                    static_cast<std::size_t>(count));
    }
  }
  if (input.bad()) {
    throw std::runtime_error("Failed while reading material map for SHA-256: " + path.string());
  }
  return hexDigest(hasher.digest());
}

inline std::string normalizeSha256(std::string value) {
  if (value.size() != 64U) {
    throw std::invalid_argument("Expected material-map SHA-256 must contain 64 hexadecimal digits");
  }
  for (char& ch : value) {
    const auto uch = static_cast<unsigned char>(ch);
    if (!std::isxdigit(uch)) {
      throw std::invalid_argument("Expected material-map SHA-256 contains a non-hexadecimal digit");
    }
    ch = static_cast<char>(std::tolower(uch));
  }
  return value;
}

struct MaterialMapSelection {
  std::string geometryDeclaredPath;
  std::string selectedPath;
  std::string expectedSha256;
  bool geometryDeclared{false};
  bool usedFallback{false};
  bool requireGeometryMaterialMap{false};
};

struct MaterialMapValidation {
  std::string declaredSha256;
  std::string selectedSha256;
  std::string expectedSha256;
  bool comparedToGeometry{false};
};

inline MaterialMapValidation validateMaterialMapSelection(const MaterialMapSelection& selection) {
  MaterialMapValidation result;
  if (selection.requireGeometryMaterialMap &&
      (!selection.geometryDeclared || selection.geometryDeclaredPath.empty())) {
    throw std::runtime_error(
        "ACTS strict material-map policy requires the loaded DD4hep geometry to declare 'material-map'");
  }
  if ((selection.requireGeometryMaterialMap || !selection.expectedSha256.empty()) &&
      selection.selectedPath.empty()) {
    throw std::runtime_error("ACTS strict material-map policy selected an empty material-map path");
  }

  if (!selection.expectedSha256.empty()) {
    result.expectedSha256 = normalizeSha256(selection.expectedSha256);
    result.selectedSha256 = sha256File(selection.selectedPath);
    if (result.selectedSha256 != result.expectedSha256) {
      throw std::runtime_error("Selected ACTS material map does not match expected SHA-256");
    }
  }

  if (selection.requireGeometryMaterialMap) {
    result.declaredSha256 = sha256File(selection.geometryDeclaredPath);
    if (result.selectedSha256.empty()) {
      result.selectedSha256 = sha256File(selection.selectedPath);
    }
    if (!result.expectedSha256.empty() && result.declaredSha256 != result.expectedSha256) {
      throw std::runtime_error(
          "Geometry-declared ACTS material map does not match expected SHA-256");
    }
    if (result.selectedSha256 != result.declaredSha256) {
      throw std::runtime_error(
          "Selected ACTS material-map content differs from the geometry-declared map");
    }
    result.comparedToGeometry = true;
  }
  return result;
}

} // namespace eicrecon::acts_material
