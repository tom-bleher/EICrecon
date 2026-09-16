// SPDX-License-Identifier: LGPL-3.0-or-later
// Copyright (C) 2026 ePIC Collaboration

#include "MaterialMapContract.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {

void require(bool condition, std::string_view message) {
  if (!condition) {
    throw std::runtime_error(std::string(message));
  }
}

template <typename Exception, typename Callable>
void requireThrows(Callable&& callable, std::string_view message) {
  bool caught = false;
  try {
    callable();
  } catch (const Exception&) {
    caught = true;
  }
  require(caught, message);
}

void writeFile(const std::filesystem::path& path, std::string_view content) {
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) {
    throw std::runtime_error("Cannot create test file: " + path.string());
  }
  out.write(content.data(), static_cast<std::streamsize>(content.size()));
  if (!out) {
    throw std::runtime_error("Cannot write test file: " + path.string());
  }
}

std::string sha256Text(std::string_view text) {
  eicrecon::acts_material::Sha256 hasher;
  hasher.update(reinterpret_cast<const std::uint8_t*>(text.data()), text.size());
  return eicrecon::acts_material::hexDigest(hasher.digest());
}

} // namespace

int main() {
  using namespace eicrecon::acts_material;

  require(sha256Text("") ==
              "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
          "SHA-256 empty-string vector failed");
  require(sha256Text("abc") ==
              "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
          "SHA-256 abc vector failed");

  const auto directory = std::filesystem::temp_directory_path() / "eicrecon_material_map_contract_test";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);
  const auto declared = directory / "declared.cbor";
  const auto relocated = directory / "relocated.cbor";
  const auto mismatch = directory / "mismatch.cbor";
  writeFile(declared, "same validated material map bytes\n");
  writeFile(relocated, "same validated material map bytes\n");
  writeFile(mismatch, "different material map bytes\n");

  const auto expected = sha256File(declared);
  require(expected == sha256File(relocated), "identical files should hash identically");
  require(expected != sha256File(mismatch), "different files should not hash identically");

  {
    MaterialMapSelection selection{
        .geometryDeclaredPath = declared.string(),
        .selectedPath = relocated.string(),
        .expectedSha256 = {},
        .geometryDeclared = true,
        .usedFallback = false,
        .requireGeometryMaterialMap = true,
    };
    const auto result = validateMaterialMapSelection(selection);
    require(result.comparedToGeometry, "strict geometry comparison was not recorded");
    require(result.declaredSha256 == expected, "declared map hash mismatch");
    require(result.selectedSha256 == expected, "relocated identical map should be accepted");
  }

  requireThrows<std::runtime_error>(
      [&] {
        MaterialMapSelection selection{
            .geometryDeclaredPath = declared.string(),
            .selectedPath = mismatch.string(),
            .expectedSha256 = {},
            .geometryDeclared = true,
            .usedFallback = false,
            .requireGeometryMaterialMap = true,
        };
        (void)validateMaterialMapSelection(selection);
      },
      "strict policy should reject different selected content");

  requireThrows<std::runtime_error>(
      [&] {
        MaterialMapSelection selection{
            .geometryDeclaredPath = {},
            .selectedPath = declared.string(),
            .expectedSha256 = {},
            .geometryDeclared = false,
            .usedFallback = true,
            .requireGeometryMaterialMap = true,
        };
        (void)validateMaterialMapSelection(selection);
      },
      "strict policy should require a geometry declaration");

  {
    std::string uppercase = expected;
    for (char& ch : uppercase) {
      if (ch >= 'a' && ch <= 'f') {
        ch = static_cast<char>(ch - 'a' + 'A');
      }
    }
    MaterialMapSelection selection{
        .geometryDeclaredPath = {},
        .selectedPath = relocated.string(),
        .expectedSha256 = uppercase,
        .geometryDeclared = false,
        .usedFallback = false,
        .requireGeometryMaterialMap = false,
    };
    const auto result = validateMaterialMapSelection(selection);
    require(result.expectedSha256 == expected, "expected SHA-256 should be normalized to lowercase");
    require(result.selectedSha256 == expected, "expected SHA-256 should validate selected map");
    require(!result.comparedToGeometry, "hash-only policy should not claim geometry comparison");
  }

  requireThrows<std::invalid_argument>(
      [&] {
        MaterialMapSelection selection{
            .selectedPath = declared.string(),
            .expectedSha256 = "not-a-sha256",
        };
        (void)validateMaterialMapSelection(selection);
      },
      "malformed expected SHA-256 should be rejected");

  requireThrows<std::runtime_error>(
      [&] {
        MaterialMapSelection selection{
            .selectedPath = mismatch.string(),
            .expectedSha256 = expected,
        };
        (void)validateMaterialMapSelection(selection);
      },
      "wrong but well-formed expected SHA-256 should reject selected content");

  std::filesystem::remove_all(directory);
  std::cout << "material-map contract tests passed\n";
  return 0;
}
