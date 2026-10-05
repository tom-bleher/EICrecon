# Validates the B0 tracker DD4hep/ACTS geometry: every sensitive sensor needs
# an ACTS surface inside its tracking volume, and readout cell IDs must
# roundtrip through both DD4hep and ACTS. With a material map, every B0
# layer must resolve mapped material: at least one concrete mapped approach
# surface per layer and no leftover prototype material (individual
# null-material approaches are allowed).
# Scope: geometry/mapping consistency, not a material-budget validation.
add_executable(b0_geometry_test b0_geometry_test.cc)

target_link_libraries(
  b0_geometry_test
  PRIVATE algorithms_tracking_library
          DD4hep::DDCore
          DD4hep::DDRec
          podio::podio
          EDM4EIC::edm4eic
          EDM4HEP::edm4hep)

# Add Acts version definitions to test executable
target_compile_definitions(
  b0_geometry_test PRIVATE "Acts_VERSION_MAJOR=${Acts_VERSION_MAJOR}"
                           "Acts_VERSION_MINOR=${Acts_VERSION_MINOR}")

install(TARGETS b0_geometry_test DESTINATION bin)

# Without a compact file the test exits 77; CTest reports that as skipped
# rather than failed, so plain `ctest` stays green.
add_test(NAME t_b0_geometry_test COMMAND b0_geometry_test)
set_tests_properties(t_b0_geometry_test PROPERTIES SKIP_RETURN_CODE 77)
