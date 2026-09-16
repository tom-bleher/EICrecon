cmake_minimum_required(VERSION 3.16)
project(material_map_contract_test LANGUAGES CXX)

enable_testing()

add_executable(material_map_contract_test test_material_map_contract.cpp)
target_compile_features(material_map_contract_test PRIVATE cxx_std_20)
target_include_directories(material_map_contract_test PRIVATE
  "${CMAKE_CURRENT_LIST_DIR}/../../services/geometry/acts")

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(material_map_contract_test PRIVATE
    -Wall -Wextra -Wpedantic -Werror)
endif()

add_test(NAME material_map_contract COMMAND material_map_contract_test)
