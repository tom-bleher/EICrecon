set(PLUGIN_NAME "algorithms_tracking")

# Function creates ${PLUGIN_NAME}_plugin and ${PLUGIN_NAME}_library targets
# Setting default includes, libraries and installation paths
plugin_add(${PLUGIN_NAME} WITH_SHARED_LIBRARY WITHOUT_PLUGIN)

# The macro grabs sources as *.cc *.cpp *.c and headers as *.h *.hh *.hpp Then
# correctly sets sources for ${_name}_plugin and ${_name}_library targets Adds
# headers to the correct installation directory
plugin_glob_all(${PLUGIN_NAME})

# Find dependencies
plugin_add_dd4hep(${PLUGIN_NAME})
plugin_add_acts(${PLUGIN_NAME})
plugin_add_cern_root(${PLUGIN_NAME})
plugin_add_eigen3(${PLUGIN_NAME})
plugin_add_event_model(${PLUGIN_NAME})

# The macro grabs sources as *.cc *.cpp *.c and headers as *.h *.hh *.hpp Then
# correctly sets sources for ${_name}_plugin and ${_name}_library targets Adds
# headers to the correct installation directory
plugin_glob_all(${PLUGIN_NAME})

plugin_include_directories(${PLUGIN_NAME} PUBLIC
                           $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

# Add libraries (same as target_include_directories but for both plugin and
# library)
plugin_link_libraries(${PLUGIN_NAME} Eigen3::Eigen particle_service_library)

# Use the stock multipoint API when available; retain ACTS 44 compatibility.
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
cmake_push_check_state(RESET)
set(CMAKE_REQUIRED_LIBRARIES Acts::Core)
set(_multipoint_probe
    "EICRECON_ACTS_${Acts_VERSION_MAJOR}_${Acts_VERSION_MINOR}_${Acts_VERSION_PATCH}_HAS_MULTIPOINT"
)
check_cxx_source_compiles(
  "#include <Acts/Seeding/EstimateTrackParamsFromSeed.hpp>
   int main() {
     return Acts::estimateTrackParamsFromSpacePoints(
       std::span<const Acts::Vector3>{}, Acts::Vector3::Zero(), 0., 0, {}, 0).ok();
   }"
  ${_multipoint_probe})
cmake_pop_check_state()
target_compile_definitions(
  ${PLUGIN_NAME}_library
  PRIVATE EICRECON_HAS_ACTS_MULTIPOINT=$<BOOL:${${_multipoint_probe}}>)
unset(_multipoint_probe)
