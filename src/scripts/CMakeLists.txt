configure_file(eicrecon-this.sh.in eicrecon-this.sh @ONLY)

install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/eicrecon-this.sh DESTINATION bin)
install(PROGRAMS eicrecon-client.py DESTINATION bin)

option(BUILD_B0_MATERIAL_MAP "Build the experimental B0 window material mapper" OFF)
if(BUILD_B0_MATERIAL_MAP)
  add_subdirectory(b0_material_map)
endif()
