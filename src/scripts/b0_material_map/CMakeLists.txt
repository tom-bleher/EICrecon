find_package(Acts REQUIRED COMPONENTS PluginRoot)

add_executable(eicrecon-b0-material-map b0_material_map.cc)
target_link_libraries(eicrecon-b0-material-map PRIVATE algorithms_tracking_library Acts::PluginRoot)
install(TARGETS eicrecon-b0-material-map DESTINATION ${CMAKE_INSTALL_BINDIR})
