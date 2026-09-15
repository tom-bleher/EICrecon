# Offline B0 tracking analysis (no JANA plugin).
add_test(
  NAME t_b0_tracking_analysis_selftest
  COMMAND ${CMAKE_COMMAND} -E env python3
          ${CMAKE_CURRENT_SOURCE_DIR}/analyze_b0_tracking.py --self-test)
