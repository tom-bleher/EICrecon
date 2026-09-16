# Dependency-light B0 tests

Run without the EICrecon simulation stack:

```sh
cmake -S src/tests/b0_standalone -B build-b0 -DCMAKE_BUILD_TYPE=Release
cmake --build build-b0
ctest --test-dir build-b0 --output-on-failure
```

These tests use exceptions rather than `assert`, so Release/NDEBUG builds do not silently skip checks. GCC/Clang warnings are errors. The station map is independent of ACTS and DD4hep; full-stack tests of the ACTS station selector remain necessary. Later telescope-math tests require Eigen3; use `-DEIGEN3_INCLUDE_DIR=/path/to/eigen3` when it is not found automatically.

The new map and the legacy seeder use the same 50 mm single-linkage station definition. The legacy seeder's event fallback is retained as a historical baseline; the new local-telescope chain uses geometry IDs throughout.
