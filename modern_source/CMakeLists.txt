cmake_minimum_required(VERSION 3.14)
project(mpass_modern CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Main executable ──────────────────────────────────────────────────────────
add_executable(pass_modern pass_modern.cpp)
target_include_directories(pass_modern PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# ── Test executable ──────────────────────────────────────────────────────────
add_executable(run_tests tests_modern.cpp)
target_include_directories(run_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

enable_testing()
add_test(NAME mpass_tests COMMAND run_tests)
