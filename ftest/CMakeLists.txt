# utest/CMakeLists.txt

add_executable(arcos ftest.cpp)

target_include_directories(arcos PRIVATE ${CMAKE_SOURCE_DIR}/_deps/googletest-src/googletest/include)
target_link_libraries(arcos gtest_main)

