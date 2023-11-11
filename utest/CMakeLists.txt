# utest/CMakeLists.txt

add_executable(Arcos test.cpp)

target_include_directories(Arcos PRIVATE ${CMAKE_SOURCE_DIR}/_deps/googletest-src/googletest/include)
target_link_libraries(Arcos gtest_main)

