cmake_minimum_required(VERSION 2.6)
project(h2asymmerrors)

add_executable(h2asymmerrors main.cpp)

install(TARGETS h2asymmerrors RUNTIME DESTINATION bin)
