cmake_minimum_required(VERSION 3.10)
project(polynomial_Ramiisaac)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(polynomial_Ramiisaac
    .gitignore
    main.cpp
    makefile
    polynomial.cpp
    polynomial.h
    README.md
    unit-test.cpp)
