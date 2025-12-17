option(CONDA_BUILD "Building inside conda-build" OFF)

cmake_minimum_required(VERSION 4.0)
project(fastder)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Threads REQUIRED)

include_directories(cpp)

add_library(fastder_lib
        cpp/Averager.cpp
        cpp/BedGraphRow.cpp
        cpp/Parser.cpp
        cpp/SJRow.cpp
        cpp/StitchedER.cpp
        cpp/Integrator.cpp
        cpp/GTFRow.cpp
)


add_executable(fastder
        cpp/main.cpp
)
target_link_libraries(fastder PRIVATE fastder_lib Threads::Threads)

enable_testing()

if(NOT CONDA_BUILD)
       add_subdirectory(tests)
endif()

# packaging
set(CPACK_PACKAGE_NAME "fastder")
set(CPACK_PACKAGE_VENDOR "Martina Lehmann")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "fastder: A Fast, Agnostic, Splice-Junction-Based Tool for the Detection of Expressed Regions")
set(CPACK_PACKAGE_VERSION "0.0.1")
set(CPACK_PACKAGE_CONTACT "martina.lavanya@gmail.com")

# Debian-specific packaging settings
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Martina Lehmann")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libstdc++6")

include(CPack)
