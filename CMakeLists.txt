cmake_minimum_required(VERSION 2.6)
project(Sound)

add_definitions(-Wall -std=c++0x)

find_package(OpenAL REQUIRED)
include_directories(${OPENAL_INCLUDE_DIR})

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

find_package(FFTW REQUIRED)
include_directories(${FFTW_INCLUDE_DIRS})

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
#add_subdirectory(src)

add_executable(sound main.cpp)
target_link_libraries(sound ${OPENAL_LIBRARY} ${OpenCV_LIBS} ${FFTW_LIBS})
