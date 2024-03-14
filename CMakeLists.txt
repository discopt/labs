cmake_minimum_required(VERSION 3.3)

project(LABS)
find_package(SCIP REQUIRED)

include_directories(${SCIP_INCLUDE_DIRS})

add_executable(labs
  main.cpp
  cons_labsnogood.cpp
  )

target_link_libraries(labs ${SCIP_LIBRARIES})

