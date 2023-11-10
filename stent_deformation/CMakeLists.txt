# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.

cmake_minimum_required(VERSION 3.1...3.23)
project(Surface_mesh_deformation_Examples)

find_package(CGAL REQUIRED)

find_package(Eigen3 3.1.91 QUIET) #(requires 3.1.91 or greater)
include(CGAL_Eigen3_support)
if(TARGET CGAL::Eigen3_support)
  create_single_source_cgal_program("stent_deformation.cpp")

  foreach(
    target
    stent_deformation)
    target_link_libraries(${target} PUBLIC CGAL::Eigen3_support)
  endforeach()

  find_package(OpenMesh QUIET)
  if(OpenMesh_FOUND)
    include(UseOpenMesh)
    message(STATUS "Found OpenMesh")
  else()
    message(STATUS "NOTICE: Examples that use OpenMesh will not be compiled.")
  endif()

else()
  message("NOTICE: These examples require Eigen 3.1.91 (or greater) or later and will not be compiled.")
endif()
