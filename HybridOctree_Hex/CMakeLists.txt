cmake_minimum_required(VERSION 3.16)
project(HexGen VERSION 1.0)

# Mesh generation is compute intensive. Preserve an explicitly selected build
# type, and let multi-configuration generators select it at build time.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        Debug Release RelWithDebInfo MinSizeRel)
endif()

set(SOURCES
    Main.cpp
    HexGen.cpp
    HexGen.h
    Initialization.h
    Mesh.cpp
    Mesh.h
    StaticVars.h
    TriangleBVH.h
    ProjectionControl.h
)

add_executable(HexGen ${SOURCES})
target_compile_features(HexGen PRIVATE cxx_std_11)

if (CMAKE_GENERATOR MATCHES "Visual Studio")
    set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT HexGen)
endif()

include(CTest)
if (BUILD_TESTING)
    add_executable(hexgen_regression tests/Regression.cpp Mesh.cpp)
    target_compile_features(hexgen_regression PRIVATE cxx_std_11)
    add_test(NAME geometry_queries COMMAND hexgen_regression geometry
        "${CMAKE_CURRENT_SOURCE_DIR}/../input boundaries/bone_tri.raw"
        "${CMAKE_CURRENT_SOURCE_DIR}/../input boundaries/sphereInShell_tri.raw")
    add_test(NAME sparse_balance COMMAND hexgen_regression balance)
    add_test(NAME projection_stopping COMMAND hexgen_regression projection)
    add_test(NAME mesh_pipeline COMMAND hexgen_regression pipeline)
    add_test(NAME mesh_pipeline_bone COMMAND hexgen_regression pipeline
        "${CMAKE_CURRENT_SOURCE_DIR}/../input boundaries/bone_tri.raw")
    foreach(test_name geometry_queries sparse_balance projection_stopping mesh_pipeline mesh_pipeline_bone)
        file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-runs/${test_name}")
        set_tests_properties(${test_name} PROPERTIES TIMEOUT 120
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/test-runs/${test_name}")
    endforeach()
endif()
