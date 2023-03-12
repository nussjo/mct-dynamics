cmake_minimum_required(VERSION 3.25.0)
project(MCTDynamics)

add_library(${PROJECT_NAME}Core)

target_sources(${PROJECT_NAME}Core
    PUBLIC 
        ${CMAKE_CURRENT_SOURCE_DIR}/constants.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Correlator.h
        ${CMAKE_CURRENT_SOURCE_DIR}/mdefs.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Mesh.h
        ${CMAKE_CURRENT_SOURCE_DIR}/Utils.h
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/Correlator.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Mesh.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp
)

target_include_directories(${PROJECT_NAME}Core
    PUBLIC 
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/.>
        $<INSTALL_INTERFACE:.>
)


add_executable(${PROJECT_NAME} MCTdynamics.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE ${PROJECT_NAME}Core)
