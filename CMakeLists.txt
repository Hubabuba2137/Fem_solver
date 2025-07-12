cmake_minimum_required(VERSION 3.24)
project(FEM) 

set(CMAKE_CXX_STANDARD 11)

include(FetchContent)

set(SOURCE_FILES
    main.cpp)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})


# Copy the data files
add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/Data/
    $<TARGET_FILE_DIR:${PROJECT_NAME}>/Data/)