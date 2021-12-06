# - Example project for accessing SuperNEMO RTD data (Raw Trigger Data)
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "The project requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "  CMakeCache.txt")
    message(STATUS "  CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there, pointing it to:")
    message(STATUS "  ${CMAKE_SOURCE_DIR}")
    message(FATAL_ERROR "in-source build detected")
endif()

# Project setup
cmake_minimum_required(VERSION 3.3 FATAL_ERROR)
project(snfee_examples_read_red VERSION 0.1.0)

# - Load Builtin/Custom Modules
include(GNUInstallDirs)
list(INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)
message(STATUS "[info] CMAKE_MODULE_PATH = '${CMAKE_MODULE_PATH}'")

find_package(SNFrontEndElectronics REQUIRED)
include_directories(${SNFrontEndElectronics_INCLUDE_DIRS})

# - Executable:
add_executable(read_red read_red.cxx)

target_link_libraries(read_red SNFrontEndElectronics::snfee)
target_link_libraries(read_red Boost::program_options)
target_link_libraries(read_red ROOT::Core)
target_link_libraries(read_red ROOT::Tree)
target_link_libraries(read_red ROOT::MathCore)
target_link_libraries(read_red ROOT::Hist)
target_link_libraries(read_red ROOT::Gpad)

set_target_properties(read_red
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
        )

install(FILES
        ${CMAKE_BINARY_DIR}/read_red
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_BINDIR}
        PERMISSIONS
        OWNER_READ OWNER_WRITE OWNER_EXECUTE
        GROUP_READ GROUP_EXECUTE
        WORLD_READ WORLD_EXECUTE
        )
# end
