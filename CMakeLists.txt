cmake_minimum_required(VERSION 3.20)
project(Base_Project_Deluxe)

set(CMAKE_CXX_STANDARD 20)

# find the libraries:
find_library(GMSH_LIBRARY gmsh REQUIRED)
find_package(MFEM REQUIRED)

# include libraries:
include_directories(${MFEM_INCLUDE_DIRS})
set(MFEM_LIBS ${MFEM_LIBRARIES})
include_directories(/usr/local/include)
include_directories("/usr/include/eigen3")
# geometry-central is in /usr/local/include

add_executable(Base_Project_Deluxe main.cpp
        include/Agent.h
        src/Agent.cpp
        include/Space.h
        src/Space.cpp
        include/utils.h
        src/utils.cpp
        include/CollisionManager.h
        src/CollisionManager.cpp
        src/Field.cpp
        include/Field.h)

target_link_libraries(Base_Project_Deluxe PRIVATE ${MFEM_LIBS} ${GMSH_LIBRARY} PRIVATE
        geometry-central
        cholmod
        umfpack
        camd
        amd
        colamd
        ccolamd
        suitesparseconfig)