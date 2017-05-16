cmake_minimum_required(VERSION 3.7)
project(FypSim)

include_directories(C:/Users/M/Documents/R/win-library/3.4/Rcpp/include)
include_directories(C:/PROGRA~1/R/R-3.4.0/include)
include_directories(C:/Users/M/Documents/R/win-library/3.4/RcppArmadillo/include)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/procedure.cpp
        src/Calibration.cpp
        src/Calibration.h
        src/Design.cpp
        src/Design.h
        src/RcppExports.cpp
        src/Cell.h
        src/Statistics.cpp
        src/Statistics.h
        src/StudyPhase.cpp
        src/StudyPhase.h
        src/EnvironmentCallsR.cpp
        src/EnvironmentCallsR.h
        src/DesignProcedure.cpp
        src/DesignProcedure.h)

add_executable(FypSim ${SOURCE_FILES})