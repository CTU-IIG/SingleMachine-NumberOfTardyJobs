cmake_minimum_required(VERSION 3.5)
project(ImprovedBaptiste)

set(CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3 -mavx -ffast-math -march=native -pthread")

link_directories(/opt/lemon/lib)
include_directories(/opt/lemon/include)

link_directories(/opt/gurobi811/linux64/lib)
include_directories(/opt/gurobi811/linux64/include)
link_libraries(gurobi_c++5.5)
link_libraries(gurobi81)

include_directories(./libs/mcf_class/MCFClass)
include_directories(./libs/mcf_class/RelaxIV)
include_directories(./libs/mcf_class/MCFSimplex)
include_directories(./libs/mcf_class/OPTUtils)

add_executable(generator generator/generator.cpp)

add_executable(algorithm
        main_baptiste.cpp
        algorithm/Baptiste.cpp
        heuristic/lowerbound/ImprovedLowerBound.cpp
        heuristic/lowerbound/OriginalLowerBound.cpp
        heuristic/upperbound/ImprovedUpperBound.cpp
        heuristic/upperbound/ImprovedUpperBoundOpt.cpp
        libs/mcf_class/MCFSimplex/MCFSimplex.C
        libs/mcf_class/RelaxIV/RelaxIV.C
        reduction/JobReduce.cpp
        reduction/ReducedSchedulingProblem.cpp
        solver/ilp/BinaryILPSolution.cpp
        solver/ilp/ILPGridSolver.cpp
        solver/ilp/ILPSolver.cpp
        solver/ilp/ILPSolver.h
        solver/ilp/OriginILPSolver.cpp
        scheduling/SchedulingProblem.cpp
        solver/lp/GurobiLPSolver.cpp
        solver/lp/LPSolution.cpp
        solver/lp/LPSolutionWithBranchPenalty.cpp
        solver/lp/LPSolutionWithEarlyRange.cpp
        solver/lp/LPSolver.cpp
        solver/lp/lemon/LemonLPSolver.cpp
        solver/lp/mcf_class/MCFClassLPSolver.cpp
        solver/relaxation/GurobiWeightRelaxedILPSolver.cpp)

add_executable(lp_solver
        main_lp_solver.cpp
        heuristic/lowerbound/ImprovedLowerBound.cpp
        heuristic/lowerbound/OriginalLowerBound.cpp
        heuristic/upperbound/ImprovedUpperBoundOpt.cpp
        libs/mcf_class/MCFSimplex/MCFSimplex.C
        libs/mcf_class/RelaxIV/RelaxIV.C
        reduction/JobReduce.cpp
        reduction/ReducedSchedulingProblem.cpp
        solver/ilp/BinaryILPSolution.cpp
        solver/ilp/ILPGridSolver.cpp
        solver/ilp/OriginILPSolver.cpp
        scheduling/SchedulingProblem.cpp
        solver/lp/GurobiLPSolver.cpp
        solver/lp/LPSolution.cpp
        solver/lp/LPSolutionWithBranchPenalty.cpp
        solver/lp/LPSolutionWithEarlyRange.cpp
        solver/lp/LPSolver.cpp
        solver/lp/lemon/LemonLPSolver.cpp
        solver/lp/mcf_class/MCFClassLPSolver.cpp
        solver/relaxation/GurobiWeightRelaxedILPSolver.cpp)

add_executable(ilp_solver
        main_ilp_solver.cpp
        solver/ilp/BinaryILPSolution.cpp
        solver/ilp/ILPSolver.h
        solver/ilp/OriginILPSolver.cpp
        scheduling/SchedulingProblem.cpp
        solver/lp/LPSolution.cpp)

add_executable(ilp_solver_grid
        main_ilp_solver_grid.cpp
        solver/ilp/BinaryILPSolution.cpp
        solver/ilp/ILPGridSolver.cpp
        solver/ilp/ILPSolver.h
        solver/ilp/OriginILPSolver.cpp
        scheduling/SchedulingProblem.cpp
        solver/lp/GurobiLPSolver.cpp
        solver/lp/LPSolution.cpp
        solver/lp/LPSolutionWithEarlyRange.cpp
        solver/lp/LPSolver.cpp
        solver/relaxation/GurobiWeightRelaxedILPSolver.cpp)