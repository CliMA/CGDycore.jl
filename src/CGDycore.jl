module CGDycore

using LinearAlgebra
using SparseArrays
using UnPack
using NCDatasets
using Dierckx
using StructArrays
using StaticArrays
using MPI
using WriteVTK
using ArgParse
#using StaticArrays
using NetCDF
using MuladdMacro
using Statistics
using StrideArraysCore: @gc_preserve, StrideArray, StaticInt
using RootSolvers

include("Grid/Geometry.jl")
include("Grid/Node.jl")
include("Grid/Edge.jl")
include("Grid/Face.jl")
include("Grid/GridStruct.jl")
include("Grid/AddVerticalGrid.jl")
include("Grid/CartGrid.jl")
include("Grid/CubedGrid.jl")
include("Grid/FacesInNodes.jl")
include("Grid/JacobiDG3.jl")
include("Grid/JacobiSphere3.jl")
include("Grid/OrientFaceCart.jl")
include("Grid/OrientFaceSphere.jl")
include("Grid/Orientation.jl")
include("Grid/Renumbering.jl")
include("Grid/Topo.jl")
include("Grid/TopoNeu.jl")
include("Grid/Trans.jl")
include("Grid/NearestNeighbour.jl")
include("Grid/vtkWriteHex.jl")
include("Grid/Connectivity.jl")
include("Grid/SubGrid.jl")
include("Grid/InputGrid.jl")

include("DG/DLagrange.jl")
include("DG/DerivativeMatrixSingle.jl")
include("DG/GaussLobattoQuad.jl")
include("DG/Lagrange.jl")
include("DG/Tools.jl")

include("DyCore/Average.jl")
include("DyCore/AverageFB.jl")
include("DyCore/Damping.jl")
include("DyCore/DiscretizationCG.jl")
include("DyCore/Operator.jl")
include("DyCore/Fcn.jl")
include("DyCore/FcnTracer.jl")
include("DyCore/FcnTestOperator.jl")
include("DyCore/HorLimiter.jl")
include("DyCore/JacSchur.jl")
include("DyCore/MassCG.jl")
include("DyCore/NumberingFemCG.jl")
include("DyCore/Project.jl")
include("DyCore/ProjectW.jl")
include("DyCore/ProjectVec.jl")
include("DyCore/Source.jl")
include("DyCore/simpson.jl")
include("DyCore/Diffusion.jl")
include("DyCore/vtkCG.jl")
include("DyCore/vtkCGGrid.jl")
include("DyCore/vtkSphere.jl")
include("DyCore/ThreadCache.jl")
include("DyCore/TopographySmoothing.jl")
include("DyCore/parse_commandline.jl")

include("IntegrationMethods/JacStruc.jl")
include("IntegrationMethods/LinIMEXMethod.jl")
include("IntegrationMethods/LinIMEXSchur.jl")
include("IntegrationMethods/IMEXMethod.jl")
include("IntegrationMethods/IMEXSchur.jl")
include("IntegrationMethods/MISMethod.jl")
include("IntegrationMethods/MISSchur.jl")
include("IntegrationMethods/SSPRungeKuttaMethod.jl")
include("IntegrationMethods/RosenbrockMethod.jl")
include("IntegrationMethods/RungeKuttaMethod.jl")
include("IntegrationMethods/RosenbrockSchur.jl")
include("IntegrationMethods/RungeKuttaExplicit.jl")
include("IntegrationMethods/SchurSolve.jl")
include("IntegrationMethods/SSPRungeKutta.jl")
include("IntegrationMethods/TimeStepper.jl")

include("Parallel/Exchange.jl")
include("Parallel/Hilbert.jl")
include("Parallel/EqualAreaPartitioner.jl")

include("Model/GlobalVariables.jl")
include("Model/Thermodynamics.jl")
include("Model/Pressure.jl")
include("Model/fRho.jl")
include("Model/fQv.jl")
include("Model/fQc.jl")
include("Model/fTSurf.jl")
include("Model/fRhoBGrd.jl")
include("Model/fT.jl")
include("Model/fTBGrd.jl")
include("Model/fTheta.jl")
include("Model/fIntEn.jl")
include("Model/fTotEn.jl")
include("Model/fThetaBGrd.jl")
include("Model/fTr.jl")
include("Model/fVel.jl")
include("Model/fVelGeo.jl")
include("Model/fPsi.jl")
include("Model/fVelW.jl")
include("Model/fpBGrd.jl")
include("Model/fTest.jl")
include("Model/KineticEnergy.jl")
include("Model/InitialConditions.jl")
include("Model/Parameters.jl")
include("Model/Fun.jl")
include("Model/InitProfileBryanFritsch.jl")

include("Statistics/Averages.jl")

include("DyCore/InitDriver.jl")
#include("GPU/CGGPU.jl")

OOP = 5

end
