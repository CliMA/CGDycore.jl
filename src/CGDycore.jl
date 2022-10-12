module CGDycore

using LinearAlgebra
using SparseArrays
using UnPack
using NCDatasets
using DelimitedFiles
using Dierckx
using MPI
using StaticArrays
using WriteVTK
using NetCDF

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
include("Grid/TransCart.jl")
include("Grid/TransSphere.jl")
include("Grid/cart2sphere.jl")
include("Grid/hS.jl")
include("Grid/sphere2cart.jl")
include("Grid/vtkWriteHex.jl")
include("Grid/Connectivity.jl")
include("Grid/Geometry.jl")
include("Grid/SubGrid.jl")
include("Grid/InputGrid.jl")

include("DG/DLagrange.jl")
include("DG/DerivativeMatrixSingle.jl")
include("DG/GaussLobattoQuad.jl")
include("DG/Lagrange.jl")

include("DyCore/Average.jl")
include("DyCore/AverageFB.jl")
include("DyCore/BoundaryW.jl")
include("DyCore/BoundaryWOutput.jl")
include("DyCore/Damping.jl")
include("DyCore/Discretization.jl")
# include("DyCore/Energy.jl")
include("DyCore/FCurlNon3Vec.jl")
include("DyCore/FDiv3Vec.jl")
include("DyCore/FDivRhoGrad2Vec.jl")
include("DyCore/FGrad3Vec.jl")
include("DyCore/FGradDiv2Vec.jl")
include("DyCore/FRotCurl2Vec.jl")
include("DyCore/FVort2VecDSS.jl")
include("DyCore/FcnNHCurlVec.jl")
include("DyCore/FcnTracer.jl")
include("DyCore/HorLimiter.jl")
include("DyCore/JacSchur.jl")
include("DyCore/MassCG.jl")
include("DyCore/NumberingFemCG.jl")
include("DyCore/Project.jl")
include("DyCore/ProjectW.jl")
include("DyCore/ProjectVec.jl")
include("DyCore/Source.jl")
include("DyCore/simpson.jl")
include("DyCore/VerticalDiffusionScalar.jl")
include("DyCore/vtkCG.jl")
include("DyCore/vtkCGGrid.jl")
include("DyCore/vtkOutput.jl")
include("DyCore/vtkSphere.jl")
include("DyCore/ThreadCache.jl")
include("DyCore/TopographySmoothing.jl")

include("IntegrationMethods/JacStruc.jl")
include("IntegrationMethods/LinIMEXMethod.jl")
include("IntegrationMethods/LinIMEXSchur.jl")
include("IntegrationMethods/IMEXMethod.jl")
include("IntegrationMethods/IMEXSchur.jl")
include("IntegrationMethods/MISMethod.jl")
include("IntegrationMethods/MISSchur.jl")
include("IntegrationMethods/RosenbrockMethod.jl")
include("IntegrationMethods/RosenbrockSchur.jl")
include("IntegrationMethods/RungeKuttaExplicit.jl")
include("IntegrationMethods/RungeKuttaMethod.jl")
include("IntegrationMethods/SchurSolve.jl")
include("IntegrationMethods/SSPRungeKutta.jl")

include("Parallel/Exchange.jl")
include("Parallel/Hilbert.jl")

include("Model/PhysParameters.jl")
include("Model/Pressure.jl")
include("Model/fRho.jl")
include("Model/fQv.jl")
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
include("Model/fVelW.jl")
include("Model/fpBGrd.jl")
include("Model/Energy.jl")
include("Model/KineticEnergy.jl")

include("Statistics/Averages.jl")

OOP = 5

end
