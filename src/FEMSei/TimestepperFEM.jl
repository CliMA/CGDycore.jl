#TimeStepper for Vectorinvariant Linear Shallow Water Equations
function TimeStepper(backend,FTB,U,dtau,Fcn,Model,Grid,nQuadM,nQuadS,Jacobi,nAdveVel,FileNameOutput,Proc,ProcNumber,cName,nPrint,ref)

  pPosS = Model.pPosS
  pPosE = Model.pPosE
  uPosS = Model.uPosS
  uPosE = Model.uPosE
  @views Up = U[pPosS:pPosE]
  @views Uu = U[uPosS:uPosE]
  Flat = false
  vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat;Refine=ref)
  
  # Output of the initial values
  FileNumber = 0
  NumRefine = size(vtkSkeletonMesh.RefineMidPoints,1)
  VelSp = zeros(Grid.NumFaces*NumRefine,2)
  hout = zeros(Grid.NumFaces*NumRefine)
  Vort = zeros(Grid.NumFaces*NumRefine)
  UCache = similar(U)
 
  uCurl = zeros(Model.DG.NumG)
  ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  ConvertVelocitySp!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  Vorticity!(backend,FTB,Vort,Model.DG,Uu,Model.RT,Model.ND,Model.Curl,Grid,Grid.Type,nQuadS,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[hout Vort VelSp],
    FileNumber,cName)

  time = 0.0
  UNew = similar(U)
  F = similar(U)
  @views Fp = F[pPosS:pPosE]
  @views Fu = F[uPosS:uPosE]
  time = 0.0
  for i = 1 : nAdveVel
    @show i,(i-1)*dtau/3600   
    Fcn(backend,FTB,F,U,Model,Grid,nQuadM,nQuadS,Jacobi;UCache)
    @show sum(abs.(F))
    @. UNew = U + 1/3 * dtau * F
    Fcn(backend,FTB,F,UNew,Model,Grid,nQuadM,nQuadS,Jacobi;UCache)
    @. UNew = U + 1/2 * dtau * F
    Fcn(backend,FTB,F,UNew,Model,Grid,nQuadM,nQuadS,Jacobi;UCache)
    @. U = U + dtau * F
    if mod(i,nPrint) == 0
      @show "Druck ",i  
      ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      ConvertVelocitySp!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      Vorticity!(backend,FTB,Vort,Model.DG,Uu,Model.RT,Model.ND,Model.Curl,Grid,Grid.Type,nQuadS,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      FileNumber += 1
      Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[hout Vort VelSp],
        FileNumber,cName)
    end  
    time += dtau
  end
  ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  ConvertVelocitySp!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  Vorticity!(backend,FTB,Vort,Model.DG,Uu,Model.RT,Model.ND,Model.Curl,Grid,Grid.Type,nQuadS,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  FileNumber += 1
  Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[hout Vort VelSp],
    FileNumber,cName)
end

function TimeStepperEul(backend,FTB,U,dtau,Fcn,Model,Grid,nQuadM,nQuadS,Jacobi,nAdveVel,GridType,Proc,ProcNumber)

  pPosS = Model.pPosS
  pPosE = Model.pPosE
  uPosS = Model.uPosS
  uPosE = Model.uPosE
  @views Up = U[pPosS:pPosE]
  if uPosS > 0
    @views Uu = U[uPosS:uPosE]
  end  
  vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid)
  
  FileNumber = 0
  if uPosS > 0
    VelCa = zeros(Grid.NumFaces,Grid.Dim)
    VelSp = zeros(Grid.NumFaces,2)
  end  
  pC = zeros(Grid.NumFaces,1)
  ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi!)
  if uPosS > 0
    ConvertVelocityCart!(backend,FTB,VelCa,Uu,Model.RT,Grid,Jacobi!)
    ConvertVelocitySp!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi!)
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [hout VelCa VelSp], FileNumber)
  else
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [hout;], FileNumber)  
  end  

  time = 0.0
  UCache = similar(U)
  F = similar(U)
  time = 0.0
  for i = 1 : nAdveVel
    @show i,time  
    Fcn(backend,FTB,F,U,Model,Grid,nQuadS,nQuadM,Jacobi!;UCache)
    @. U = U + dtau * F
    if mod(i,1) == 0
      ConvertScalar!(backend,FTB,pC,Up,Model.DG,Grid,Jacobi)
      if uPosS > 0
        ConvertVelocityCart!(backend,FTB,VelCa,Uu,Model.RT,Grid,Jacobi!)
        ConvertVelocitySp!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi!)
        FileNumber += 1
        Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [pC VelCa VelSp], FileNumber)
      else
        FileNumber += 1
        Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [pC;], FileNumber)  
      end  
    end  
    time += dtau
  end
  ConvertScalar!(backend,FTB,pC,Up,Model.DG,Grid,Jacobi!)
  if uPosS > 0
    ConvertVelocityCart!(backend,FTB,VelCa,Uu,Model.RT,Grid,Jacobi!)
    ConvertVelocitySp!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi!)
    FileNumber += 1
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [pC VelCa VelSp], FileNumber)
  else
    FileNumber += 1
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [pC;], FileNumber)  
  end  
end

