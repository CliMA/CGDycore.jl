#TimeStepper for Linear Shallow Water Equations
function TimeStepperLin(backend,FTB,U,dtau,Fcn,Model,Grid,nQuadM,nQuadS,Jacobi,nAdveVel,FileNameOutput,Proc,ProcNumber,nPrint,Flat,ref)

  pPosS = Model.pPosS
  pPosE = Model.pPosE
  uPosS = Model.uPosS
  uPosE = Model.uPosE
  @views Up = U[pPosS:pPosE]
  @views Uu = U[uPosS:uPosE]
  UCache = zeros(Model.DG.NumG+Model.RT.NumG)

  vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat;Refine=ref)
  
  # Output of the initial values
  FileNumber = 0
  NumRefine = size(vtkSkeletonMesh.RefineMidPoints,1)
  if Grid.Form == Grids.SphericalGrid()
    VelOut = zeros(Grid.NumFaces*NumRefine,2)
    hout = zeros(Grid.NumFaces*NumRefine)
    Vort = zeros(Grid.NumFaces*NumRefine)
    cName = ["Vort";"h";"uS";"vS"]
  else  
    VelOut = zeros(Grid.NumFaces*NumRefine,3)
    hout = zeros(Grid.NumFaces*NumRefine)
    Vort = zeros(Grid.NumFaces*NumRefine)
    cName = ["Vort";"h";"uC";"vC";"wC"]
  end  
 
  ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  ConvertVelocity!(backend,FTB,VelOut,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints,
    Grid.Form)
  Vorticity!(backend,FTB,Vort,Model.CG,Uu,Model.RT,Grid,Grid.Type,nQuadM,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[Vort hout VelOut],
    FileNumber,cName)

  time = 0.0
  UNew = similar(U)
  F = similar(U)
  time = 0.0
  for i = 1 : nAdveVel
    @show i,(i-1)*dtau/3600   
    Fcn(backend,FTB,F,U,Model,Grid,nQuadM,nQuadS,Jacobi;UCache)
    @. UNew = U + 1/3 * dtau * F
    Fcn(backend,FTB,F,UNew,Model,Grid,nQuadM,nQuadS,Jacobi;UCache)
    @. UNew = U + 1/2 * dtau * F
    Fcn(backend,FTB,F,UNew,Model,Grid,nQuadM,nQuadS,Jacobi;UCache)
    @. U = U + dtau * F

    # Output
    if mod(i,nPrint) == 0
      @show "Druck ",i  
      ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      ConvertVelocity!(backend,FTB,VelOut,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints,
        Grid.Form)
      Vorticity!(backend,FTB,Vort,Model.CG,Uu,Model.RT,Grid,Grid.Type,nQuadM,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      FileNumber += 1
      Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[Vort hout VelOut],
        FileNumber,cName)
    end  
    time += dtau
  end
  ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  ConvertVelocity!(backend,FTB,VelOut,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints,
    Grid.Form)
  Vorticity!(backend,FTB,Vort,Model.CG,Uu,Model.RT,Grid,Grid.Type,nQuadM,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  FileNumber += 1
  Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[Vort hout VelOut],
    FileNumber,cName)
end

#TimeStepper for Vectorinvariant Nonlinear Shallow Water Equations
function TimeStepperVecI(backend,FTB,U,dtau,Fcn,Model,ExchangeDG,ExchangeCG,Grid,
  nQuadM,nQuadS,Jacobi,nAdveVel,FileNameOutput,Proc,ProcNumber,nPrint,Flat,ref)

  pPosS = Model.pPosS
  pPosE = Model.pPosE
  uPosS = Model.uPosS
  uPosE = Model.uPosE
  @views Up = U[pPosS:pPosE]
  @views Uu = U[uPosS:uPosE]
  UCache = zeros(Model.DG.NumG+Model.RT.NumG+Model.CG.NumG)

  vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat;Refine=ref)
  
  # Output of the initial values
  FileNumber = 0
  NumRefine = size(vtkSkeletonMesh.RefineMidPoints,1)
  if Grid.Form == Grids.SphericalGrid()
    VelOut = zeros(Grid.NumFaces*NumRefine,2)
    hout = zeros(Grid.NumFaces*NumRefine)
    Vort = zeros(Grid.NumFaces*NumRefine)
    cName = ["h";"Vort";"uS";"vS"]
  else  
    VelOut = zeros(Grid.NumFaces*NumRefine,3)
    hout = zeros(Grid.NumFaces*NumRefine)
    Vort = zeros(Grid.NumFaces*NumRefine)
    cName = ["h";"Vort";"uC";"vC";"wC"]
  end  
 
  ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  ConvertVelocity!(backend,FTB,VelOut,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints,
    Grid.Form)
  Vorticity!(backend,FTB,Vort,Model.CG,Uu,Model.RT,Grid,Grid.Type,nQuadM,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[hout Vort VelOut],
    FileNumber,cName)

  time = 0.0
  UNew = similar(U)
  F = zeros(FTB,Model.DG.NumG+Model.RT.NumG)
  time = 0.0
  for i = 1 : nAdveVel
    if Proc == 1  
      @show i,(i-1)*dtau/3600   
    end  
    Fcn(backend,FTB,F,U,Model,Grid,nQuadM,nQuadS,Jacobi,ExchangeDG,ExchangeCG;UCache)
    @. UNew = U + 1/3 * dtau * F
    Fcn(backend,FTB,F,UNew,Model,Grid,nQuadM,nQuadS,Jacobi,ExchangeDG,ExchangeCG;UCache)
    @. UNew = U + 1/2 * dtau * F
    Fcn(backend,FTB,F,UNew,Model,Grid,nQuadM,nQuadS,Jacobi,ExchangeDG,ExchangeCG;UCache)
    @. U = U + dtau * F

    # Output
    if mod(i,nPrint) == 0
      if Proc == 1  
        @show "Druck",i
      end  
      ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      ConvertVelocity!(backend,FTB,VelOut,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints,
        Grid.Form)
      Vorticity!(backend,FTB,Vort,Model.CG,Uu,Model.RT,Grid,Grid.Type,nQuadM,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      FileNumber += 1
      Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[hout Vort VelOut],
        FileNumber,cName)
      MPI.Barrier(MPI.COMM_WORLD)
    end  
    time += dtau
  end
  ConvertScalar!(backend,FTB,hout,Up,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  ConvertVelocity!(backend,FTB,VelOut,Uu,Model.RT,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints,
    Grid.Form)
  Vorticity!(backend,FTB,Vort,Model.CG,Uu,Model.RT,Grid,Grid.Type,nQuadM,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  FileNumber += 1
  Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[hout Vort VelOut],
    FileNumber,cName)
end

#TimeStepper for Conservative Nonlinear Shallow Water Equations
function TimeStepperCons(backend,FTB,U,dtau,Fcn,Model,Grid,nQuad,nQuadM,nQuadS,Jacobi,nAdveVel,FileNameOutput,Proc,ProcNumber,nPrint,Flat,ref)
time = 0.0

  hPosS = Model.hPosS # hPosS
  hPosE = Model.hPosE # hPosE
  huPosS = Model.huPosS # huPosS
  huPosE = Model.huPosE # huPosE
  uPosVec = Model.uPosVec
  @views Uh = U[hPosS:hPosE]
  @views Uhu = U[huPosS:huPosE]

  UNew = zeros(FTB,huPosE)
  @views UNewh = UNew[hPosS:hPosE]  
  @views UNewhu = UNew[huPosS:huPosE]  
  @views UNewu = UNew[huPosS:huPosE]  
  F =  zeros(FTB,huPosE)
  @views Fh = F[hPosS:hPosE]  
  @views Fhu = F[huPosS:huPosE]  
  uRec = zeros(FTB,uPosVec)

  vtkSkeletonMesh = Outputs.vtkStruct{Float64}(backend,Grid,Grid.NumFaces,Flat;Refine=ref)
  
  # Output of the initial values
  FileNumber = 0
  NumRefine = size(vtkSkeletonMesh.RefineMidPoints,1)
  if Grid.Form == Grids.SphericalGrid()
    VelOut = zeros(Grid.NumFaces*NumRefine,2)
    hout = zeros(Grid.NumFaces*NumRefine)
    Vort = zeros(Grid.NumFaces*NumRefine)
    cName = ["h";"Vort";"uS";"vS"]
  else  
    VelOut = zeros(Grid.NumFaces*NumRefine,3)
    hout = zeros(Grid.NumFaces*NumRefine)
    Vort = zeros(Grid.NumFaces*NumRefine)
    cName = ["h";"Vort";"uC";"vC";"wC"]
  end  

  ConvertScalar!(backend,FTB,hout,Uh,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  ConvertScalarVelocity!(backend,FTB,VelOut,Uhu,Model.RT,Uh,Model.DG,Grid,Jacobi,
    vtkSkeletonMesh.RefineMidPoints,Grid.Form)
  ProjectHDivScalarHDiv(backend,FTB,UNewu,Model.RT,Uh,Model.DG,Uhu,Model.RT,
  Grid,nQuadM,Jacobi)
  Vorticity!(backend,FTB,Vort,Model.CG,UNewu,Model.RT,Grid,Grid.Type,nQuadM,Jacobi,vtkSkeletonMesh.RefineMidPoints)
  Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[hout Vort VelOut],FileNumber,cName)

  for i = 1 : nAdveVel
    @show i,(i-1)*dtau/3600 
    @. F = 0  
    FcnConsNonLinShallow!(backend, FTB, F, U, Model, Grid, nQuad, Jacobi)
    @. UNew = U + 1 / 3 * dtau * F

    @. F = 0  
    FcnConsNonLinShallow!(backend, FTB, F, UNew, Model, Grid, nQuad, Jacobi)
    @. UNew = U + 0.5 * dtau * F

    @. F = 0  
    FcnConsNonLinShallow!(backend, FTB, F, UNew, Model, Grid, nQuad, Jacobi)
    @. U = U + dtau * F
      
    # Output
    if mod(i,nPrint) == 0 
      FileNumber += 1
      @show "Print",i
      ConvertScalar!(backend,FTB,hout,Uh,Model.DG,Grid,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      ConvertScalarVelocity!(backend,FTB,VelOut,Uhu,Model.RT,Uh,Model.DG,Grid,Jacobi,
        vtkSkeletonMesh.RefineMidPoints,Grid.Form)
      ProjectHDivScalarHDiv(backend,FTB,UNewu,Model.RT,Uh,Model.DG,Uhu,Model.RT,
      Grid,nQuadM,Jacobi)
      Vorticity!(backend,FTB,Vort,Model.CG,UNewu,Model.RT,Grid,Grid.Type,nQuadM,Jacobi,vtkSkeletonMesh.RefineMidPoints)
      Outputs.vtkSkeleton!(vtkSkeletonMesh,FileNameOutput,Proc,ProcNumber,[hout Vort VelOut],FileNumber,cName)
    end
  end
  @show "finished"
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
    ConvertScalarVelocityCart!(backend,FTB,VelCa,Uu,Model.RT,Grid,Jacobi!)
    ConvertScalarVelocitySp!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi!)
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
        ConvertScalarVelocityCart!(backend,FTB,VelCa,Uu,Model.RT,Grid,Jacobi!)
        ConvertScalarVelocity!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi!,Grid.Form)
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
    ConvertScalarVelocityCart!(backend,FTB,VelCa,Uu,Model.RT,Grid,Jacobi!)
    ConvertScalarVelocity!(backend,FTB,VelSp,Uu,Model.RT,Grid,Jacobi!,Grid.Form)
    FileNumber += 1
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [pC VelCa VelSp], FileNumber)
  else
    FileNumber += 1
    Outputs.vtkSkeleton!(vtkSkeletonMesh, GridType, Proc, ProcNumber, [pC;], FileNumber)  
  end  
end

