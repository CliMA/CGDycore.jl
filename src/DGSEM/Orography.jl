function OrographyDG(backend,FT,Grid, CellToProc, ParallelCom, OrdPoly)

  Proc = ParallelCom.Proc
  ProcNumber = ParallelCom.ProcNumber
  Trans = Grids.TransSphereX!
  CG = FiniteElements.CGQuad{FT}(backend,OrdPoly,1,0,Grid)
  Exchange = Parallels.ExchangeStruct{FT}(backend,Grid,CG,CellToProc,Proc,ProcNumber,
      false;Discretization="CG")
  zS = Grids.Orography4(backend,FT,CG,Exchange,Grid,ParallelCom)
end


