function GeoPot(Aux,DG,Metric,Exchange,Global)
  NumberThreadGPU = Global.ParallelCom.NumberThreadGPU
  @views Sources.GeoPotential!(Global.Model.GeoPotential,Aux[:,:,1:DG.NumI,:],DG.Glob,Metric.X,NumberThreadGPU)
  GPAuxPos = Global.Model.GPAuxPos
  @views Parallels.ExchangeData3DSendGPU(Aux[:,:,:,GPAuxPos:GPAuxPos],Exchange)
  @views Parallels.ExchangeData3DRecvSetGPU!(Aux[:,:,:,GPAuxPos:GPAuxPos],Exchange)
end
