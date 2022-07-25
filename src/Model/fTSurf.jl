function fTSurf(x,time,Global,Param)
  Model=Global.Model
  Phys=Global.Phys
  str = lowercase(Model.Problem)
  if str == "heldsuarezmoistsphere"
    (Lon,Lat,R)=cart2sphere(x[1],x[2],x[3])  
    TS = Param.DeltaTS * exp(-0.5 * Lat^2 / Param.DeltaLat^2) + Param.TSMin 
  end
  return TS
end    
