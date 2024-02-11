struct SurfaceData{FT<:AbstractFloat}
         TS::FT
         RhoVS::FT
         z0M::FT
         z0H::FT
         uStar::FT
         CT::FT
         CH::FT
       end
FT=Float32       
a=SurfaceData{FT}(0.0,0.0,0.0,0.0,0.0,0.0,0.0)       
SurfaceData1=zeros(Float32,4,200)
