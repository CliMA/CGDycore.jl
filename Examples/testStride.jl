using StrideArraysCore: @gc_preserve, StrideArray, StaticInt
using Base

function fVelT(x,time)
  uS = x[1]
  vS = x[2]
end

function FcnTracerConvA!(time)
  x = StrideArray{Float64}(undef, StaticInt(3))
  x[1] = 1.0
  x[2] = 1.0
  x[3] = 1.0
  @gc_preserve fVelT(x,time)
  @gc_preserve fVelT(x,time)
  @gc_preserve fVelT(x,time)
  @gc_preserve fVelT(x,time)
end

  time1 = 10.0
  for i = 1 : 20
    @time FcnTracerConvA!(time1)
  end


