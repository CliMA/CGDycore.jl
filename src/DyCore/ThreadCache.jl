function CreateCache(OP,nz,NumV,NumTr)
  D2 = Array{Float64, 2}
  D3 = Array{Float64, 3}
  D4 = Array{Float64, 4}
  TCacheC1 = D2[D2(undef, OP, OP)
    for _ in 1:Threads.nthreads()]
  TCacheC2 = D2[D2(undef, OP, OP)
    for _ in 1:Threads.nthreads()]
  TCacheC3 = D2[D2(undef, OP, OP)
    for _ in 1:Threads.nthreads()]
  TCacheC4 = D2[D2(undef, OP, OP)
    for _ in 1:Threads.nthreads()]
  TCacheC5 = D2[D2(undef, OP, OP)
    for _ in 1:Threads.nthreads()]
  TRhoCG = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  Tv1CG = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  Tv2CG = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  TwCG = D3[D3(undef, OP, OP, nz+1)
    for _ in 1:Threads.nthreads()]
  TwCCG = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  TThCG = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  TCacheCC1 = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  TCacheCC2 = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  TCacheCC3 = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  TCacheCC4 = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  TCacheCC5 = D3[D3(undef, OP, OP, nz)
    for _ in 1:Threads.nthreads()]
  TFCG = D4[D4(undef, OP, OP, nz, NumV)
    for _ in 1:Threads.nthreads()]
  TFTrCG = D4[D4(undef, OP, OP, nz, NumTr)
    for _ in 1:Threads.nthreads()]
  TTrCG = D4[D4(undef, OP, OP, nz, NumTr)
    for _ in 1:Threads.nthreads()]
  TDivTrCG = D4[D4(undef, OP, OP, nz, NumTr)
    for _ in 1:Threads.nthreads()]
  TCacheCF1 = D3[D3(undef, OP, OP, nz+1)
    for _ in 1:Threads.nthreads()]
  TCacheCF2 = D3[D3(undef, OP, OP, nz+1)
    for _ in 1:Threads.nthreads()]
  TCacheCF3 = D3[D3(undef, OP, OP, nz+1)
    for _ in 1:Threads.nthreads()]
  TCacheCF4 = D3[D3(undef, OP, OP, nz+1)
    for _ in 1:Threads.nthreads()]
  return(; TCacheC1, TCacheC2, TCacheC3, TCacheC4, TCacheC5,
           TRhoCG, Tv1CG, Tv2CG, TwCG, TwCCG, TThCG, TFCG,
           TCacheCC1, TCacheCC2, TCacheCC3, TCacheCC4, TCacheCC5, 
           TCacheCF1, TCacheCF2, TCacheCF3, TCacheCF4)
end  
        
