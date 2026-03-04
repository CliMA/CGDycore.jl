mutable struct IMEXDirkMethod{FT<:AbstractFloat}
  name::String
  nStage::Int
  AE::Array{FT, 2}
  bE::Array{FT, 1}
  cE::Array{FT, 1}
  AI::Array{FT, 2}
  bI::Array{FT, 1}
  cI::Array{FT, 1}
end

function IMEXDirkMethod{FT}() where FT<:AbstractFloat
  nStage=0
  AE=zeros(FT,0,0)
  bE=zeros(FT,0)
  cE=zeros(FT,0)
  AI=zeros(FT,0,0)
  bI=zeros(FT,0)
  cI=zeros(FT,0)
  name = ""
  return IMEXDirkMethod{FT}(
    name,
    nStage,
    AE,
    bE,
    cE,
    AI,
    bI,
    cI,
  )
end

function IMEXDirkMethod{FT}(Method) where FT<:AbstractFloat
  str = Method
   
  if str == "ExImpEul"
    nStage = 2
    AI = [0 0 ;
         0 1.0]
    bI = [0, 1.0]
    cI = [0, 1.0]
    AE = [0 0 ;
         1.0 0]
    bE = [0, 1.0]
    cE = [0, 1.0]
  elseif str == "ARS343"
    nStage = 4
    cE=zeros(FT,nStage)
    cI=zeros(FT,nStage)

    gamma = 0.4358665215084590
    a42 = 0.5529291480359398
    a43 = a42
    b1 = -3 / 2 * gamma^2 + 4 * gamma - 1 / 4
    b2 = 3 / 2 * gamma^2 - 5 * gamma + 5 / 4
    a31 =
        (1 - 9 / 2 * gamma + 3 / 2 * gamma^2) * a42 +
        (11 / 4 - 21 / 2 * gamma + 15 / 4 * gamma^2) * a43 - 7 / 2 + 13 * gamma - 9 / 2 * gamma^2
    a32 =
        (-1 + 9 / 2 * gamma - 3 / 2 * gamma^2) * a42 +
        (-11 / 4 + 21 / 2 * gamma - 15 / 4 * gamma^2) * a43 + 4 - 25 / 2 * gamma +
        9 / 2 * gamma^2
    a41 = 1 - a42 - a43
    AE = [0 0 0 0
            gamma 0 0 0
            a31 a32 0 0
            a41 a42 a43 0]
    bE = [0, b1, b2, gamma]
    AI = [0 0 0 0
         0 gamma 0 0
         0 (1 - gamma)/2 gamma 0
         0 b1 b2 gamma]
    bI = [0, b1, b2, gamma]  
elseif str == "ARS222"
  nStage = 3
  gammaD = (2 − sqrt(2))/2
  delta = 1 - 1 / (2 * gammaD)
  AE = [0 0 0
        gammaD 0 0
        delta 1 - delta 0]
  bE = [delta; 1 - delta; 0]      
  cE = zeros(FT,nStage)
  AI = [0 0 0
        0 gammaD 0
        0 1 - gammaD gammaD]
  bI = [0; 1 - gammaD; gammaD]      
  cI = zeros(FT,nStage)
elseif str == "BHR(5,5,3)"
  nStage = 5
  # Koeffizienten für BHR(5,5,3)
  gamma = 0.435866521508459

# Implizites Tableau (At)
  AI = [0 0 0 0 0;
        0.4358665215 gamma 0 0 0;
        0.1190425715 0.0354994875 gamma 0 0;
        0.1265147814 0.0463934375 0.0985534810 gamma 0;
        0.1553066127 0.1340156906 0.2443425946 0.0304685799 gamma]
  bI = AI[5, :] # Stiffly Accurate
  cI = [0.0, 0.8717330430, 0.5904085806, 0.7073282215, 1.0]

# Explizites Tableau (A)
  AE = [0 0 0 0 0;
     0.8717330430 0 0 0 0;
     0.4137255570 0.1766830235 0 0 0;
     0.4184618797 0.1691230198 0.1197433219 0 0;
     0.2789178303 0.1565576757 0.2843477169 0.2801767769 0]
  bE = [0.2789178303, 0.1565576757, 0.2843477169, 0.2801767769, 0.0]
  cE = cI
elseif str == "Boscarino"
  nStage = 5
  AE = [    0    0   0    0 0
          1/2    0   0    0 0
        11/18 1/18   0    0 0
         5/6  -5/6 1/2    0 0
         1/4   7/4 3/4 -7/4 0]
  bE = [1/4;   7/4; 3/4; -7/4; 0]       
  cE = zeros(FT,nStage)
  AI = [0   0     0   0   0
        0  1/2    0   0   0
        0  1/6  1/2   0   0
        0 -1/2  1/2 1/2   0
        0  3/2 -3/2 1/2 1/2]  
  bI = [0;  3/2; -3/2; 1/2; 1/2]      
  cI = zeros(FT,nStage)

  elseif str == "KenCar_3_2003"
    nStage = 4
    AE=zeros(FT,nStage,nStage)
    bE=zeros(FT,nStage)
    cE=zeros(FT,nStage)
    AI=zeros(FT,nStage,nStage)
    bI=zeros(FT,nStage)
    cI=zeros(FT,nStage)

    cE[2] = 1767732205903 / 2027836641118
    cE[3] = 3 / 5
    cE[4] = 1

    AE[2,1] = 1767732205903 / 2027836641118
    AE[3,1] = 5535828885825 / 10492691773637
    AE[3,2] = 788022342437 / 10882634858940
    AE[4,1] = 6485989280629 / 16251701735622 
    AE[4,2] = -4246266847089 / 9704473918619 
    AE[4,3] = -10755448449292 / 10357097424841 

    bE[1] = 1471266399579 / 7840856788654 
    bE[2] = -4482444167858 / 7529755066697 
    bE[3] = 11266239266428 / 11593286722821 
    bE[4] = 1767732205903 / 4055673282236

    cI .= cE

    AI[2,1] = 1767732205903 / 2027836641118
    AI[2,2] = 1767732205903 / 2027836641118
    AI[3,1] = 2746238789719 / 10658868560708
    AI[3,2] = -640167445237 / 6845629431997
    AI[3,3] = 1767732205903 / 2027836641118
    AI[4,1] = 1471266399579 / 7840856788654 
    AI[4,2] = -4482444167858 / 7529755066697 
    AI[4,3] = 11266239266428 / 11593286722821 
    AI[4,4] = 1767732205903 / 2027836641118

    bI .= bE
  end
  return IMEXDirkMethod{FT}(
    str,
    nStage,
    AE,
    bE,
    cE,
    AI,
    bI,
    cI,
  )
end
