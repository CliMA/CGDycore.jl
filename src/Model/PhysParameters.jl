struct Cache{A3,A4,A5,A6,A7}
    lat::A3
    JC::A4
    J::A5
    X::A6
    dXdxIF::A6
    dXdxIC::A6
    dXdx::A7
    dXdxI::A7
end
function Cache(OrdPoly, OrdPolyZ, nz, NumFaces)
    # 3
    lat    = zeros(OrdPoly+1,OrdPoly+1,NumFaces)
    # 4
    JC     = zeros(OrdPoly+1,OrdPoly+1,NumFaces,nz)
    # 5
    J      = zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,NumFaces,nz)
    # 6
    X      = zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,3,NumFaces,nz)
    dXdxIF = zeros(OrdPoly+1,OrdPoly+1,NumFaces,nz+1,3,3)
    dXdxIC = zeros(OrdPoly+1,OrdPoly+1,NumFaces,nz,3,3)
    # 7
    dXdx   = zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,NumFaces,nz,3,3)
    dXdxI  = zeros(OrdPoly+1,OrdPoly+1,OrdPolyZ+1,NumFaces,nz,3,3)
    # dXdx   = nothing
    # dXdxI  = nothing

    A3 = typeof(lat)
    A4 = typeof(JC)
    A5 = typeof(J)
    A6 = typeof(dXdxIF)
    A7 = typeof(dXdx)
    return Cache{A3,A4,A5,A6,A7}(
        lat,
        JC,
        J,
        X,
        dXdxIF,
        dXdxIC,
        dXdx,
        dXdxI,
    )
end

Base.@kwdef mutable struct PhysParametersStruct{C}
    cache::C
    Lx = nothing
    Ly = nothing
    H = nothing
    OrdPoly = nothing
    hS = nothing
    Grid = nothing
    ProfRho = nothing
    ProfRhoBGrd = nothing
    ProfpBGrd = nothing
    ProfTheta = nothing
    ProfVel = nothing
    Damping = nothing
    Coriolis = nothing
    Buoyancy = nothing
    Source = nothing
    xc = nothing
    d = nothing
    a = nothing
    T0 = nothing
    DeltaT = nothing
    uMax = nothing
    vMax = nothing
    NBr = nothing
    Equation = nothing
    TopoS = nothing
    NumV = nothing
    RhoPos = nothing
    uPos = nothing
    vPos = nothing
    wPos = nothing
    ThPos = nothing
    Thermo = nothing
    HyperVisc = nothing
    HyperDCurl = nothing
    HyperDGrad = nothing
    HyperDDiv = nothing
    Upwind = nothing
    Flat = nothing
    vtkFileName = nothing
    vtk = nothing
    cNames = nothing
    TBGrd = nothing
    T = nothing
    RhoBGrdF = nothing
    pBGrd = nothing
    EndTime = nothing
    RK = nothing
    ROS = nothing
    Grav = nothing
    cS = nothing
    Cpd = nothing
    Cvd = nothing
    Rd = nothing
    p0 = nothing
    Cpv = nothing
    Gamma = nothing
    kappa = nothing
    Omega = nothing
    RadEarth = nothing
    JF = nothing
    dXdxIC = nothing
    dXdxIC11 = nothing
    dXdxIC12 = nothing
    dXdxIC21 = nothing
    dXdxIC22 = nothing
    dXdxIC31 = nothing
    dXdxIC32 = nothing
    dXdxIC33 = nothing
    dXdxIF = nothing
    dXdxIF11 = nothing
    dXdxIF12 = nothing
    dXdxIF21 = nothing
    dXdxIF22 = nothing
    dXdxIF31 = nothing
    dXdxIF32 = nothing
    dXdxIF33 = nothing
    latN = nothing
    nPanel = nothing
    ModelType = nothing
    Deep = nothing
    HeightLimit = nothing
    T0E = nothing
    T0P = nothing
    B = nothing
    K = nothing
    LapseRate = nothing
    U0 = nothing
    PertR = nothing
    Up = nothing
    PertExpR = nothing
    PertLon = nothing
    PertLat = nothing
    PertZ = nothing
    StrideDamp = nothing
    Relax = nothing
    CoriolisType = nothing
    Th0 = nothing
    ExpDist = nothing
    lat0 = nothing
    lon0 = nothing
    day = nothing
    k_a = nothing
    k_f = nothing
    k_s = nothing
    DeltaT_y = nothing
    DeltaTh_z = nothing
    T_equator = nothing
    T_min = nothing
    sigma_b = nothing
    z_D = nothing
    level = nothing
    fig = nothing
    SliceXY = nothing
    RadPrint = nothing
    RefProfile = nothing
    Pres = nothing
    KE = nothing
    xC0 = nothing
    zC0 = nothing
    rC0 = nothing
    DeltaTh = nothing
    hC = nothing
    x0C = nothing
    aC = nothing
    alphaG = nothing
    betaG = nothing
    hH = nothing
    H0G = nothing
    uM = nothing
    lat0G = nothing
    lat1G = nothing
    eN = nothing
    CacheC1 = nothing
    CacheC2 = nothing
    CacheC3 = nothing
    CacheC4 = nothing
    CacheC5 = nothing
    CacheC6 = nothing
    CacheC7 = nothing
    CacheF1 = nothing
    CacheF2 = nothing
    CacheF3 = nothing
    CacheF4 = nothing
    CacheF5 = nothing
    CacheF6 = nothing
    CacheF7 = nothing
    Cache1 = nothing
    Cache2 = nothing
    Cache3 = nothing
    Cache4 = nothing
    FCG = nothing
    fV = nothing
    k = nothing
    Vn = nothing
    RhoCG = nothing
    v1CG = nothing
    v2CG = nothing
    wCG = nothing
    wCCG = nothing
    ThCG = nothing
    J = nothing
end

function PhysParameters(cache)
# Physical parameters
Grav=9.81;
cS=360;
Grav=9.81e0;
Cpd=1004.0e0;
Cvd=717.0e0;
Rd=Cpd-Cvd;
p0=1.0e5;
Cpv=1885.0e0;
Gamma=Cpd/Cvd;
kappa=Rd/Cpd;

# Sphere
Omega=2*pi/(24*3600);
RadEarth=6.37122e+6;

return PhysParametersStruct{typeof(cache)}(;
cache,
Grav,
cS,
Cpd,
Cvd,
Rd,
p0,
Cpv,
Gamma,
kappa,
Omega,
RadEarth)
end
