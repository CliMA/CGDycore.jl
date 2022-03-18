Base.@kwdef mutable struct PhysParametersStruct
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
    J = nothing
    lat = nothing
    X = nothing
    dXdx = nothing
    dXdxI = nothing
    JC = nothing
    JF = nothing
    dXdxIC = nothing
    dXdxIF = nothing
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
end

function PhysParameters()
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

return PhysParametersStruct(;Grav,
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
