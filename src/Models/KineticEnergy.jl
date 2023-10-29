function KineticEnergy!(KE,v1CG,v2CG,wCG,CG,Global,iF)    
#   Kinetic energy for a column 
    nz = Global.Grid.nz
#   Boundary value for vertical velocity and cell center
    BoundaryW!(view(wCG,:,:,1),v1CG,v2CG,CG,Global,iF)
#   Kinetic energy
    @views @. KE = 0.5*(v1CG*v1CG + v2CG*v2CG + 0.5 * (wCG[:,:,1:nz]*wCG[:,:,1:nz] + wCG[:,:,2:nz+1]*wCG[:,:,2:nz+1]))
end
