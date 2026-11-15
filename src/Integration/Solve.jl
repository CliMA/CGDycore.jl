# DGSEM works with Cartesian velocity components only and does not take the
# velocity form, CGSEM gets it as the last argument.
CallFcn!(Fcn,F,U,FE::FiniteElements.DGElement,Metric,Phys,Aux,Exchange,Global,VelForm) =
  Fcn(F,U,FE,Metric,Phys,Aux,Exchange,Global)
CallFcn!(Fcn,F,U,FE::FiniteElements.CGElement,Metric,Phys,Aux,Exchange,Global,VelForm) =
  Fcn(F,U,FE,Metric,Phys,Aux,Exchange,Global,VelForm)

# Cartesian velocity frame of CGSEM: after every step the tangential velocity of the state
# is projected onto the tangent plane. The tendencies are tangential, but the rounding of the
# state update (each Cartesian component separately) has the same sign step by step in a
# nearly steady flow, so that n . uh grows linearly (2e-3 m/s in two days in Float32).
ProjectState!(U,FE::FiniteElements.CGElement,Metric,Global) =
  CGSEM.ProjectTangential!(U,Metric,Global,Global.Model.VelocityFrame)
ProjectState!(U,FE,Metric,Global) = nothing

CallJac!(Jac,U,fac,FE::FiniteElements.DGElement,Metric,Phys,Aux,JCache,Global,VelForm) =
  Jac(U,fac,FE,Metric,Phys,Aux,JCache,Global)
CallJac!(Jac,U,fac,FE::FiniteElements.CGElement,Metric,Phys,Aux,JCache,Global,VelForm) =
  Jac(U,fac,FE,Metric,Phys,Aux,JCache,Global,VelForm)

function Solve!(k,v,Jac,fac,DG::FiniteElements.DGElement,Metric,Global,VelForm)
  DGSEM.Solve!(k,v,Jac,fac,DG,Metric,Global)
end
function Solve!(k,v,Jac,fac,DG::FiniteElements.CGElement,Metric,Global,VelForm)
  CGSEM.Solve!(k,v,Jac,fac,DG,Metric,Global)
end
