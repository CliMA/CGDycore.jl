# DGSEM works with Cartesian velocity components only and does not take the
# velocity form, CGSEM gets it as the last argument.
CallFcn!(Fcn,F,U,FE::FiniteElements.DGElement,Metric,Phys,Aux,Exchange,Global,VelForm) =
  Fcn(F,U,FE,Metric,Phys,Aux,Exchange,Global)
CallFcn!(Fcn,F,U,FE::FiniteElements.CGElement,Metric,Phys,Aux,Exchange,Global,VelForm) =
  Fcn(F,U,FE,Metric,Phys,Aux,Exchange,Global,VelForm)

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
