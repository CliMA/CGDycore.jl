function ROS=RosenbrockMethod(Method)
switch Method
  case 'RK3_H'
    ROS.transformed=true;
    ROS.nStage=3;
    ROS.alpha=[0 0 0
      1/3 0 0 
      0 1/2 0];
    g=(3.0+sqrt(3.0))/6.0;
    ROS.Gamma=[g 0 0
               (1-12*g^2)/(-9+36*g) g 0
               -1/4+2*g 1/4-3*g g];
    ROS.b=[0,0,1];
    ROS.a=ROS.alpha/ROS.Gamma;
    ROS.c=-inv(ROS.Gamma);
    ROS.m=ROS.b/ROS.Gamma;
    ROS.d=ROS.Gamma(1,1);
  case 'RODAS_N'
    ROS.transformed=false;
    ROS.nStage=4;
    ROS.alpha=[0 0 0 0
      0 0 0 0
      1. 0 0 0
      0.75 -0.25 0.5 0];
    ROS.Gamma=[0.5 0 0 0
      1. 0.5 0 0
      -0.25 -0.25 0.5 0
      1./12 1./12 -2./3 0.5];
    ROS.b  = [5./6 -1./6 -1./6 0.5];
  case 'RODAS'
    ROS.transformed=true;
    ROS.nStage=4;
    ROS.alpha=[0 0 0 0
      0 0 0 0
      1. 0 0 0
      0.75 -0.25 0.5 0];
    ROS.Gamma=[0.5 0 0 0
      1. 0.5 0 0
      -0.25 -0.25 0.5 0
      1./12 1./12 -2./3 0.5];
    ROS.b  = [5./6 -1./6 -1./6 0.5];
    ROS.a=ROS.alpha/ROS.Gamma;
    ROS.c=-inv(ROS.Gamma);
    ROS.m=ROS.b/ROS.Gamma;
    ROS.d=ROS.Gamma(1,1);
  case 'TSROSWSANDU3_N'
    ROS.transformed=false;
    ROS.nStage=3;
    ROS.alpha=[0 0 0
      0.43586652150845899941601945119356 0 0
      0.43586652150845899941601945119356 0 0];
    ROS.Gamma=[0.43586652150845899941601945119356 0 0
      -0.19294655696029095575009695436041 0.43586652150845899941601945119356 0
      0 1.74927148125794685173529749738960 0.43586652150845899941601945119356];
    ROS.b=[-0.75457412385404315829818998646589,1.94100407061964420292840123379419...
      ,-0.18642994676560104463021124732829];
  case 'TSROSWSANDU3'
    ROS.transformed=true;
    ROS.nStage=3;
    ROS.alpha=[0 0 0
      0.43586652150845899941601945119356 0 0
      0.43586652150845899941601945119356 0 0];
    ROS.Gamma=[0.43586652150845899941601945119356 0 0
      -0.19294655696029095575009695436041 0.43586652150845899941601945119356 0
      0 1.74927148125794685173529749738960 0.43586652150845899941601945119356];
    ROS.b=[-0.75457412385404315829818998646589,1.94100407061964420292840123379419...
      ,-0.18642994676560104463021124732829];
    ROS.a=ROS.alpha/ROS.Gamma;
    ROS.c=-inv(ROS.Gamma);
    ROS.m=ROS.b/ROS.Gamma;
    ROS.d=ROS.Gamma(1,1);
  case 'TROSWLASSP3P4S2C'
    r=5;
    ROS.transformed=false;
    ROS.nStage=4;
    ROS.alpha = [0 0 0 0 
      1./2. 0 0 0 
      1./2. 1./2. 0 0 
      1./6. 1./6. 1./6. 0 ];
    ROS.Gamma= [1./2. 0 0 0
      0.0 3./4. 0 0
      -2./3. -23./9. 2./9. 0
      1./18. 65./108. -2./27 0];
    ROS.b = [1./6.,1./6.,1./6.,1./2.];
    ROS.b2= [3./16.,10./16.,3./16.,0];
    
    %ROS.m
  case 'ROSEul'
    ROS.nStage=1;
    
  case 'ROS3Pw'
    ROS.nStage=3;
    ROS.beta0=7.88675134594812865529e-01;
    ROS.alpha(1,1)=2.0e0;
    ROS.alpha(2,1)=6.33974596215561403412e-01;
    ROS.alpha(3,1)=1.63397459621556140341e+00;
    ROS.alpha(3,2)=2.94228634059947813384e-01;
    ROS.alpha(3,3)=1.07179676972449078320e+00;
    
    ROS.beta(2,1)=-2.53589838486224561365e+00;
    ROS.beta(3,1)=-1.62740473580835520728e+00;
    ROS.beta(3,2)=-2.74519052838329002952e-01;
    ROS.alpha=ROS.alpha*ROS.beta0;
    ROS.beta=ROS.beta*ROS.beta0;
    ROS.d=ROS.beta0;
  case 'ROS-AMF'
    ROS.nStage=2;
    ROS.transformed=true;
    d=1/3+1/6*sqrt(3);

    ROS.alpha=[0 0 
               2/3 0];
    ROS.Gamma=[d 0
               -4/3*d d];
    ROS.b=[1/4,3/4];
    ROS.a=ROS.alpha/ROS.Gamma;
    ROS.c=-inv(ROS.Gamma);
    ROS.m=ROS.b/ROS.Gamma;
    ROS.d=ROS.Gamma(1,1);
  case 'ROS-AMF_N'
    ROS.nStage=2;
    ROS.transformed=false;
    d=1/3+1/6*sqrt(3);

    ROS.alpha=[0 0 
               2/3 0];
    ROS.Gamma=[d 0
               -4/3*d d];
    ROS.b=[1/4,3/4];
    ROS.a=ROS.alpha/ROS.Gamma;
    ROS.c=-inv(ROS.Gamma);
    ROS.m=ROS.b/ROS.Gamma;
    ROS.d=ROS.Gamma(1,1);  
  case 'ROS2'
    ROS.nStage=2;
    % Matrix form aftrer Wolfbrandt
    %     d=.5;
    %     a(2,1)=2./3.;
    %     d(2,1)=.5;
    %     b(1)=.25;
    %     b(2)=.75
    alpha=zeros(ROS.nStage,ROS.nStage);
    gamma=zeros(ROS.nStage,ROS.nStage);
    ROS.d=0.5+sqrt(3)/6;
    d=ROS.d;
    alpha(2,1)=2/3;
    b(1,1)=0.25;
    b(1,2)=0.75;
    gamma(1,1)=d;
    gamma(2,1)=-d/b(1,2);
    gamma(2,2)=d;
    ROS.a=alpha/gamma;
    ROS.m=b/gamma;
    ROS.c=inv(diag(diag(gamma)))-inv(gamma);
  case 'ROSRK3'
    ROS.transformed=true;
    ROS.nStage=3;
    ROS.a=zeros(ROS.nStage,ROS.nStage);
    ROS.c=zeros(ROS.nStage,ROS.nStage);
    beta0=1;
    %beta0=1.0/2.0
    %beta0=(3.0+sqrt(3.0))/6.0;
    alpha(1,1)=1.0/(3.0*beta0);
    alpha(2,1)=(-1.0 + 12.0*beta0^2.0)/(18.0*beta0^2.0*(-1.0 + 4.0*beta0));
    alpha(2,2)=1.0/(2.0*beta0);
    alpha(3,1)=   (1.0 - 3.0*beta0*(7.0 + 16.0*beta0*(-2.0 + 3.0*beta0)))/...
      (36.0*beta0^3.0*(-1.0 + 4.0*beta0));
    alpha(3,2)=   (-1.0 + 12.0*beta0)/(4.0*beta0^2.0);
    alpha(3,3)=   1.0/beta0;
    beta(2,1)= (1.0 - 12.0*beta0^2.0)/(beta0^2.0*(-9.0 + 36.0*beta0));
    beta(3,1)= (-1.0 + 3.0*beta0*(7.0 + 16.0*beta0*(-2.0 + 3.0*beta0)))/...
      (36.*beta0^3.0*(-1.0 + 4.0*beta0));
    beta(3,2)= (0.250 - 3.0*beta0)/beta0^2.0;
    alpha=alpha*beta0;
    beta=beta*beta0;
    ROS.a(2:3,1:2)=alpha(1:2,1:2);
    ROS.m=alpha(3,:);
    ROS.c=beta;
    ROS.d=1;
end
end
