function RK=MethodRKIMEX(name)
switch name
  case 'ARK2'
    RK.nStage=3;
    RK.aE=zeros(RK.nStage,RK.nStage);
    RK.aI=zeros(RK.nStage,RK.nStage);
    RK.cE=zeros(RK.nStage,1);
    RK.cI=zeros(RK.nStage,1);
    RK.bE=zeros(1,RK.nStage);
    RK.bI=zeros(1,RK.nStage);
    RK.gamma=1-1/sqrt(2);
    alpha=(3+2*sqrt(2))/6;
    delta=1/(2*sqrt(2));
    gamma=1-1/sqrt(2);
    
    RK.gamma=gamma;
    RK.aE(2,1)=2-sqrt(2);
    RK.aE(3,1)=1-alpha;
    RK.aE(3,2)=alpha;
    RK.cE(2)=2-sqrt(2);
    RK.cE(3)=1;
    RK.bE(1)=delta;
    RK.bE(2)=delta;
    RK.bE(3)=gamma;
    
    RK.aI(2,1)=gamma;
    RK.aI(2,2)=gamma;
    RK.aI(3,1)=delta;
    RK.aI(3,2)=delta;
    RK.aI(3,3)=gamma;
    RK.cI=RK.cE;
    RK.bI=RK.bE;
    
    case 'HSDIRK2(2,2,2)'
    RK.nStage=2;
    RK.aE=zeros(RK.nStage,RK.nStage);
    RK.aI=zeros(RK.nStage,RK.nStage);
    RK.cE=zeros(RK.nStage,1);
    RK.cI=zeros(RK.nStage,1);
    RK.bE=zeros(1,RK.nStage);
    RK.bI=zeros(1,RK.nStage);
    RK.gamma=0.5;
    
    RK.aE(2,1)=1;
    RK.cE(2)=1;
    RK.bE(1)=.5;
    RK.bE(2)=.5;

    
    RK.aI(1,1)=.5;
    RK.aI(2,2)=.5;
    RK.cI(1)=.5;
    RK.cI(2)=.5;
    RK.bI=RK.bE;
  case default
    a='Keine Methode'
end
end
