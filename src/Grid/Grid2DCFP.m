function [Points,Faces,Cell]=Grid2DCFP(nx,ny,GridType,lenx,leny)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creation of the grid, with point coordinates, point numeration, face
% numeration, face definition, cell numeration and cell definition.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
option=GridType;
switch option
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rectangular grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'quad'
    EpsRel=0.0;
    Points((nx+1)*(ny+1)).P.x=0;
    iP=0;
    dx=1/nx;
    dy=1/ny;
    for iy=1:ny+1
      for ix=1:nx+1
        iP=iP+1;
        Points(iP).iP=iP;
        Points(iP).x=(ix-1)/nx;
        Points(iP).y=(iy-1)/ny;
        Points(iP).ix=ix;
        Points(iP).iy=iy;
        if ix==1 || ix==nx+1 || iy==1 || iy==ny+1
          Points(iP).Boundary=1;
        else
          Points(iP).Boundary=0;
        end
      end
    end
    for iP=1:size(Points,2)
      if Points(iP).Boundary==0
        Pert=rand()-0.5;
        Points(iP).x=Points(iP).x+Pert*dx*EpsRel;
        Pert=rand()-0.5;
        Points(iP).y=Points(iP).y+Pert*dy*EpsRel;
      end
    end
    
    Faces(nx*(ny+1)+(nx+1)*ny).P(2)=0;
    iF=0;
    for iy=1:ny
      for ix=1:nx+1
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='x';
        Faces(iF).P(1)=ix+(nx+1)*(iy-1);
        Faces(iF).P(2)=ix+(nx+1)*iy;
        if ix==1
          Faces(iF).Boundary=1;
        elseif ix==nx+1
          Faces(iF).Boundary=1;
        else
          Faces(iF).Boundary=0;
        end
        Faces(iF).uE=zeros(2,2);
      end
    end
    for iy=1:ny+1
      for ix=1:nx
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='y';
        Faces(iF).P(1)=ix+(nx+1)*(iy-1);
        Faces(iF).P(2)=ix+1+(nx+1)*(iy-1);
        if iy==1
          Faces(iF).Boundary=1;
        elseif iy==ny+1
          Faces(iF).Boundary=1;
        else
          Faces(iF).Boundary=0;
        end
        Faces(iF).uE=zeros(2,2);
      end
    end
    
    Cell(nx*ny).P(4)=0;
    iC=0;
    for iy=1:ny
      for ix=1:nx
        iC=iC+1;
        Cell(iC).iC=iC;
        Cell(iC).nP=4;
        Cell(iC).P(1)=ix+(nx+1)*(iy-1);
        Cell(iC).P(2)=ix+1+(nx+1)*(iy-1);
        Cell(iC).P(3)=ix+1+(nx+1)*iy;
        Cell(iC).P(4)=ix+(nx+1)*iy;
        Cell(iC).P(5)=ix+(nx+1)*(iy-1);
        Cell(iC).F(1)=(nx+1)*ny+ix+nx*(iy-1);
        Cell(iC).F(2)=ix+1+(nx+1)*(iy-1);
        Cell(iC).F(3)=(nx+1)*ny+ix+nx*iy;
        Cell(iC).F(4)=ix+(nx+1)*(iy-1);
      end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% triangular grid on base of the rectangular grid with one diagonal put in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'quadtri'
    EpsRel=0.0;
    Points((nx+1)*(ny+1)).P.x=0;
    iP=0;
    dx=1/nx;
    dy=1/ny;
    for iy=1:ny+1
      for ix=1:nx+1
        iP=iP+1;
        Points(iP).iP=iP;
        Points(iP).x=(ix-1)/nx;
        Points(iP).y=(iy-1)/ny;
        Points(iP).ix=ix;
        Points(iP).iy=iy;
        if ix==1 || ix==nx+1 || iy==1 || iy==ny+1
          Points(iP).Boundary=1;
        else
          Points(iP).Boundary=0;
        end
      end
    end
    for iP=1:size(Points,2)
      if Points(iP).Boundary==0
        Pert=rand()-0.5;
        Points(iP).x=Points(iP).x+Pert*dx*EpsRel;
        Pert=rand()-0.5;
        Points(iP).y=Points(iP).y+Pert*dy*EpsRel;
      end
    end
    
    Faces(nx*(ny+1)+(nx+1)*ny).P(2)=0;
    iF=0;
    for iy=1:ny
      for ix=1:nx+1
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='x';
        Faces(iF).P(1)=ix+(nx+1)*(iy-1);
        Faces(iF).P(2)=ix+(nx+1)*iy;
        if ix==1
          Faces(iF).Boundary=1;
        elseif ix==nx+1
          Faces(iF).Boundary=1;
        else
          Faces(iF).Boundary=0;
        end
        Faces(iF).uE=zeros(2,2);
      end
    end
    for iy=1:ny+1
      for ix=1:nx
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='y';
        Faces(iF).P(1)=ix+(nx+1)*(iy-1);
        Faces(iF).P(2)=ix+1+(nx+1)*(iy-1);
        if iy==1
          Faces(iF).Boundary=1;
        elseif iy==ny+1
          Faces(iF).Boundary=1;
        else
          Faces(iF).Boundary=0;
        end
        Faces(iF).uE=zeros(2,2);
      end
    end
    for iy=1:ny
      for ix=1:nx
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='diag';
        Faces(iF).P(1)=ix+(nx+1)*(iy-1);
        Faces(iF).P(2)=ix+1+(nx+1)*(iy);
        Faces(iF).Boundary=0;
        Faces(iF).uE=zeros(2,2);
      end
    end
    
    Cell(nx*ny).P(4)=0;
    iC=0;
    for iy=1:ny
      for ix=1:nx
        iC=iC+1;
        Cell(iC).iC=iC;
        Cell(iC).nP=3;
        Cell(iC).P(1)=ix+(nx+1)*(iy-1);
        Cell(iC).P(2)=ix+1+(nx+1)*(iy-1);
        Cell(iC).P(3)=ix+1+(nx+1)*iy;
        Cell(iC).P(4)=ix+(nx+1)*(iy-1);
        Cell(iC).F(1)=(nx+1)*ny+ix+nx*(iy-1);
        Cell(iC).F(2)=ix+1+(nx+1)*(iy-1);
        Cell(iC).F(3)=(nx+1)*ny+nx*(ny+1)+ceil(iC/2);
        iC=iC+1;
        Cell(iC).iC=iC;
        Cell(iC).nP=3;
        Cell(iC).P(1)=ix+(nx+1)*(iy-1);
        Cell(iC).P(2)=ix+1+(nx+1)*iy;
        Cell(iC).P(3)=ix+(nx+1)*iy;
        Cell(iC).P(4)=ix+(nx+1)*(iy-1);
        Cell(iC).F(1)=(nx+1)*ny+nx*(ny+1)+ceil(iC/2);
        Cell(iC).F(2)=ix+(nx+1)*(iy-1);
        Cell(iC).F(3)=(nx+1)*ny+ix+nx*iy;
      end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% triangular grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'tri'
    EpsRel=0.0;
    Points(nx*ny+nx+ny+2).P.x=0;
    iP=0;
    dx=1/nx;
    dy=1/ny;
    for iy=1:ny+1
      for ix=1:nx+1
        iP=iP+1;
        Points(iP).iP=iP;
        if mod((iy-1),2)==0
          Points(iP).x=(ix-1)/nx;
          Points(iP).y=(iy-1)/ny;
          Points(iP).ix=ix;
          Points(iP).iy=iy;
        else
          if ix==1
            Points(iP).x=(ix-1)/nx;
            Points(iP).y=(iy-1)/ny;
            Points(iP).ix=ix;
            Points(iP).iy=iy;
            iP=iP+1;
            Points(iP).iP=iP;
            Points(iP).x=(2*ix-1)/(2*nx);
            Points(iP).y=(iy-1)/ny;
            Points(iP).ix=ix;
            Points(iP).iy=iy;
          elseif ix==nx+1
            Points(iP).x=(ix-1)/nx;
            Points(iP).y=(iy-1)/ny;
            Points(iP).ix=ix;
            Points(iP).iy=iy;
          else
            Points(iP).x=(2*ix-1)/(2*nx);
            Points(iP).y=(iy-1)/ny;
            Points(iP).ix=ix;
            Points(iP).iy=iy;
          end
        end
        if ix==1 || ix==nx+1 || iy==1 || iy==ny+1
          Points(iP).Boundary=1;
        else
          Points(iP).Boundary=0;
        end
      end
    end
    for iP=1:size(Points,2)
      if Points(iP).Boundary==0
        Pert=rand()-0.5;
        Points(iP).x=Points(iP).x+Pert*dx*EpsRel;
        Pert=rand()-0.5;
        Points(iP).y=Points(iP).y+Pert*dy*EpsRel;
      end
    end
    
    Faces(2*nx*ny+2*ny+ceil((ny+1)/2)*nx+floor((ny+1)/2)*(nx+1)).P(2)=0;
    iF=0;
    for iP=1:size(Points,2)-nx
      if and(Points(iP).x==0,Points(iP).y<1)
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='x';
        if mod(Points(iP).iy-1,2)==0
          Faces(iF).P(1)=iP;
          Faces(iF).P(2)=iP+nx+1;
        else
          Faces(iF).P(1)=iP;
          Faces(iF).P(2)=iP+nx+2;
        end
        Faces(iF).Boundary=1;
        Faces(iF).uE=zeros(2,2);
      elseif and(Points(iP).x==1,Points(iP).y<1)
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='x';
        if mod(Points(iP).iy-1,2)==0
          Faces(iF).P(1)=iP;
          Faces(iF).P(2)=iP+nx+2;
        else
          Faces(iF).P(1)=iP;
          Faces(iF).P(2)=iP+nx+1;
        end
        Faces(iF).Boundary=1;
        Faces(iF).uE=zeros(2,2);
      end
    end
    iF=iF+1;
    for iP=1:size(Points,2)-1
      Faces(iF).iF=iF;
      Faces(iF).Type='y';
      if Points(iP).x~=1
        Faces(iF).P(1)=iP;
        Faces(iF).P(2)=iP+1;
        if or(Points(iP).y==0,Points(iP).y==1)
          Faces(iF).Boundary=1;
        else
          Faces(iF).Boundary=0;
        end
        Faces(iF).uE=zeros(2,2);
        iF=iF+1;
      end
    end
    iF=iF-1;
    for iP=1:size(Points,2)
      if and(and(mod(Points(iP).iy-1,2)==1,Points(iP).x>0),...
          Points(iP).x<1)
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='diag';
        Faces(iF).P(1)=iP-nx-2;
        Faces(iF).P(2)=iP;
        Faces(iF).Boundary=0;
        Faces(iF).uE=zeros(2,2);
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='diag';
        Faces(iF).P(1)=iP-nx-1;
        Faces(iF).P(2)=iP;
        Faces(iF).Boundary=0;
        Faces(iF).uE=zeros(2,2);
      end
    end
    for iP=1:size(Points,2)
      if and(and(and(mod(Points(iP).iy-1,2)==1,Points(iP).x>0),...
          Points(iP).x<1),Points(iP).y<1)
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='diag';
        Faces(iF).P(1)=iP;
        Faces(iF).P(2)=iP+nx+1;
        Faces(iF).Boundary=0;
        Faces(iF).uE=zeros(2,2);
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='diag';
        Faces(iF).P(1)=iP;
        Faces(iF).P(2)=iP+nx+2;
        Faces(iF).Boundary=0;
        Faces(iF).uE=zeros(2,2);
      end
    end
    
    Cell((2*nx+1)*ny).P(4)=0;
    iC=0;
    for iP=1:size(Points,2)
      if and(mod(Points(iP).iy-1,2)==0,Points(iP).y<1)
        if Points(iP).x==1
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=3;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP+nx+2;
          Cell(iC).P(3)=iP+nx+1;
          Cell(iC).P(4)=iP;
        else
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=3;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP+nx+2;
          Cell(iC).P(3)=iP+nx+1;
          Cell(iC).P(4)=iP;
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=3;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP+1;
          Cell(iC).P(3)=iP+nx+2;
          Cell(iC).P(4)=iP;
        end
      elseif and(mod(Points(iP).iy-1,2)==1,Points(iP).y<1)
        if Points(iP).x==0
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=3;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP+1;
          Cell(iC).P(3)=iP+nx+2;
          Cell(iC).P(4)=iP;
        elseif and(Points(iP).x<1,Points(iP).x>0)
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=3;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP+nx+2;
          Cell(iC).P(3)=iP+nx+1;
          Cell(iC).P(4)=iP;
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=3;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP+1;
          Cell(iC).P(3)=iP+nx+2;
          Cell(iC).P(4)=iP;
        end
      end
    end
    
    for iC=1:size(Cell,2)
      for iF=1:size(Faces,2)
        for i=1:Cell(iC).nP
          if and(Cell(iC).P(i)==Faces(iF).P(1),Cell(iC).P(i+1)==Faces(iF).P(2))
            Cell(iC).F(1,i)=Faces(iF).iF;
          elseif and(Cell(iC).P(i)==Faces(iF).P(2),Cell(iC).P(i+1)==Faces(iF).P(1))
            Cell(iC).F(1,i)=Faces(iF).iF;
          end
        end
      end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hexagonal Grid in Progress
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  case 'hex'
    ny=2*ny;
    EpsRel=0.0;
    if mod(nx,2)==1
      Points((ny/2+1)*(nx+1)+(ny/2-1)*nx+nx+2).P.x=0;
    else
      Points((ny/2)*(nx+1)+(ny/2-1)*nx+2*(nx+2)).P.x=0;
    end
    iP=0;
    dx=lenx/nx;
    dy=leny/ny;
    for iy=1:ny+1
      for ix=1:nx+1
        iP=iP+1;
        Points(iP).iP=iP;
        if or(iy==1,iy==ny+1)
          if or(mod((iy-1),4)==1,mod(iy-1,4)==2)
            Points(iP).x=(ix-1)/nx;
            Points(iP).y=(iy-1)/(ny);
            Points(iP).ix=ix;
            Points(iP).iy=iy;
          else
            if ix==1
              Points(iP).x=(ix-1)/nx;
              Points(iP).y=(iy-1)/ny;
              Points(iP).ix=ix;
              Points(iP).iy=iy;
              iP=iP+1;
              Points(iP).iP=iP;
              Points(iP).x=(2*ix-1)/(2*nx);
              Points(iP).y=(iy-1)/ny;
              Points(iP).ix=ix;
              Points(iP).iy=iy;
            elseif ix==nx+1
              Points(iP).x=(ix-1)/nx;
              Points(iP).y=(iy-1)/ny;
              Points(iP).ix=ix;
              Points(iP).iy=iy;
            else
              Points(iP).x=(2*ix-1)/(2*nx);
              Points(iP).y=(iy-1)/ny;
              Points(iP).ix=ix;
              Points(iP).iy=iy;
            end
          end
        else
          if or(mod((iy-1),4)==1,mod(iy-1,4)==2)
            Points(iP).x=(ix-1)/nx;
            Points(iP).y=(iy-1)/(ny);
            Points(iP).ix=ix;
            Points(iP).iy=iy;
          else
            if (2*ix-1)/(2*nx)>1
              iP=iP-1;
              continue
            else
              Points(iP).iP=iP;
              Points(iP).x=(2*ix-1)/(2*nx);
              Points(iP).y=(iy-1)/ny;
              Points(iP).ix=ix;
              Points(iP).iy=iy;
            end
          end
        end
        if ix==1 || ix==nx+1 || iy==1 || iy==ny+1
          Points(iP).Boundary=1;
        else
          Points(iP).Boundary=0;
        end
      end
    end
    for iP=1:size(Points,2)
      if Points(iP).Boundary==0
        Pert=rand()-0.5;
        Points(iP).x=Points(iP).x+Pert*dx*EpsRel;
        Pert=rand()-0.5;
        Points(iP).y=Points(iP).y+Pert*dy*EpsRel;
      end
    end
    
    Faces(nx*ny).P(2)=0;
    iF=1;
    for iP=1:size(Points,2)-1
      for iPP=iP+1:size(Points,2)
        if Points(iP).x==Points(iPP).x
          if and(and(Points(iPP).y-Points(iP).y>=2/ny,Points(iP).x>0)...
              ,Points(iP).x<1)
            break
          else
            Faces(iF).iF=iF;
            Faces(iF).Type='x';
            Faces(iF).P(1)=iP;
            Faces(iF).P(2)=iPP;
            if and(Points(iP).x==0,Points(iPP).x==0)
              Faces(iF).Boundary=1;
            elseif and(Points(iP).x==1,Points(iPP).x==1)
              Faces(iF).Boundary=1;
            else
              Faces(iF).Boundary=0;
            end
            iF=iF+1;
            break;
          end
        end
      end
    end
    for iP=1:size(Points,2)-1
      if Points(iP).y==Points(iP+1).y
        if or(Points(iP).y==1,Points(iP).y==0)
          Faces(iF).iF=iF;
          Faces(iF).Type='y';
          Faces(iF).P(1)=iP;
          Faces(iF).P(2)=iP+1;
          Faces(iF).Boundary=1;
          iF=iF+1;
        end
      end
    end
    iF=iF-1;
    for iP=1:size(Points,2)
      if mod(Points(iP).iy-1,4)==1
        if Points(iP).iy==2
          if Points(iP).x==0
            iF=iF+1;
            Faces(iF).iF=iF;
            Faces(iF).Type='diag';
            Faces(iF).P(1)=iP-nx-1;
            Faces(iF).P(2)=iP;
            Faces(iF).Boundary=0;
            Faces(iF).uE=zeros(2,2);
          elseif Points(iP).x==1
            iF=iF+1;
            Faces(iF).iF=iF;
            Faces(iF).Type='diag';
            Faces(iF).P(1)=iP-nx-2;
            Faces(iF).P(2)=iP;
            Faces(iF).Boundary=0;
            Faces(iF).uE=zeros(2,2);
          else
            iF=iF+1;
            Faces(iF).iF=iF;
            Faces(iF).Type='diag';
            Faces(iF).P(1)=iP-nx-2;
            Faces(iF).P(2)=iP;
            Faces(iF).Boundary=0;
            Faces(iF).uE=zeros(2,2);
            iF=iF+1;
            Faces(iF).iF=iF;
            Faces(iF).Type='diag';
            Faces(iF).P(1)=iP-nx-1;
            Faces(iF).P(2)=iP;
            Faces(iF).Boundary=0;
            Faces(iF).uE=zeros(2,2);
          end
        else
          if Points(iP).x==0
            iF=iF+1;
            Faces(iF).iF=iF;
            Faces(iF).Type='diag';
            Faces(iF).P(1)=iP-nx;
            Faces(iF).P(2)=iP;
            Faces(iF).Boundary=0;
            Faces(iF).uE=zeros(2,2);
          elseif Points(iP).x==1
            iF=iF+1;
            Faces(iF).iF=iF;
            Faces(iF).Type='diag';
            Faces(iF).P(1)=iP-nx-1;
            Faces(iF).P(2)=iP;
            Faces(iF).Boundary=0;
            Faces(iF).uE=zeros(2,2);
          else
            iF=iF+1;
            Faces(iF).iF=iF;
            Faces(iF).Type='diag';
            Faces(iF).P(1)=iP-nx-1;
            Faces(iF).P(2)=iP;
            Faces(iF).Boundary=0;
            Faces(iF).uE=zeros(2,2);
            iF=iF+1;
            Faces(iF).iF=iF;
            Faces(iF).Type='diag';
            Faces(iF).P(1)=iP-nx;
            Faces(iF).P(2)=iP;
            Faces(iF).Boundary=0;
            Faces(iF).uE=zeros(2,2);
          end
        end
      end
    end
    for iP=1:size(Points,2)
      if mod(Points(iP).iy-1,4)==3
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='diag';
        Faces(iF).P(1)=iP;
        Faces(iF).P(2)=iP-nx-1;
        Faces(iF).Boundary=0;
        Faces(iF).uE=zeros(2,2);
        iF=iF+1;
        Faces(iF).iF=iF;
        Faces(iF).Type='diag';
        Faces(iF).P(1)=iP;
        Faces(iF).P(2)=iP-nx;
        Faces(iF).Boundary=0;
        Faces(iF).uE=zeros(2,2);
      end
    end
    
    Cell(nx*ny/2).P(4)=0;
    iC=0;
    for iP=1:size(Points,2)-1
      if and(Points(iP).y==0,Points(iP+1).y==0)
        iC=iC+1;
        Cell(iC).iC=iC;
        Cell(iC).nP=3;
        Cell(iC).P(1)=iP;
        Cell(iC).P(2)=iP+1;
        Cell(iC).P(3)=iP+nx+2;
        Cell(iC).P(4)=iP;
      elseif and(Points(iP).y==1,Points(iP+1).y==1)
        if and(Points(iP).x==0,mod(ny,4)==0)
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=4;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP-2*nx-1;
          Cell(iC).P(3)=iP-nx;
          Cell(iC).P(4)=iP+1;
          Cell(iC).P(5)=iP;
        elseif and(Points(iP+1).x==1,mod(ny,4)==0)
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=4;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP-nx-1;
          Cell(iC).P(3)=iP-2*nx-1;
          Cell(iC).P(4)=iP+1;
          Cell(iC).P(5)=iP;
        else 
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=5;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP-nx-1;
          Cell(iC).P(3)=iP-2*nx-1;
          Cell(iC).P(4)=iP-nx;
          Cell(iC).P(5)=iP+1;
          Cell(iC).P(6)=iP;
        end
      end
      if iP<size(Points,2)-nx
        if and(and(Points(iP).x==Points(iP+nx).x,Points(iP).y~=1),...
            Points(iP).x<1-1/(2*nx))
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=6;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP-nx;
          Cell(iC).P(3)=iP+1;
          Cell(iC).P(4)=iP+nx+1;
          Cell(iC).P(5)=iP+2*nx+1;
          Cell(iC).P(6)=iP+nx;
          Cell(iC).P(7)=iP;
        end
      end
      if iP<size(Points,2)-nx-1
        if and(and(and(and(Points(iP).x==Points(iP+nx+1).x,...
            Points(iP).y~=1),Points(iP).y~=0),Points(iP).x~=1),...
            Points(iP).iy<ny-1)
          if Points(iP).iy==2
            iC=iC+1;
            Cell(iC).iC=iC;
            Cell(iC).nP=6;
            Cell(iC).P(1)=iP;
            Cell(iC).P(2)=iP-nx-1;
            Cell(iC).P(3)=iP+1;
            Cell(iC).P(4)=iP+nx+2;
            Cell(iC).P(5)=iP+2*nx+2;
            Cell(iC).P(6)=iP+nx+1;
            Cell(iC).P(7)=iP;
          else
            iC=iC+1;
            Cell(iC).iC=iC;
            Cell(iC).nP=6;
            Cell(iC).P(1)=iP;
            Cell(iC).P(2)=iP-nx;
            Cell(iC).P(3)=iP+1;
            Cell(iC).P(4)=iP+nx+2;
            Cell(iC).P(5)=iP+2*nx+2;
            Cell(iC).P(6)=iP+nx+1;
            Cell(iC).P(7)=iP;
          end
        end
      end
      if and(iP<size(Points,2)-3*nx-1,ny~=4)
        if and(and(Points(iP).x==0,Points(iP+3*nx+1).x==0),Points(iP).iy<ny-1)
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=4;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP+nx+1;
          Cell(iC).P(3)=iP+2*nx+1;
          Cell(iC).P(4)=iP+3*nx+1;
          Cell(iC).P(5)=iP;
        elseif and(Points(iP).x==1,Points(iP+3*nx+1).x==1)
          iC=iC+1;
          Cell(iC).iC=iC;
          Cell(iC).nP=4;
          Cell(iC).P(1)=iP;
          Cell(iC).P(2)=iP+3*nx+1;
          Cell(iC).P(3)=iP+2*nx;
          Cell(iC).P(4)=iP+nx;
          Cell(iC).P(5)=iP;
        end
      end
    end
    
    for iC=1:size(Cell,2)
      for iF=1:size(Faces,2)
        for i=1:Cell(iC).nP
          if and(Cell(iC).P(i)==Faces(iF).P(1),Cell(iC).P(i+1)==Faces(iF).P(2))
            Cell(iC).F(1,i)=Faces(iF).iF;
          elseif and(Cell(iC).P(i)==Faces(iF).P(2),Cell(iC).P(i+1)==Faces(iF).P(1))
            Cell(iC).F(1,i)=Faces(iF).iF;
          end
        end
      end
    end
end
end