syms ult urt
syms ul ur
syms ulb urb
syms wlt wt wrt
syms wlb wb wrb
syms vortlb vortrb
syms vortlt vortrt

syms dz1 dz2 dz4 dz4




%           -----------
%           |         |
%           |         | 
%           ult       urt   
%           |         | 
%           |         | 
% ----wlt--------wt-------wrt----
% |         |         |         |
% |         |         |         |
% |         ul        ur        |
% |         |         |         |
% |         |         |         |
% ----wlb--------wb-------wrb----
%           |         |
%           |         |
%           ulb       urb
%           |         |
%           |         |
%           -----------

% vortlb = ulb + wb - ul - wlb;
% vortrb = urb + wrb - ur - wb;
% vortrt = ur + wrt - urt - wt;
% vortlt = ul + wt - ult - wlt;

urdot = 0.25 * (wb + wrb) * vortrb + ...
        0.25 * (wt + wrt) * vortrt;
      
uldot = 0.25 * (wlb + wb) * vortlb + ...
        0.25 * (wlt + wt) * vortlt;
      
wbdot = -0.25 * (ulb + ul) * vortlb - ...
         0.25 * (urb + ur) * vortrb;
       
wtdot = -0.5 * (ul + ult) * vortlt - ...
         0.5 * (ur + urt) * vortrt;
kin = 0.5 * ul * uldot + 0.5 * ur * urdot + ...
      0.5 * wb * wbdot + 0.5 * wt * wtdot;
kin = simplify(kin);    
aa=3;    
       

%          ult-------urt
%           |         |         
%           |         |        
%           |         |
%           |         |         
%           |         |         
%          ulb-------urb
%
%          vlt-------vrt
%           |         |         
%           |         |        
%           |         |
%           |         |         
%           |         |         
%          vlb-------vrb


%
%          z4 w4
%           |
%           | dz4 
%           |
%          z3 w3
%           | 
%           | dz3 
%           |
%          z2 w2
%           | 
%           | dz2
%           |
%          z1 w1
%           | 
%           | dz1
%           |
%          z0  w0

%  dw/dt =  - w dw/dz = -1/2 dw_i^2/dz
%
%  dw_i/dt = - w_i * (w_{i+1}-w_{i-1})/(dz_{i+1}+dz_i)
%
%  Kinetic energy
%
%  sum w_i dw_i/dt 1/2(dz_{i+1}+dz_i) = sum -1/2 w_i^2 (w_{i+1}-w_{i-1})
%
%  dc/dt = -w dc/dz = -d(wc)/dz + c * dw/dz 



           
           

 
 
 
 
 
