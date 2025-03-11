function FluxVolumeNonLin!(F,V,DG,dXdxI,Grid,Phys)
  NX = DG.OrdPoly + 1
  nV = size(V,4)
  FLoc = zeros(3,nV,NX,NX)
  ConX = zeros(NX,NX)
  ConY = zeros(NX,NX)
  iDG = 0
  @inbounds for iF = 1 : Grid.NumFaces
    iDG = (iF - 1) * NX *NX  
    @inbounds for j = 1 : NX
      @inbounds for i = 1 : NX
        iDG += 1  
        @views FluxNonLin(FLoc[:,:,i,j],V[1,1,iDG,:],Phys)
      end
    end  
    @inbounds for iv = 1 : nV
      iD = 0  
      @inbounds for j = 1 : NX  
        @inbounds for i = 1 : NX
          iD += 1  
          ConX[i,j] = dXdxI[1,1,1,iD,1,iF] * FLoc[1,iv,i,j] +
          + dXdxI[1,2,1,iD,1,iF] * FLoc[2,iv,i,j] +
          + dXdxI[1,3,1,iD,1,iF] * FLoc[3,iv,i,j]
          ConY[i,j] = dXdxI[2,1,1,iD,1,iF] * FLoc[1,iv,i,j] +
          + dXdxI[2,2,1,iD,1,iF] * FLoc[2,iv,i,j] +
          + dXdxI[2,3,1,iD,1,iF] * FLoc[3,iv,i,j]
        end
      end
      ConX = DG.DW * ConX
      ConY = ConY * DG.DW'
      iD = 0  
      iDG = (iF - 1) * NX *NX  
      @inbounds for j = 1 : NX  
        @inbounds for i = 1 : NX
          iDG += 1  
          F[1,1,iDG,iv] = ConX[i,j] + ConY[i,j]
        end
      end
    end   
  end  
end


@inline function FluxNonLin(flux,V,Phys)
  uPos = 2
  vPos = 3
  wPos = 4
  hPos = 1
  p = PresSh(V,Phys)
  u = V[uPos] / V[hPos]
  v = V[vPos] / V[hPos]
  w = V[wPos] / V[hPos]
  flux[1,hPos] = -V[uPos]
  flux[1,uPos] = -V[uPos] * u - p
  flux[1,vPos] = -V[uPos] * v
  flux[1,wPos] = -V[uPos] * w

  flux[2,hPos] = -V[vPos]
  flux[2,uPos] = -V[vPos] * u
  flux[2,vPos] = -V[vPos] * v - p
  flux[2,wPos] = -V[vPos] * w


  flux[3,hPos] = -V[wPos]
  flux[3,uPos] = -V[wPos] * u
  flux[3,vPos] = -V[wPos] * v
  flux[3,wPos] = -V[wPos] * w - p
end

@inline function FluxNonLinAver(flux,VL,VR,m_L,m_R,Phys)
  uPos = 2
  vPos = 3
  wPos = 4
  hPos = 1
  pL = PresSh(VL,Phys)
  pR = PresSh(VR,Phys)
  hL = VL[hPos]
  hR = VR[hPos]
  uL = VL[uPos] / hL
  vL = VL[vPos] / hL
  wL = VL[wPos] / hL
  uR = VR[uPos] / hR
  vR = VR[vPos] / hR
  wR = VR[wPos] / hR

  pAv = 0.5 * (pL + pR)
  uAv = 0.5 * (uL + uR)
  vAv = 0.5 * (vL + vR)
  wAv = 0.5 * (wL + wR)
  hAv = 0.5 * (hL + hR)
  mAv1 = 0.5 * (m_L[1] + m_R[1])
  mAv2 = 0.5 * (m_L[2] + m_R[2])
  mAv3 = 0.5 * (m_L[3] + m_R[3])
  qHat = mAv1 * uAv + mAv2 * vAv + mAv3 * wAv
  flux[1] = hAv * qHat
  flux[2] = flux[1] * uAv + mAv1 * pAv
  flux[3] = flux[1] * vAv + mAv2 * pAv
  flux[4] = flux[1] * wAv + mAv3 * pAv
end    

function VolumeFluxAver(FluxAver,fTilde,gTilde,V,dXdxI,Phys)
  nV = size(fTilde,1);
  NX = size(fTilde,2);
  @inbounds for j = 1 : NX
    @inbounds for i = 1 : NX
      @inbounds for l = i + 1 : NX
        @views FluxAver(fTilde[:,l,i,j],V[i,j,:],V[l,j,:],dXdxI[1,:,i,j],
          dXdxI[1,:,l,j],Phys)  
        @. fTilde[:,i,l,j] = fTilde[:,l,i,j]
      end
      @. fTilde[:,i,i,j] = 0.0
      @inbounds for l = j + 1 : NX
        @views FluxAver(gTilde[:,l,i,j],V[i,j,:],V[i,l,:],dXdxI[2,:,i,j],
          dXdxI[2,:,i,l],Phys)  
        @. gTilde[:,j,i,l] = gTilde[:,l,i,j]
      end
      @. gTilde[:,j,i,j] = 0.0
    end
  end
end

function FluxSplitVolumeNonLin!(FluxAver,F,V,DG,dXdxI,Grid,Phys)
  NX = DG.OrdPoly + 1
  NV = 4
  fTilde = zeros(NV,NX,NX,NX)
  gTilde = zeros(NV,NX,NX,NX)
  @inbounds for iF = 1 : Grid.NumFaces
    iDG1 = (iF - 1) * NX *NX + 1  
    iDG2 = iF * NX * NX   
    @views VolumeFluxAver(FluxAver,fTilde,gTilde,reshape(V[1,1,iDG1:iDG2,:],NX,NX,NV),
      reshape(dXdxI[:,:,1,:,1,iF],3,3,NX,NX),Phys)
    @views FLoc = reshape(F[1,1,iDG1:iDG2,:],NX,NX,NV)
    @inbounds for j = 1 : NX
      @inbounds for i = 1 : NX
        @inbounds for l = 1 : NX
          @. FLoc[i,j,:] -= DG.DVT[l,i] * fTilde[:,l,i,j] +
           DG.DVT[l,j] * gTilde[:,l,i,j]
        end
      end
    end
  end
end
