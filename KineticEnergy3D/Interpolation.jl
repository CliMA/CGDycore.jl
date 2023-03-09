function ChangeBasis!(XOut,XIn,Fe)

  nxOut = size(XOut,1)
  nyOut = size(XOut,2)
  nzOut = size(XOut,3)
  nxIn = size(XIn,1)
  nyIn = size(XIn,2)
  nzIn = size(XIn,3)

  Buf1 = zeros(nxOut,nyIn,nzIn,3)
  Buf2 = zeros(nxOut,nyOut,nzIn,3)

  for kIn = 1 : nzIn
    for jIn = 1 : nyIn
      for iIn = 1 : nxIn
        for iOut = 1 : nxOut
          @views @.  Buf1[iOut,jIn,kIn,:] = Buf1[iOut,jIn,kIn,:] + 
            Fe.IntXE2F[iOut,iIn] * XIn[iIn,jIn,kIn,:]
        end
      end
    end
  end  
  for kIn = 1 : nzIn
    for jIn = 1 : nyIn
      for jOut = 1 : nyOut
        for iOut = 1 : nxOut
          @views @.  Buf2[iOut,jOut,kIn,:] = Buf2[iOut,jOut,kIn,:] + 
            Fe.IntYE2F[jOut,jIn] * Buf1[iOut,jIn,kIn,:]
        end
      end
    end
  end  
  for kIn = 1 : nzIn
    for kOut = 1 : nzOut
      for jOut = 1 : nyOut
        for iOut = 1 : nxOut
          @views @. XOut[iOut,jOut,kOut,:] = XOut[iOut,jOut,kOut,:] + 
            Fe.IntZE2F[kOut,kIn] * Buf2[iOut,jOut,kIn,:]
        end
      end
    end
  end  
end

function ComputeOutputC!(cOut,cIn,Fe)

  nxOut = size(cOut,1)
  nyOut = size(cOut,2)
  nzOut = size(cOut,3)
  nxIn = size(cIn,1)
  nyIn = size(cIn,2)

  Buf1 = zeros(nxOut,nyIn)

  for jIn = 1 : nyIn
    for iIn = 1 : nxIn
      for iOut = 1 : nxOut
         Buf1[iOut,jIn] = Buf1[iOut,jIn] + 
          Fe.IntXF2cE[iOut,iIn] * cIn[iIn,jIn]
      end
    end
  end
  @. cOut = 0.0 
  for jIn = 1 : nyIn
    for jOut = 1 : nyOut
      for iOut = 1 : nxOut
        cOut[iOut,jOut] = cOut[iOut,jOut] + 
          Fe.IntYF2cE[jOut,jIn] * Buf1[iOut,jIn]
      end
    end
  end
end

