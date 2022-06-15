  mutable struct Section 
    pStart::Int
    pEnd::Int
    atlasDof::Array{Int,1}
    atlasOff::Array{Int,1}
    maxDof::Int
    pointMajor::Bool
    bc::Section
    perm::Bool
    setup::Bool
#   bcIndices::Array{Int,1}
#   setup::Bool
    numFields::Int
#   fieldNames::Char
#   numFieldComponents::Array{Int,1}
#   field::Array{Section,1}
    Section() = (x = new(); x.bc = x)
    Section(pStart,
            pEnd, 
            atlasDof, 
            atlasOff, 
            maxDof, 
            pointMajor, 
            bc, 
            perm, 
            setup,
            numFields) = new(
      pStart,
      pEnd,
      atlasDof,
      atlasOff,
      maxDof,
      pointMajor,
      bc,
      perm,
      setup,
      numFields,
    )

  end

  function SectionOut()
    pStart = 0
    pEnd = 0
    atlasDof = zeros(Int,0)
    atlasOff = zeros(Int,0)
    maxDof = -1
    pointMajor = true
    bc = Section()
    perm = false
    setup = false
    numFields = 0
    return Section(
      pStart,
      pEnd,
      atlasDof,
      atlasOff,
      maxDof,
      pointMajor,
      bc,
      perm,
      setup,
      numFields,
    )
  end  

  function SectionSetChart!(s::Section, pStart::Int, pEnd::Int)

    if pStart == s.pStart && pEnd == s.pEnd
    else  
      sizehint!(s.atlasDof, 0)
      sizehint!(s.atlasOff, 0)
      s.pStart = pStart
      s.pEnd   = pEnd
      s.atlasDof = zeros(Int,pEnd - pStart)
      s.atlasOff = zeros(Int,pEnd - pStart)
    end
  end  

  function SectionSetDof!(s::Section, point::Int, numDof::Int)
    s.atlasDof[point-s.pStart] = numDof  
    s.maxDof = -1
  end

  function SectionGetMaxDof!(s::Section)
    if s.maxDof == -1
     s.maxDof = 0;
     for p = 1 :  s.pEnd - s.pStart
       s.maxDof = max(s.maxDof, s.atlasDof[p])
     end
    end
    maxDof = s.maxDof
  end

  function SectionGetChart(s::Section)
    return (s.pStart,s.pEnd)
  end  

  function SectionGetDof(s::Section, point::Int)
    return s.atlasDof[point - s.pStart]
  end  

  function SectionGetOffset(s::Section, point::Int)
    return s.atlasOff[point - s.pStart]
  end

  function SectionSetUp!(s::Section)
    offset = 0
    if s.setup
      return  
    end
    s.setup = true
    # Set offsets and field offsets for all points */
    # Assume that all fields have the same chart */
    if s.perm
      pind = ISGetIndices(s.perm)
    end  
    if s.pointMajor
      for p = 1 : s.pEnd - s.pStart
        if s.perm
          q = pind[p]  
        else
          q = p
        end  
        # Set point offset */
        s.atlasOff[q] = offset
        offset += s.atlasDof[q]
        #  Set field offset */
        foff = s.atlasOff[q]
        for f = 1 : s.numFields
          sf = s.field[f]
          sf.atlasOff[q] = foff
          foff += sf.atlasDof[q]
        end
      end
    end
  end  

  function SectionGetStorageSize(s::Section)
    size = 0
    for p = 1 : s.pEnd - s.pStart
      size += s.atlasDof[p] > 0 ? s.atlasDof[p] : 0
    end  
    return size
  end  

  function SectionAddDof!(s::Section, point::Int , numDof::Int)
    s.atlasDof[point - s.pStart] += numDof
  end  
