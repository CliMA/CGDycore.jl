  mutable struct DM_Plex 
    ConeSection::Section   
    Cones::Array{Int,1}
    ConeOrientations::Array{Int,1}
    SupportSection::Section
    Supports::Array{Int,1}
  end    

  function DM_Plex()
    ConeSection = SectionOut()
    Cones = zeros(Int,0)
    ConeOrientations = zeros(Int,0)
    SupportSection = SectionOut()
    Supports = zeros(Int,0)
    return DM_Plex(
      ConeSection,
      Cones,
      ConeOrientations,
      SupportSection,
      Supports,
    )
  end  

#=  
struct _p_DM {
  PETSCHEADER(struct _DMOps);
  Vec                     localin[DM_MAX_WORK_VECTORS],localout[DM_MAX_WORK_VECTORS];
  Vec                     globalin[DM_MAX_WORK_VECTORS],globalout[DM_MAX_WORK_VECTORS];
  DMNamedVecLink          namedglobal;
  DMNamedVecLink          namedlocal;
  DMWorkLink              workin,workout;
  DMLabelLink             labels;            /* Linked list of labels */
  DMLabel                 depthLabel;        /* Optimized access to depth label */
  DMLabel                 celltypeLabel;     /* Optimized access to celltype label */
  void                    *ctx;    /* a user context */
  PetscErrorCode          (*ctxdestroy)(void**);
  ISColoringType          coloringtype;
  MatFDColoring           fd;
  VecType                 vectype;  /* type of vector created with DMCreateLocalVector() and DMCreateGlobalVector() */
  MatType                 mattype;  /* type of matrix created with DMCreateMatrix() */
  PetscInt                bind_below; /* Local size threshold (in entries/rows) below which Vec/Mat objects are bound to CPU */
  PetscInt                bs;
  ISLocalToGlobalMapping  ltogmap;
  PetscBool               prealloc_skip; // Flag indicating the DMCreateMatrix() should not preallocate (only set sizes and local-to-global)
  PetscBool               prealloc_only; /* Flag indicating the DMCreateMatrix() should only preallocate, not fill the matrix */
  PetscBool               structure_only; /* Flag indicating the DMCreateMatrix() create matrix structure without values */
  PetscInt                levelup,leveldown;  /* if the DM has been obtained by refining (or coarsening) this indicates how many times that process has been used to generate this DM */
  PetscBool               setupcalled;        /* Indicates that the DM has been set up, methods that modify a DM such that a fresh setup is required should reset this flag */
  PetscBool               setfromoptionscalled;
  void                    *data;
  /* Hierarchy / Submeshes */
  DM                      coarseMesh;
  DM                      fineMesh;
  DMCoarsenHookLink       coarsenhook; /* For transfering auxiliary problem data to coarser grids */
  DMRefineHookLink        refinehook;
  DMSubDomainHookLink     subdomainhook;
  DMGlobalToLocalHookLink gtolhook;
  DMLocalToGlobalHookLink ltoghook;
  /* Topology */
  PetscInt                dim;                  /* The topological dimension */
  /* Auxiliary data */
  PetscHMapAux            auxData;              /* Auxiliary DM and Vec for region denoted by the key */
  /* Flexible communication */
  PetscSF                 sfMigration;          /* SF for point distribution created during distribution */
  PetscSF                 sf;                   /* SF for parallel point overlap */
  PetscSF                 sectionSF;            /* SF for parallel dof overlap using section */
  PetscSF                 sfNatural;            /* SF mapping to the "natural" ordering */
  PetscBool               useNatural;           /* Create the natural SF */
  /* Allows a non-standard data layout */
  PetscBool               adjacency[2];         /* [use cone() or support() first, use the transitive closure] */
  PetscSection            localSection;         /* Layout for local vectors */
  PetscSection            globalSection;        /* Layout for global vectors */
  PetscLayout             map;
  /* Constraints */
  struct {
    PetscSection section;
    Mat          mat;
    Vec          bias;
  } defaultConstraint;
  /* Basis transformation */
  DM                      transformDM;          /* Layout for basis transformation */
  Vec                     transform;            /* Basis transformation matrices */
  void                   *transformCtx;         /* Basis transformation context */
  PetscErrorCode        (*transformSetUp)(DM, void *);
  PetscErrorCode        (*transformDestroy)(DM, void *);
  PetscErrorCode        (*transformGetMatrix)(DM, const PetscReal[], PetscBool, const PetscScalar **, void *);
  /* Coordinates */
  PetscInt                dimEmbed;             /* The dimension of the embedding space */
  DM                      coordinateDM;         /* Layout for coordinates */
  Vec                     coordinates;          /* Coordinate values in global vector */
  Vec                     coordinatesLocal;     /* Coordinate values in local  vector */
  PetscBool               periodic;             /* Is the DM periodic? */
  DMField                 coordinateField;      /* Coordinates as an abstract field */
  PetscReal              *L, *maxCell;          /* Size of periodic box and max cell size for determining periodicity */
  DMBoundaryType         *bdtype;               /* Indicates type of topological boundary */
  /* Null spaces -- of course I should make this have a variable number of fields */
  NullSpaceFunc           nullspaceConstructors[10];
  NullSpaceFunc           nearnullspaceConstructors[10];
  /* Fields are represented by objects */
  PetscInt                Nf;                   /* Number of fields defined on the total domain */
  RegionField            *fields;               /* Array of discretization fields with regions of validity */
  DMBoundary              boundary;             /* List of boundary conditions */
  PetscInt                Nds;                  /* Number of discrete systems defined on the total domain */
  DMSpace                *probs;                /* Array of discrete systems */
  /* Output structures */
  DM                      dmBC;                 /* The DM with boundary conditions in the global DM */
  PetscInt                outputSequenceNum;    /* The current sequence number for output */
  PetscReal               outputSequenceVal;    /* The current sequence value for output */
  PetscErrorCode        (*monitor[MAXDMMONITORS])(DM, void *);
  PetscErrorCode        (*monitordestroy[MAXDMMONITORS])(void **);
  void                   *monitorcontext[MAXDMMONITORS];
  PetscInt                numbermonitors;

  PetscObject             dmksp,dmsnes,dmts;
#ifdef PETSC_HAVE_LIBCEED
  Ceed                    ceed;                 /* LibCEED context */
  CeedElemRestriction     ceedERestrict;        /* Map from the local vector (Lvector) to the cells (Evector) */
#endif
};
=#  

  mutable struct DM
    data::DM_Plex
    label::DMLabelLink
  end

  function DM()
    data = DM_Plex()
    label = DMLabelLink()
    return DM(
      data,
      label,
    )
  end  

  function SetChart!(dm::DM, pStart::Int, pEnd::Int)
    mesh = dm.data
    SectionSetChart!(mesh.ConeSection, pStart, pEnd)
    SectionSetChart!(mesh.SupportSection, pStart, pEnd)
  end

  function SetConeSize!(dm::DM, p::Int, size::Int)
    mesh = dm.data
    SectionSetDof!(mesh.ConeSection, p, size)
  end 

  function SetCone!(dm::DM, p::Int, cone::Array{Int,1})
    mesh = dm.data
    (pStart,pEnd) = SectionGetChart(mesh.ConeSection)
    dof = SectionGetDof(mesh.ConeSection, p)
    off = SectionGetOffset(mesh.ConeSection, p)
    for c = 1 : dof
       mesh.Cones[off+c] = cone[c]  
    end  
  end

  function GetChart(dm::DM)
    mesh = dm.data
    (pStart, pEnd) = SectionGetChart(mesh.ConeSection)
  end

  function GetConeSize(dm::DM, p::Int)
    mesh = dm.data
    SectionGetDof(mesh.ConeSection, p)
  end

  function GetSupportSize(dm::DM, p::Int)
    mesh = dm.data
    SectionGetDof(mesh.SupportSection, p)
  end

  function SetConeOrientation!(dm::DM, p::Int, coneOrientation::Array{Int,1})
    mesh = dm.data
    (pStart,pEnd) = SectionGetChart(mesh.ConeSection)
    dof = SectionGetDof(mesh.ConeSection, p)
    off = SectionGetOffset(mesh.ConeSection, p)
    for c = 1 : dof
       o = coneOrientation[c] 
       cdof = SectionGetDof(mesh.ConeSection, mesh.Cones[off+c])
       mesh.ConeOrientations[off+c] = o
    end  
  end  

#=
 function CreateLabel(dm::DM, name::Char)
 end
  flg = HasLabel(dm, name)
  if (!flg) {
    PetscCall(DMLabelCreate(PETSC_COMM_SELF, name, &label));
    PetscCall(DMAddLabel(dm, label));
    PetscCall(DMLabelDestroy(&label));
  }
  PetscFunctionReturn(0);
}
  function HasLabel(dm::DM, char name, hasLabel::Bool)
  next = dm.labels
  const char    *lname;

  hasLabel = false
  while (next) {
    PetscCall(PetscObjectGetName((PetscObject) next->label, &lname));
    PetscCall(PetscStrcmp(name, lname, hasLabel));
    if (*hasLabel) break;
    next = next->next;
  }
  PetscFunctionReturn(0);
}
=#

  function DMSetUp_Plex!(dm::DM)
    mesh = dm.data
    SectionSetUp!(mesh.ConeSection)
    size = SectionGetStorageSize(mesh.ConeSection)
    mesh.Cones=zeros(Int,size)
    mesh.ConeOrientations=zeros(Int,size)
    maxSupportSize = SectionGetMaxDof!(mesh.SupportSection)
    if maxSupportSize > 0
      SectionSetUp!(mesh.SupportSection)
      size = SectionGetStorageSize(mesh.SupportSection)
      mesh.Supports = zeros(Int, size)
    end  
  end


  function Symmetrize(dm::DM)
    mesh = dm.data
    # Calculate support sizes */
    (pStart,pEnd) = SectionGetChart(mesh.ConeSection)
    for p = pStart + 1 : pEnd 
      dof = SectionGetDof(mesh.ConeSection, p)
      off = SectionGetOffset(mesh.ConeSection, p)
      for c = off + 1 : off+dof
        SectionAddDof!(mesh.SupportSection, mesh.Cones[c], 1)
      end  
    end  
    SectionSetUp!(mesh.SupportSection)
    #  Calculate supports */
    supportSize = SectionGetStorageSize(mesh.SupportSection)
    mesh.Supports = zeros(Int, supportSize)
    offsets = zeros(Int, pEnd - pStart)
    for p = pStart + 1 : pEnd
      dof = SectionGetDof(mesh.ConeSection, p)  
      off = SectionGetOffset(mesh.ConeSection, p)
      for c = off + 1 : off+dof
        q = mesh.Cones[c]
        offS = SectionGetOffset(mesh.SupportSection, q)
        mesh.Supports[offS + offsets[q] + 1] = p
        offsets[q] += 1
      end
    end 
  end  


  function Stratify(dm::DM)
    mesh = dm.data

    numRoots = 0
    numLeaves = 0
    #  Create depth label */
    (pStart,pEnd) = GetChart(dm)
    # Initialize roots and count leaves */
    sMin = 10000 #MAX_INT
    sMax = -10000 #MIN_INT;
    for p = pStart + 1 : pEnd
      coneSize = GetConeSize(dm, p)
      supportSize = GetSupportSize(dm, p)
      if coneSize == 0 && supportSize > 0
        sMin = min(p, sMin)
        sMax = max(p, sMax)
        numRoots += 1
      elseif supportSize == 0 && coneSize > 0
        numLeaves += 1
      elseif supportSize == 0 && coneSize == 0
        # Isolated points */
        sMin = min(p, sMin)
        sMax = max(p, sMax)
      end
    end
    @show sMin
    @show sMax
#   PetscCall(DMPlexCreateDepthStratum(dm, label, 0, sMin, sMax+1));

  #=
  if (numRoots + numLeaves == (pEnd - pStart)) {
    PetscInt sMin = PETSC_MAX_INT;
    PetscInt sMax = PETSC_MIN_INT;
    PetscInt coneSize, supportSize;

    for (p = pStart; p < pEnd; ++p) {
      PetscCall(DMPlexGetConeSize(dm, p, &coneSize));
      PetscCall(DMPlexGetSupportSize(dm, p, &supportSize));
      if (!supportSize && coneSize) {
        sMin = PetscMin(p, sMin);
        sMax = PetscMax(p, sMax);
      }
    }
    PetscCall(DMPlexCreateDepthStratum(dm, label, 1, sMin, sMax+1));
  } else {
    PetscInt level = 0;
    PetscInt qStart, qEnd, q;

    PetscCall(DMLabelGetStratumBounds(label, level, &qStart, &qEnd));
    while (qEnd > qStart) {
      PetscInt sMin = PETSC_MAX_INT;
      PetscInt sMax = PETSC_MIN_INT;

      for (q = qStart; q < qEnd; ++q) {
        const PetscInt *support;
        PetscInt        supportSize, s;

        PetscCall(DMPlexGetSupportSize(dm, q, &supportSize));
        PetscCall(DMPlexGetSupport(dm, q, &support));
        for (s = 0; s < supportSize; ++s) {
          sMin = PetscMin(support[s], sMin);
          sMax = PetscMax(support[s], sMax);
        }
      }
      PetscCall(DMLabelGetNumValues(label, &level));
      PetscCall(DMPlexCreateDepthStratum(dm, label, level, sMin, sMax+1));
      PetscCall(DMLabelGetStratumBounds(label, level, &qStart, &qEnd));
    }
  }
  { /* just in case there is an empty process */
    PetscInt numValues, maxValues = 0, v;

    PetscCall(DMLabelGetNumValues(label, &numValues));
    PetscCallMPI(MPI_Allreduce(&numValues,&maxValues,1,MPIU_INT,MPI_MAX,PetscObjectComm((PetscObject)dm)));
    for (v = numValues; v < maxValues; v++) {
      PetscCall(DMLabelAddStratum(label, v));
    }
  }
  PetscCall(PetscObjectStateGet((PetscObject) label, &mesh->depthState));
  PetscCall(PetscLogEventEnd(DMPLEX_Stratify,dm,0,0,0));
  PetscFunctionReturn(0);
  =#
  end
