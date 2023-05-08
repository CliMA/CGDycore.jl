#=
  struct _p_IS {
  PETSCHEADER(struct _ISOps);
  PetscLayout  map;
  PetscInt     max,min;         /* range of possible values */
  void         *data;
  PetscInt     *total, *nonlocal;   /* local representation of ALL indices across the comm as well as the nonlocal part. */
  PetscInt     local_offset;        /* offset to the local part within the total index set */
  IS           complement;          /* IS wrapping nonlocal indices. */
  PetscBool    info_permanent[2][IS_INFO_MAX]; /* whether local / global properties are permanent */
  ISInfoBool   info[2][IS_INFO_MAX];         /* local / global properties */
};
=#
  mutable struct DMLabel 
    #PETSCHEADER(int);
    numStrata::Int      # Number of integer values */
    defaultValue::Int   # Background value when no value explicitly given */
    stratumValues::Array{Int,1}  # Value of each stratum */
    # Basic IS storage */
    validIS::Array{Bool,1}        # The IS is valid (no additions need to be merged in) */
    stratumSizes::Array{Int,1}   # Size of each stratum */
    points::Array{IS,1}         # Points for each stratum, always sorted */  IS
    # Hash tables for fast search and insertion */
    #PetscHMapI  hmap           # Hash map for fast strata search */
    #PetscHSetI *ht             # Hash set for fast insertion */
    # Index for fast search */
    pStart::Int
    pEnd::Int                     # Bounds for index lookup */
    ba::Bool                       # A bit-wise index */
    # Propagation */
    propArray::Array{Int,1}        # Array of values for propagation */
  end

  mutable struct DMLabelLink 
    label::DMLabel
    output::Bool
    next::DMLabelLink
    DMLabelLink() = (x = new(); x.next = x)
    output = false
    Section(label,
            output,
            next,) = new(
      label,
      output,
      next)
  end          
