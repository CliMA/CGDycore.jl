# ADIOS2 output (BP5, VTX data model) of the same cells and cell data as unstructured_vtkSphere.
# One file Global.Output.vtkFileName * ".bp" for the whole run, one step per output, one block per
# MPI rank. ParaView reads it with the ADIOS2 VTX reader.

const MEMORY_SPACE_HOST = Cint(1)

mutable struct ADIOSWriter
  adios::Adios
  io::AIO
  engine::Engine
  variables::Dict{String, Variable}
  geometry::Matrix{Float64}
  connectivity::Matrix{Int64}
  celltype::UInt32
  open::Bool
end

function adios_schema(names)
  arrays = join(("        <DataArray Name=\"$name\" />\n" for name in names))
  return """
  <VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
    <UnstructuredGrid>
      <Piece NumberOfPoints="NumberOfNodes" NumberOfCells="NumberOfCells">
        <Points>
          <DataArray Name="geometry" />
        </Points>
        <Cells>
          <DataArray Name="connectivity" />
          <DataArray Name="types" />
        </Cells>
        <CellData>
  $(arrays)        <DataArray Name="TIME">step</DataArray>
        </CellData>
      </Piece>
    </UnstructuredGrid>
  </VTKFile>
  """
end

function ADIOSWriter(filename, pts, cells, names, nparts)
  nodes = length(cells[1].connectivity)
  connectivity = Matrix{Int64}(undef, nodes + 1, length(cells))
  for (k, cell) in enumerate(cells)
    connectivity[1,k] = nodes
    @views connectivity[2:end,k] .= cell.connectivity .- 1
  end
  adios = nparts > 1 ? adios_init_mpi(MPI.COMM_WORLD) : adios_init_serial()
  io = declare_io(adios, "CGDycore")
  set_engine(io, "BP5")
  define_attribute(io, "vtk.xml", adios_schema(names))
  variables = Dict{String, Variable}(
    "NumberOfNodes" => define_variable(io, "NumberOfNodes", UInt32, LocalValue()),
    "NumberOfCells" => define_variable(io, "NumberOfCells", UInt32, LocalValue()),
    "types" => define_variable(io, "types", UInt32),
    "step" => define_variable(io, "step", Float64),
    "geometry" => define_variable(io, "geometry", Float64, nothing, nothing, size(pts)),
    "connectivity" => define_variable(io, "connectivity", Int64, nothing, nothing, size(connectivity)))
  for name in names
    variables[name] = define_variable(io, name, Float64, nothing, nothing, (length(cells),))
  end
  engine = open(io, filename, mode_write)
  return ADIOSWriter(adios, io, engine, variables, Matrix{Float64}(pts), connectivity,
    UInt32(cells[1].ctype.vtk_id), true)
end

function adios_vtkSphere(U,Trans,FE,Metric,Phys,Global,VelForm,GridForm,part,nparts;
  Thermo=zeros(0,0,0,0),KV=zeros(0,0,0))

  fields = Tuple{String,Vector{Float64}}[]
  CellFieldsSphere!((name, data) -> push!(fields, (name, copy(data))),
    U,Trans,FE,Metric,Phys,Global,VelForm,GridForm;Thermo=Thermo,KV=KV)
  if Global.Output.Writer === nothing
    Global.Output.Writer = ADIOSWriter(Global.Output.vtkFileName * ".bp",
      Global.vtkCache.pts, Global.vtkCache.cells, first.(fields), nparts)
  end
  w = Global.Output.Writer
  begin_step(w.engine)
  put!(w.engine, w.variables["step"], Float64(Global.Output.vtk))
  put!(w.engine, w.variables["types"], w.celltype)
  put!(w.engine, w.variables["NumberOfNodes"], UInt32(size(w.geometry,2)))
  put!(w.engine, w.variables["NumberOfCells"], UInt32(size(w.connectivity,2)))
  put!(w.engine, w.variables["geometry"], w.geometry)
  put!(w.engine, w.variables["connectivity"], w.connectivity)
  for (name, data) in fields
    put!(w.engine, w.variables[name], data)
  end
  end_step(w.engine)
  Global.Output.vtk = Global.Output.vtk + 1
  return String[]
end

function close_output!(Global)
  w = Global.Output.Writer
  if w !== nothing && w.open
    close(w.engine)
    adios_finalize(w.adios)
    w.open = false
  end
  return nothing
end
