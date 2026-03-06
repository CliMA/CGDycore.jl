sing MPI
using HDF5
using CUDA

# MPI initialisieren
MPI.Init()
comm = MPI.COMM_WORLD
info = MPI.Info()

# Prüfen, ob MPI CUDA-aware ist (optional, abhängig von MPI-Implementierung)
# @show MPI.has_cuda() 

# 1. Datei mit MPI-IO-Treiber öffnen
fapl = H5P.create(H5P_FILE_ACCESS)
H5P.set_fapl_mpio(fapl, comm, info)
fid = H5F.create("cuda_aware_julia.h5", H5F_ACC_TRUNC, H5P_DEFAULT, fapl)

# 2. Dataset erstellen
dims = (1024,)
space = H5S.create_simple(length(dims), reverse(dims))
dset = H5D.create(fid, "gpu_data", H5T_NATIVE_FLOAT, space)

# 3. Daten auf der GPU vorbereiten
d_data = CUDA.rand(Float32, 1024)

# 4. Direkt vom GPU-Pointer schreiben
# HDF5.jl erkennt das CuArray und nutzt den Device-Pointer
H5D.write(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, d_data)

# Cleanup
H5D.close(dset)
H5S.close(space)
H5F.close(fid)
MPI.Finalize()
