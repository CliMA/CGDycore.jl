h5open("$file_name.hdf", "w") do file
        VTKHDF = create_group(file, "VTKHDF")
        HDF5.attributes(VTKHDF)["Version"] = [2,0]
        
        type = "UnstructuredGrid"
        dspace = HDF5.dataspace(type)
        dtype = HDF5.datatype(type)
        HDF5.h5t_set_cset(dtype, HDF5.H5T_CSET_ASCII)
        attr = create_attribute(VTKHDF, "Type", dtype, dspace)
        write_attribute(attr, dtype, type)

        VTKHDF["NumberOfConnectivityIds"] = [length(connectivity)]
        VTKHDF["NumberOfPoints"] = [size(points, 2)]
        VTKHDF["NumberOfCells"] = [size(connectivity, 2)]
        VTKHDF["Points"] = points
        VTKHDF["Types"] = types
        VTKHDF["Connectivity"] = connectivity[:] .- 1
        VTKHDF["Offsets"] = offsets # starting with 0 (Ncells+1)

        CellData = create_group(VTKHDF, "CellData")
        
        for i in eachindex(solution_coeffs)
                interpolation = Interpolator * solution_coeffs[i]
                name = @sprintf "%s_%03d" array_name i
                CellData[name] = interpolation  
        end
    end
