using GLMakie
using GeometryBasics

# 1. Define a Sphere primitive
using GLMakie
using GeometryBasics

# 1. Define a Sphere primitive
# Sphere(center, radius)
s = Sphere(Point3f(0), 1.0f0)

# 2. Tesselate it into a triangular mesh
# '64' here defines the level of detail
m = geometry_mesh(s, 64)

# 3. Plotting
fig = Figure(resolution = (800, 800))
ax = LScene(fig[1, 1], show_axis = false)

# Plot the surface with a colormap
mesh!(ax, m, color = [p[3] for p in coordinates(m)], colormap = :magma)

# Overlay the "unstructured" wireframe to see the grid lines
wireframe!(ax, m, color = (:white, 0.2), linewidth = 1)

display(fig)
