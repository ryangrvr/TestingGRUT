# X5 · Refinement and spatial control

The registered 720-to-1,440-bin comparison changes the N3 terminal
count by more than 5%. The registered cutoff change also moves the N1
energy by more than 1%. The pressure-work test passes and the momentum
test requires further resolution.

Extend the frozen H55 comparison with 2,880 and 5,760 bins at the
45,000 comoving cutoff, with 4,320 bins at the 67,500 cutoff, and
with the half timestep at 2,880 bins. Test adjacent fine grids,
matched spacing for the cutoff comparison, and halved timestep.
Keep the same potential, initial data, adiabatic prescription and
output interval.

Build a separate flat periodic one-dimensional spatial control. Use
a two-component Dirac field with every negative-energy mode of the
uniform initial Hamiltonian filled, and a scalar mass phase. Evolve
all occupied spinors and the scalar field together. Match its force
and stress at the finite regulator. Use a spatially homogeneous
history and a zero-mean scalar wave, comparing the latter with a run
started at the same spatial mean. Check the spatial and temporal
refinements and a direct finite-matrix ODE. This control tests the
effect of spatial modes on production in a resolved example. It
does not substitute for a spatially resolved H55 cosmology.
