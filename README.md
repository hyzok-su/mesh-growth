![UUV Mapping](./MeshGrowth.gif) 
# Mesh Growth

This project is inspired by Nervous System’s Floraform framework, where form is generated through differential growth, local interaction rules, and adaptive remeshing.

In a similar spirit, this project treats geometry as a dynamic system that evolves through:
- local growth rates
- remeshing
- feedback between environment and geometry

Reference: https://n-e-r-v-o-u-s.com/blog/?p=6721

## Dependencies

This project uses:

- Platform: Rhino, Grasshopper

- Kangaroo Physics by Daniel Piker
  https://github.com/Dan-Piker

- Plankton (C# half-edge mesh library)
  https://github.com/meshmash/Plankton

## Acknowledgements

Special thanks to the developers of Kangaroo and Plankton for providing the simulation and mesh data structures used in this project.
