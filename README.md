# AstroGen 2.0

Welcome to AstroGen 2.0, a modern n-body simulation software designed to model the dynamics of galaxies and clusters. Developed in C++, AstroGen2 leverages optimized algorithms and data structures to deliver highly accurate simulations of cosmic phenomena.

## Features
- **n-body Simulation**: Simulate complex gravitational interactions between a large number of objects with high precision.
- **tree-based Force Calculation**: Efficient force calculation by dividing space into a tree structure (Octree).
- **Smoothed Particle Hydrodynamics (SPH)**: Models fluid-like behavior of interstellar gas.
- **Radiative Cooling**: Integrates radiative cooling to realistically model the thermodynamics of interstellar gas and promote star formation.
- **Star Formation**: Simulates the process of star formation from dense, cold gas.
- **Supernova feedback**: Integrates the feedback processes from Supernovas including models for both thermal and kinetic feedback
  
## Installation
See the GetttingStarted.md file for further information. 
  
## Simulation
The Simulation Engine is responsible for calculating the n-body interactions.
See /docs folder for the Users-guide, code paper and information about Dataformats etc.

Example use case: The simulation of the merger of the Andromeda and Milky Way galaxies: 

![Unbenannt](https://github.com/user-attachments/assets/b0d95ab1-4729-4994-ae02-8854382334e0)

The Plots were created with AGRender wich is compatible with the AstroGen2 Formats: https://github.com/Philip-Spaeth/AGRender

## Contributing
We welcome contributions from the community! Whether you are fixing bugs, adding new features, or improving documentation, your help is greatly appreciated. Please submit pull requests to the main branch and ensure your code follows the project's coding standards.

## License
AstroGen 2.0 is released under the GNU 3.0 License. See 'LICENSE' for more information.
