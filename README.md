# Simulation of Radiation Power on an Orbiting Satellite 

This repository contains a numerical code for computing the radiative power absorbed by a satellite in orbit around the Earth. The simulation accounts for radiation fluxes from the Sun, terrestrial albedo, and terrestrial thermal radiation. The code is designed to work for a generic orbit with user-specified orbital elements, although a simplified version with zero eccentricity is also acceptable.

## Features
- **Orbital Elements:** Specify the orbital elements of the satellite, including its position, velocity, and orientation.
- **Radiation Contributions:** Compute absorbed power from solar radiation, terrestrial albedo, and terrestrial thermal radiation.
- **Zonal Model:** Utilize the zonal model for terrestrial albedo and infrared radiation as illustrated in the provided class material.
- **Flexible Resolution:** The code generates output tables and/or plots of absorbed power over time, allowing for a customizable resolution (e.g., 100 points over one orbit).
- **Earth Surface Discretization:** Achieve accurate results by discretizing the Earth's surface into sufficiently small elements and summing contributions from individual elements.
- **Nadir Pointing Satellite:** Simulate a thin, Nadir-pointing satellite with zero thickness, moving under the gravitational attraction of a spherical, homogeneous Earth (ignoring perturbations).

## Requirements
- [Specify any external libraries or dependencies required to run the code]

## Usage
[Provide instructions on how to use the code, input requirements, and any other relevant details.]

## Example Output
[Include an example of the output generated by the code, such as tables or plots.]

## Algorithm Description
[Include a detailed description of the algorithm, potentially with a flowchart to illustrate the key steps.]

## Physical Constants
[Specify the set of physical constants used in the simulation, ensuring they are obtained from a reputable source and internally consistent.]

## Language and Platform
[Specify the programming language used (e.g., Fortran) and the platform on which the code was developed.]

## Contributors
[List the contributors and their respective roles in the project.]

## License
[Specify the license under which the code is distributed.]

Feel free to reach out for any questions or issues. Happy coding!