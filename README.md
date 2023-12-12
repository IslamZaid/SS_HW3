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
![albedo_elements](https://github.com/IslamZaid/SS_HW3/assets/114555306/80907101-8609-480d-a132-e6fa0a075560) 
![QQQ_time](https://github.com/IslamZaid/SS_HW3/assets/114555306/b7d9d381-5950-4ccc-8e11-3574cb4abb2c)

## Algorithm Description
[Include a detailed description of the algorithm, potentially with a flowchart to illustrate the key steps.]

## Language and Platform
This code is developed mainly using MATLAB. 
## Contributors
Islam Mohamed Zaid, and Mohamed Sherif
Supervision: Dr. Elena Fantino
## License
CC-BY-4.0
Feel free to reach out for any questions or issues. Happy coding!
