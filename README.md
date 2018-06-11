# Beluga-TELEMAC-MASCARET
The Beluga branch of the source code is being developed to take these meteorological conditions into account through NetCDF files

In the present century, the evolution of hydrodynamics models along with the development of atmosphere prediction models and Earth observing satellites, provide an appropriate platform for accurate real-time simulation of complex phenomena.
Storm Surges and associated extreme flooding, which cause significant loss of life and damage to the properties with large economic costs in low-lying and densely populated areas, is such an example of those complicated events.
The main objectives here in terms of long-term coastal protection and sustainable development, are dealing with possible inundation due to the storm surges and determining the potential impacts on Belgian coastal area.
To this end, among existing numerical modelling systems, TELEMAC 2D module of TELEMAC-MASCARET suite is chosen for calculating the hydrodynamics and TOMAWAC module for the wave propagation over the entire North Sea.
Since wind plays an important role in the actual threat of surges in the coastal areas (although low pressure also contributes a much smaller amount), the importance of Wind and Atmospheric Pressure Module of TELEMAC 2D is vital. 
The development of the Beluga branch of the TELEMAC source code aims to impose meteorological forcing through NetCDF files. Provided they fulfil the Climate and Forecast (CF) convention, doing so ensures that the forcing files are self-describing. These changes allow the users to take the wind and pressure into account, which vary not only in time (as the standard version) but also in space.
These modifications include adding the functionality of reading a Binary Wind and Pressure data in netCDF (Network Common Data Form) format and interpolating the data from regular grid to TELEMAC computational mesh in Mercator projection coordinate through METEO subroutine. In practical terms, the user defines the NetCDF file in the steering file, and the data are interpolated in space and time by TELEMAC to provide meteorological conditions at each grid node and at each time step.

During stormy condition as sea foam layers develop at the air-sea interface, momentum transfer decreases at the sea surface leading to a reduction in drag coefficient. So, to better parameterize the surface roughness with increasing wind speed, the constant sea surface drag coefficient has been replaced with sea-state-dependent coefficients within PROSOU subroutine using three different formulations, Flather (1976), Fairall (2003), and Makin (2005).
Keywords: TELEMAC 2D, netCDF, Wind and Atmospheric Pressure, Sea Surface Drag coefficient


https://www.dropbox.com/sh/ifp8nbtglqn21wj/AABuEPBl5MhG15ZrKHGRzVwxa?dl=0
