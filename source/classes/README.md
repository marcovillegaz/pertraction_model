# src/class
## perstractionModel Class
The compound that are in CompoundLibrary struct are the ones which represent the extraction system, including the polymer that constitutes the membrane. This data is use to feed the the mass trasnfer model, which is defined as a class.

This class contains as properties all the variables that characterized the system: concentrations, area, volumetric flows, properties of each compound, etc. 


# src/calculations
## GCVOL60.m 
Function tha takes the critical properties of a compound and return a function handle that could estimate the density as function of temperature. 

*Modify to return a function handle and make it more general*