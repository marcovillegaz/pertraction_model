# src/classes

## Libraries Classes
This class stores and make operation over the library. A library is an object than can store compounds data or thermodinamyc data. Compounds data correspond to all the properties like, density, visicosity, critical properties, and some of them are function handle Matlab objects. Thermodynamic data is inormation about thermodynamic models that are need depending in the model. For example, unifac models employs group contributions parameters that are unique for each system. At this level of development, the user should do the task of decompose the molecula in groups. 

For more detail, see "Determination of Activity Coefficients" chapter 8.10 from Poiling et al. 

### ThermoLibrary 
ThermoLibrary is the parent class of every library of thermodynamic models. It construct and store the complete path of the file that contains the model information. 

### CompoundsLibrary (TASK)
*Combine the methods extractPropertyAsArray() and extractPropertyAsFunc() in only one function that return and array. If propertyPath correspond to a functionHandle it returns a functionHandle vector*

*Move this extractPropertyAs function into @CompoundsLibrary*






## perstractionModel Class
The compound that are in CompoundLibrary struct are the ones which represent the extraction system, including the polymer that constitutes the membrane. This data is use to feed the the mass trasnfer model, which is defined as a class.

This class contains as properties all the variables that characterized the system: concentrations, area, volumetric flows, properties of each compound, etc. 


# src/calculations
## GCVOL60.m 
Function tha takes the critical properties of a compound and return a function handle that could estimate the density as function of temperature. 

*Modify to return a function handle and make it more general*

# References
1. Poling, B. E., Prausnitz, J. M., & O’Connell, J. P. (2001). The properties of gases and liquids (5. ed). McGraw-Hill.