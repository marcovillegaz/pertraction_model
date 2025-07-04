#### extractCompoundProperty()
Extracts a specified property from a list of compounsd in the 
compoundLibrary and return a cell or array depending on the data type. This function is useful because other function like
siddiqiLuca or maxwell stefan, use the properties of the 
components as arrays intead of passing all the compound library, making them more scalable. 

#### loadCompoundData(filePath):
This function is important if you modify the structure in compounds spreadsheets. 

This function take the property information stored in different spread sheet and perform the next actions: 

From OtherProps sheet, read the cell and rearrange them by grouop to store them as struct field. 

Perform a correlation of density and viscosity eperimental data

Estiamte viscosity and density from grouop contribution models information. 

## initCompoundsLibrary.m
Iimport al the properties from data/input and rearranged them in a struct in a master struct named
"compoundsStruct" which is used by other scripts and functions. 
### loadCOmpoundData.m
This function load the data for one compounds. IF there is experimetnal density or temeprature, the experimental data is fitted and save as function handle. On the other hand, if there is no experimental data, the density and viscosity are estimatated and the save as function handle. 