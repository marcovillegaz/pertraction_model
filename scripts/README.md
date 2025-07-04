# \scripts
## processExpeirments.m
This scripts opens an excel spreadsheet that containd all the experimenal data obtain from the perstarction experience and the partitions constants experiments. 

In the case of the perstraction experiments, the data is rearranged ans stored as struct.

In the case of the partitions constants, that represent the equilibrium between a membrane and a liquid phase, the data is extracted, and the partition constant are computed. Then the results are stored in a table for further use in the massTransferModel.m. 

*Modifiy line 75 to  80*

## freeVolumeFit.m
Obtain teh free volume parameters from density and viscosity temeprature data for compounds stored in the main struct. 