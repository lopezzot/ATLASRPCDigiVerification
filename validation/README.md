Tools to create plots using the output file of the ATLAS RPC digitization

### Usage
```sh
git clone https://github.com/lopezzot/ATLASRPCDigiVerification.git
mkdir build && cd build
cmake ../
make
./DigiValidation path-to-legacydigitizationfile path-to-run4digitizationfile
```
alternatively, it is possible to set the plot labels for Run3 and Run4 digitization with
```sh
export LegacyLabel="AString"
export Run4Label="AString"
./DigiValidation path-to-legacydigitizationfile path-to-run4digitizationfile
```
