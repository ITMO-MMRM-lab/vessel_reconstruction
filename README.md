# Vessel reconstruction (in-progress)
Creation of reconstruction of blood vessel based on clinical data.
This project was created to generate input data for the [ISR3D model](https://github.com/ISR3D/ISR3D). 
For generation, we use the reconstruction of the vessel after stenting and information about the initial measurements.

## Quick start
1. Clone a repository:
```console
$ git clone https://github.com/ITMO-MMRM-lab/vessel_reconstruction.git
```

2. Create and activate conda environment:
```console
$ conda create -n <name_env> --file requirements.txt
$ conda activate <name_env> 
```

>[!IMPORTANT]
> The VMTK library is available only for conda,  http://www.vmtk.org/download/.


6. Run: 
```console
$ python vessel_reconstruction/__main__.py
```
## [Requirements](requirements.txt)
- numpy=1.25.2
- pandas=2.1.0
- pygalmesh=0.10.7
- python=3.10.12
- vmtk=1.5.0
- vtk=9.1.0
- openpyxl=3.1.2
