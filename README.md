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
> The [vmtk](http://www.vmtk.org/download/) library is available only for conda.

6. Run: 
```console
$ python vessel_reconstruction/__main__.py
```
## [Requirements](requirements.txt)
- [numpy](https://numpy.org/)=1.25.2
- [pandas](https://pandas.pydata.org/)=2.1.0
- [progress](https://pypi.org/project/progress/)=1.6
- [pygalmesh](https://github.com/meshpro/pygalmesh)=0.10.7
- [python](https://www.python.org/)=3.10.12
- [scipy](https://scipy.org/)=1.11.2
- [vmtk](http://www.vmtk.org/)=1.5.0
- [vtk](https://vtk.org/)=9.1.0
- [openpyxl](https://openpyxl.readthedocs.io/en/stable/)=3.1.2
