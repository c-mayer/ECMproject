# ECMproject - Elementary conversion mode enumeration via mplrs project

ECMproject is an open source program which fastenes up the enumeration of Elementary Conversion Modes (ECMs) on more complex metabolic models compared to [ecmtool](https://github.com/SystemsBioinformatics/ecmtool). It can be used as a stand-alone program to enumerate all ECMs and therefore get a broad overview about all metabolic possibilities an organism can perform. Via altering the input model, different environments can be simulated.

## Background

The idea of enumerating ECMs dervies from the enumeration of EFMs (Elementary Flux Modes). ECMs provide a better alternative to EFMs, if just overall conversions of metabolites are in question, due to less combinatorial explosion on more complex models. Instead of getting reaction rates with all possible elementary fluxes a metabolic model can perform, the enumeration of ECMs provides uptake (negative value) and emission (positive value) rates of metabolites. Therefore, the enumeration of all ECMs provides a framework of all possible conversions of external metabolites. ECMproject uses lexicographic reverse search via 'mplrs project' to project the whole model on only external reactions and uses mplrs to solve the vertex enumeration problem in a parallel manner. Therefore, ECMproject is a new tool with its individual solution of solving the enumeration of ECMs.

## Introduction

ECMproject is a standalone command line program for Linux systems written in Python. It got developed in the course of two Master Theses to get a faster alternative to [ecmtool](https://github.com/SystemsBioinformatics/ecmtool) on more complex models. Performance of ECMproject strongly depends on the structure of the metabolic model. All metabolic models which were used to develop ECMproject are avaiable in this GitHub repository. Up to the model 'mmsyn_sm09' with 278 reactions and 269 metabolites, [ecmtool](https://github.com/SystemsBioinformatics/ecmtool) is faster than ECMproject. For more complex models, ECMproject is the faster tool, despite both tools reach very high computational times on complex models quite fast.

## Requirements

* Python 3.11 or higher

### Required python packages

* cobra 0.26.3 or higher
* lrslib 70.a (python library for mplrs)
* numpy 1.23.5 or higher
* tcsh (just for '*--mfel option*')

Note, that just all important packages needed are shown in this list. You may require more python packages for environmental reasons or to run other scripts (e.g. for plotting) in this GitHub repository. For easy installation of all these packages in a functioning environment, we provide .yml files. The 'ECMproject_env.yml' is the recommended environment for using ECMproject. For completion, we also provide an 'ecmtool_env.yml' file to run ecmtool which got used in the Master Thesis.

#### Help with conda

```
conda env create -f ECMproject.yml
```

For a more comprehensive environment which also can run plotting scripts etc.

or:

```
conda create -n env python

conda activate env

pip install cobra

conda install -c conda-forge lrslib=70.a
```

`conda install -c conda-forge tcsh` (just for '*--mfel*' option required)

For a lean environment which can run ECMproject with all of its options.


### Tool for parallelized lexicographic reverse search

* mplrs with mplrs project

Download latest version from http://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/.

Note, that mplrs uses MPI for parallelization and therefore requires a MPI library. Used MPI was 'mpirun (Open MPI) 4.1.0'.

Note, if you installed mplrs into another directory than your PATH, you have to set the '*--mplrs*' option in the command line and give the path to the installed mplrs file (see Usage examples). Of course, you also can add the mplrs directory to your PATH instead or create a symbolic link. If you have multiple mplrs versions installed, taking the mplrs version you want with the '*--mplrs*' option is recommended.

#### Help

NOTE: mplrs uses MPI for parallelization and therefore requires the MPI library. For this example, Open MPI is used.

```
apt install libopenmpi-dev
```

Now we can install mplrs.

```
wget http://cgm.cs.mcgill.ca/~avis/C/lrslib/archive/lrslib-072.tar.gz

tar -xzf lrslib-072.tar.gz

cd lrslib-072

make && make mplrs && make install
```

## Installation

* Install mplrs.
* Install all required packages in a new Python 3.11 (or higher) environment. (recommended: install the ECMproject_env.yaml file in a new environment)
* Download the latest version of ECMproject via git clone, or as a zip file from https://github.com/c-mayer/ECMproject.

Navigate to the directory and run ECMproject via `./ECMproject.py -h`

## Description

Metabolic models as input are given via the SBML-format. Cobrapy is used to parse the input file into python objects, which get processed and written into an equation system like format, the H-representation (like it is done in [EFMlrs](https://github.com/BeeAnka/EFMlrs)). This H-representation gets projected onto external reactions via 'mplrs project' and afterwards converted to a V-representation via mplrs. By combining the stochiometric matrix with this V-representation, we get our Conversion Cone with all ECMs.

## Usage

ECMproject is a standalone command line tool for Linux systems. All avaiable options are shown in the help page (`./ECMproject.py -h`). Two options are required at least. The input SBML-File is given with the '*--file*' or '*-f*' option. Secondly, a model name has to be given with the '*--model_name*' or '*-m*' option which will name all further produced files. It is recommended to use the '*--parallel*' or '*-p*' option for parallelizing the postprocessing step and therefore fasten up the analysis. You can enter the number of used cores for mplrs and postprocessing with the '*n_processes*' or '*-n*' option. It is recommended to use higher core numbers since it can speed up the analysis quite strong on more complex models. Be aware, that you always need enough free CPUs if high numbers are chosen. With the '*--gzipped*' or '*-gz*' option, result files get output in a compressed manner. Only use this option, if disc space is a problem, because it slows down the postprocessing step significantly.

### Usage examples

Be aware, that all examples are written with relative paths constructed on the repository structure with the directory of ECMproject as working directory. Of course, you can replace these paths with your own paths.

```
./ECMproject.py -f ./metabolic_models/mmsyn_sm05.xml -m mmsyn_sm05 -n 10
```

This example takes the 'mmsyn_sm05.xml' SBML-model and saves the result file named 'mmsyn_sm05.csv' with all ECMs into the current working directory. Ten cores are used to perform **mplrs**.

<br/>

```
./ECMproject.py -f ./metabolic_models/mmsyn_sm05.xml -m mmsyn_sm05 -o ./results/ -n 10 -p
```

In this example, the result file gets saved into a result directory. Ten cores are used to perform **mplrs and postprocessing**.

<br/>

```
./ECMproject.py -f ./metabolic_models/mmsyn_sm05.xml -m mmsyn_sm05 -o ./results/ -n 10 -p -gz -t
```

In this example, the result file gets saved into the results directory in compressed format, together with a file which shows the needed time for each part of the analysis in seconds.

<br/>

```
./ECMproject.py -f ./metabolic_models/mmsyn_sm05.xml -m mmsyn_sm05 -o ./results/ -n 10 -p -mp ~/lrslib/v072/mplrs
```

In this example, mplrs got installed into the home directory. Therefore, the '*-mp*' option has to be given, if the home directory is not in the PATH.

### Advanced usage examples

```
./ECMproject.py -f ./metabolic_models/mmsyn_sm05.xml -m mmsyn_sm05 -o ./results/ -n 10 -p -dv -tmp ./new_tmp/
```

In this example, additionally to the result file, all intermediate results are saved into the 'new_tmp' directory. The '*-dv*' option asssures, that all intermediate results are not deleted at the end of the analysis.

<br/>

```
./ECMproject.py -f ./metabolic_models/mmsyn_sm05.xml -m mmsyn_sm05 -o ./results/ -n 10 -p -po ./tmp/mmsyn_sm05.projected
```

In this example, only the postprocessing step gets performed. Be aware you have a saved temporary V-representation (xxx.projected) beforehand (with the '*-dv*' option). The '*-po*' option got implemented to run the analysis fast if only postprocessing throws an error. You save time by not doing the projection and the conversion again (which need the most time in the analysis).

## Bug reports

Unusual high computation times for one run during the projection step can happen sometimes.

Note: Parsing of input SBML-file should be no problem if input file is written well.

## Authors

- Christian Mayer
- Marcus Holzer

## Acknowledgements

The source code was written by Christian Mayer and Marcus Holzer with the help of Bianca Allegra Buchner. Some code snippets of [EFMlrs](https://github.com/BeeAnka/EFMlrs) where used as inspiration. The tool was developed in the working group of Univ.-Prof. Dipl.-Ing. Dr. Jürgen Zanghellini.

## License

Free software: GPLv3 license
