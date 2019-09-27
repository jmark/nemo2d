```
                                                 .                 
                            ____     _          ":"                
  _ __   ___ _ __ ___   ___|___ \ __| |       ___:____     |"\/"|  
 | '_ \ / _ \ '_ ` _ \ / _ \ __) / _` |     ,'        `.    \  /   
 | | | |  __/ | | | | | (_) / __/ (_| |     |  O        \___/  |   
 |_| |_|\___|_| |_| |_|\___/_____\__,_|   ~^~^~^~^~^~^~^~^~^~^~^~^~
                                            ~^~^ ~^~ ^~^~^ ~^~^    
```

# nemo2d (beta v0.9) #

It is a lightweight, easy to understand 2D hydro code which illustrates the usage
of the 'Discontinuous Galerkin Spectral Element Method' (DGSEM).

For reference and quick insights into the datastructures a simple 'Finite-Volume'
scheme (FV) is provided as well.

The midterm goal is to grow this code into an 'easy-to-comprehend',
'copy-and-paste-ready' snippet library for integration into bigger projects or as
starting point for small scale projects like teaching/workshops/thesis-writing, etc.

The manageable code size invites for easy fiddling and experimentation.

Happy coding!

## Installation ##

```bash
git clone https://github.com/jmark/nemo2d.git
cd nemo2d # ... and modify Makefile to your needs
make
```
## Running ##

```bash
./build/nemo2d
```

## Visualization ##

In regular time intervalls 'vtk'-files (ASCII, legacy format) are written out.
It is a widely adopted file format and can be opened, for example, by 'Paraview'
or 'VisIt'.

## Parallelization ##

Parallelization is provided via OpenMP. Enable '-fopenmp' flag in the Makefile.

```bash
export OMP_NUM_THREADS=${NUMBER_OF_THREADS}
./build/nemo2d
```

## Workflow ##

It is recommended to use a little script which cleans up old files, compiles
and runs the program in one go. Following example also allows different source
and build directories.

```bash
#!/bin/bash

WORKDIR="$(pwd)" # current directory
BUILDDIR="${WORKDIR}/build"
SOURCEDIR="${HOME}/nemo2d"

NUMBER_OF_THREADS=4 # quadcore cpu

rm -vf "${WORKDIR}/"*.vtk # remove all plot files

make BUILDDIR="${BUILDDIR}" -C "${SOURCEDIR}" || exit $?

export OMP_NUM_THREADS="${NUMBER_OF_THREADS}"
./build/nemo2d
```

## Tips ##
After major changes in the code it is advisable to delete the build directory
and compile from a clean state.
```bash
rm -r build
make
```
