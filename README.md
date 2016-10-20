README
=======

*Version 1.1.0*

ChainAMPS is a set of libraries and driver routines, designed to implement MPS algorithms for large systems of coupled chains (i.e. large local or *physical* dimension).
However it also works perfectly well for canonical 1D MPS algorithms, though it is not optimised for this task.

Currently there are distinct drivers for DMRG; both infinite (iDMRG) and finite (fDMRG); and for TEBD, including the infinite case (iTEBD).
The driver naming convention for algorithm X is X_DRV.bin.

Several example models are implemented, via example files in the examples sub directory. In addition user defined models can be studied.

Further documentation is provided via doc/html/index.html

Quick start / How do I get set up?
-------------

* Dependencies: Requires blas and lapack ("framework accelerate" is used on mac), arpack (some issues with certain versions), and SuiteSparse's cxsparse. All available through e.g. macports.
* Optionally can use Intel's freely available Threaded Building Blocks (TBB) to perform threaded sparse multiplication, by specifying the "USETBB" environment variable (-DUSETBB is included by default in the makefiles).
* A rudimentary CMakeLists.txt is included that works well on mac systems with macports.
* Try 'mkdir BUILD; cd BUILD; cmake ..' to set up the build environment
* Then 'make'.
* For usage info on a particular driver run 'X_DRV.bin' where X is replaced by the name of the driver, e.g. iTEBD.
* An example is "./iDMRG_DRV.bin -B 100 -N 20 ../examples/spin_half_example_input" for 20 iDMRG iterations on a system of coupled spin half objects, with bond dimension 100.

