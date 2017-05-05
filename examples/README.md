#HOW TO USE THESE EXAMPLES

Copy the contents of the chosen example directory *making sure to preserve any subdirectory structure*  into the directory containing the compiled binaries.

For example if the compiled binaries are in the directory `BUILD` :
~~~~
cp -r examples/user_defined/spin_half_example_1/* BUILD/
~~~~
then run a chosen driver on the example input file.
Input files describe a vertex model and specify the inter-vertex coupling parameters.

The simplest driver is `STORE_MODEL.bin`, which takes no optional arguments, but requires one input file defining the vertex model.
It generates and stores the vertex basis and operators for information and use in post-processing data.
~~~~
cd BUILD
./STORE_MODEL MyInputFile
~~~~
This will produce a human readable .BASIS file containing the ordered basis states of the vertex described in `MyInputFile`, including the values of their various Abelian charges (if conserved by the model).
It also produces a binary format .SPARSEMATRIX file for every vertex operator that is defined by the model.
~~~~
MyInputFile.BASIS
MyInputFile-OperatorName1.SPARSEMATRIX
MyInputFile-OperatorName2.SPARSEMATRIX
MyInputFile-OperatorName3.SPARSEMATRIX
 etc.
~~~~
In general the most pertinent output files for all drivers are prefaced by the input filename, `MyInputFile`.

Input files can either call on a built in (predefined in the source code) model, or specify a user defined model, by pointing to user created files that define the vertex basis, vertex hamiltonian and any other operators.

The current set of built in vertex models includes:
* 1D Ising field theory (continuum limit of the quantum Ising chain) `CONTINUUM ISING`
* 1D free fermion field theory `CONTINUUM FREE FERMION`
* Luttinger liquid `LL`
* xxx lattice chains `OLD_XXX`
In each case the vertex itself is a single instance of the named model.

Further example documentation, for both built in and user defined models, is included in the individual subdirectories.



