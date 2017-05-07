# Lattice Quantum Ising Quench #

This example takes the form of a quench of the (lattice) quantum Ising chain, with an additional longitudinal field.
To use it first copy the files and folder into your run (build) directory e.g. BUILD/
The two files, `initial_input` and `quench_input` are input files for the drivers.
The directories `initial_defs/` and `quench_defs/` contain the defintions required to build the pre and post quench quantum Ising chain.
The quench is performed by first using DMRG on a finite size system to create the prequench state, then time evolving this with a different hamiltonian using TEBD.

The initial hamiltonian is defined as a sum of local field terms, with the field applied along the z direction:
~~~~
H = \sum_i h_z sigma^z_i
~~~~
with h_z=-10.
This value is controlled by the file `initial_defs/vertex_hamiltonian.in` because it is a single vertex term.
The coupling between vertices is set to zero in `initial_input`.

First run the fDMRG driver in your build directory:
~~~~
./fDMRG_DRV.bin -F 1 initial_input 20
~~~~
This should create a (trivial, product state) ground state for a chain of 20 sites (and open boundary conditions) with all spins aligned with the +z direction.
Because the state is so simple it converges immediately, but one finite size sweep, -F 1, is required to create all the left canonical matrices that are fed to the TEBD routine.
The file `initial_input_Energies.dat`, should convince you that the state is correct.
Note that the operators are all in terms of Pauli matrices, rather than spin operators, so the energy may differ from your expectations by a factor of 1/2.

We then perform a quench by altering the hamiltonian to
~~~~
H = \sum_i J sigma^x_i sigma^x_{i+1} - h_x sigma^x_i + h_z sigma^z_i
~~~~
with J=1.0, h_z=0.2, h_x=-1.0.
We time evolve our state by running:
~~~~
./TEBD_DRV.bin -i GroundState -B 50 -n 1000 -s 0.01 -m 10 -M Sx,5 quench_input 20
~~~~
Which will run with bond dimension 50 to time t=10, measuring every 10 steps, with step size 0.01 (2nd order Trotter by default).
It measures Sx (sigma_x really) on the 5th chain. You can add more measurements using more -M options (or you should be able to).
Including correlations of the form -M Sx,6,Sx,7 which measures the spin-spin correlator between chains 6 and 7.
The chains are numbered starting at 1 (and ending at 20 in this example). 

Note that the operator names (such as Sx etc.) correspond to operators defined in the `vertex_operators.in` files.
You can edit J in `quench_input`, while h_x and h_z, being local, can be changed by altering the matrix defined in `quench_defs/vertex_hamiltonian.in`.
