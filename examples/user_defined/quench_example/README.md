# Lattice Quantum Ising Quench #

In this example we perform a quench of the quantum Ising chain (on a lattice), in the presence of a longitudinal field.
First copy the files and subfolders into your run (build) directory.
The two files, `initial_input` and `quench_input` are input files for the drivers.
The directories `initial_defs/` and `quench_defs/` contain the defintions required to build the pre and post quench quantum Ising chain.

We generate the ground state of a finite system, using fMDRG, then time evolve this with a different hamiltonian using TEBD.

The initial hamiltonian is defined as a sum of local field terms, with the field applied along the z direction:
~~~~
H = \sum_i h_z sigma^z_i
~~~~
with `h_z=-10`.
For the ground state of this model all the spins point in the +z direction.

Let's examine the input files to see how we set up the initial model.
In `initial_input` we have
~~~~
USER DEFINED
./initial_defs/vertex_basis.in ./initial_defs/vertex_hamiltonian.in ./initial_defs/vertex_operators.in
hz,(Sz):-10.0
~~~~
The first line tells the driver that the model is user defined.
The first and second entries on the next line must be the file that defines the vertex (local) basis and the file containing the vertex hamiltonian respectively.
The third line defines our only "coupling" a local (single site) term.

First run the fDMRG driver in your build directory:
~~~~
./fDMRG_DRV.bin -F 1 initial_input 20
~~~~
This should create a (trivial, product state) ground state for a chain of 20 sites (and open boundary conditions) with all spins aligned with the +z direction.
Because the state is so simple it converges immediately, but one finite size sweep, -F 1, is required to create all the left canonical matrices that are fed to the TEBD routine.
The contents of the file `initial_input_Energies.dat`, should convince you that the state has the correct energy:
~~~~
# Index, Energy, Energy/vertices, Entropy, Truncation, Fidelity
0 -20 -10 0 0 0
1 -40 -10 0 0 0
2 -60 -10 0 0 0
3 -80 -10 0 0 0
4 -100 -10 0 0 0
5 -120 -10 0 0 0
6 -140 -10 0 0 0
7 -160 -10 0 0 0
8 -180 -10 0 0 0
9 -200 -10 0 0 0
10 -200 -10 0 0 1
11 -200 -10 0 0 1
~~~~
The zeroth step is just the energy of the two vertex system.
Each iteration up to the 9th increases the system size by two vertices (Ising spins in this example).
The first line gives column titles, and we see that the energy of the 20 vertex (20 lattice sites) system (reached at the 9th iteration) is -200.
The 10th and 11th iterations occur during the finite size sweep.
We can also see that we have a product state, as the entropy is zero for all iterations.
Currently the fidelity is not calculated during the finite size sweeps, and is just reported as 1.

Note that the operators are all in terms of Pauli matrices, rather than spin operators, so the energy may differ from your expectations by a factor of 1/2.

We now perform a quench by altering the hamiltonian to
~~~~
H = \sum_i J sigma^x_i sigma^x_{i+1} + h_x sigma^x_i + h_z sigma^z_i
~~~~
with` J=1.0`, `h_x=-0.2`, `h_z=1.0`.
We time evolve our state by running:
~~~~
./TEBD_DRV.bin -i GroundState -B 50 -n 1000 -s 0.01 -m 10 -M Sx,5 quench_input 20
~~~~
Which will perform 1000 time steps of size 0.01, and bond dimension 50, performing a measurement every 10 steps.
The Trotter decomposition is 2nd order by default.
The measurement requested is specified by `-M Sx,5` which measures Sx (sigma_x really) on the 5th chain.
Operator names (such as Sx etc.) correspond to operators defined in the `vertex_operators.in` files.
The chains are numbered starting at 1 (and ending at 20 in this example). 
You can add more measurements using more -M options.
So to additionally measure correlations between two vertices (say 6 and 7) we could have used:
~~~~
./TEBD_DRV.bin -i GroundState -B 50 -n 1000 -s 0.01 -m 10 -M Sx,5 -M Sx,6,Sx,7  quench_input 20
~~~~
You can edit the post quench  J in `quench_input`, while `h_x` and `h_z`, being local, can be changed by altering the matrix defined in`quench_defs/vertex_hamiltonian.in` (if we only want to alter them at t=0), or in `quench_input` if desired.

