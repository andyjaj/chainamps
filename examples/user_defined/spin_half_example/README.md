# EXCITED STATES EXAMPLE #
In this example we generate some excited states for a (finite length) spin chain, and take advantage of symmetries to use good quantum numbers.

Our hamilonian is 
`H=-sum_i (sigma^x_i sigma^x_{i+1}+sigma^y_i sigma^y_{i+1}+sigma^z_i sigma^z_{i+1}+hz sigma^z_i)`
The ground state of this ferromagnet Heisenberg model has all spins pointing in the same direction, aligned along z because of the hz term.
However as the code implictly targets the total charge 0 sector, and we will label the basis states by their z spin component, we will instead find a state with zero magentisation along the z direction.
The actual target charges can be controlled spearately.
The input file is `example_input` and contains:
~~~~
USER DEFINED
./vertex_basis.in ./vertex_hamiltonian.in ./Sx.in ./Sy.in ./Sz.in
Jx,(Sx,Sx):-1.0 Jy,(Sy,Sy):-1.0 Jz,(Sz,Sz):-1.0 hz,(Sz):0.00001
~~~~
The model is user defined, and the 2nd line tells us the names of the files defining the local (vertex) basis, the local hamiltonian (we leave this file empty) and any operators.
The basis file contents are
~~~~
Defs 0
0 1
1 -1
~~~~
The `Defs 0` line tells the driver that one quantum number is defined, of type 0, which is interpreted as `Z` i.e.  the conserved charge takes integer (positive or negative) values.
The next two lines index the states (1st integer) and give the value of the quantum number (2nd integer), sigma^z=+1,-1 in this case.
The files `Sx.in` define the Pauli matrices in a simple sparse triplet format.
The last line defines the nearest neighbour couplings, the local field and all parameter values.

We run the finite DMRG driver for 10 sites, calculating the ground state and two excited states, with 2 finite size sweeps:
`./fDMRG_DRV.bin -B 50 -F 2 -X 2 example_input 10`
The output for the ground state energy is in `example_input_Energies.dat`:
~~~~
# Index, Energy, Energy/vertices, Entropy, Truncation, Fidelity
0 -1 -0.5 0.6931471805599454 0 0
1 -3.000000000000001 -0.7500000000000002 0.8675632284814612 0 2.220446049250313e-16
2 -4.999999999999997 -0.8333333333333329 1.018230153951394 0 0.002396096296213623
3 -6.999999999999997 -0.8749999999999997 1.138073514962318 0 0.002167073412880804
4 -8.999999999999993 -0.8999999999999992 1.235866276237509 0 0.001673310200597466
5 -8.999999999999995 -0.8999999999999995 1.235866276237508 0 1
6 -8.999999999999998 -0.8999999999999998 1.235866276237506 0 1
7 -8.999999999999998 -0.8999999999999998 1.235866276237503 0 1
8 -8.999999999999998 -0.8999999999999998 1.235866276237501 0 1
~~~~
Energies calculated at the centre of the system are listed for the growth steps (iterations 0 to 4) and finite size sweeps (iterations 5 to 8).
This the energy we expect, as we have open boundary conditions.
Currently there is no specific driver for measurements on a finite chain.
Instead we perform a dummy time evolution step:
`./TEBD_DRV.bin -O 1 -s 0.00001 -M Sx,5 -M Sz,5 -M Sx,5,Sx,6 -i GroundState example_input 10`
which is a first order Trotter evolution, using the same hamiltonian (supposedly we have an eigenstate so nothing should change) and a tiny time step.
The output in `example_input_TEBD_10_GroundState_Evolution.dat` is
~~~~
# Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap), Re(Sx,5), Im(Sx,5) , Re(Sz,5), Im(Sz,5) , Re(Sx,5:Sx,6), Im(Sx,5:Sx,6) 
0 1e-05 0 1.235866276237501 1 -0.99999999595 -8.999999987850001e-05  0 0  -3.791238156747312e-15 0  0.5555555555555551 6.924660822426996e-22 
~~~~
The magnetisation at site 5 is zero, as expected by symmetry, but there are finite correlations between chains (second to last column).

The excited state energies are in `example_input_Excited_Energies.dat`:
~~~~
# Sweep Index, State Index, Energy, Entropy, Truncation, 1
0 1 -8.771632395987295 1.595649670719961 0 1
1 1 -8.790765503376415 1.610730405930463 0 1
2 1 -8.804226065180616 1.625094540056334 0 1
3 1 -8.804226065180616 1.625094540056335 0 1
4 1 -8.804226065180609 1.625094540056337 0 1
5 2 -8.590550223937081 1.876600417364971 0 1
6 2 -8.590561788655107 1.87658610439692 0 1
7 2 -8.590573354892323 1.876571800059066 0 1
8 2 -8.590573354892323 1.876571800059066 0 1
9 2 -8.590573354892319 1.876571800059069 0 1
~~~~
Here the second column gives the index of the excited state.

