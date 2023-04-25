# 1D LATTICE QUANTUM ISING MODEL QUENCH #

In this example we perform a quench of a single quantum Ising (lattice) chain, in the presence of a longitudinal field.
First copy the files and subfolders into your run (build) directory.
The two files, `initial_input` and `quench_input` are input files for the drivers.
The directories `initial_defs/` and `quench_defs/` contain the defintions required to build the pre and post quench quantum Ising chain.

We generate the ground state of a finite system, using fMDRG, then time evolve this with a different hamiltonian using TEBD.

## GENERATING THE INITIAL STATE ##

### Explanation of input files ###

The initial hamiltonian is defined as a sum of local field terms, with the field applied along the z direction:
~~~~
H = \sum_i h_z sigma^z_i
~~~~
with `h_z=-10`.
For the ground state of this model all the spins point in the +z direction.
As this is such a simple state, we could tell the TEBD_DRV driver directly to start with it.
However here we will demonstrate using fDMRG to find the state, then feed the answer to TEBD_DRV.

Let's examine the input files to see how we set up the initial model.

In `initial_input` we have
~~~~
USER DEFINED
./initial_defs/vertex_basis.in ./initial_defs/vertex_hamiltonian.in ./initial_defs/vertex_operators.in
hz,(Sz):-10.0
~~~~
The first line tells the driver that the model is user defined, and therefore that it will need to load information from several user specified files.
On the next line, the first and second entries are, respectively, the file that defines the vertex (local) basis and the file containing the vertex hamiltonian. These two files must appear in this order.
Any files after these two are interpreted as files defining operators.
The third line defines our only "coupling": a local (single site) term.
The operator involved in this term, `Sz`, is named and defined in `vertex_operators.in`.

The contents of `./initial_defs/vertex_basis.in` are
~~~~
Defs 1
0	 0
1	 0
~~~~
The first line defines the conserved Abelian (`Z_n`) charges in the model.
Because in this simple example we won't use conservation laws, we define only one charge and set n=1 for safety.
Note to achieve a `Z` charge (i.e. the integers) we would set n=0 or to a negative number.
The next two lines define the two local basis states.
The first integer on each line is just a dummy index, useful when organising the states.
The second integer is the value of the single `Z_n` charge, which is zero in both cases.
If we wanted to switch to conserved magnetisation (total Sz) we could use:
~~~~
Defs 0
0	 1
1	 -1
~~~~
This would allow for a charge with `Z` symmetry and associate the 0 index state with spin up, and the 1 index state with spin down.

As our system is very simple and we have included the necessary local term in initial_input, the file
`./initial_defs/vertex_hamiltonian.in` could be left empty.
However for clarity we define all the elements of this two by two matrix to be zero=(0.0,0.0) in complex form.
~~~~
0 0 (0.0,0.0)
0 1 (0.0,0.0)
1 0 (0.0,0.0)
1 1 (0.0,0.0)
~~~~

We choose to define all of our operators in a single file, `./initial_defs/vertex_operators.in` because they are simple:
~~~~
       Sz              Sx              Sy
        NM              NM              NM
0 0     (1.0)           (0.0)           (0.0)
1 1     (-1.0)          (0.0)           (0.0)
0 1     (0.0)           (1.0)           (0.0,-1.0)
1 0     (0.0)           (1.0)           (0.0,1.0)

~~~~
Here the first line gives the names of the operators in order.
The second line is a dummy line that can be used for comments if you wish.
The lines afterwards define the row, column and then values for the elements.
Omitted elements are assumed equal to zero.

### Running fDMRG_DRV ###

We now run the fDRMG driver.
In your build directory, execute
~~~~
./fDMRG_DRV.bin -F 1 initial_input 20
~~~~
This should create a (trivial, product state) ground state for a chain of 20 sites (and open boundary conditions) with all spins aligned with the +z direction.
Because the state is so simple it converges immediately, but one finite size sweep, -F 1, is required to create all the left canonical matrices that are fed to the TEBD routine.
In general, for a more sophisticated initial hamiltonian one should set the number of sweeps higher and check for convergence of (at the very least) the energy.
Also we have not set a maximum bond dimension here using -B (the default is no truncation), because we know that the state is simple.
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

The finite DMRG driver will produce output of the form `GroundState_X_N.TYPE` where X=Left or Right, for left or right canonical MPS matrices or sweep direction, N is the vertex number, measured from the left or right end of the system, depending on the canonisation, and TYPE=MPS_matrix or BLOCK.
The former are the stored MPS matrices that define the state.
The BLOCK files are intermediate storage files used by the fDRMG algorithm, and can be deleted if you wish.

If we had requested excited states with `-X 2` for example, we would also have files prefixed by `Excited_1` and `Excited_2`, and an extra file with the excited energies (and entropy etc.) in it.

## THE QUENCH ##

We now perform a quench by altering the hamiltonian to
~~~~
H = \sum_i J_x sigma^x_i sigma^x_{i+1} + h_x sigma^x_i + h_z sigma^z_i
~~~~
with` J=1.0`, `h_x=-0.2`, `h_z=1.0`.

The `quench_input` file contains
~~~~
USER DEFINED
2 ./quench_defs/vertex_hamiltonian.in ./quench_defs/vertex_operators.in
Jx,(Sx,Sx):1.0 hx,(Sx):-0.2 hz,(Sz):1.0 
~~~~
Again, the model is user defined, and we have to let the driver know where the definition files are.
However instead of specifying a basis file we have simply given a positive integer.
This tells the driver to simply create a basis of size 2 with no conserved quantum numbers.
The files `./quench_defs/vertex_hamiltonian.in` and `./quench_defs/vertex_operators.in` are the same as before in this example.
The last line defines the post quench couplings and parameters.

We time evolve our state by running:
~~~~
./TEBD_DRV.bin -i GroundState -B 50 -n 100 -s 0.01 -m 10 -M Sx,4 -M Sz,5 quench_input 20
~~~~
Which will perform 100 time steps of size 0.01, and bond dimension 50, performing a measurement every 10 steps.
The Trotter decomposition is 2nd order by default.
The measurement requested is specified by `-M Sx,5` which measures Sx (sigma_x really) on the 5th chain.
Operator names (such as Sx etc.) correspond to operators defined in the `vertex_operators.in` files.
The chains are numbered starting at 1 (and ending at 20 in this example). 
You can add more measurements using more -M options.
So to additionally measure correlations between two vertices (say 6 and 7) we could have used  `-M Sx,6,Sx,7`.

For this example the contents of the out file `quench_input_TEBD_20_GroundState_Evolution.dat` should look something like
~~~~
# Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap), Re(Sx,4), Im(Sx,4) , Re(Sz,5), Im(Sz,5) 
0 0.09999999999999999 2.056800830138675e-31 0.05524791340753726 0.906877637508299 -0.3033123572169272 -0.8546511939804127  -0.003933755180302386 0  0.9602657645177681 0 
1 0.2 2.887791217784476e-31 0.1597882325141478 0.6852772104859244 -0.5003227277432383 0.468275584794825  -0.01497332214868944 0  0.8531897109163494 0 
2 0.3000000000000001 3.05528761691764e-31 0.2735126239201029 0.4499018213158437 0.4332951088065349 0.1211073800709372  -0.03108668318749463 -1.852960857822178e-25  0.710520580271787 0 
3 0.4000000000000002 3.706541937330327e-31 0.3759220914665607 0.2763905354583518 -0.1207142819890908 -0.2486358586664698  -0.04964382924503302 1.214995697926664e-17  0.5716984968244793 0 
4 0.5000000000000002 3.84981617479902e-31 0.4587468171664387 0.1750493331599814 -0.03766603971225978 0.1709489353349375  -0.06826280306531989 1.187856672441692e-18  0.4696439858017911 0 
5 0.6000000000000003 4.347566194517904e-31 0.522836699221925 0.1241748744185748 0.08646088037588485 -0.08912864635736237  -0.08549260104635911 -2.962338190050741e-18  0.4201106585244942 0 
6 0.7000000000000004 3.489444957524124e-31 0.5749102271803508 0.1011715198574469 -0.0976387200131585 0.02650201472450365  -0.1010441980498919 5.451822975316556e-21  0.4186729083777906 0 
7 0.8000000000000005 2.518715624135342e-31 0.6239223249834378 0.09120611774931656 0.08618723343718702 0.02983817533539922  -0.1154881144630473 6.898065147620218e-21  0.4458948103470515 0 
8 0.9000000000000006 3.125407870699134e-31 0.6772930100419619 0.08534445894527673 -0.04419151833824479 -0.07301223445164925  -0.1296097646972653 -1.12223476569221e-29  0.4776400316367463 0 
9 1.000000000000001 4.061372012275019e-31 0.73817093544951 0.07863622906292486 -0.01816263681471499 0.07650996761973924  -0.1437873911014145 -1.854906670661934e-17  0.4955330688229499 0 
~~~~
We can see that the entropy starts near zero at t=0.1 (or 0.09999... here) as expected for an initial product state, and that the value of Sx on chain 4 (8th column) starts near zero as expected, while the value of Sz on chain 5 (10th column, 2nd from last) starts near +1 before decreasing.
