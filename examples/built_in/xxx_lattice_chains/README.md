# iDMRG EXAMPLE #
In this example we perform iDMRG on a system of coupled lattice xxx (Heisenberg) spin chains.
The chains are 4 sites long each, and we build an infinitely long cylinder of them.
Copy the file and folder to your run directory.
As this is a built in model the input file `xxx_input` is somewhat specific and fragile. The contents are:
~~~~
OLD_XXX
N 4.0
J 1.0
~~~~
The first line tells the driver to use the built in xxx model (the necessary files for this model are in the folder `xxx_files`).
The second line tells the driver to use 4 site long chains (only the 4 site data is included in the distribution).
The last line defines the *interchain* coupling parameter value.

To perform 10 iDMRG steps with bond dimension 50 we execute:
`./iDMRG_DRV.bin -B50 -N 10 xxx_input`
The results in `xxx_input_Energies.dat` are the raw DMRG energies, from growing the system by two vertices at each step.
~~~~
# Index, Energy, Energy/vertex, Entropy, Truncation, Fidelity
0 -4.820089374374787 -2.410044687187394 1.213702157008032 0 1.110223024625157e-16
1 -10.26428962097885 -2.566072405244713 0.9902948375147761 2.235632608666451e-05 0.0008040856330850943
2 -15.72808737031446 -2.62134789505241 1.063136610955919 2.362659753254866e-05 0.0001227665059162497
3 -21.19338287800255 -2.649172859750318 1.070924785801976 3.250689073942141e-05 1.576444779949338e-05
4 -26.65891718125858 -2.665891718125858 1.080850833767454 3.769018756618419e-05 6.082171550403359e-06
5 -32.12451133143336 -2.67704261095278 1.085358330031867 4.134034504301785e-05 1.856778964470784e-06
6 -37.59011717188822 -2.685008369420587 1.088087929509814 4.357133871693812e-05 7.646076590761552e-07
7 -43.05572367990966 -2.690982729994354 1.089603552383166 4.49265027057455e-05 4.784700518012386e-07
8 -48.52132853571689 -2.695629363095383 1.090473958002445 4.572699996835082e-05 4.405426354825437e-07
9 -53.98693181233662 -2.699346590616831 1.090970317247042 4.619526273561347e-05 4.593480406045458e-07
~~~~
Clearly these aren't converging very quickly.
On the other hand, the translationally invariant unit cell energies are in `xxx_input_iDMRGEnergies.dat`:
~~~~
# Index, Energy/vertex, Entropy
0 -2.732699938427228 1.070940173777507
1 -2.732783281508454 1.082287901099962
2 -2.732794934780465 1.085936253278948
3 -2.732797253003403 1.088019499900724
4 -2.732797157952852 1.089055766595033
5 -2.732796682018099 1.089628173051084
6 -2.732796277980149 1.089943744541338
7 -2.732795997051221 1.090121597833738
~~~~
These look a lot better, though the bond dimension is small, so they are not hugely accurate.
The unit cell at each iteration is stored in a UNITCELL file (in binary format).
Note that these translationally invariant cells are not generated in the very first two DMRG steps.

To perform post-processing measurements on the generated state we must first store the operators and basis using `STORE_MODEL.bin`.
`./STORE_MODEL.bin xxx_input`.
This produces the required BASIS and SPARSEMATRIX files.
We can now make a measurement, for example:
`/UNITCELL_MEASURE.bin -O xxx_input-Sz.SPARSEMATRIX -S 1 *.UNITCELL`
This measures correlations between Sz on nearest neighbour chains on all the generated UNITCELL files.
The output is in `UNITCELL_Results_xxx_input-Sz_xxx_input-Sz_1.dat`
~~~~
# Index,Re(,xxx_input-Sz(i),xxx_input-Sz(i+1)),Im(,xxx_input-Sz(i),xxx_input-Sz(i+1))
0 -0.1063686728140055 2.202373606458316e-19 
1 -0.1064560519487928 9.988263640493825e-19 
2 -0.1065214571355158 -1.747358944081729e-18 
3 -0.1065342297321957 1.615486827028107e-18 
4 -0.1065420124674465 5.518261156145134e-19 
5 -0.1065449811681246 2.810878673389773e-18 
6 -0.1065465386788641 -2.754507206699945e-18 
7 -0.106547321160436 4.066724044631399e-19 
~~~~

We can also consider correlations between different operators at different separation:
`./UNITCELL_MEASURE.bin -O xxx_input-S+.SPARSEMATRIX -O xxx_input-S-.SPARSEMATRIX -S 5 *.UNITCELL`
The results appear in the file `UNITCELL_Results_xxx_input-S+_xxx_input-S-_5.dat`.
Other possible measurements such as the entanglement, and transfer matrix eignevalues can be found.
For help type `./UNITCELL_MEASURE.bin`.

