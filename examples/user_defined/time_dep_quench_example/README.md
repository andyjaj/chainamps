# TIME DEPENDENT QUENCH PROTOCOL EXAMPLE #

In this example we consider evolving a state using a time dependent hamiltonian.
We will take a simple lattice spin chain, initially with spins on alternating sites pointing in the +x and -x directions.
We then quench by turning on a field in the x direction (which should do nothing) from time 0 to 0.015.
From time 0.015 to 0.1 we switch off the x field and turn on a field along the z direction.
From time 0.1 onwards we switch off the z field and turn on an S^x S^x type coupling between neighbouring sites.

The input file, `time_dep_input` contains the following:
~~~~
USER DEFINED
2 ./vertex_hamiltonian.in ./Sx.in ./Sy.in ./Sz.in
0.0 hx,(Sx):10.0
0.015 hz,(Sz):1.0
0.1 Jx,(Sx,Sx):10.0
~~~~
The model is user defined, the 2 tells us that the basis consists of only two states and that we are not taking advantage of symmetries, (so no quantum number conservation information is provided).
To take advantage of conservation laws we would need to specify the basis in a separate file (replacing 2 with its filename); see the `spin_half_example` directory for this.
The `vertex_hamiltonian.in` file is empty, because for this example we can write local contributions simply in terms of the operators we define.
The next three files separately define the operators for our model.
For example `Sy.in`:
~~~~
Sy
0 1	(0.0,-1.0)
1 0	(0.0,1.0)
~~~~
Which defines the operator's name and the non zero elements, for the Pauli sigma^y matrix.
The last three lines of `time_dep_input` define the times at which we apply certain hamiltonians.
The lines are *not cumulative, the hx term is only applied from time 0.0 to 0.015* and is not part of the evolution hamiltonian at later times.

## Using a c-specifier file ##

For this quench we use a c-number specification file to setup the initial (time<=0.0) state.
The contents of this file, `x-state` are:
~~~~
0 1.0 1 1.0
0 1.0 1 -1.0
~~~~
The first line produces a superposition on the first vertex, with equal coefficients for the 0 and 1 basis states (spin pointing in the +x direction).
The second line produces a superposition with equal but opposite coeffs on the second vertex (spin pointing in the -x direction).
These vertex states will be normalised by the driver, only the relative weights of basis states for a single vertex matter, i.e.
~~~~
0 30.0 1 30.0
0 5.0 1 -5.0
~~~~
should produce exactly the same results.
Also, as we have only specified the vertex states for the first two vertices, if we choose to model a system with more vertices, the values will cycle.
So for 4 vertices, vertex 1 and vertex 3 will be in the same initial state, and 2 will be the same as 4.

## The quench ##
Copy the files in this directory to your run directory (containing the binaries).
We quench by executing
`./TEBD_DRV.bin -c x-state -B 50 -n 20 -s 0.01 -m 2 -M Sx,5 -M Sx,6 -M Sz,5 time_dep_input 10`
which specifies a chain of 10 sites initially setup according to our definitions, does 20 time steps (of varying size, due to the time dependent hamiltonian, but with a maximum size of 0.01), and every two steps measures Sx on sites 5 and 6 and Sz on site 5.
The bond dimension is set at 50, but can grow slightly if there are degenerate singular values.
The main output file is `time_dep_input_TEBD_10_x-state_Evolution.dat` and should look like:
~~~~
# Index, Time, Truncation, Entropy, abs(Overlap), Real(Overlap), Im(Overlap), Re(Sx,5), Im(Sx,5) , Re(Sx,6), Im(Sx,6) , Re(Sz,5), Im(Sz,5) 
0 0.015 0 0 0.9999999999999997 0.9887710779360419 -0.1494381324735997  1 0  -0.9999999999999998 0  1.447857316597016e-15 0 
1 0.03388888888888889 0 0 0.9983956373289458 0.9420715830325234 -0.33059791451258  0.999286504615177 -1.387778780781446e-17  -0.9992865046151772 0  1.877404481875899e-15 0 
2 0.05277777777777777 0 0 0.9935968377823938 0.9375435056687534 -0.3290088920219585  0.9971470366120379 -3.469446951953614e-18  -0.997147036612038 0  2.019163915928379e-15 0 
3 0.07166666666666666 0 0 0.9856462712099729 0.930041466840815 -0.326376227545311  0.993584648991674 3.469446951953614e-18  -0.9935846489916742 0  2.003822455187709e-15 0 
4 0.09055555555555556 0 0 0.974614405636922 0.9196319591510475 -0.3227232551010492  0.988604425248338 0  -0.988604425248338 0  1.81452075587174e-15 0 
5 0.11 0 0.0001135324469132679 0.9659727189964782 -0.8171399885879703 -0.5151558336036784  0.9839430018294039 0  -0.9855847669095603 0  -0.06850261659298423 0 
6 0.13 0 0.0007652891622311049 0.9659727189964792 0.6873392051526953 -0.6787253574941133  0.9839430018294048 0  -0.985584766909561 0  -0.1639550817693264 0 
7 0.15 0 0.001570958532181881 0.9659727189964796 0.5048101742372266 0.8235714794916302  0.9839430018294042 -8.961326462605135e-21  -0.985584766909561 0  -0.1599545944075233 0 
8 0.17 0 0.002084709635840938 0.9659727189964789 -0.9167270631708528 0.304491026954151  0.983943001829404 0  -0.9855847669095605 0  -0.05892779645992556 0 
9 0.1900000000000001 0 0.002041023245540212 0.9659727189964792 -0.08824555620871076 -0.9619334777700931  0.9839430018294043 5.549914795405862e-17  -0.9855847669095608 0  0.07784381208617534 0 
~~~~
Note that at time 0.015, the value of Sx on site 5 is +1, while on 6 it is -1 (columns 8 and 10).
The entanglement is zero up to time 0.1, which we expect, because the sites were not coupled before that time, and the system started in a pure state.
Also note that the time dependency of our hamiltonian means that the time step size have to change.
Because the hamiltonian changes at time 0.015, and our maximum time step size is 0.01, the step size is decreased to fit an integer number of steps (2) of the same size into this interval.
A similar effect occurs from 0.015 to 0.1.
Remember we measure on every 2nd step.
