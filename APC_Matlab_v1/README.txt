===========================================================================

                  Adaptive-Picard-Chebyshev Iteration 
                     (perturbed two-body problem)

                              VERSION 1
             
                        Texas A&M University
                  Department of Aerospace Engineering
                               May 2017

R. Woollands and J. Junkins, Nonlinear Differential Equation Solvers via 
Adaptive Picard-Chebyshev Iteration: Applications in Astrodynamics, Journal
of Guidance, Control and Dynamics, Jan 2019.

         contact: junkins@tamu.edu OR robyn.woollands@gmail.com

===========================================================================

1. run_example.m

This file allows the user to set up a test. The user should only make changes
between "BEGIN USER INPUT" and "END USER INPUT". The Adaptive-Picard-Chebyshev
integrator/propagator is called from this file. After integration the Hamiltonian
can also be computed as a check of solution accuracy.

2. adaptive_picard_chebyshev.m

This file performs the orbit propagation using variable Chebyshev order (N) 
and segment size.

3. jacobi_integral.m

This file computes the Hamiltonian at each point along the trajectory for the
purpose of checking the solution accuracy.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                           ADAPTIVE_PICARD_CHEBYSHEV.M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PART 1 - ASSIGN VARIABLES

PART 2 - DETERMINE SEGMENTATION SCHEME

The first step in is to find the Keplerian perigee. This is required 
since the segmentation scheme for the perturbed two-body problem is based on
symmetry about perigee. More details can be found in the paper. 

Using the state at perigee we start with 3 segments per orbit (equal in true 
anomaly) and compute the Keplerian solution for the first segment using Battin's 
analytical F&G solution (FnG.m). This gives us a trajectory that experiences 
a very similar gravity field to that of the final solution.

Using N=10, thus 11 discrete sample points along the Keplerian trajectory, we 
compute the acceleration at each point considering the full spherical harmonic 
gravity series for the user specified gravity degree. Chebyshev polynomials of 
order N-1 are using to do a least squares fit of these discrete acceleration
points. The magnitudes last 3 coefficients for each component (x,y,z) of the 
acceleration fit are computed. An adequate fit is achieved once the magnitude
of the last 3 coefficients of each component of the acceleration falls below
tol/100. If this does not occur using N=10 then the number of nodes
is double to 2*N, and then 4*N. (See the discussion on node doubling in the 
paper). If three segments per orbit and N=40 nodes still does not satisfy tol/100
then 2 more segments are added to the orbit (now a total of 5) and the same
node doubling procedure is repeated. For low Earth orbits considering a high
fidelity gravity model, 11 segments per orbit with 40 nodes per segment is common.

Once the number of segments per orbit and the number of nodes per segment has 
been computed (using the above techniques), one final check is performed. That
is, the acceleration fit is computed using least squares and then a difference 
is calculated between the least squares fitted acceleration and the actual 
acceleration computed using the spherical harmonic series. The number of nodes
is adjusted one final time if necessary. That is, if the acceleration fit is a 
lot more accurate than tol/100 then N is reduced incrementally until the minimum value
of N is obtained that still satisfies tol/100. 

PART 3 - Adaptive Picard-Chebyshev Iteration Initial Setup

The segment start and end times are computed for the optimal true anomaly 
segmentation scheme discussed in the previous section. Note that if the user
specified intial conditions are not at perigee then a short first/last segment 
is needed prior to beginning the optimal segmentation scheme. The times for
these short segments are also computed in this section. Finally the Clenshaw-Curtis
quadrature constant matrices are computed in this section. This computation may
be done a priori and all the matrices can be stored and called when needed. This
is not done in the MATLAB code but is done in the C code as the C code is where
optimiaztion for computational efficiency is done. The MATLAB code is a proof of
concept and an optimization is done with regard to function evaluations, as timing in 
MATLAB can be unreliable.

PART 4 - ADAPTIVE PICARD-CHEBYSHEV ITERATION (loop through all segments)

Picard iteration is done sequentially for each of the segments until the final
time is reached. In PART 4A (see code) the radially adaptive gravity, variable
fidelity gravity model, hot start and quasilinearization are applied and computed
is the user has specified this in run_example.m. The paper has more details on 
these features of the algorithm.

In PART 4B (see code) "Perigee Passaged" is computed and re-oculated at the end
of each orbit. If left for several orbits the perturbed trajectory perigee would 
start to shift from the Keplerian perigee and the precomputed segmentation times
would not synchronize with the Keplerian true anomaly breaks. This leads to a redeuction 
in the accuracy of the final solution. To prevent this we simply re-osculate perigee
after each orbit propagation.

PART 5 - INTERPOLATE SOLUTION

Once the solution for each segment has been obtained using Picard iteration, the
coefficients are stored in an array. These coeffients are used for interpolating
the solution onto the user specified output times.

Finally, the interpolated solution (position and velocity), the corresponding 
times, the number and segments per orbit, the number of nodes per segment and
the total number of equivalent function evaluations are output back to run_example.m
The Hamiltonian may be computed if desired.









