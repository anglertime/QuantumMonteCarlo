Quick description of the files
------------------------------


[A] Basic files (modified frequently):

* startup.m - initialization file. Add the path of 'BoostGraphLibrary' into MATLAB directory.
* vmc.m - Starting file. Contains most of the metropolis algorithm.
* read_configuration.m - Contains configuration information like K (currently set to 1) and grid size.
                         Also generates the box for the system.
* fgmres_reorder_local.m - Iterative methods to find the approximate ratio. We decompose the domain
                    into two parts. With reordering, we pick paired particles and orbitals locally.
* inexact_gmres.m - Solve Schur system (A11 - A12*A22^(-1)*A21)z1 = e1 with GMRES. When we apply the 
                    optimal reordering, we use diagonal preconditioner.
                         
[B] Files computing observables (no description as name says it all):

* Energy.m
* CalcPairCorrelation.m
* CalcStructureFactor.m

[C] Reordering files

* BoostGraphLibrary - Contains the main function max_flow.m, which implement the max flow method to find a
                      perfect/maximum matching.
* Munkres_Hungarian.m - Munkres (Hungarian) Algorithm for Linear Assignment Problem. This is matrix 
                       interpretation of Hungarian algorithm for max weighted matching of bipartite graph.
* Hungarian.m - A function for finding a minimum edge weight matching given a MxN Edge weight matrix 
                WEIGHTS using the Hungarian Algorithm. This algorithm is slower than Munkres_Hungarian.m in MATLAB.
* max_min_diagonal.m - Maximize the smallest element in the diagonal of a matrix. 
* max_weighted_match.m - Find Maximum weighted matching with threshold of diagonal elements.

[D] Infrastructure files (hardly modified):

* mygmres.m - GMRES algorithm.
* rotmat.m - Compute the Givens rotation matrix parameters.
* dense_detRatio.m - Computes the dense determinant ratio.
* compute_matrix.m - Builds the Slater matrix.
* convert_dist.m - Compute the distance taking account of the periodicity.
* pos_to_string.m - Performs formatting (of particle and orbital positions) for display.
* put_in_box.m - Ensures that the domain is the box (e.g. particles remain inside the box).
* random_move.m - Randomly identifies a particle to move, and also the amount by which this particle is to be moved.
* move_particle.m - Finds the particle's configuration after a random move. Also computes 'u' to find 
                   the determinant ratio (as 1+mhat'*u).
* save_config.m - Saves the configuration to a file.


[E] Data and text files:

* particle30K - Particle configuration after 30,000 Monte Carlo steps for a 1024 system. Notice this 
                configuration has already been reordered.
* matrix30K - The slater matrix in 30,000 Monte Carlo step corresponding to the orbital1024 and 
               particle30K configuration.
* particle1024 - Particle configuration for a 1024 system. All particles are about 0.3h far from its 
                 coupled orbital in a random direction.
* orbital1024 - Orbital configuration for a 1024 system.



Execution
---------
>> startup
>> vmc
