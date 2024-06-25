var documenterSearchIndex = {"docs":
[{"location":"api/","page":"API","title":"API","text":"Modules = [RGRDMT]","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = RGRDMT","category":"page"},{"location":"#Renormalization-Group-Reduced-Density-Matrix","page":"Home","title":"Renormalization Group Reduced Density Matrix","text":"","category":"section"},{"location":"#Motivation","page":"Home","title":"Motivation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"We want to find the ground state of a local Hamiltonian for a spin-system on some lattice.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A spin-system is a collection spin each with Hilbert space mathbbC^d. For simplicity, we will consider spins residing on a chain with uniform distance a.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The local Hamiltonian is denoted as H = sum_i^N h_i where h_i only acts non-trivially on a few spins that are physically close together. For simplicity, we will consider a 2-local Hamiltonian H = sum_i^N h_i^(2). More concretely, we could think of the Hamiltonian having the form:","category":"page"},{"location":"","page":"Home","title":"Home","text":"H = sum_i^N X_iX_i+1 + Y_iY_i+1 + Z_iZ_i+1","category":"page"},{"location":"","page":"Home","title":"Home","text":"The ground state energy may diverge in case N rightarrow infty. Therefore, we rephrase our goal into finding the ground state energy per-site for a Hamiltonian.","category":"page"},{"location":"","page":"Home","title":"Home","text":"e_0 = min_ketpsi frac1N brapsi H ketpsi","category":"page"},{"location":"","page":"Home","title":"Home","text":"Due to translational invariance of the Hamiltonian, we would expect the ground state energy to also be reached with the following equation","category":"page"},{"location":"","page":"Home","title":"Home","text":"e_0 = min_ketpsi brapsi h_i^(2) ketpsi","category":"page"},{"location":"","page":"Home","title":"Home","text":"Since h_i^(2) has no support on spins other than the two at location i and i+1, we could safely simplify the minimization. Instread of trying to minimize e_0 over an exponentially large state ketpsi, similar result is achieved when you do the minimization over reduced density matrix on spin i and i+1, we denote it as rho^(2).","category":"page"},{"location":"","page":"Home","title":"Home","text":"e_0 = min_rho^(2) tr( rho^(2) h_i^(2))","category":"page"},{"location":"","page":"Home","title":"Home","text":"However, this naive reduction will bring problems because rho^(2) may not correspond to a physical state. Therefore, we need to add the constraint that rho^(2) be the reduced density matrix on two spins obtained from a reduced density matrix on three spins. This constraint alone will make the solution more physical. However, we would still improve by posting a series of constraint for our 1D example:","category":"page"},{"location":"","page":"Home","title":"Home","text":"rho^(i) = tr_L(rho^(i+1)) = tr_R(rho^(i+1))","category":"page"},{"location":"","page":"Home","title":"Home","text":"where tr_LR denotes the partial trace of the left/right most spin's Hilbert space. Using a graphical tensor network representation, we could visualize the constraints as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"When the above constraints are taken into account up to n spins' reduced density matrix, we denote the all such two-spin reduce density matrix rho^(2) with set mathcalS.","category":"page"},{"location":"","page":"Home","title":"Home","text":"A good news is that mathcalS is a convex set. Furthermore, the e = tr(rho^(2) h_i^(2)) is a linear function and the constrains are also linear. Therefore, we could use convex optimization to find the ground state energy per-site.","category":"page"},{"location":"","page":"Home","title":"Home","text":"rho^(2)","category":"page"},{"location":"","page":"Home","title":"Home","text":"that gives the same e = tr(rho^(2) h_i^(2)) forms a hyperplane. Finding e_0 is equivalent to finding the hyperplane with lowest energy. This is visualized in the following diagram","category":"page"},{"location":"","page":"Home","title":"Home","text":"<img src=\"pics/space-shrinking.png\" width=\"400\">","category":"page"},{"location":"","page":"Home","title":"Home","text":"When we add more constraints, this amounts to cutting mathcalS. When we follow the lines of previous reasoning and adds more restrictions in the form of","category":"page"},{"location":"","page":"Home","title":"Home","text":"rho^(i) = tr_L(rho^(i+1)) = tr_R(rho^(i+1))","category":"page"},{"location":"","page":"Home","title":"Home","text":"we pay a price of adding exponentially more constraints since the dimension of rho^(i) is d^2i. Another observation is that not all constraints are essential to obtain the hyperplane with e = e_0. This promotes simplification to the algorithm. Rather than strictly requiring all the constraints be satisfied","category":"page"},{"location":"","page":"Home","title":"Home","text":"rho^(i) = tr_L(rho^(i+1)) = tr_R(rho^(i+1))","category":"page"},{"location":"","page":"Home","title":"Home","text":"we could instead do isometries and pick out part of the reduced density matrix that are closestly related to ground states. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"<img src=\"pics/isometry_mapping.png\" width=\"400\">","category":"page"},{"location":"#Implementation","page":"Home","title":"Implementation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"As a result, we obtain the following plot for how close the lower bound was for the ground state energy. We obtain the ground state energy for the Heisenberg XXX model H = sum_i^N X_iX_i+1 + Y_iY_i+1 + Z_iZ_i+1. The theoretical value for the ground state energy per site is e_0 = 14 - ln2. Let the ground state energy per site obtained from convex optimization with relaxed constrains be e_rlx, we define the gap of groun state energy as Delta E_relax = e_0 - e_rlx","category":"page"},{"location":"","page":"Home","title":"Home","text":"<img src=\"pics/test_xxx.png\" width=\"400\">","category":"page"},{"location":"","page":"Home","title":"Home","text":"This result should be compared with the result from paper ","category":"page"},{"location":"","page":"Home","title":"Home","text":"<img src=\"pics/paper_xxx.png\" width=\"400\">","category":"page"},{"location":"#Reference","page":"Home","title":"Reference","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"[1] I. Kull et al., Phys. Rev. X 14, 021008 (2024).","category":"page"},{"location":"","page":"Home","title":"Home","text":"[2] V. Zauner-Stauber et al., Phys. Rev. B 97, 045145 (2018).","category":"page"},{"location":"","page":"Home","title":"Home","text":"Documentation for RGRDMT.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"}]
}
