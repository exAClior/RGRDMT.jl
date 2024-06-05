# Variational Uniform Matrix Product State

---

A **variational** algorithm for finding the **ground state** of a **translational invariant Hamiltonian**. 

## Key Ingredients of Variational Method 

1. translational invariant Hamiltonian: $H$  
2. Matrix Product State (MPS) $\Psi(A)$ parameterized by tensor $A$
3. Methods to take gradient and update ansatz.  
    ![MPS Parameterized by A](image-6.png){ width=90% }

--- 

### How to update

- Gradient g: 
![Gradient definition](image-13.png)
![Simplification Condition](image-14.png)
![Simplified Gradient Definition](image-15.png)
![Unsimplified Tensor for gradient](image-16.png)
![Simplified Version](image-19.png)
![Partial Contraction Definition](image-17.png)
![Transfer Matrix Definition](image-18.png)

$A' = A .- g$

- 


### Hello

- Changes that keeps MPS on the manifold representable by $| \Psi(A) \rangle$
![Tangent Space](image-7.png){ width=70%}
- Tangent vector: $| \Phi(B;A) \rangle := B^i \frac{\partial}{\partial A_{i}} | \Psi(A) \rangle := B^i |\partial_{i}\Psi(A)\rangle$
![Tangent Space graph representation](image-8.png)

#### Choice of B
- Gauge freedom of $A$ gives gauge freedom of $B$
- There are some gauge that makes computation of $\langle \Phi(B;A) | \Psi(B;A) \rangle$ easier
![Definition of Good B](image-9.png)
- E.g for easiness
![Overlap is easy to compute](image-10.png)
- Construction of B
![Construct V_l](image-11.png)
![Construction from X](image-12.png)
---