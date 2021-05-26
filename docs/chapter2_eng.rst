Special methods of modeling quantum spin systems
=================================================

This section presents selected numerical methods that are commonly used to solve problems of quantum many-body systems. The code and descriptions below serve as an example of how these methods can be implemented. For most tasks, this code needs some additional work. I would recommend using the more advanced libraries that can be found on github.com and other platforms.


Numerical Renormalization Group (NRG) method
----------------------------------------------

.. seealso:: The used method  is described in detail in section 2.2 of the book `Strongly Correlated Systems. Numerical Methods.`_


.. _`Strongly Correlated Systems. Numerical Methods.`: https://www.springer.com/gp/book/9783642351051

The goal of the NRG method is to find the ground basis state of a large system when direct diagonalization is impossible due to the enormous computational time.

Consider an example with the Heisenberg Hamiltonian, although this example is not optimal, since The NRG method does not solve the problem well for systems in which there is a lot of entanglement. For entangled systems, it is optimal to use the DMRG method, which will be discussed below.

1. We form the Heisenberg Hamiltonian of size :math:`2^m\times2^m` in an iterative way. 

2. We repeat the following operations :math:`n-m-1` times (find the ground base state of the system :math:`2^n\times2^n`)

- **find the eigenvectors of the Hamiltonian**

- **we rewrite the Hamiltonian in the basis of eigenvectors**. In this case, the first state should correspond to the state with the lowest energy.

- **truncate the Hamiltonian** from the size :math:`2^{m+1}\times2^{m+1}` to :math:`2^m\times2^m`, deleting states with energy far from the ground state (at the first step, we skip this point, since the Hamiltonian of the initial size :math:`2^m\times2^m`).

- **truncate the basis of eigenvectors** from size :math:`2^{m+1}\times2^{m+1}`  to size :math:`2^m\times2^m`  (at the first step we skip this point, because the basis of the original size :math:`2^m\times2^m`

- **rewrite operators** :math:`\hat{S}_x`, :math:`\hat{S}_y` and :math:`\hat{S}_z` that are necessary for the iterative formation of the Hamiltonian of the larger system, in the basis of the truncated eigenvectors of the Hamiltonian.

- **increase the Hamiltonian by adding one spin**, as in the iterative method, using the previous Hamiltonian and :math:`\hat{S}_x`, :math:`\hat{S}_y` and :math:`\hat{S}_z`

3. After the end of the cycle, we obtain a Hamiltonian whose eigenvalues are "close" to the ground state of the Hamiltonian of size math:`2^n\times2^n`

.. literalinclude:: _static/nrg.py
   :language: python
   :linenos:

   
Density matrix renormalization group (DMRG)
----------------------------------------------

.. seealso:: The used method  is described in detail in sections 2.3 - 2.5 of the book `Strongly Correlated Systems. Numerical Methods.`_ 

.. _`Strongly Correlated Systems. Numerical Methods. Editors: Adolfo Avella, Ferdinando Mancini.`: https://www.springer.com/gp/book/9783642351051

The goal of the DMRG method is to find the ground state of a large system when direct diagonalization is impossible due to the enormous computational time. Unlike the NRG method, this method is more accurate because it accounts for entangled states. 

Consider an example with the Heisenberg Hamiltonian, where we need to find the eigenvalues of the Hamiltonian of size :math:`2^n\times2^n`.

1. We form the left Heisenberg Hamiltonian of size :math:`2^m\times2^m` in an iterative way. 

2. We form the right Heisenberg Hamiltonian of size :math:`2^m\times2^m` in an iterative way. 

3. Combining the left and right Hamiltonian, we get the superblock Hamiltonian :math:`2^{2m}\times2^{2m}`

4. Find the ground state of the superblock system.

5. Get the truncated density matrices of the left and right subsystems from the ground state of the superblock system, each size :math:`2^m\times2^m`.

6. We repeat the following operations :math:`n-m-1` times (where n – the size of the system whose ground state we find):

- **find the eigenvectors** of each truncated density matrix.

- **change the order of the basis vectors**, the first one should become the last, because the eig () function sorts vectors starting with the vectors with the lowest eigenvalue. The eigenvalues have the meaning of the probability of a given state, therefore, it is necessary to rearrange the order of the eigenvectors, the first should be the vector with the maximum eigenvalue. 

- keep only the first :math:`2^m` eigenvectors, **truncate** the rest

- **rewrite the current left and right Hamiltonian** size :math:`2^m\times2^m` each in basis of :math:`2^m` eigenvectors of the corresponding density matrix

- **rewrite operators** :math:`\hat{S}_x`, :math:`\hat{S}_y` and :math:`\hat{S}_z`, necessary for the iterative formation of the Hamiltonian of the larger system, in the basis of the eigenvectors of the corresponding density matrix

- **form the left and right Hamiltonians in a new basis**, increased by 1 site size :math:`2^{m+1}\times2^{m+1}`

- creating new operators :math:`\hat{S}_x`, :math:`\hat{S}_y` and :math:`\hat{S}_z`, increased by 1 site. This is necessary to combine the left and right blocks.

- we combine the left and right Hamiltonian, we get the superblock Hamiltonian :math:`2^{2m+2}\times2^{2m+2}`

- find the ground state of the superblock system.

- get truncated density matrices of the left and right subsystems from the ground state of the superblock system, each of size :math:`2^m\times2^m`.

The results of calculating the ground state energy will be equal to the ground state energy values of the left or right matrix.

.. literalinclude:: _static/dmrg.py
   :language: python
   :linenos:
   

Bethe ansatz
--------------

The method used is described in detail in the article `Introduction to the Bethe ansatz I`_

.. seealso:: M.Karbach and G.Muller, “Introduction to the Bethe ansatz i,” (1998), arXiv:cond-mat/9809162 [cond-mat.stat-mech].

.. _`Introduction to the Bethe ansatz I`: https://arxiv.org/pdf/cond-mat/9809162.pdf

Here is a description of a method for obtaining and classifying two-magnon states. The one-magnon states are plane waves with a wave vector :math:`k_m`:

.. math::

   |\psi\rangle=\sum_{x=1}^n\frac{1}{\sqrt{n}}\exp{(ik_max)}|x\rangle,

where :math:`n` - number of spins, :math:`a` - lattice constant. The dispersion law for one-magnon states:

.. math::

   E_m=4J(1-\cos{k_ma})-E_0.

where :math:`E_0` is the ground state energy of a ferromagnet.

Two-magnon states are characterized by two wave numbers :math:`k_1` and :math:`k_1`, which can be both real and complex, and satisfy :math:`k_1+k_2\equiv k_{2m}=2\pi/n(\lambda_1+\lambda_2)`, where :math:`\lambda_i` Bethe integers. The dispersion law for two-magnon states:

.. math::

   E_{2m}(k)=4J\sum_{i}(1-\cos{k_{i}a})-E_0,

where :math:`i=1,2`. For 10 spins, there are 45 two-magnon states in the chain.

All states can be divided into 3 classes - C1, C2 and C3. The first class represents waves with the following values of the Bethe numbers :math:`\lambda_1=0`, :math:`\lambda_2=0,1...,n-1`. The number of such states is :math:`n`.

The second class satisfies the following conditions :math:`\lambda_2 - \lambda_1 > 1`. The number of such states is :math:`\frac{n(n-5)}{2}+3`. 

The third class satisfies the following conditions:math:`\lambda_1 - \lambda_2 < 2`.  The number of such states is :math:`2n-3`. 

We get the Hamiltonian block for 10 spins with 2 magnons:


Below is the code of the program, which successively finds three classes of two-magnon eigenstates using the example of the Hamiltonian for 10 spins.

.. literalinclude:: _static/bethe_get_10_spins.py
   :language: python
   :linenos:

Getting and checking the found states can be done as follows:

.. literalinclude:: _static/bethe_ansatz.py
   :language: python
   :linenos:
   
The code of the functions used is presented below:

.. literalinclude:: _static/bethe_ansatz_check.py
   :language: python
   :linenos:

Code for helper functions:

.. literalinclude:: _static/bethe_helpers.py
   :language: python
   :linenos:

