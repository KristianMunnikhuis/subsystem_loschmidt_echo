# subsystem_loschmidt_echo

kmunnik@BU.edu

We are currently doing research on the subsystem Loschmidt Echo, a quantity closely related to the loschmidt echo in a many bodied quantum system. We are leveraging DMRG, Exact Diagonalization, and analytical methods to study this quantity.

This project uses both julia and python.

#Packages
Primarily we use ITensor, ITensorMPS for our julia coding.

We utilize QuSpin in python for exact diagonalization.


/src - function files. The model we use is the TFIM. Currently only for Analytic and Numerical functions for Julia.

/scripts - where a lot of the work is done:
##It should be noted that "mutual info" really refers to the connected correlation function of which mutual info is an upper bound.
      /mutual_info/exact_diag/ -includes all methods that use Exact Diagonalization to solve the TFIM (including thermal and adiabatic methods)
      /mutual_info/DMRG/  - Includes calculations of connected correlation functions using DMRG       
      /mutual_info/analytical/ - Includes analytical results for projector averages, thermal averages, mutual info, and functiosn to compare DMRG to the analytics.
                          - Analytics can be found summarized in "Quench dynamics in the Transverse Field Ising Model" lecture notes, but eventually I'll link my own summary.

      /correlations/ - Legacy folder for an older experiment. Will eventually house connected and unconnected correlation functions.

      /cluster_jobs/ - Messy folder for preparing jobs for the cluster. 

/results/ - Several *messy* folders. Needs clean up. Ideally data should be stored here.

Helpful Resources:

The quantum Ising chain for beginners:
-Incredible (but dense) covering of TFIM methods. Opened my eyes. Specifically, we are able to use ED with very large system sizes for integrable chains through methods used in this report.
https://arxiv.org /abs/2009.09208
Introduction to the TFIM:
-Contains analytic calculations of post-quench parameters and notation used in analytical calculations 
https://www.pks.mpg.de/fileadmin/user_upload/MPIPKS/group_pages/DQI/Introduction_TFIM.pdf
Quench dynamics and relaxation in isolated integrable quantum spin chains:
-Summary of many main results in quench dynamics for the TFIM.
https://arxiv.org/pdf/1603.06452

Late-time critical behavior of local stringlike observables under quantum quenches
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.064105
Observing Dynamical Quantum Phase Transitions through Quasilocal String Operators
https://link.aps.org/accepted/10.1103/PhysRevLett.126.200602
Probing quantum many-body dynamics using subsystem Loschmidt echos
https://arxiv.org/pdf/2501.16995



