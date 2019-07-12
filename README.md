# Numerical-Unconstrained Optimization Solver

In this project, I implemented Gradient Descent, Newton, Modified Newton, DFP and BFGS numerical optimization algorithms in C programming language. All of these algorithms are boosted with Armijo step size rule which helps them to decide an appropriate step size at each iteration.

Note 1: The implementation requires that the user has to provide the gradient (and the Hessian if requested by the algorithm) of the objective function!

Note 2: The uploaded versions of the algorithms have been provided with the gradient and the Hessian of the extended Rosenbrock function which constitutes a benchmark in numerical optimization. The algorithms have been also tested for several other benchmark functions as shown in the file "Unconstrained Presentation.pdf". In this file, one can also see a comparison made between my solver and AMPL's MINOS solver (commercial solver) where it is shown that in some cases, our algorithms outperform AMPL's MINOS solver in the number of iterations required until convergence.  
