# Numerical Unconstrained Optimization Solver

In this project, I implemented Gradient Descent, Newton, Modified Newton, DFP and BFGS numerical optimization algorithms in C programming language. All of these algorithms are boosted with Armijo step size rule which helps them to decide an appropriate step size at each iteration.

Several "tricks" were used in order to implement the aforementioned algorithms in C. For example, Newton's method require the knowledge of the Hessian matrix and someone needs to take its inverse in order to solve the direction finding problem. In order to solve this problem in C, we applied Cholesky Decomposition of the inverse Hessian and combined the Cholesky factors in order to solve the direction finding problem.

Last, we made a comparison between the implemented algorithms and AMPL's MINOS solver (commercial solver) in 10 benchmark functions. Interestingly, we found that our implementation of the Modified Newton algorithm converges in less steps than the MINOS solver and with much less gradient evaluations. 

Note 1: The implementation requires that the user has to provide the gradient (and the Hessian if requested by the algorithm) of the objective function!

Note 2: The uploaded versions of the algorithms have been provided with the gradient and the Hessian of the extended Rosenbrock function which constitutes a benchmark in numerical optimization. The algorithms have been also tested for several other benchmark functions as shown in the file "Unconstrained Presentation.pdf".  
