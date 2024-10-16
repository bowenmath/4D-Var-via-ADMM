# ADMM for 4D-Var

This is the repository for the source code of the paper

> [Numerical Solution for Nonlinear 4D Variational Data Assimilation (4D-Var) via ADMM](https://arxiv.org/abs/2410.04471), Bowen Li and Bin Shi.

The codes are implementations of linearized multi-block ADMM for solving four-dimension variational data assimilation problem (4D-Var) of Lorenz63 model, Burgers equation, and 2D vorticity equation.

## Requirements

- Python 3
- Numpy
- Scipy
- Matplotlib
- Matlab

## File structure

The file structure is as following.

- lorenz63  
    > ADMMlorenz63 - rk4 - noise.py  
    *This is the main program with ADMM implementation for Lorenz 63 model.*  

    > objective_function - landscape.py  
    > objective_function - X.py  
    > objective_function - Y.py  
    > objective_function - Z.py  
    *These are the programs generating the landscape figures of the objective function in the paper.*  

- burgers  
    > ADMMburgers - difference - noise.py  
    > ADMMburgers - FEM - noise.py  
    > ADMMburgers - spectral - noise.py  
    *These are the main programs with ADMM implementation for Burgers equation with finite difference method, finite element method, spectral method, respectively.*  

    > difference - Dirichlet  - simulation.py  
    > FEM - Dirichlet  - simulation.py  
    > spectral - Dirichlet  - simulation.py  
    *These are the numerical simulation for Burgers equation by finite difference method, finite element method, spectral method, respectively.*  

    > verify_gradient.py  
    *This is the verification of the tangent linear model.*  

The above codes are implemented by Python. The below uses Matlab for the implementation.

- qgmodel  
    > ADMMqgmodel.m  
    *This is the main program that implements ADMM for 4D-Var of 2D vorticity equation.*  

    > arakawa.m  
    > arakawa_g.m  
    > bdcondition.m  
    > f.m  
    > f_adj.m  
    > inversepoisson.m  
    > laplacian.m  
    > subproblem.m  
    *These are the functions used in the implementation.*  

    > draw.m  
    *Use this program to draw the figure after running the main program.*

    > f_g.m  
    > verify_gradient.m  
    *These are the programs to verify the tangent linear model and adjoint model.*

    > bluered.m  
    *This is the customized display color.*
