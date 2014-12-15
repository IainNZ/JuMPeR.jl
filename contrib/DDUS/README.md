Data-Driven Uncertainty Sets (DDUS)
==========

The package DDUS contains implementations of many of the constructions for uncertainty sets from the paper "Data-Driven Robust Optimization" by D. Bertsimas, V. Gupta and N. Kallus, available from [here](http://arxiv.org/abs/1401.0212).   
Specifically, we have implemented oracles for each of the following sets (Eq. numbers refer to previous paper):
- UM - (Eq. 28) 
- UI - (Eq. 18) 
- UFB - (Eq. 23)
- UCS - (Eq. 35)
- ULCX - (Eq. 31)

More sets and additional features may be added going forward based on interest.  

For the most part, our implementations closely follow the descriptions in the paper.  In a few places, we have opted for simpler, approximate formulae for improved efficiency where we felt the difference in practice was neglible.   

##Usage

Each of our constructions are suitable to be used with JuMPeR as *cutting plane* oracles.  Reformulation may be supported in the future based on need.  A typical invocation might be:

```julia
dd_oracle = UCSOracle(data, epsilon, alpha)

m = RobustModel()
# ... Build up the model #

setDefaultOracle!(m, dd_oracle)  # to uses this oracle for all constraints
```

or 
``` julia
addConstraint(m, x[1] * us[1] + xs[2] * us[2] <= 5, dd_oracle)  #only for this one constraint
```

Most oracles support a simple constructor as above, taking in the data and two parameters, epsilon and alpha.  Some oracles require additional information, such as the support of the uncertainty. (When in doubt, check the source file for the interface marked "preferred interface.")  

All oracles assume that the data are given with each example in a row, and each column representing one component of the uncertainty.  **The ordering of the columns is important** and is assumed to correspond to the index of the uncertainties in the optimization model.  (That is, u[1] is the uncertainty whose data is given by column 1.)  The parameters epsilon and alpha are described in detail the above paper, and roughly control the probability of infeasibility and the decisionmaker's tolerance for ambiguity, respectively.  See also below on tuning these parameters.

Although fairly robust (*punny*), the preferred constructors for oracles can sometimes be slow because they perform all of the data analysis required to construct the set.  When possible, one can reuse the same oracle for multiple constraints.  When solving different optimizaiton problems in a loop, one can also used the specialized constructors for the oracles to customize the data analysis step.  (See the comments in the source code.)

## Examples
To be added.

## Choosing the "Right" set and Tuning Epsilon and Alpha in Practice
The cited paper proves that under certain conditions, each of the above sets satisfy a strong probabilistic guarantee.  In applications where it is important to have a provable guarantee on feasibilty, those results can help guide the choice of set. 

Many applicaitons, however, do not require provably good performance, just *practically* good performance.  In these cases, we suggest following the suggestions in Section 10 of the paper, and choosing the set, epsilon and alpha via cross-validation.  Some generic functionality to do this will (hopefully) be added soon.  In the meantime, ????? in the examples folder illustrates one possible cross-validation scheme for a particular example.  

## Citation
If you find this package useful, please consider citing the above paper as:

> @ARTICLE{2014arXiv1401.0212B,
   author = {Bertsimas, D. and Gupta, V. and Kallus, N.},
    title = "{Data-Driven Robust Optimization}",
    journal = {ArXiv e-prints},
   eprint = {1401.0212},
   keywords = {Mathematics - Optimization and Control},
   year = 2014,
  month = dec,
  url = {http://arxiv.org/abs/1401.0212},
}




