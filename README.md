# gelation

This is the repository of my first year of graduate studies lab project I am broadly studying percolation theory in gel systems and the conformations of low-functionnality polymers.

## Repository structure 

```
* simple-models/
  ** bethe.py  : Describes gelation using mean-field theory approximation
  Gelation transition occurs on a Bethe lattice, a special type of bond-
  percolation model (oftenly called Bernouilli percolation)
```


## DLCA Notes / Ideas 

Criterion to stop the simulation earlier for the computation of $P(\phi_0, L)$ ; 

- If $N = \phi_0 L < L $ the system cannot percolate $\rightarrow$ we continue to the next run.
- The system exhibits low convergence for low densities $\phi_0$, due to two main reasons. First, the configuration of the system may be too large to be explored in finite time. And oftenly, the number of configurations that leads to percolation by far fewer. Second, the more the particles aggregate, the more the dynamics become slow as probability acceptance is inversely proportionnal to the mass of the cluster. So, after a finite number of steps, if no aggregation occurs, it may be wise to stop the simulation there. May introduce some bias.

## References

- Polymer physics, Rubinstein & Colby 
