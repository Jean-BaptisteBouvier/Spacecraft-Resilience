# Spacecraft-Resilience

This project gathers the MATLAB codes used in the paper "Resilient Trajectory Tracking to Partial Loss of Control Authority over Actuators with Actuation Delay" by Jean-Baptiste Bouvier, Himmat Panag, Robyn Woollands and Melkior Ornik.
This work is not published and not available yet.

The fuel-optimal trajectory used as reference trajectory is designed with the convex optimization method described in ["Autonomous optimal trajectory planning for orbital rendezvous, satellite inspection, and final approach based on convex optimization"](https://link.springer.com/article/10.1007/s40295-021-00260-5) from Nicholas Ortolano, see reference below.


**Prerequisite**
---
The MATLAB toolboxes [mpt](https://www.mpt3.org/), [cvx](http://cvxr.com/cvx/) and [CORA](https://tumcps.github.io/CORA/) are required to run these codes.


**File Structure**
---

- Folder `functions` contains all the functions that will be called by the scripts of the repository.
- Folder `data` contains the precomputed reference trajectories and undesirable thrust signals. These are computed by the `main` file if not found in the `data` folder. The reference trajectory can be long to compute because of the optimization problem to solve. We store the undesirable thrust signals to be able to reuse them and compare different scenarios.
- The file `spacecraft_resilience.m` verifies the resilience of the spacecraft dynamics and compute its reachable set.
- The file `feedback_control.m` computes the parameters $\alpha$, $\beta$, $\gamma$, $\varepsilon$ and $\rho_{max}$ for a choice of matrices $K$, $P$ and $Q$ for the resilient feedback control Theorem 4.
- The file `main.m` computes the resilient trajectory tracking simulation.







**Citations**
---
Our work can be cited with:
```
@article{bouvier2023tracking,  
  title = {Resilient Trajectory Tracking to Partial Loss of Control Authority over Actuators with Actuation Delay},   
  author = {Jean-Baptiste Bouvier, Himmat Panag, Robyn Woollands and Melkior Ornik},    
  journal = {ArXiv},    
  year = {2023},   
  volume = {},  
  number = {},  
  pages = {}  
}
```

The optimal trajectory design can be cited with:
```
@article{Ortolano,
  title = {Autonomous optimal trajectory planning for orbital rendezvous, satellite inspection, and final approach based on convex optimization},
  author = {Ortolano, Nicholas and Geller, David K and Avery, Aaron},
  journal = {Journal of the Astronautical Sciences},
  pages = {444 -- 479},
  year = {2021},
  volume = {68},
  publisher = {Springer}
}
```



**Contributors**
---
- [Jean-Baptiste Bouvier](https://github.com/Jean-BaptisteBouvier)
- [Himmat Panag](https://www.linkedin.com/in/himmatpanag/?originalSubdomain=au)
- [Robyn Woollands](https://woollands.web.illinois.edu/index.html)
- [Melkior Ornik](https://mornik.web.illinois.edu/)



