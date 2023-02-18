# Spacecraft-Resilience

This project gathers the MATLAB codes used in the paper "Resilient Trajectory Tracking to Partial Loss of Control Authority over Actuators with Actuation Delay" by Jean-Baptiste Bouvier, Himmat Panag, Robyn Woollands and Melkior Ornik.
This work is not published and not available yet.


**Prerequisite**
---
The MATLAB toolboxes [mpt](https://www.mpt3.org/), [cvx](http://cvxr.com/cvx/) and [CORA](https://tumcps.github.io/CORA/) are required to run these codes.


**File Structure**
---

- The file `spacecraft_resilience.m` verifies the resilience of the spacecraft dynamics and compute its reachable set.
- The file `feedback_control.m` computes the parameters $\alpha$, $\beta$, $\gamma$, $\varepsilon$ and $\rho_{max}$ for a choice of matrices $K$, $P$ and $Q$ for the resilient feedback control Theorem 4.
- The file `main.m` computes the resilient trajectory tracking simulation.
- All the functions called by these aforementioned files are gathered in the folder `functions`.




**Running**
---




**Citation**
---
```
@article{bouvier2023tracking,  
  title = {Resilient Trajectory Tracking to Partial Loss of Control Authority over Actuators with Actuation Delay},   
  author = {Jean-Baptiste Bouvier, Himmat Panag, Robyn Woollands and Melkior Ornik},    
  journal = {Journal of Guidance Dynamics and Control},    
  year = {2023},   
  volume = {},  
  number = {},  
  pages = {}  
}
```

**Contributors**
---
- [Jean-Baptiste Bouvier](https://github.com/Jean-BaptisteBouvier)
- [Himmat Panag](https://www.linkedin.com/in/himmatpanag/?originalSubdomain=au)
- [Robyn Woollands](https://woollands.web.illinois.edu/index.html)
- [Melkior Ornik](https://mornik.web.illinois.edu/)



