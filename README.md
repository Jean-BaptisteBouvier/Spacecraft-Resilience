# Spacecraft-Resilience

This project gathers the MATLAB codes used in the paper "Resilient Trajectory Tracking to Partial Loss of Control Authority over Actuators with Actuation Delay" by Jean-Baptiste Bouvier, Himmat Panag, Robyn Woollands and Melkior Ornik.
This work is not published and not available yet.


**Prerequisite**
---
The MATLAB toolboxes [mpt](https://www.mpt3.org/), [cvx](http://cvxr.com/cvx/) and [CORA](https://tumcps.github.io/CORA/) are required to run these codes.


**File Structure**
---

1- The file `spacecraft_resilience.m` verifies the resilience of the spacecraft dynamics and compute its reachable set.
2- The file `feedback_control.m` computes the parameters $\varepsilon$, `rho_max` depending on matrices $K$, $P$ and $Q$ for the resilient feedback control Theorem 4.
3-
4- The file `main.m` computes the resilient trajectory tracking.




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



