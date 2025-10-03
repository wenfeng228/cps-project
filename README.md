# Simulating Policy Priority Inference

We develop a streamlined, network-based variant of Policy Priority Inference (PPI) to study how budget allocations propagate across interdependent policy domains under structural spillovers. Our formulation isolates core computational mechanisms, diffusion on a directed weighted graph and heuristic allocation, while retaining interpretability and scalability for iterative simulation and sensitivity analysis.


## Context
Policy leaders must allocate scarce budgets across many goals with incomplete information and cross-policy spillovers. Estimating how a change in one domain affects others and doing so fast is difficult.


## Approach
We model policy domains (or SDG indicators) as nodes in a directed weighted graph. Positive edge weights represent synergies; negative weights represent trade-offs. At each iteration, a budget allocator prioritizes nodes with large target gaps and strong structural influence; updates then diffuse through the network. We evaluate convergence, target accuracy, and allocation efficiency, and compare outcomes with a no-spillover baseline.


## Our Contributions
- **Governance-light design:** Focus on structural spillovers and budget dynamics, omitting governance/behavioral submodules to keep the mechanism interpretable.  
- **Network-aware allocator:** A simple, transparent priority score couples target gaps with out-degree (or strength), yielding an efficient greedy-style allocation.  

