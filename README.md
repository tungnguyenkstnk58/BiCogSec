# BiCogSec
- MATLAB simulation code for Bilateral Cognitive Security in Networked Control Systems under Stealthy Injection Attacks
- Paper: [https://arxiv.org/abs/2411.15319](https://arxiv.org/abs/2505.01232)

### Requirements
- MATLAB R2023 or later
- YALMIP toolbox (https://yalmip.github.io/tutorial/installation/)
- mosek solver (https://www.bing.com/search?q=mosek&qs=n&form=QBRE&sp=-1&lq=0&pq=mosek&sc=13-5&sk=&cvid=91562ECDC09B48479E2A648C52A8ED5C)

### How to Run
- Download ConnectedGraph.m, defender_chk.m, adversary_chk.m, and CogSec_01.m
- adversary_chk.m is for Lemma 1
- defender_chk.m is or Theorem 1
- To get an experiment, run CogSec_01.m, which will call adversary_chk.m and defender_chk.m in the computation.
