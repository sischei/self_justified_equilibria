# Self-Justified Equilibria: Existence and Computation

<p align="center">
<img src="screens/SJE.png" width="800px"/>
</p>

## Description of programs and datasets used

### Organization of the repository

This Julia-based code repository supplements the work of [Felix Kuebler](https://sites.google.com/site/fkubler/) and [Simon Scheidegger](https://sites.google.com/site/simonscheidegger), titled _[Self-Justified Equilibria: Existence and Computation](#citation)_ (Kubler and Scheidegger; 2025).

* This repository contains three distinct folders:
  1. ["Replication Codes"](code): Replication codes for *Section  - Application of SJE: A High-Dimensional Stochastic Production Economy*, where non-stationary integrated assessment models (IAMs) are solved by adopting ["Deep Equilibrium Nets (DEQN)"](https://onlinelibrary.wiley.com/doi/epdf/10.1111/iere.12575) to the context of climate economic models. Notice that the codes provided here complement the Online Appendix D of our article, where the formal underpinnings of the code are outlined.
      - How to compute the individual model calibrations is detailed in the various readmes that provided under the following three links: [optimal results](DEQN_for_IAMs/dice_generic/README.md), [business as usual results](DEQN_for_IAMs/gdice_baseline/README.md), and [results from the appendix](DEQN_for_IAMs/dice_generic_FEX).
      - The content and usage of the generic Deep Equilibrium Nets for Integrated Assessment Models framework are outlined in the corresponding [README](DEQN_for_IAMs/README.md).
      
  2. ["Replication of Figures"](figures_replication): Replication routines for plotting all the figures that are presented in the paper.

  
### Replication of the numerical results

* To replicate the results of the article step-by-step, a detailed set of instructions is provided _[here](#Replication)_.
  

### Datasets used

* The detailed references to the datasets used in our work are provided _[here](#Datasets)_.
    
    
## Computational requirements

### Software requirements

* We provide implementations that use Julia 1.9.

* For the  The basic dependencies are [Tensorflow==2.x](https://www.tensorflow.org/), [hydra](https://hydra.cc/) and [Tensorboard](https://www.tensorflow.org/tensorboard) (for monitoring).

* The file ``Project.toml `` lists the detailed dependencies for Julia. See [here](https://pkgdocs.julialang.org/v1/toml-files/) for detailed instructions how on how to set up and use the ``Project.toml`` file in a Julia project.


### Controlled randomness Project.toml

The random seed for our computations in *Section 6 - The social cost of carbon and optimal abatement in the DICE economy* is set at ``Climate_in_Climate_Economics/DEQN_for_IAMS/config/config.yaml``, line 10.


### Memory and runtime requirements

* To solve one overlapping generations model solved with SJE as discussed in *Section 6 - The social cost of carbon and optimal abatement in the DICE economy* until full convergence, it requires about 15 min on an ordinary laptop. All those models presented in the paper were solved using our [DEQN library](DEQN_for_IAMs), which we ran on an 8-core Intel compute node on [https://nuvolos.cloud](https://nuvolos.cloud) with 64GB of RAM, and 100Gb of fast local storage (SSD).

* All the postprocessing codes (to produce the summary statistics, plots, and so forth) were run on an ordinary 4-core Intel-based laptop with Ubuntu version 18.04.5 LTS and consume typically few seconds to run.

* The approximate time needed to reproduce all the analyses for this paper on a standard (year 2023) desktop machine is 1-3 days of human time.

## Authors

* [Felix Kuebler](https://sites.google.com/site/fkubler/) (the University of Zuerich, Department for Banking and Finance, and Swiss Finance Institute)
* [Simon Scheidegger](https://sites.google.com/site/simonscheidegger) (the University of Lausanne, Department of Economics)


## Citation

Please cite [Self-Justified Equilibria: Existence and Computation](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3494876) in your publications if it helps your research:

```
@article{kubler2019self,
  title={Self-justified equilibria: Existence and computation},
  author={Kubler, Felix and Scheidegger, Simon},
  journal={Available at SSRN 3494876},
  year={2019}
}

```


## Support

This work was generously supported by grants from the [Swiss National Science Foundation](https://www.snf.ch) under project IDs "Can economic policy mitigate climate change," "New methods for asset pricing with frictions,” and the [Enterprise for Society (E4S)](https://e4s.center).

