# Replication of the figures of the manuscript "Self-Justified Equilibria: Existence and Computation"

* Each Julia script in this folder allows to replicate the respective figure
in the paper. 

* To run the Julia files (fig_8.jl, fig_9a.jl, fig_9b.jl, fig_10.jl), make sure you have the packages "XLSX" and "Plots" installed. You can do this by launching Julia in a terminal, and then type the following commands:

```
julia> using Pkg

julia> Pkg.add("XLSX")
julia> Pkg.add("Plots")
```

* The data for the figures 1-7 is located and explained in the folder
[calibration_data](../calibration_data).

* Figures 8-21: The produced figures will be stored in the folder [figs](figs).


* The data for the figures 11-21 is being produced from the solutions of
the models. Replication code for the solutions can be found in the folder
[DEQN_for_IAMs](../DEQN_for_IAMs).

* To replicate the figures please make sure you are in the folder /figures_replication and do following:

**Fig.1**
Run the script fig_1.py. This file requires the two additional files ``ClimDICE.py`` and ``TestDefs.py``, which are also located in the working directory.

```
python fig_1.py
```
