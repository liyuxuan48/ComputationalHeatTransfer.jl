# ComputationalHeatTransfer.jl
tools for numerical simulation of conductive and convective heat transfer

[![Dev](https://img.shields.io/badge/docs-stable-blue.svg)](https://liyuxuan48.github.io/ComputationalHeatTransfer.jl)

This readme file assumes user already installed the **IJulia** package and got the jupyter notebook running.

# installation
To install this package, open julia REPL, and then run:

```julia
using Pkg
Pkg.add(url="https://github.com/liyuxuan48/ComputationalHeatTransfer.jl.git")
```

You may also need to install some other packages to make the notebook run smoothly, such as:

```julia
Pkg.add("Plots") # plotting functions
Pkg.add("Interact") # interactive features
Pkg.add("XLSX") # reading the experiment file in .xlsx format
Pkg.add("ProgressMeter") # having a progressbar when running the simulation
```

After installing the packages, we need to run

```julia
Pkg.add("Conda")
using Conda
Conda.pip_interop(true)
Conda.pip("install", "webio_jupyter_extension") # link the jupyter notebook with the Interact package
```

this will make the interactive features for **Interact** package in the postprocessing notebook avaliable.

# Setup
After installing the packages, we should perform one more step to migrate the example notebooks (hidden deeply among the julia source code files) to your working directory. I assume you already created a new directory.

```julia
using ComputationalHeatTransfer
cd("/path/of/working/directory") # navigate to your working directory
setup_examples(pwd())
```

This will create three folders in your working directory: examples, expdata, and numedata. Then you can enjoy the notebooks in examples. To open jupyter notebook, you can run:

```julia
using IJulia
notebook()
```

# Update

Because this package is still work in progress, I may make small modifications and fix some bugs. So I would suggest you to update the package frequently (before loading the package on your first run everyday:) )
```julia
using Pkg
Pkg.update("ComputationalHeatTransfer")
```

And then, open your jupyter notebook as usual

```julia
using IJulia
notebook()
```
