[![Dev](https://img.shields.io/badge/docs-stable-blue.svg)](https://liyuxuan48.github.io/ComputationalHeatTransfer.jl)

This readme file assumes user already installed the **IJulia** package and jupyter notebook.

# installation
To install this package, open julia REPL, and then run:

```julia
] add https://github.com/liyuxuan48/ComputationalHeatTransfer.jl.git
```

You may also need to install some other packages to make the notebook run smoothly, such as:

```julia
] add Plots
] add Interact
] add WebIO
] add ProgressMeter
```

After installing the packages, we need to run

```julia
# within a Julia REPL
using Conda
Conda.pip_interop(true)
Conda.pip("install", "webio_jupyter_extension")
```

this will make some interactive features in postprocessing notebook avaliable, but it is not the end of world if you could not make it work:)

# Setup
After installing the packages, we should perform one more step to migrate the example notebooks (hidden deeply among the julia source code files) to your working directory. Firstly, navigate to your working directory and open julia REPL there. And then run:

```julia
using ComputationalHeatTransfer
```

```julia
setup_examples(pwd())
```

This will create three folders in your working directory: examples, expdata, and numedata. Then you can enjoy the notebooks in examples once you figure out how to use jupyter notebook:)