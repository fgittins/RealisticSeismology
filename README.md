# RealisticSeismology

This project presents an implementation of realistic, finite-temperature nuclear-matter models in neutron-star seismology.

## Installation

The software is developed using the Julia programming language. To use it:

1. Install [Julia](https://julialang.org/downloads/)

2. Download this repository

3. Run Julia in the repository directory

4. Type `]` to enter Julia's package manager (Pkg.jl) REPL,

```julia-repl
(@v1.10) pkg>
```

5. `activate` the project environment with

```julia-repl
(@v1.10) pkg> activate .
```

6. `instantiate` the project,

```julia-repl
(RealisticSeismology) pkg> instantiate
```

For more information on Julia packages and environments, see the [Pkg.jl documentation](https://pkgdocs.julialang.org/v1/).

## Getting started

General use of this software is demonstrated in the `scripts` and `notebooks` directories. The notebooks are written in Julia Markdown and may be compiled using [Weave.jl](https://weavejl.mpastell.com/stable/). For example,

```julia-repl
using Weave
weave("mode_demo.jmd")
```
