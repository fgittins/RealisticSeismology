# RealisticSeismology

This project presents an implementation of realistic, finite-temperature nuclear-matter models in neutron-star seismology. It was developed to support [Gittins and Andersson (Phys. Rev. D **111**, 083024, 2025)](https://doi.org/10.1103/PhysRevD.111.083024) and [Gittins _et al._ (Phys. Rev. D **111**, 023049, 2025)](https://doi.org/10.1103/PhysRevD.111.023049).

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
julia> using Weave

julia> weave("mode_demo.jmd")
```

## Citation

If you found this project to be useful in academic work, please cite it using the following references:

```bibtex
@article{gittins2024neutronstar,
        title="{Neutron-star seismology with realistic, finite-temperature nuclear matter}", 
       author={{Gittins}, F. and {Andersson}, N.},
      journal={Phys.\ Rev.\ D},
         year={2025},
        month=apr,
       volume={111},
        issue={8},
        pages={083024},
          doi={10.1103/PhysRevD.111.083024},
       eprint={2406.05177},
archivePrefix={arXiv},
 primaryClass={gr-qc}
}

@article{gittins2025problematicsystematics,
        title="{Problematic systematics in neutron-star merger simulations}", 
       author={{Gittins}, F. and {Matur}, R. and {Andersson}, N. and {Hawke}, I.},
      journal={Phys.\ Rev.\ D},
         year={2025},
        month=jan,
       volume={111},
        issue={2},
        pages={023049},
          doi={10.1103/PhysRevD.111.023049},
       eprint={2409.13468},
archivePrefix={arXiv},
 primaryClass={gr-qc}
}
```
