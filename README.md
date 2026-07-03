# Liftover.jl

[![CI](https://github.com/mashu/Liftover.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/mashu/Liftover.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/mashu/Liftover.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mashu/Liftover.jl)
[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://mashu.github.io/Liftover.jl/stable/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Julia liftover for point coordinates between genome assemblies, using UCSC chain files.
Motivated by [jeremymcrae/liftover](https://github.com/jeremymcrae/liftover); a Julia rewrite, not a direct port of that codebase.

```julia
using Pkg; Pkg.add(url="https://github.com/mashu/Liftover.jl")
using Liftover

converter = get_lifter("hg19", "hg38"; one_based=true)
converter["1"][103786442]
```
