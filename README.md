This repository contains code, data, and a written report comparing three convex relaxations of the **Maximum Cut (Max-Cut)** problem on standard benchmark instances.
## Running `maxcut.jl`

`maxcut.jl` solves a single Max-Cut instance with one of the three relaxations and prints the resulting bound.

### Prerequisites

- Julia (≥ 1.9 recommended).
- Julia packages: `JuMP`, `Gurobi`, `SCS`, `LinearAlgebra` (the last is in the standard library).
- A working Gurobi license (free academic licenses are available at gurobi.com). `SCS` is open-source and installs without a license.

Install the Julia packages once:

```julia
using Pkg
Pkg.add(["JuMP", "Gurobi", "SCS"])
```

### Instance file format

A plain-text file. The first line contains `n m` (number of vertices and edges). Each of the following `m` lines contains an edge `i j w` (1-indexed endpoints and a weight). Both the `.rud` and `.mc` extensions used in `set1/` and `set2/` follow this format.

### Command

```
julia maxcut.jl <instance_file> <formulation>
```

where `<formulation>` is:

- `1` — SDP relaxation (`C^SDP`)
- `2` — McCormick + Triangle LP relaxation (`C^{MC+Tri}`)
- `3` — Combined SDP + McCormick + Triangle (`C^{MC+Tri+SDP}`)
