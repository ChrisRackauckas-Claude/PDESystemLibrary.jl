module PDESystemLibrary
using ModelingToolkit, DomainSets, FunctionMaps
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK
using Interpolations

import SciMLBase

using Markdown
using Random

Random.seed!(100)

all_systems = []

include("../lib/burgers.jl")
include("../lib/linear_diffusion.jl")
include("../lib/linear_convection.jl")
include("../lib/nonlinear_diffusion.jl")
include("../lib/general_linear_system.jl")
include("../lib/brusselator.jl")

"""
    get_pdesys_with_tags(withtags; without = [], f = all)

Return the PDE systems whose `metadata` tags match `withtags`.

`withtags` is an iterable of tags to include. By default all requested tags must be
present. Pass `f = any` to select systems that contain at least one requested tag.
Use `without` to exclude systems containing any of those tags.

# Examples

```julia
diffusion_systems = get_pdesys_with_tags(["Diffusion"])
one_dimensional_or_heat = get_pdesys_with_tags(["1D", "Heat"]; f = any)
smooth_diffusion = get_pdesys_with_tags(["Diffusion"]; without = ["Discontinuous"])
```
"""
function get_pdesys_with_tags(withtags; without = [], f = all)
    return filter(all_systems) do ex
        b = f(t -> t ∈ ex.metadata, withtags)
        b && all(t -> t ∉ ex.metadata, without)
    end
end

export get_pdesys_with_tags
end # module PDESystemLibrary
