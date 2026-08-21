using PrecompileTools: @compile_workload

@compile_workload begin
    get_pdesys_with_tags(["Diffusion"])
    get_pdesys_with_tags(["1D", "Linear"]; without = ["Discontinuous"])
    get_pdesys_with_tags(["Burgers", "Reaction"]; f = any)
end
