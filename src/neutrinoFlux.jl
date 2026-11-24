module neutrinoFlux
    #  (c) Yael Deniz and Isabel Goos 2026

    using CSV,DataFrames,Interpolations
    using Downloads,CodecZlib,DelimitedFiles
    using CairoMakie
    using ..commonBatchs: myProduceOrLoad
    include("neutrinoFlux/config.jl")
    include("neutrinoFlux/generate_neutrino_flux.jl")
    export read_neutrino_flux_table, produce_neutrino_flux
    export produce_or_load_neutrino_flux
    export prefetch_daemonflux_data!
    export plot_neutrino_flux_spectra, plot_neutrino_flux_heatmap
    export NeutrinoFluxParameters, NeutrinoFluxTable
    export configure_neutrino_flux, set_default_neutrino_flux!
    export DEFAULT_NUFLUX_MODEL, DEFAULT_NUFLUX_PARAMS

end
