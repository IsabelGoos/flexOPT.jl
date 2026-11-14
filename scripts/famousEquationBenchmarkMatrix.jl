module FamousEquationBenchmarkMatrix

export SCHEMES, EQUATIONS, MATERIAL_SCENARIOS, material_scenarios
export space_sizes, benchmark_schema_version

const benchmark_schema_version = "famous_equations_hierarchical_taylor_v2"

"Symmetric interpolation geometry around the physical centre of a p-point stencil."
function interpolation_geometry(points::Int, mode::Symbol)
    centre = (points - 1) / 2
    if mode === :center
        return (field_points=1, field_offset=centre, field_order=-1,
                material_points=1, material_offset=centre, material_order=-1)
    end

    # Odd nodal stencils have an integer centre and their staggered grid has
    # p-1 half-grid points.  Even stencils have a half-grid centre and the
    # opposite-parity samples are the p nodal points.
    same_points = isodd(points) ? points : points - 1
    same_offset = isodd(points) ? 0.0 : 0.5
    opposite_points = isodd(points) ? points - 1 : points
    opposite_offset = isodd(points) ? 0.5 : 0.0
    if mode === :material_staggered
        return (field_points=1, field_offset=centre, field_order=-1,
                material_points=opposite_points,
                material_offset=opposite_offset, material_order=1)
    elseif mode === :multipoint
        return (field_points=same_points, field_offset=same_offset, field_order=1,
                material_points=same_points, material_offset=same_offset,
                material_order=1)
    end
    throw(ArgumentError("unknown interpolation geometry: $mode"))
end

function scheme(name, family, points; geometry=:center,
                supplementary_order=family in (:explicitFD, :convFD) ? 0 : 2,
                role=family in (:explicitFD, :convFD) ? :reference : :candidate)
    itpl = interpolation_geometry(points, geometry)
    return (; name, family, points_space=points, points_time=3,
        order_b_space=family in (:explicitFD, :convFD) ? -1 : 1,
        order_b_time=family in (:explicitFD, :convFD) ? -1 : 1,
        supplementary_order,
        interpolation_geometry=geometry,
        taylor_inverse_mode=:hierarchical_constrained,
        itpl..., recommended_role=role)
end

const SCHEMES = [
    scheme("explicitFD", :explicitFD, 3),
    scheme("convFD", :convFD, 3),
    scheme("OPT3", :OPT, 3; geometry=:center, role=:robust_candidate),
    # For an even stencil, :center places the single Taylor centre at the
    # geometric midpoint (offset 1.5 for four points).
    scheme("OPT4", :OPT, 4; geometry=:center),
    scheme("OPT5", :OPT, 5; geometry=:center,
        role=:high_accuracy_candidate),
]

const EQUATIONS = [
    (
        label="Poisson 1-D",
        equation="1DpoissonHetero",
        physics=:poisson,
        space_dimension=1,
        has_time=false,
        fields=1,
        branches=(:elliptic,),
        material_variables=(:kappa,),
    ),
    (
        label="Poisson 2-D",
        equation="2DpoissonHetero",
        physics=:poisson,
        space_dimension=2,
        has_time=false,
        fields=1,
        branches=(:elliptic,),
        material_variables=(:kappa,),
    ),
    (
        label="SH frequency 1-D",
        equation="1DsismoFreqHetero",
        physics=:sh_frequency,
        space_dimension=1,
        has_time=false,
        fields=1,
        branches=(:SH,),
        material_variables=(:rho, :mu),
    ),
    (
        label="SH time 1-D",
        equation="1DsismoTime",
        physics=:sh_time,
        space_dimension=1,
        has_time=true,
        fields=1,
        branches=(:SH,),
        material_variables=(:rho, :mu),
    ),
    (
        label="Acoustic time 2-D",
        equation="2DacousticTime",
        physics=:acoustic,
        space_dimension=2,
        has_time=true,
        fields=1,
        branches=(:P,),
        material_variables=(:velocity,),
    ),
    (
        label="Elastic time 2-D",
        equation="2DsismoTimeIsoHeteroForce",
        physics=:elastic,
        space_dimension=2,
        has_time=true,
        fields=2,
        branches=(:P, :S),
        material_variables=(:rho, :lambda, :mu),
    ),
    (
        label="Elastic time 3-D",
        equation="3DsismoTimeIsoHeteroForce",
        physics=:elastic,
        space_dimension=3,
        has_time=true,
        fields=3,
        branches=(:P, :S1, :S2),
        material_variables=(:rho, :lambda, :mu),
    ),
]

# Every PDE receives the subset relevant to its declared material variables.
# Frequencies are integer vectors on the 2π-periodic domain.  The second
# component is ignored for one-dimensional equations.
const MATERIAL_SCENARIOS = [
    (
        name="homogeneous",
        active=(),
        material_wave=(0, 0, 0),
        phase=0.0,
        amplitude_fraction=0.0,
        relation=:homogeneous,
    ),
    (
        name="same_wave_phase0",
        active=:all,
        material_wave=(2, 1, 1),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:same_wavelength,
    ),
    (
        name="same_wave_phase_pi2",
        active=:all,
        material_wave=(2, 1, 1),
        phase=pi / 2,
        amplitude_fraction=0.15,
        relation=:same_wavelength,
    ),
    (
        name="long_material_phase0",
        active=:all,
        material_wave=(1, 1, 1),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:material_longer,
    ),
    (
        name="short_material_phase0",
        active=:all,
        material_wave=(6, 3, 3),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:material_three_times_shorter,
    ),
    (
        name="short_material_phase_pi2",
        active=:all,
        material_wave=(6, 3, 3),
        phase=pi / 2,
        amplitude_fraction=0.15,
        relation=:material_three_times_shorter,
    ),
    (
        name="density_only_phase0",
        active=(:rho,),
        material_wave=(2, 1, 1),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:single_variable,
    ),
    (
        name="stiffness_only_phase0",
        active=(:kappa, :mu, :lambda),
        material_wave=(2, 1, 1),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:single_variable_group,
    ),
    (
        name="rho_mu_opposite_phase",
        active=(:rho, :mu),
        material_wave=(2, 1, 1),
        phase=(rho=0.0, mu=pi),
        amplitude_fraction=0.15,
        relation=:opposite_phase,
    ),
    (
        name="lambda_mu_quadrature",
        active=(:lambda, :mu),
        material_wave=(2, 1, 1),
        phase=(lambda=0.0, mu=pi / 2),
        amplitude_fraction=0.15,
        relation=:quadrature,
    ),
    (
        name="constant_varying_impedance",
        active=(:rho, :mu),
        material_wave=(2, 1, 1),
        phase=(rho=0.0, mu=0.0),
        amplitude_fraction=(rho=0.15, mu=0.15),
        relation=:approximately_constant_wave_speed,
    ),
]

function material_scenarios(equation)
    variables = Set(equation.material_variables)
    return filter(MATERIAL_SCENARIOS) do scenario
        scenario.active === :all && return true
        isempty(scenario.active) && return true
        any(variable -> variable in variables, scenario.active)
    end
end

space_sizes(equation) = equation.space_dimension == 1 ?
    [16, 24, 32, 48, 64, 96, 128, 192] :
    equation.space_dimension == 2 ? [8, 12, 16, 24, 32, 48, 64] :
    [6, 8, 10, 12, 16, 20, 24]

end
