module FamousEquationBenchmarkMatrix

export SCHEMES, EQUATIONS, MATERIAL_SCENARIOS, material_scenarios
export space_sizes, benchmark_schema_version

const benchmark_schema_version = "famous_equations_fd_opt_v1"

const SCHEMES = [
    (
        name="FD3",
        family=:FD,
        points_space=3,
        points_time=3,
        order_b_space=-1,
        order_b_time=-1,
        supplementary_order=0,
        interpolation_order=-1,
        recommended_role=:reference,
    ),
    (
        name="OPT3",
        family=:OPT,
        points_space=3,
        points_time=3,
        order_b_space=1,
        order_b_time=1,
        supplementary_order=2,
        interpolation_order=-1,
        recommended_role=:robust_candidate,
    ),
    (
        name="FD5space-FD3time",
        family=:FD,
        points_space=5,
        points_time=3,
        order_b_space=-1,
        order_b_time=-1,
        supplementary_order=0,
        interpolation_order=-1,
        recommended_role=:reference,
    ),
    (
        name="OPT5space-OPT3time-hat-supp0",
        family=:OPT,
        points_space=5,
        points_time=3,
        order_b_space=1,
        order_b_time=1,
        supplementary_order=0,
        interpolation_order=-1,
        recommended_role=:high_accuracy_candidate,
    ),
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
]

# Every PDE receives the subset relevant to its declared material variables.
# Frequencies are integer vectors on the 2π-periodic domain.  The second
# component is ignored for one-dimensional equations.
const MATERIAL_SCENARIOS = [
    (
        name="homogeneous",
        active=(),
        material_wave=(0, 0),
        phase=0.0,
        amplitude_fraction=0.0,
        relation=:homogeneous,
    ),
    (
        name="same_wave_phase0",
        active=:all,
        material_wave=(2, 1),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:same_wavelength,
    ),
    (
        name="same_wave_phase_pi2",
        active=:all,
        material_wave=(2, 1),
        phase=pi / 2,
        amplitude_fraction=0.15,
        relation=:same_wavelength,
    ),
    (
        name="short_material_phase0",
        active=:all,
        material_wave=(6, 3),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:material_three_times_shorter,
    ),
    (
        name="short_material_phase_pi2",
        active=:all,
        material_wave=(6, 3),
        phase=pi / 2,
        amplitude_fraction=0.15,
        relation=:material_three_times_shorter,
    ),
    (
        name="density_only_phase0",
        active=(:rho,),
        material_wave=(2, 1),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:single_variable,
    ),
    (
        name="stiffness_only_phase0",
        active=(:kappa, :mu, :lambda),
        material_wave=(2, 1),
        phase=0.0,
        amplitude_fraction=0.15,
        relation=:single_variable_group,
    ),
    (
        name="rho_mu_opposite_phase",
        active=(:rho, :mu),
        material_wave=(2, 1),
        phase=(rho=0.0, mu=pi),
        amplitude_fraction=0.15,
        relation=:opposite_phase,
    ),
    (
        name="lambda_mu_quadrature",
        active=(:lambda, :mu),
        material_wave=(2, 1),
        phase=(lambda=0.0, mu=pi / 2),
        amplitude_fraction=0.15,
        relation=:quadrature,
    ),
    (
        name="constant_varying_impedance",
        active=(:rho, :mu),
        material_wave=(2, 1),
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

space_sizes(equation) =
    equation.space_dimension == 1 ?
    [16, 24, 32, 48, 64, 96, 128, 192] :
    [8, 12, 16, 24, 32, 48, 64]

end
