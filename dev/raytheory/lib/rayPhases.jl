# dev/raytheory/lib/rayPhases.jl
#
# A phase is NOT a separate equation: it is one Hamiltonian branch (P or S)
# plus an ordered list of interface events. Everything that happens between
# two events is a free characteristic of H_mode; at every discontinuity that
# is not the next event the ray is transmitted in the same mode (except at a
# fluid-solid boundary, which needs an explicit event, as in TauP: "P" never
# enters the core, "PKP" does).
#
#     P, S, PP, SS, PS, SP, pP, sP, sS, PcP, ScS, PcS, PKP, SKS, PKKP, PKiKP,
#     PKIKP, PmP, Pms, P410s, S410p, Pv410P, P^410P (underside), and for
#     Cartesian interfaces with names: Pv{i1}P, Pv{i1}S, P{i1}s ...
#
# Conventions (TauP-like):
#   upper-case leg  : the leg starts by going DOWN (first leg) / after an event
#   lower-case p, s : up-going leg from the source (depth phases) or an
#                     up-going converted leg (P410s)
#   c / i           : topside reflection at the CMB / ICB
#   K, I (J)        : P in the outer core, P (S) in the inner core
#   m, <number>, {name}
#                   : the discontinuity; bare between X and lower-case y:
#                     up-going transmission/conversion (P410s, Pms); bare
#                     between X and upper-case Y: down-going transmission
#                     (P410S), except "m": topside reflection (PmP, SmS)
#   v<d>, ^<d>      : topside / underside reflection at d (Pv410P, P^410P)

"""
    PhaseEvent(boundary, action, mode, direction)

At the next crossing of `boundary` (a name, see boundaryMatches) by a ray
travelling in `direction` (:down, :up, :any), do `action` (:transmit or
:reflect) and continue in `mode` (:P or :S).
"""
struct PhaseEvent
    boundary::Symbol
    action::Symbol
    mode::Symbol
    direction::Symbol
end

struct Phase
    name::String
    mode::Symbol              # initial mode
    startDirection::Symbol    # :down, :up, :any
    events::Vector{PhaseEvent}
end
Phase(name, mode, dir, ev::PhaseEvent...) = Phase(String(name), mode, dir, collect(ev))

Base.show(io::IO, ph::Phase) = print(io, "Phase(\"", ph.name, "\": start ", ph.mode, " ", ph.startDirection,
    isempty(ph.events) ? "" : ", " * join(["$(e.boundary):$(e.action)→$(e.mode)($(e.direction))" for e in ph.events], ", "), ")")

legMode(c::Char) = c in ('P', 'p', 'K', 'I') ? :P : c in ('S', 's', 'J') ? :S : error("not a leg: $c")
isLeg(c::Char) = c in ('P', 'S', 'p', 's', 'K', 'I', 'J')

"""
    parsePhase(name) -> Phase

Small TauP-like parser (see the header of rayPhases.jl). Throws on what it
does not understand rather than guessing.
"""
function parsePhase(name::AbstractString)
    s = collect(name)
    legs = Char[]
    markers = String[]   # marker between legs[i] and legs[i+1]
    i = 1
    cur = ""
    while i <= length(s)
        c = s[i]
        if isLeg(c)
            if !isempty(legs) || !isempty(cur)
                isempty(legs) && error("phase $name: cannot start with a marker")
                push!(markers, cur)
            end
            push!(legs, c); cur = ""; i += 1
        elseif c == '{'
            j = findnext(==('}'), s, i)
            j === nothing && error("phase $name: unclosed {")
            cur *= String(s[i:j]); i = j + 1
        else
            cur *= string(c); i += 1
        end
    end
    isempty(cur) || error("phase $name: trailing marker '$cur'")
    isempty(legs) && error("phase $name: no leg")
    ev = PhaseEvent[]
    for k in 1:length(legs)-1
        a, b, m = legs[k], legs[k+1], markers[k]
        mb = legMode(b)
        bname(x) = startswith(x, "{") ? Symbol(x[2:end-1]) : x == "m" ? :moho : Symbol("d", x)
        if m == ""
            if a in ('P', 'S', 'p', 's') && b in ('P', 'S', 'p', 's')
                push!(ev, PhaseEvent(:surface, :reflect, mb, :up))
            elseif a in ('P', 'S', 'p', 's') && b == 'K'
                push!(ev, PhaseEvent(:cmb, :transmit, :P, :down))
            elseif a == 'K' && b in ('P', 'S', 'p', 's')
                push!(ev, PhaseEvent(:cmb, :transmit, mb, :up))
            elseif a == 'K' && b == 'K'
                push!(ev, PhaseEvent(:cmb, :reflect, :P, :up))
            elseif a == 'K' && b in ('I', 'J')
                push!(ev, PhaseEvent(:icb, :transmit, mb, :down))
            elseif a in ('I', 'J') && b == 'K'
                push!(ev, PhaseEvent(:icb, :transmit, :P, :up))
            elseif a in ('I', 'J') && b in ('I', 'J')
                error("phase $name: inner-core underside reflections not supported yet")
            else
                error("phase $name: cannot connect $a and $b")
            end
        elseif m == "c"
            push!(ev, PhaseEvent(:cmb, :reflect, mb, :down))
        elseif m == "i"
            push!(ev, PhaseEvent(:icb, :reflect, :P, :down))
        elseif startswith(m, "v")
            push!(ev, PhaseEvent(bname(m[2:end]), :reflect, mb, :down))
        elseif startswith(m, "^")
            push!(ev, PhaseEvent(bname(m[2:end]), :reflect, mb, :up))
        elseif m == "m" && isuppercase(b)
            push!(ev, PhaseEvent(:moho, :reflect, mb, :down))
        else
            push!(ev, PhaseEvent(bname(m), :transmit, mb, islowercase(b) ? :up : :down))
        end
    end
    first = legs[1]
    first in ('P', 'S', 'p', 's') || error("phase $name must start with P, S, p or s")
    startDir = islowercase(first) ? :up : (isempty(ev) ? :any : :down)
    return Phase(String(name), legMode(first), startDir, ev)
end

"Phases given as names or Phase objects"
asPhase(p::Phase) = p
asPhase(p::AbstractString) = parsePhase(p)
