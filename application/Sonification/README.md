# Sonification — listen to a trace

Turns any time series (a flexOPT synthetic, a recorded seismogram, a
normal-mode spectrum…) into a WAV file, using `src/sonification.jl` (ported
from SonifSismo.jl). It is an *output*: any application can hand its traces to it.

```sh
julia --project=. application/Sonification/run.jl                                   # illustrative Earth "hum" (6 normal modes, 4 days)
julia --project=. application/Sonification/run.jl days=2 acceleration=50000
julia --project=. application/Sonification/run.jl file=mytrace.csv dt=0.01 acceleration=100
```

`acceleration` compresses time and raises the pitch by the same factor:
a 0.8 mHz mode played ×10⁵ is heard at 80 Hz, and 4 days last about 3.5 s.
The default signal is a sum of damped cosines at the frequencies of a few
fundamental spheroidal modes; it is a listening example, not a computed synthetic.

Outputs: `outputs/<name>.wav` and `outputs/<name>.png` (waveform).
If `WAV` is missing: `julia --project=. -e 'using Pkg; Pkg.add("WAV")'`.
The full demo with real synthetics is `notebooks/sonificationDemo.ipynb`.
