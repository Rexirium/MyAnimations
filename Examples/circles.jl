using LinearAlgebra
using CairoMakie

let 
    R, r = 2.0, 0.5
    Ω = 1.0
    ω = (R - r) / r * Ω
    T = 2π / Ω
    tf = 4T

    nsteps = 1000
    time_len = 20
    fps = 25

    ts = range(0, tf, nsteps+1)

    xs = (R - r) * cos.(Ω .* ts) + r * cos.(ω .* ts)
    ys = (R - r) * sin.(Ω .* ts) - r * sin.(ω .* ts)

    xc = (R - r) * cos.(Ω .* ts)
    yc = (R - r) * sin.(Ω .* ts)

    nframes = fps * time_len
    spf = nsteps ÷ nframes

    set_theme!()

end