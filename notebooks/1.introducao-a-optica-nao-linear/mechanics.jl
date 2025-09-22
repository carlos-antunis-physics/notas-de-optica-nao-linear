module mechanics

    export harmonicOscillator, dampedHarmonicOscillator, drivenHarmonicOscillator, anharmonicOscillator, evolve;

    mutable struct oscillator_t
        t::Float64
        x::Float64
        v::Float64
        ω::Float64
        ζ::Float64
        drivingField::Function
        non_linearity::Function
    end

    harmonicOscillator(x, v, ω) = oscillator_t(0, x, v, ω, 0, _ -> 0.0, _ -> 0.0);
    dampedHarmonicOscillator(x, v, ω, ζ) = oscillator_t(0, x, v, ω, ζ, _ -> 0.0, _ -> 0.0);
    drivenHarmonicOscillator(x, v, ω, ζ, drivingField) = oscillator_t(0, x, v, ω, ζ, drivingField, _ -> 0.0);
    anharmonicOscillator(x, v, ω, ζ, drivingField, non_linearity) = oscillator_t(0, x, v, ω, ζ, drivingField, non_linearity);

    evolve(state::oscillator_t, dt::Float64) = begin
        kx = [0., 0., 0., 0.];
        kv = [0., 0., 0., 0.];
        ω² = state.ω * state.ω;

        _t = state.t;
        _x = state.x;
        _v = state.v;

        kx[1] = _v;
        kv[1] =
            state.drivingField(_t) -
            ω² * _x -
            state.ζ * _v -
            state.non_linearity(_x);

        _t = state.t + 0.5dt;
        _x = state.x + 0.5dt * kx[1];
        _v = state.v + 0.5dt * kv[1];

        kx[2] = _v;
        kv[2] =
            state.drivingField(_t) -
            ω² * _x -
            state.ζ * _v -
            state.non_linearity(_x);

        _t = state.t + 0.5dt;
        _x = state.x + 0.5dt * kx[2];
        _v = state.v + 0.5dt * kv[2];

        kx[3] = _v;
        kv[3] =
            state.drivingField(_t) -
            ω² * _x -
            state.ζ * _v -
            state.non_linearity(_x);

        _t = state.t + dt;
        _x = state.x + dt * kx[3];
        _v = state.v + dt * kv[3];

        kx[4] = _v;
        kv[4] =
            state.drivingField(_t) -
            ω² * _x -
            state.ζ * _v -
            state.non_linearity(_x);

        state.t += dt;
        state.x += (dt / 6.) * (kx[1]+ 2kx[2] + 2kx[3] + kx[4]);
        state.v += (dt / 6.) * (kv[1]+ 2kv[2] + 2kv[3] + kv[4]);
    end

end