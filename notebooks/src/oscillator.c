typedef double drivingForce_t(double);
typedef double nonLinearity_t(double);

typedef struct oscillator_t {
    double t, x, v;
    double omega;
    double drag;
    drivingForce_t *impulse;
    nonLinearity_t *non_linearity;
} oscillator_t;

double __none__(double __u__) {
    return 0.0;
}

static inline oscillator_t harmonicOscillator(
    double x, double v,
    double omega
) {
    return (oscillator_t) {
        .t = 0.,
        .x = x, .v = v,
        .omega = omega,
        .drag = 0.0,
        .impulse = &__none__,
        .non_linearity = &__none__
    };
}
static inline oscillator_t dampedHarmonicOscillator(
    double x, double v,
    double omega, double drag
) {
    return (oscillator_t) {
        .t = 0.,
        .x = x, .v = v,
        .omega = omega,
        .drag = drag,
        .impulse = &__none__,
        .non_linearity = &__none__
    };
}
static inline oscillator_t drivenHarmonicOscillator(
    double x, double v,
    double omega, double drag,
    drivingForce_t *impulse
) {
    return (oscillator_t) {
        .t = 0.,
        .x = x, .v = v,
        .omega = omega,
        .drag = drag,
        .impulse = impulse,
        .non_linearity = &__none__
    };
}

static inline oscillator_t anharmonicOscillator(
    double x, double v,
    double omega,
    nonLinearity_t *non_linearity
) {
    return (oscillator_t) {
        .t = 0.,
        .x = x, .v = v,
        .omega = omega,
        .drag = 0.0,
        .impulse = &__none__,
        .non_linearity = non_linearity
    };
}
static inline oscillator_t dampedAnharmonicOscillator(
    double x, double v,
    double omega, double drag,
    nonLinearity_t *non_linearity
) {
    return (oscillator_t) {
        .t = 0.,
        .x = x, .v = v,
        .omega = omega,
        .drag = drag,
        .impulse = &__none__,
        .non_linearity = non_linearity
    };
}
static inline oscillator_t drivenAnharmonicOscillator(
    double x, double v,
    double omega, double drag,
    drivingForce_t *impulse,
    nonLinearity_t *non_linearity
) {
    return (oscillator_t) {
        .t = 0.,
        .x = x, .v = v,
        .omega = omega,
        .drag = drag,
        .impulse = impulse,
        .non_linearity = non_linearity
    };
}

static inline oscillator_t evolve(
    oscillator_t state,
    double dt
) {
    double kx[4], kv[4];
    double _t, _x, _v;

    _t = state.t;
    _x = state.x;
    _v = state.v;

    kx[0] = _v;
    kv[0] = 
        state.impulse(_t) -
        state.omega * state.omega * _x -
        state.drag * _v -
        state.non_linearity(_x);
    
    _t = state.t + 0.5 * dt;
    _x = state.x + 0.5 * dt * kx[0];
    _v = state.v + 0.5 * dt * kv[0];

    kx[1] = _v;
    kv[1] =
        state.impulse(_t) -
        state.omega * state.omega * _x -
        state.drag * _v -
        state.non_linearity(_x);

    _t = state.t + 0.5 * dt;
    _x = state.x + 0.5 * dt * kx[1];
    _v = state.v + 0.5 * dt * kv[1];

    kx[2] = _v;
    kv[2] =
        state.impulse(_t) -
        state.omega * state.omega * _x -
        state.drag * _v -
        state.non_linearity(_x);

    _t = state.t + 0.5 * dt;
    _x = state.x + 0.5 * dt * kx[2];
    _v = state.v + 0.5 * dt * kv[2];

    kx[3] = _v;
    kv[3] =
        state.impulse(_t) -
        state.omega * state.omega * _x -
        state.drag * _v -
        state.non_linearity(_x);

    state.t += dt;
    state.x += (dt / 6.) * (
        kx[0] + 2 * kx[1] + 2 * kx[2] + kx[3]
    );
    state.v += (dt / 6.) * (
        kv[0] + 2 * kv[1] + 2 * kv[2] + kv[3]
    );

    return state;
}