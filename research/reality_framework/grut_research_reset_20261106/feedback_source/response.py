"""Retarded and finite-grid linear response for one Majorana mass phase."""
import numpy as np
from scipy.fft import irfft, next_fast_len, rfft, rfftfreq
from scipy.integrate import quad
from scipy.signal import fftconvolve

COEFFICIENT = 1/(32*np.pi**2)
ENVELOPE_WEIGHTS = np.array([1/16, 1/4, 3/8, 1/4, 1/16])


def carriers(pulse):
    return pulse.frequency + np.arange(-2, 3)*np.pi/pulse.duration


def profile(times, pulse):
    """Amplitude, speed and acceleration of the compact cosine pulse."""
    frequencies = carriers(pulse)
    arguments = np.asarray(times)[..., None] * frequencies
    cosine, sine = np.cos(arguments), np.sin(arguments)
    support = np.abs(times) < pulse.duration
    weights = pulse.amplitude*ENVELOPE_WEIGHTS
    return (support * (cosine @ weights),
            -support.astype(float) * (sine @ (frequencies*weights)),
            -support.astype(float) * (cosine @ (frequencies**2*weights)))


def transform(nu, pulse):
    frequency = carriers(pulse)
    width = pulse.duration
    near = np.sinc((nu-frequency)*width/np.pi)
    far = np.sinc((nu+frequency)*width/np.pi)
    return pulse.amplitude*width*((near+far) @ ENVELOPE_WEIGHTS)


def kernel(s):
    """K_R(s) at positive frequency, with K_R(0) and K_R'(0) subtracted."""
    s = np.asarray(s, dtype=float)
    log_integral = np.empty(s.shape, dtype=complex)
    near = np.abs(s) < .01
    z = s[near]
    log_integral[near] = -z*(1/6+z*(1/60+z*(1/420+z*(1/2520+z/13860))))

    spacelike = s <= -.01
    beta = np.sqrt(1-4/s[spacelike])
    log_integral[spacelike] = -2+2*beta*np.arctanh(1/beta)

    below = (s >= .01) & (s < 4)
    gamma = np.sqrt(4/s[below]-1)
    log_integral[below] = -2+2*gamma*np.arctan(1/gamma)

    above = s >= 4
    beta = np.sqrt(1-4/s[above])
    log_integral[above] = -2+2*beta*np.arctanh(beta)-1j*np.pi*beta
    return COEFFICIENT*s*log_integral


def parameter_integral(s):
    """Independent branch of the Feynman parameter representation."""
    def real_integrand(x):
        return np.log(abs(1-s*x*(1-x)))

    if s > 4:
        width = np.sqrt(1-4/s)
        zeros = [(1-width)/2, (1+width)/2]
    else:
        width, zeros = 0., None
    real = quad(real_integrand, 0, 1, points=zeros, epsabs=1e-12, epsrel=1e-12)[0]
    return COEFFICIENT*s*(real-1j*np.pi*width)


def static_reference(wave):
    s = -wave**2

    def integrand(x):
        z = x*(1-x)
        return np.log1p(-s*z)-s*z/(1-s*z)

    return float(kernel(s).real), COEFFICIENT*quad(integrand, 0, 1)[0]


def continuum_force(amplitude, dt, wave, padding=512):
    """Pad a compact input before evaluating the causal kernel by FFT."""
    length = next_fast_len(len(amplitude)+2*round(padding/dt))
    start = (length-len(amplitude))//2
    stop = start+len(amplitude)
    padded = np.zeros(length)
    padded[start:stop] = amplitude
    omega = 2*np.pi*rfftfreq(length, dt)
    # Numpy's positive-frequency inverse has exp(+i*omega*t).
    signal = irfft(rfft(padded)*kernel(omega**2-wave**2).conj(), n=length)
    probe = min(round(4/dt), start)
    return signal[start:stop], signal[start-probe:start]


def spectral_work(pulse, wave):
    """Positive-frequency pair work, evaluated without a time transform."""
    threshold = np.sqrt(wave**2+4)

    def integrand(u):
        nu = threshold+u*u
        excess = u*u*(nu+threshold)
        return 2*u*nu*transform(nu, pulse)**2*np.sqrt((4+excess)*excess)

    # The compact cosine-fourth envelope has a rapidly decaying tail.
    panels = np.sqrt(np.linspace(threshold, 80, 81)-threshold)
    return COEFFICIENT*sum(
        quad(integrand, a, b, epsabs=1e-12, epsrel=1e-10)[0]
        for a, b in zip(panels, panels[1:]))


def split_tangent(amplitude, dt, frequencies, strengths, contact):
    """Linearized midpoint kicks, with a separate second-vertex contact."""
    lag = np.arange(len(amplitude))*dt
    sine = np.zeros(len(lag))
    for start in range(0, len(frequencies), 256):
        batch = slice(start, start+256)
        sine += strengths[batch] @ np.sin(frequencies[batch, None]*lag)
    return contact*amplitude-dt*fftconvolve(sine, amplitude)[:len(amplitude)]
