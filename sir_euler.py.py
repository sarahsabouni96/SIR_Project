from math import log
import matplotlib.pyplot as plt

# SIR parameters and initial populations
N = 1000
beta = 0.3
gamma = 0.1
S0, I0, R0 = 990.0, 10.0, 0.0
days = 160


def euler_sir(h):
    """Compute the SIR populations over 160 days with step size h."""
    if h <= 0:
        raise ValueError("The step size must be positive.")

    steps = round(days / h)
    if abs(steps * h - days) > 1e-10:
        raise ValueError("Choose a step size that divides the time interval.")

    time = [0.0]
    S = [S0]
    I = [I0]
    R = [R0]

    for n in range(1, steps + 1):
        # Use the populations at the start of the step for both transfers.
        s, i, r = S[-1], I[-1], R[-1]
        infections = h * beta * s * i / N
        removals = h * gamma * i

        S.append(s - infections)
        I.append(i + infections - removals)
        R.append(r + removals)
        time.append(n * h)

        assert min(S[-1], I[-1], R[-1]) >= 0
        assert abs(S[-1] + I[-1] + R[-1] - N) < 1e-8

    return time, S, I, R


# Analytical peak for the chosen initial conditions
S_peak = gamma * N / beta
exact_peak = I0 + S0 - S_peak + S_peak * log(S_peak / S0)

# Compare step sizes
print(f"Exact peak I = {exact_peak:.6f}")
print("h (days)   peak day    peak I   abs. error    S(160)")

plot_step = 0.125
for h in (2, 1, 0.5, 0.25, plot_step):
    time, S, I, R = euler_sir(h)

    peak_I = max(I)
    peak_day = time[I.index(peak_I)]
    error = abs(peak_I - exact_peak)
    print(f"{h:8.3f} {peak_day:10.2f} {peak_I:9.3f} {error:12.3f} {S[-1]:9.3f}")

    # Draw the figure using the h = 0.125 day run.
    if h == plot_step:
        plt.figure(figsize=(7, 4))
        plt.plot(time, S, "-",  color="#0072B2", label="Susceptible")
        plt.plot(time, I, "--", color="#D55E00", label="Infectious")
        plt.plot(time, R, "-.", color="#333333", label="Removed")
        plt.xlabel("Time (days)")
        plt.ylabel("Number of individuals")
        plt.xlim(0, days)
        plt.ylim(0, N)
        plt.legend()
        plt.tight_layout()
        plt.show()