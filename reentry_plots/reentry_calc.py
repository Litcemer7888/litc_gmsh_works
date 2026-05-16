import numpy as np
from scipy.integrate import solve_ivp
from gost_atmosphere import density_gost81
from apollo_aero import get_aero_coefficients
import pandas as pd
import matplotlib.pyplot as plt

# physical constants and parameters
R_earth = 6.371e6
GM = 3.986004418e14
g0 = 9.80665

mass = 5500.0             # kg
nose_radius = 3.5
frontal_area = np.pi * nose_radius**2
alpha_deg = 27.0
alpha_rad = np.radians(alpha_deg)

# simulation parameters
t_span = (0.0, 2000.0)
rtol = 1e-8
atol = 1e-10

h0 = 120_000.0
v0 = 7_830.0             # 11.03 km/s - lunar return
gamma_deg0 = -5.6
gamma0 = np.radians(gamma_deg0)
s0 = 0.0

y0 = [h0, v0, gamma0, s0]

def aerodynamic_forces(velocity_mps, altitude_m, alpha_deg):
    """Returns (drag, lift) in Newtons."""
    rho = density_gost81(altitude_m)
    q = 0.5 * rho * velocity_mps**2
    
    CL, CD, mach, LD = get_aero_coefficients(velocity_mps, altitude_m, alpha_deg)
    
    drag = q * frontal_area * CD
    lift = q * frontal_area * CL
    
    return drag, lift#, mach, LD

def reentry_odes(t, y):
    h, v, gamma, s = y
    r = R_earth + h
    g = GM / r**2
    
    D, L, = aerodynamic_forces(v, h, alpha_deg)
    a_drag = D / mass
    a_lift = L / mass
    
    dv_dt = -a_drag - g * np.sin(gamma)
    dgamma_dt = (v / r - g / v) * np.cos(gamma) + (a_lift / v)
    dh_dt = v * np.sin(gamma)
    ds_dt = v * np.cos(gamma)
    
    return [dh_dt, dv_dt, dgamma_dt, ds_dt]


def velocity_event(t, y):
    return y[1] - 500.0
velocity_event.terminal = True
velocity_event.direction = -1

sol = solve_ivp(reentry_odes, t_span, y0,
                method='DOP853',
                events=velocity_event,
                rtol=1e-8, atol=1e-10,
                max_step=5.0)

t = sol.t
h = sol.y[0]
v = sol.y[1]
gamma = sol.y[2]
s = sol.y[3]

h_km = h / 1000.0
v_kms = v / 1000.0
gamma_deg = np.degrees(gamma)
s_km = s / 1000.0

mach = np.zeros_like(v)
LD_history = np.zeros_like(v)
for i, (vi, hi) in enumerate(zip(v, h)):
    _, _, mach[i], LD_history[i] = get_aero_coefficients(vi, hi, alpha_deg)

g0 = 9.80665
accel_mps2 = np.gradient(v, t)
accel_g = -accel_mps2 / g0

# print("\n" + "=" * 60)
print("=" * 60)
# print(f"Mass = {mass} kg, Nose radius = {nose_radius} m")
# print(f"Reference area = {frontal_area:.2f} m2")
print(f"AoA = const = {alpha_deg}")
print(f"initial: h = {h0/1000:.1f} km, v = {v0/1000:.2f} km/s, γ = {np.degrees(gamma0):.2f}")
print(f"final:   t = {t[-1]:.1f} s, h = {h_km[-1]:.1f} km, v = {v_kms[-1]:.2f} km/s")
print(f"max deceleration: {np.max(accel_g):.1f} g")
print(f"max L/D ratio: {np.max(LD_history):.3f}")
print(f"downrange distance: {s_km[-1]:.1f} km")

df = pd.DataFrame({
    'time (s)': t,
    'altitude (km)': h_km,
    'velocity (km/s)': v_kms,
    'mach': mach,
    'flight_path_angle (deg)': gamma_deg,
    'downrange (km)': s_km,
    'acceleration (g)': accel_g,
    'lift/drag Ratio': LD_history,
    'angle_of_attack': alpha_deg
})
df.to_csv(f'reentry_v_gamma_{v0}{gamma_deg0}.csv', index=False)
print("\saved to" f'reentry_v_gamma_{v0}{gamma_deg0}.csv')



plt.figure(figsize=(16, 10))

plt.subplot(2,3,1)
plt.plot(t, h_km, 'b-', linewidth=2)
plt.xlabel('time, s')
plt.ylabel('altitude, km')
plt.grid(True, alpha=0.3)
plt.title('Altitude | Time')

plt.subplot(2,3,2)
plt.plot(t, v_kms, 'r-', linewidth=2)
plt.xlabel('time, s')
plt.ylabel('Velocity, km/s')
plt.grid(True, alpha=0.3)
plt.title('velocity | time')

plt.subplot(2,3,3)
plt.plot(t, accel_g, 'g-', linewidth=2)
plt.xlabel('time, s')
plt.ylabel('acceleration, g')
plt.grid(True, alpha=0.3)
plt.title('g-force | time')

plt.subplot(2,3,4)
plt.plot(t, gamma_deg, 'm-', linewidth=2)
plt.xlabel('time, s')
plt.ylabel('angle, deg')
plt.grid(True, alpha=0.3)
plt.title('flight path angle | time')

plt.subplot(2,3,5)
plt.plot(s_km, h_km, 'c-', linewidth=2)
plt.xlabel('downrange, km')
plt.ylabel('altitude, km')
plt.grid(True, alpha=0.3)
plt.title('horizontally scewed trajectory')
#plt.gca().invert_yaxis()

plt.subplot(2,3,6)
plt.plot(t, LD_history, 'orange', linewidth=2, label='L\D')
plt.plot(t, mach, 'purple', linewidth=2, label='mach')
plt.xlabel('time, s')
plt.ylabel('LD / Mach')
plt.grid(True, alpha=0.3)
plt.title('L/D Ratio and mach | time')
plt.legend()

plt.tight_layout()
plt.savefig(f'reentry_v_gamma_{v0}{gamma_deg0}.png', dpi=150)
plt.show()

# print("\nPlots saved to", f'reentry_v_gamma_{v0}{gamma_deg0}.png')
