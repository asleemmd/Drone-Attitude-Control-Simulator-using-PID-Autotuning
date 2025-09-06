# drone_visualize_metrics.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ---------------- Load CSV ----------------
df = pd.read_csv("sim_results.csv")

time = df['time'].values
roll = df['roll_deg'].values
pitch = df['pitch_deg'].values
yaw = df['yaw_deg'].values
p = df['p_deg'].values
q = df['q_deg'].values
r = df['r_deg'].values
tau_x = df['tau_x'].values
tau_y = df['tau_y'].values
tau_z = df['tau_z'].values

# ---------------- Step setpoints ----------------
roll_sp = 10.0
pitch_sp = -6.0
yaw_sp = 30.0

# ---------------- Metric functions ----------------
def compute_metrics(signal, sp, time, tol=0.02):
    """
    Compute Rise Time, Settling Time, Overshoot, Steady-State Error
    tol: settling error band fraction
    """
    steady_val = signal[-1]
    sse = abs(sp - steady_val)

    # Rise Time (10% -> 90%)
    try:
        if sp >= 0:
            t_r_start = time[np.where(signal >= 0.1*sp)[0][0]]
            t_r_end   = time[np.where(signal >= 0.9*sp)[0][0]]
        else:  # negative setpoint
            t_r_start = time[np.where(signal <= 0.1*sp)[0][0]]
            t_r_end   = time[np.where(signal <= 0.9*sp)[0][0]]
        rise_time = t_r_end - t_r_start
    except:
        rise_time = np.nan

    # Overshoot
    if sp >= 0:
        peak = np.max(signal)
        overshoot = (peak - sp)/abs(sp) * 100.0
    else:
        peak = np.min(signal)
        overshoot = (sp - peak)/abs(sp) * 100.0 

    # Settling Time
    if sp >= 0:
        idx_settle = np.where(np.abs(signal - sp) <= tol*abs(sp))[0]
    else:
        idx_settle = np.where(np.abs(signal - sp) <= tol*abs(sp))[0]

    if len(idx_settle) == 0:
        settling_time = np.nan
    else:
        # first time after which it stays within tol band
        for i in range(len(idx_settle)):
            if np.all(np.abs(signal[idx_settle[i]:] - sp) <= tol*abs(sp)):
                settling_time = time[idx_settle[i]]
                break
        else:
            settling_time = np.nan

    return rise_time, settling_time, overshoot, sse

# ---------------- Compute metrics ----------------
metrics = {}
for axis, sig, sp in zip(['Roll','Pitch','Yaw'], [roll,pitch,yaw], [roll_sp,pitch_sp,yaw_sp]):
    rt, st, os, sse = compute_metrics(sig, sp, time)
    metrics[axis] = {'Rise Time (s)':rt, 'Settling Time (s)':st,
                     'Overshoot (%)':os, 'SSE (deg)':sse}

# Print metrics
print("Control Metrics for Autotuned Drone PID:")
for axis, m in metrics.items():
    print(f"{axis:5} | Rise Time: {m['Rise Time (s)']:.3f}s | "
          f"Settling Time: {m['Settling Time (s)']:.3f}s | "
          f"Overshoot: {m['Overshoot (%)']:.2f}% | SSE: {m['SSE (deg)']:.3f}")

# ---------------- Visualization ----------------
plt.figure(figsize=(14,10))

# 1. Attitude angles
plt.subplot(3,1,1)
plt.plot(time, roll, label='Roll')
plt.plot(time, pitch, label='Pitch')
plt.plot(time, yaw, label='Yaw')
plt.axhline(roll_sp, color='r', linestyle='--', alpha=0.3)
plt.axhline(pitch_sp, color='g', linestyle='--', alpha=0.3)
plt.axhline(yaw_sp, color='b', linestyle='--', alpha=0.3)
plt.ylabel('Angle (deg)')
plt.title('Drone Attitude Response')
plt.legend()
plt.grid(True)

# 2. Body rates
plt.subplot(3,1,2)
plt.plot(time, p, label='p (roll rate)')
plt.plot(time, q, label='q (pitch rate)')
plt.plot(time, r, label='r (yaw rate)')
plt.ylabel('Rate (deg/s)')
plt.title('Body Angular Rates')
plt.legend()
plt.grid(True)

# 3. Control torques
plt.subplot(3,1,3)
plt.plot(time, tau_x, label='tau_x')
plt.plot(time, tau_y, label='tau_y')
plt.plot(time, tau_z, label='tau_z')
plt.xlabel('Time (s)')
plt.ylabel('Torque (Nm)')
plt.title('Control Torques')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()
