# Drone Attitude Control Simulator

**PID Autotuning and Manual Gain Optimization in C++ & Python**

A simulator for quadcopter attitude control (roll, pitch, yaw) implemented in C++ with PID controllers. Includes autotuning (step-response based) and manual fine-tuning of gains, with performance evaluation and visualization in Python.

**Project Highlights:**

Built quaternion-based attitude dynamics model in C++.

Implemented PID autotuning and refined gains with manual tuning.

Computed control metrics: rise time, settling time, overshoot, steady-state error.

Visualized attitude response, angular rates, and control torques using Python.

**Tech Stack:**

C++ → Drone dynamics simulation & PID control

Python (NumPy, Pandas, Matplotlib) → Data analysis & visualization

CSV logging → Data exchange between simulator and analysis

**Results:**

Roll  | Rise Time: 0.180s | Settling Time: 2.280s | Overshoot: 6.53% | SSE: 0.005

Pitch | Rise Time: 0.176s | Settling Time: 2.554s | Overshoot: 4.42% | SSE: 0.003

Yaw   | Rise Time: 0.240s | Settling Time: 2.782s | Overshoot: 4.53% | SSE: 0.046


The tuned controller achieved fast response, minimal overshoot, and negligible steady-state error across all three axes.
