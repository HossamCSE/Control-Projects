# Adaptive PID Tuner (MATLAB Project)

This project implements an **adaptive PID tuning** tool using **MATLAB**. It allows users to input a transfer function via a GUI and compares the system response before and after PID optimization. The project visualizes performance improvements and displays key control metrics such as **overshoot**, **settling time**, and **steady-state error**.

---

## 🔧 Features

- GUI-based input for system transfer function (supports time delay).
- Initial PID tuning using **Skogestad method**.
- Adaptive PID optimization using gradient descent.
- Plots system response before and after tuning.
- Calculates and displays control metrics:
  - Overshoot
  - Settling Time
  - Steady-State Error

---

## 📊 Example Use

1. Run the `update_adaptive_pid_tuner.m` file.
2. Enter a transfer function like:
