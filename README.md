# 🧬 SEIHFR Ebola Model (Discrete-Time)

A discrete-time SEIHFR model for Ebola transmission, combining simulation and next generation matrix analysis. It generates plots to illustrate epidemic dynamics and validate theoretical results.

---

## 📌 Overview

This project studies the transmission dynamics of Ebola using a discrete-time SEIHFR model.

The population is divided into:

- S — Susceptible  
- E — Exposed  
- I — Infectious  
- H — Hospitalized (not contributing to transmission)  
- F — Funeral-related infectious  
- R — Removed  

---

## ⚙️ Discrete-Time Simulation

The model evolves step-by-step in time.

At each time step:

- Infection comes from I and F  
- E → I at rate α  
- I → H at rate δ  
- I → F at rate γ  
- H → R at rate η  
- F → R at rate σ  

The simulation is implemented using:

- `next_step()`  
- `simulate_model()`

---

## 📊 Infection Measure

We define:

D = E + I + F

This represents the total infection-related population.

---

## 📐 Next Generation Matrix

We consider infectious compartments:

E, I, F

H is excluded because hospitalized individuals do not generate new infections.

The reproduction value is:

R = (S/N) * (β_I σ + β_F γ) / (σ(δ + γ))

---

## 📈 Results

The script generates two plots:

1. SEIHFR population dynamics  
2. D and R over time  

These show that when R < 1, the infection level begins to decline.

---

## ▶️ How to Run

```r
source("Ebola_SEIFHR_model.R")
