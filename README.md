# Measles Agent-Based Model (NetLogo)

This repository contains a fully structured agent-based model (ABM) developed in NetLogo for simulating measles transmission dynamics and evaluating the impact of public health interventions such as vaccination and quarantine. Designed to operate under both dynamic and scenario-based input conditions, the model is suitable for both policy exploration and academic research.

The model was primarily developed by **Parthiv Patel** at the NYC Department of Health and Mental Hygiene, with structural input and guidance from **Rebecca Conley** (Data Scientist) and **Steffen Helmer** (Director, Data Unit).

---

## Model Objective

To simulate measles outbreak dynamics at an individual level, incorporating stochastic disease progression, immunity (including pre-existing seroprevalence), vaccination processes, quarantine strategies, and post-exposure prophylaxis (PEP), under flexible user-defined intervention schedules.

---

## Core Features

### 1. Disease Dynamics

- Discrete-time simulation (1 tick = 1 day)
- Daily contact rates drawn from a **Poisson distribution** (`mean-contacts`)
- Transmission occurs with fixed **probability** (`transmission-prob`) upon contact
- Latency assigned using a **log-normal distribution** (`latent-period-mean`, `latent-period-sd`)
- Fixed **infectious period** (`infectious-duration`)
- Permanent recovery and immunity after infection
- A proportion of agents begin immune based on `seroprevalence-prob`

### 2. Vaccination Module

- Targets agents who are **susceptible**, **exposed**, or **seropositive** (but unvaccinated)
- Immunity develops after a delay (`vax-immune-onset`) with a probability (`vax-dose-1-eff`)
- Daily dose counts are supplied:
  - In **scenario mode**: constant value (`vax_val`)
  - In **dynamic mode**: date-based values (`vaccination-dates`, `vaccination-data`)
- Doses used during quarantine (`quar_vax`) are subtracted from available total
- Doses are allocated proportionally across eligible agent states

### 3. Quarantine Module

- Triggered based on:
  - Specific dates (dynamic mode)
  - Daily constant value (`quar_val`) in scenario mode
- Candidate selection:
  - If infectious agents exist, their **contact lists** are prioritized
  - Otherwise, a **random selection** is made
- A **binomial draw** determines how many exposed/infectious contacts are selected:
  - `random-binomial(n, quar_eff)`
- Quarantined agents:
  - Do not interact while in quarantine (perfect isolation)
  - Are released after `quarantine-duration` days
  - Become immune if uninfected and unvaccinated at the end

### 4. Post-Exposure Prophylaxis (PEP)

- If an agent is vaccinated within 3 days of exposure, they may bypass full infection
- This outcome depends on vaccine efficacy (`vax-dose-1-eff`)

---

## Input Modes and Control Flags

### `is-scenario?`

- **TRUE**: Ignores time-based inputs; uses fixed values for interventions
- **FALSE**: Uses lists of dates and values for dynamic vaccination/quarantine

### `continue-vax?` and `continue-quar?`

- If enabled, interventions continue using the last available value even after the input lists end

---

## Additional Functional Details

- Agents track **daily and cumulative contacts** to avoid redundant interactions
- Simulation automatically ends when no exposed or infectious agents remain
- Quarantined agents are **fully isolated** and cannot transmit or become infected
- Built-in **plots and reporters** display outbreak status, intervention performance, and agent state distributions

---

## Interface Parameters

| Parameter                | Description                                   |
|--------------------------|-----------------------------------------------|
| `mean-contacts`          | Mean daily contacts (Poisson-distributed)     |
| `transmission-prob`      | Probability of infection per contact          |
| `latent-period-mean`     | Mean of latent period distribution            |
| `latent-period-sd`       | SD of latent period distribution              |
| `infectious-duration`    | Number of days agents remain infectious       |
| `seroprevalence-prob`    | Probability of initial immunity               |
| `vax-dose-1-eff`         | Probability of gaining immunity after dose    |
| `vax-immune-onset`       | Delay between vaccination and immunity onset  |
| `quar-eff`               | Contact tracing effectiveness (binomial prob) |
| `quarantine-duration`    | Number of days in quarantine                  |

---

## Output Metrics

- **Time series plots**:
  - Daily new infections
  - SEIR state transitions
  - Cumulative cases and recoveries
- **Vaccination metrics**:
  - Daily doses administered
  - Vaccine success rate
- **Quarantine metrics**:
  - Daily quarantined agents
  - Quarantine source effectiveness
- **Epidemic outcomes**:
  - Total attack rate
  - Duration of outbreak
  - Final immune population count

---

## How to Run the Model

1. Download and install [NetLogo](https://ccl.northwestern.edu/netlogo/)
2. Open `Measles_v_12.nlogo`
3. Set simulation parameters via the interface
4. Toggle `is-scenario?` based on desired input mode
5. Load dynamic input lists if using historical mode
6. Click **Setup**, then **Go**
7. Observe real-time plots and reporter windows

---

## Contributors

- **Primary Developer**: Parthiv Patel  
  NYC Department of Health and Mental Hygiene

- **Model Structuring Support**:  
  Rebecca Goldberg, Data Scientist  
  Steffen Foerster, Director, Data Unit

---

## License

This model is released under the GNU General Public License (GPL v3), consistent with the NetLogo platform. It is intended for academic, research, and public health policy use.

---

## Citation

If using or referencing this model, please cite:

**Patel, P., Conley, R., & Helmer, S. (2025).** *Agent-Based Measles Transmission and Intervention Model in NetLogo*. GitHub. https://github.com/your-username/measles-abm-netlogo
