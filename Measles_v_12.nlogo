;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                              MODEL PREAMBLE                                          ;;
;;                                                                                      ;;
;;                                                                                      ;;
;;         ================== KEY ASSUMPTIONS ===================                       ;;
;;                                                                                      ;;
;; 1. Disease Dynamics:                                                                 ;;
;;    - The model operates in discrete time steps (each tick = 1 day).                  ;;
;;    - Infected agents interact with others once per day, with their daily             ;;
;;      contact limit drawn from a Poisson distribution (mean-contacts).                ;;
;;    - Transmission is modeled by having infectious agents expose susceptible          ;;
;;      agents with a probability equal to the base transmission rate.                  ;;
;;    - Upon exposure, agents enter a latent period (assigned via a log-normal          ;;
;;      distribution) before becoming infectious.                                       ;;
;;    - Infectious agents remain so for a defined infectious period, after              ;;
;;      which they recover and become permanently immune.                               ;;
;;    - A proportion of agents start with pre-existing immunity (seroprevalence)        ;;
;;      based on a defined probability; these agents are marked as sero-positive        ;;
;;      and remain immune.                                                              ;;
;;                                                                                      ;;
;; 2. Vaccination Assumptions:                                                          ;;
;;    - Eligible agents include those in one of the following states who have           ;;
;;      not been vaccinated before:                                                     ;;
;;         • Susceptible                                                                ;;
;;         • Exposed                                                                    ;;
;;         • Seroprevalence positive (i.e., already have antibodies from prior          ;;
;;           exposure but are not vaccinated)                                           ;;
;;    - After vaccination, there is a delay (vax-immune-onset) before an agent          ;;
;;      “tries” to gain full immunity. Even then, full immunity is achieved only        ;;
;;      with a probability (vax-dose-1-eff).                                            ;;
;;    - The daily number of vaccine doses is provided via interface widgets.            ;;
;;      In dynamic mode, a list of dates and corresponding dose values is used;         ;;
;;      in scenario mode, a constant value (vax_val) is applied each day.               ;;
;;    - If any doses have been used during quarantine (tracked as quar_vax), they       ;;
;;      are subtracted from the daily total.                                            ;;
;;                                                                                      ;;
;; 3. Quarantine Assumptions:                                                           ;;
;;    - Quarantine interventions are triggered on specific dates (provided via          ;;
;;      interface inputs) or, in scenario mode, applied every day with a fixed          ;;
;;      target number (quar_val).                                                       ;;
;;    - When infectious agents are present, the model uses their contact lists to       ;;
;;      identify candidates for quarantine; otherwise, a random selection from          ;;
;;      eligible agents (susceptible, exposed, or seropositive) is made.                ;;
;;    - A binomial process (via the random-binomial function) is used to determine      ;;
;;      the number of contacts to quarantine from those already exposed or              ;;
;;      infectious. For example, if an infectious agent has 10 contacts and             ;;
;;      quar_eff is 0.3, then on average about 3 contacts are selected.                 ;;
;;    - Quarantined agents remain isolated for a fixed duration (quarantine-duration).  ;;
;;      After this period, if they are still not vaccinated or seropositive, they       ;;
;;      are released and automatically marked as immune.                                ;;
;;    - The model also includes a post-exposure prophylaxis (PEP) check where, if       ;;
;;      an agent is vaccinated shortly after exposure (within 3 days), it may           ;;
;;      recover directly (bypassing full-blown infection) if the efficacy check         ;;
;;      passes.                                                                         ;;
;;                                                                                      ;;
;;       ================== MODE SWITCHES & DATA INPUT ===================              ;;
;;                                                                                      ;;
;; - is-scenario? Flag:                                                                 ;;
;;    • When OFF (false):                                                               ;;
;;         - The model uses dynamic intervention data provided via interface            ;;
;;           widgets (lists of dates and corresponding values for vaccination           ;;
;;           and quarantine).                                                           ;;
;;         - The flags continue-vax? and continue-quar? determine whether the           ;;
;;           interventions continue with the last provided value after the              ;;
;;           historical data is exhausted.                                              ;;
;;    • When ON (true):                                                                 ;;
;;         - The model ignores the date lists and uses constant daily values            ;;
;;           (vax_val for vaccination and quar_val for quarantine) throughout           ;;
;;           the simulation.                                                            ;;
;;                                                                                      ;;
;;           ================== VACCINATION PROCESS ===================                 ;;
;;                                                                                      ;;
;; 1. Eligibility & Daily Dose Calculation:                                             ;;
;;    - Eligible agents are those in one of the following states who have not           ;;
;;      been vaccinated:                                                                ;;
;;         • Susceptible                                                                ;;
;;         • Exposed                                                                    ;;
;;         • Seroprevalence positive                                                    ;;
;;    - In dynamic mode (is-scenario? OFF), the model checks the current date           ;;
;;      against the vaccination-dates provided via interface. If the date matches,      ;;
;;      the corresponding dose (from vaccination-data) is administered. In              ;;
;;      scenario mode, the constant vax_val is used each day.                           ;;
;;    - Any doses already applied during quarantine (quar_vax) are subtracted           ;;
;;      from the daily total.                                                           ;;
;;                                                                                      ;;
;; 2. Dose Allocation Example:                                                          ;;
;;    - Suppose on a given day, 100 doses are available and there are 200               ;;
;;      eligible agents (e.g., 150 susceptible, 30 exposed, 20 sero-positive).          ;;
;;      Then the doses are roughly allocated proportionally:                            ;;
;;         • Susceptible agents receive (150/200)*100 ≈ 75 doses                        ;;
;;         • Exposed agents receive (30/200)*100 ≈ 15 doses                             ;;
;;         • Sero-positive agents receive (20/200)*100 ≈ 10 doses                       ;;
;;    - Each vaccinated agent is flagged as vaccinated, its time-vaccinated             ;;
;;      counter is reset, and after vax-immune-onset days, it may gain full             ;;
;;      immunity based on the probability (vax-dose-1-eff).                             ;;
;;                                                                                      ;;
;;            ================== QUARANTINE PROCESS ===================                 ;;
;;                                                                                      ;;
;; 1. Triggering Quarantine:                                                            ;;
;;    - On specific dates (in dynamic mode) or continuously (in scenario mode),         ;;
;;      the model attempts to quarantine a fixed number of agents as defined by         ;;
;;      the interface (quarantine-dates and quarantine-data in dynamic mode, or         ;;
;;      quar_val in scenario mode).                                                     ;;
;;    - If infectious agents are detected, their contact lists are used to              ;;
;;      identify candidates for quarantine; if not, a random selection is made          ;;
;;      from eligible agents (susceptible, exposed, or sero-positive).                  ;;
;;                                                                                      ;;
;; 2. Selection via Binomial Process:                                                   ;;
;;    - The contacts of an infectious agent are divided into two groups:                ;;
;;         • Exposed/infectious contacts (E_I_contacts)                                 ;;
;;         • Susceptible/immune contacts (S_IM_contacts)                                ;;
;;    - To determine how many of the exposed/infectious contacts to quarantine, a       ;;
;;      binomial draw is used:                                                          ;;
;;         let quar_from_exp random-binomial (length E_I_contacts) quar_eff             ;;
;;    - Example: If there are 10 contacts and quar_eff = 0.3, then on average,          ;;
;;      about 3 contacts are selected for quarantine.                                   ;;
;;    - The remaining number needed to reach the target is selected from the            ;;
;;      susceptible/immune group.                                                       ;;
;;                                                                                      ;;
;; 3. Application & Release:                                                            ;;
;;    - Selected agents are marked as quarantined (their epi-state is set to            ;;
;;      "Quarantine", they are colored orange, and their quarantine timers reset).      ;;
;;    - After a fixed quarantine-duration, agents who have not been vaccinated or       ;;
;;      become seropositive are released and automatically marked as immune.            ;;
;;                                                                                      ;;
;; ================== BINOMIAL PROCESS & VARIABILITY ===================                ;;
;;                                                                                      ;;
;; - The function random-binomial(n, p) simulates n independent trials, where           ;;
;;   each trial has a success probability p (here, p = quar_eff).                       ;;
;; - This ensures that only a fraction of contacts (on average) are selected for        ;;
;;   quarantine. For instance, with 10 trials and p = 0.3, the expected outcome is 3.   ;;
;; - Although we use a simple binomial process, a negative binomial distribution        ;;
;;   could be implemented if greater variability (over-dispersion) is needed. This      ;;
;;   would allow for occasional days with significantly more or fewer contacts          ;;
;;   being quarantined than the mean.                                                   ;;
;;                                                                                      ;;
;; ================== MODE FLAGS FOR CONTINUATION ===================                   ;;
;;                                                                                      ;;
;; - continue-vax?: If true, the vaccination intervention continues using the last      ;;
;;   provided dose value even after the scheduled (interface-defined) dates have        ;;
;;   been exhausted.                                                                    ;;
;; - continue-quar?: Similarly, if true, the quarantine intervention continues with     ;;
;;   the last provided target number once the scheduled dates have passed.              ;;
;;                                                                                      ;;
;; ================== ADDITIONAL CLARIFICATIONS ===================                     ;;
;;                                                                                      ;;
;; - All “historical” intervention data (for both vaccination and quarantine) is        ;;
;;   supplied via interface widgets (lists of dates and corresponding values), not      ;;
;;   from CSV files.                                                                    ;;
;; - Agents track their daily and cumulative contacts to prevent duplicate              ;;
;;   interactions in a single day.                                                      ;;
;; - The simulation stops automatically when there are no exposed or infectious         ;;
;;   agents remaining, indicating the outbreak has ended.                               ;;
;; - Once agents enter quarantine or isolation, they are assumed to have zero           ;;
;;   interactions with the general population (perfect isolation).                      ;;
;; - Plots and reporters continuously update to reflect daily cases, vaccination        ;;
;;   counts, quarantine counts, and other key metrics to monitor outbreak               ;;
;;   progression.                                                                       ;;
;;                                                                                      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



extensions [time csv]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                BREED AND GLOBALS                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

breed [persons person]

globals [
  ;; Global simulation variables
  date                          ; Current simulation date, created from start-date widget
  interaction-radius            ; Radius within which agents can interact (set to 2 units)

  ;; Status codes for epidemiological states
  susceptible-code              ; "Susceptible" state for uninfected agents
  exposed-code                  ; "Exposed" state (infected but not yet infectious)
  infectious-code               ; "Infectious" state (able to transmit measles)
  recovered-code                ; "Recovered" state (recovered and immune)
  immune-code                   ; "Immune" state (fully immune, possibly post-vaccination)
  quarantine-code               ; "Quarantine" state for agents who are isolated

  day                           ; Simulation day counter
  new-cases-pop-1               ; Daily new infection counter for population group 1
  new-cases-pop-2               ; Daily new infection counter for population group 2
  new-cases-pop-3               ; Daily new infection counter for population group 3
  new-cases-pop-4               ; Daily new infection counter for population group 4
  new-cases-pop-5               ; Daily new infection counter for population group 5

  exposed-cases-pop-1           ; Daily exposed cases for group 1
  exposed-cases-pop-2           ; Daily exposed cases for group 2
  exposed-cases-pop-3           ; Daily exposed cases for group 3
  exposed-cases-pop-4           ; Daily exposed cases for group 4
  exposed-cases-pop-5           ; Daily exposed cases for group 5

  susceptible-pop-1             ; Number of susceptible persons (group 1)
  susceptible-pop-2             ; Number of susceptible persons (group 2)
  susceptible-pop-3             ; Number of susceptible persons (group 3)
  susceptible-pop-4             ; Number of susceptible persons (group 4)
  susceptible-pop-5             ; Number of susceptible persons (group 5)

  recovered-pop-1               ; Daily recovered count for group 1
  recovered-pop-2               ; Daily recovered count for group 2
  recovered-pop-3               ; Daily recovered count for group 3
  recovered-pop-4               ; Daily recovered count for group 4
  recovered-pop-5               ; Daily recovered count for group 5

  pop-1-total-cases             ; Cumulative case count for group 1
  pop-2-total-cases             ; Cumulative case count for group 2
  pop-3-total-cases             ; Cumulative case count for group 3
  pop-4-total-cases             ; Cumulative case count for group 4
  pop-5-total-cases             ; Cumulative case count for group 5

  total-cases                   ; Global counter for total cases (if needed)
  daily-vax-dose                ; Daily vaccination dose available (set via interface or CSV)
  param_data                    ; List to store parameters loaded from file (if used)
  vaccination-dates             ; List of dates on which vaccination doses are specified
  vaccination-data              ; List of vaccination dose values corresponding to the dates
  quarantine-data               ; List of quarantine dose values from interface or CSV
  quarantine-dates              ; List of dates when quarantine is applied
  latent-dur-min                ; Minimum latent period duration (in days)
  latent-dur-max                ; Maximum latent period duration (in days)
  follow-up-period              ; Follow-up period for the study (in days)
  daily-quarantine-count        ; Daily count of quarantined persons
  today-vaccinate-dose          ; Number of vaccine doses to be given on a given day
  quar-vax                      ; Count of vaccinations performed during quarantine
  total-quar                    ; Total number of agents that have ever been quarantined
  inf-quar                      ; Number of infected agents that moved to quarantine
  effective-capture-count       ; Count of quarantined agents who were actually exposed/infectious
  vax-started?                  ; Boolean flag to indicate if vaccination has begun
  expo_leng                     ; (Not used in current version, placeholder for exposure length)
  mean-contacts
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                        PERSONS-OWN VARIABLES (AGENT ATTRIBUTES)            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
persons-own [
  ;; Attributes for each person/agent
  age-group                     ; Age group of the agent (e.g., "1", "2", "3", "18+", "5")
  epi-state                     ; Epidemiological state (e.g., susceptible, exposed, infectious, etc.)
  base-transmission-rate        ; Base probability of transmitting the disease upon contact
  latent-period                 ; Duration (in days) before an exposed agent becomes infectious
  infectious-period             ; Duration (in days) during which the agent is infectious
  contacts                      ; List of agents with which this agent has had contact
  immunity                      ; Immunity level (1 indicates full immunity)
  vaccinated?                   ; Boolean flag: true if the agent has been vaccinated
  isolation?                    ; Boolean flag: true if the agent is in isolation (used for quarantine logic)
  quarantine?                   ; Boolean flag: true if the agent is quarantined
  quarantine-time               ; Number of days the agent has been in quarantine
  time-exposed                  ; Counter for days since the agent was exposed
  time-infected                 ; Counter for days since the agent became infectious
  time-vaccinated               ; Counter for days since the agent was vaccinated
  daily-contacts                ; List tracking daily contacts (to limit multiple interactions per day)
  num-daily-contacts            ; Poisson-distributed number of daily contacts (simulating realistic behavior)
  serology?                     ; Boolean flag: true if the agent tested positive for antibodies
  isolation-time                ; Counter for days in isolation (used during recovery/quarantine)
]

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                             SETUP PROCEDURES                               ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; Main setup procedure: Clears the simulation, sets globals, creates agents, and
;;; initializes the infected and seroprevalence populations.
to setup
  ca  ;; clear-all the agents and graphics
  ;profiler:start
  setup-globals   ; Initialize global variables and parameters
  setup-persons   ; Create agents, assign them age groups and initial states
  reset-timer     ; Reset simulation timer (NetLogo built-in)
end

;;; setup-globals: Initializes global variables such as status codes, simulation date,
;;; vaccination and quarantine data, and other counters.
to setup-globals
  random-seed seed  ;; Ensures reproducibility if seed is set from interface
  set susceptible-code "Susceptible"
  set exposed-code "Exposed"
  set infectious-code "Infectious"
  set recovered-code "Recovered"
  set immune-code "Immune"
  set quarantine-code "Quarantine"

  ;; Set the simulation start date from an interface widget (start-date)
  set date time:create-with-format start-date "yyyy-MM-dd"
  set day 0
  set interaction-radius 2

  ;; Initialize daily new case counters for each population group
  set new-cases-pop-1 0
  set new-cases-pop-2 0
  set new-cases-pop-3 0
  set new-cases-pop-4 0
  set new-cases-pop-5 0

  ;; Initialize exposed case counters for each group
  set exposed-cases-pop-1 0
  set exposed-cases-pop-2 0
  set exposed-cases-pop-3 0
  set exposed-cases-pop-4 0
  set exposed-cases-pop-5 0

  ;; Initialize recovered counts for each group
  set recovered-pop-1 0
  set recovered-pop-2 0
  set recovered-pop-3 0
  set recovered-pop-4 0
  set recovered-pop-5 0

  ;; Set cumulative total cases for each group (using initial infection values)
  set pop-1-total-cases initial-inf-pop-1
  set pop-2-total-cases initial-inf-pop-2
  set pop-3-total-cases initial-inf-pop-3
  set pop-4-total-cases initial-inf-pop-4
  set pop-5-total-cases initial-inf-pop-5

  set daily-vax-dose 0
  set daily-quarantine-count 0
  set quar-vax 0
  set total-quar 0
  set inf-quar 0
  set effective-capture-count 0

  ;; Set the mean contact rates
  set mean-contacts R0 / (inf-dur-index * transmission-prob)
  ;; Disable vaccination if vaccine dose value is 0
  if vax_val = 0 [
    set is-vaccination? false
  ]
  set vax-started? false

  ;; Initialize vaccination data lists
  set vaccination-dates []
  set vaccination-data []

  ;; Set up vaccination data: either from scenario mode (continuous value) or
  ;; historical data loaded via CSV (using interface widgets)
  ifelse is-scenario?  [
    ; No historical data: use the provided continuous value
    set vaccination-data vax_val
  ] [
    ; Historical data exists: process the CSV rows for dates and values
    let vaccination-raw-dates csv:from-row vax-raw-dates
    foreach vaccination-raw-dates [
      row ->
      set vaccination-dates lput time:create-with-format row "yyyy-MM-dd" vaccination-dates
    ]
    let vaccination-raw-values csv:from-row vax-raw-values
    set vaccination-data map [[row] -> row ] vaccination-raw-values

    if length vaccination-dates != length vaccination-data [
      user-message "The number of vaccination dates and values do not match"
      stop
    ]
  ]

  ;; Same approach for quarantine data:
  set quarantine-dates []
  set quarantine-data []

  ifelse is-scenario? [
    ; No historical data: use a single continuous value for quarantine
    set quarantine-data quar_val
  ] [
    let quarantine-raw-dates csv:from-row quar-raw-dates
    foreach quarantine-raw-dates [
      row ->
      set quarantine-dates lput time:create-with-format row "yyyy-MM-dd" quarantine-dates
    ]
    let quarantine-raw-values csv:from-row quar-raw-values
    foreach quarantine-raw-values [
      row ->
      set quarantine-data lput row quarantine-data
    ]
    if length quarantine-dates != length quarantine-data [
      user-message "The number of quarantine dates and values do not match"
      stop
    ]
  ]
end

;;; setup-persons: Creates agents for each population group based on provided
;;; population sizes (pop-1, pop-2, etc.). Each agent is assigned an age group.
to setup-persons
  create-persons pop-1 [
    set age-group "1"
    initialize-person  ;; Set initial properties for the agent
  ]
  create-persons pop-2 [
    set age-group "2"
    initialize-person
  ]
  create-persons pop-3 [
    set age-group "3"
    initialize-person
  ]
  create-persons pop-4 [
    set age-group "18+"
    initialize-person
  ]
  create-persons pop-5 [
    set age-group "5"
    initialize-person
  ]
  initial-infected         ;; Seed the model with initial infectious persons
  initial-seroprevalence   ;; Randomly assign seroprevalence (pre-existing immunity)
end

;;; initialize-person: Sets the default attributes for a person/agent.
to initialize-person
  set color blue                ;; Default color for susceptible agents
  set size 1                    ;; Visual size of the agent
  set shape "person"            ;; Custom shape representing a person
  setxy random-xcor random-ycor ;; Place agent at a random location
  set epi-state susceptible-code  ;; All agents start as susceptible
  set contacts []               ;; Initialize contacts list (for tracking interactions)
  set vaccinated? false         ;; Not yet vaccinated
  set isolation? false          ;; Not in isolation by default
  set quarantine? false         ;; Not quarantined by default
  set quarantine-time 0         ;; No quarantine time initially
  set base-transmission-rate 0  ;; Will be set when agent becomes infectious
  set latent-period 0           ;; Latent period (will be assigned when exposed)
  set infectious-period 0       ;; Infectious period (will be assigned upon infection)
  set time-exposed 0            ;; Counter for time since exposure
  set time-infected 0           ;; Counter for time since becoming infectious
  set time-vaccinated 0         ;; Counter for time since vaccination
  set immunity 0                ;; No immunity initially (0 means susceptible)
  set daily-contacts []         ;; Daily interaction tracking list
  set serology? false           ;; Not serology positive initially
  set isolation-time 0          ;; Counter for time in isolation/quarantine
end

;;; initial-infected: Randomly selects a specified number of agents from each age group
;;; and sets them to the infectious state. It also assigns their infectious period and
;;; transmission rate.
to initial-infected
  ask n-of initial-inf-pop-1 persons with [age-group = "1"] [
    set epi-state infectious-code
    set color red
    set infectious-period inf-dur-index
    assign-transmission-rate
    set num-daily-contacts random-poisson mean-contacts
  ]
  ask n-of initial-inf-pop-2 persons with [age-group = "2"] [
    set epi-state infectious-code
    set color red
    set infectious-period inf-dur-index
    assign-transmission-rate
    set num-daily-contacts random-poisson mean-contacts
  ]
  ask n-of initial-inf-pop-3 persons with [age-group = "3"] [
    set epi-state infectious-code
    set color red
    set infectious-period inf-dur-index
    assign-transmission-rate
    set num-daily-contacts random-poisson mean-contacts
  ]
  ask n-of initial-inf-pop-4 persons with [age-group = "18+"] [
    set epi-state infectious-code
    set color red
    set infectious-period inf-dur-index
    assign-transmission-rate
    set immunity 0
    set num-daily-contacts random-poisson mean-contacts
  ]
  ask n-of initial-inf-pop-5 persons with [age-group = "5"] [
    set epi-state infectious-code
    set color red
    set infectious-period inf-dur-index
    assign-transmission-rate
    set num-daily-contacts random-poisson mean-contacts
  ]
end

;;; initial-seroprevalence: Randomly marks some susceptible agents as already immune
;;; to reflect pre-existing immunity. These agents are colored violet.
to initial-seroprevalence
  ask persons with [epi-state = susceptible-code] [
    if random-float 1 < per-seroprevalence [
      set epi-state immune-code
      set immunity sero-immunity
      set serology? true
      set color violet
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                           MAIN SIMULATION LOOP                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; go: This is the main simulation loop. Each tick (day), the following steps occur:
;;;   1. Check if there are any exposed or infectious agents. If not, the simulation stops.
;;;   2. Update the simulation clock and reset daily counters.
;;;   3. Move agents (simulate natural movement).
;;;   4. Allow disease exposure via interactions.
;;;   5. Transition exposed agents to infectious when their latent period elapses.
;;;   6. Handle recovery and quarantine of infectious agents.
;;;   7. Apply interventions (vaccination and quarantine) as needed.
;;;   8. Update plots with current data.
to go
  carefully [let t ticks] [reset-ticks reset-timer]

  ;; If there are no exposed or infectious persons, end the simulation.
  if count (persons with [epi-state = exposed-code]) = 0 and count (persons with [epi-state = infectious-code]) = 0 [
    ;print profiler:report
    ;print timer
    stop
  ]

  clock                 ; Advance the simulation day and update the date
  reset-daily-counts    ; Reset daily counters (cases, vaccinations, quarantine counts, etc.)
  interact              ; Move agents and simulate their interactions
  expose                ; Infected agents expose susceptible agents (disease transmission)
  make-infectious       ; Convert exposed agents to infectious based on latent period
  recover               ; Process recovery (or quarantine) of infectious agents
  intervention-setup    ; Apply daily interventions (vaccination and quarantine)
  ;update-susceptible-pop  ;; (Optional) update count of susceptible persons

  plot-data           ; Update all plots (daily cases, vaccination, quarantine, etc.)
  calculate-total-cases

  ; Debug print statements can be enabled for step-by-step tracking:
;  print(word "Day: " day)
;  print(word "Susceptible People: " susc-pop-4)
;  print(word "Exposed People: " count persons with [epi-state = exposed-code])
;  print(word "Recovered People: " count persons with [epi-state = recovered-code])
;  print(word "Infectious People: " count persons with [epi-state = infectious-code])
;  print(word "Immune People: " count persons with [epi-state = immune-code])
;  print(word "Quarantine People: " count persons with [epi-state = quarantine-code])
;  print(word "Total population: " total-sum)
;  print(word "Cumulative Total quarantined for this infection: " current-quarantine-count)
;  print(word "Daily Total quarantined for this infection: " daily-quarantine-count)
;  print(word "Daily cases: " daily-cases-pop-4)
;  print(word "Total cases: " total-cases-pop-4)
;  print(word "Total quarantined: " total-quar)
;  print(word "________________________________________")

  tick-advance 1   ;; Advance simulation tick (NetLogo built-in)
end

;;; total-sum: Reporter to calculate the total population across all epidemiological states.
to-report total-sum
  let total-pop susc-pop-4 + (count persons with [epi-state = exposed-code]) +
                (count persons with [epi-state = recovered-code]) +
                (count persons with [epi-state = infectious-code]) +
                (count persons with [epi-state = immune-code]) +
                (count persons with [epi-state = quarantine-code])
  report total-pop
end

;;; update-susceptible-pop: Updates the global counters for susceptible persons per age group.
to update-susceptible-pop
  set susceptible-pop-1 count persons with [age-group = "1" and epi-state = susceptible-code]
  set susceptible-pop-2 count persons with [age-group = "2" and epi-state = susceptible-code]
  set susceptible-pop-3 count persons with [age-group = "3" and epi-state = susceptible-code]
  set susceptible-pop-4 count persons with [age-group = "18+" and epi-state = susceptible-code]
  set susceptible-pop-5 count persons with [age-group = "5" and epi-state = susceptible-code]
end

;;; intervention-setup: Applies quarantine and vaccination interventions each day.
to intervention-setup
  if is-quarantine? [    ;; Check if quarantine intervention is enabled.
      check-quarantine         ;; Check for contacts to quarantine based on current data
      update-quarantine-status ;; Update status of agents currently in quarantine
    ]
  if is-vaccination? [   ;; Check if vaccination intervention is enabled.
      calculate-today-vaccinate  ;; Determine the number of vaccines available today
      vaccinate                  ;; Vaccinate eligible agents
      update-vaccine-immunity    ;; Update immunity for agents after vaccine onset period
    ]
end

;;; clock: Advances the simulation day and updates the current date.
to clock
  set day day + 1
  set date time:plus date 1 "days"
end

;;; reset-daily-counts: Resets daily counters and updates cumulative totals.
to reset-daily-counts
  ;; Reset daily new cases and exposed cases counters to zero.
  ifelse day = 1 [
    set new-cases-pop-1 initial-inf-pop-1
    set new-cases-pop-2 initial-inf-pop-2
    set new-cases-pop-3 initial-inf-pop-3
    set new-cases-pop-4 initial-inf-pop-4
    set new-cases-pop-5 initial-inf-pop-5
  ][
    set new-cases-pop-1 0
    set new-cases-pop-2 0
    set new-cases-pop-3 0
    set new-cases-pop-4 0
    set new-cases-pop-5 0
  ]

  set exposed-cases-pop-1 0
  set exposed-cases-pop-2 0
  set exposed-cases-pop-3 0
  set exposed-cases-pop-4 0
  set exposed-cases-pop-5 0
  set recovered-pop-1 0
  set recovered-pop-2 0
  set recovered-pop-3 0
  set recovered-pop-4 0
  set recovered-pop-5 0

  set daily-vax-dose 0
  set daily-quarantine-count 0
  set today-vaccinate-dose 0
  set quar-vax 0
  set inf-quar 0
  set effective-capture-count 0

  ;; Reset each infectious agent's daily contacts and assign a new contact limit (unless quarantined)
  ask persons with [epi-state = infectious-code] [
    ifelse quarantine? [
      set num-daily-contacts 0
    ] [
      set num-daily-contacts random-poisson mean-contacts
    ]
    set daily-contacts []
  ]
end

;;; plot-data: Calls individual plotting procedures to update visualizations.
to plot-data
  plot-cases           ;; Plot daily new cases per age group
  ;plot-epidata       ;; (Optional) Plot overall epidemiological data (exposed, infected, recovered)
  plot-vaccination     ;; Plot daily vaccinations given
  plot-quarantine      ;; Plot total number of quarantined persons
  plot-contacts        ;; Plot daily contact statistics
end


to calculate-total-cases
  ;; Update cumulative totals by adding today's new cases to the running total.
  if day > 1 [
    set pop-1-total-cases pop-1-total-cases + new-cases-pop-1
    set pop-2-total-cases pop-2-total-cases + new-cases-pop-2
    set pop-3-total-cases pop-3-total-cases + new-cases-pop-3
    set pop-4-total-cases pop-4-total-cases + new-cases-pop-4
    set pop-5-total-cases pop-5-total-cases + new-cases-pop-5
  ]
end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                           VACCINATION PROCEDURES                         ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; calculate-today-vaccinate: Determines the number of vaccine doses available today.
;;; It uses either continuous (scenario) mode or historical CSV data.
to calculate-today-vaccinate
  ifelse is-scenario? [
    if vax-started? [
      set today-vaccinate-dose vaccination-data
    ]
  ][
    ;; Loop through vaccination dates and check if today's date matches.
    foreach vaccination-dates [
      i ->
      if time:is-equal? date i [
        set today-vaccinate-dose item (position i vaccination-dates) vaccination-data
      ]
    ]
    let last-vax-date last vaccination-dates
    ;; After historical data is exhausted, if continue-vax? is true, use last provided value.
    if time:is-after? date last-vax-date and continue-vax? [
      set today-vaccinate-dose last vaccination-data
    ]
  ]
end

;;; vaccinate: Vaccinates eligible agents. Eligible agents include susceptible,
;;; exposed, or immune (but not vaccinated) persons who are not in quarantine.
to vaccinate
  let daily-vax-sus 0  ;; Counter for susceptible agents vaccinated today
  let daily-vax-exp 0  ;; Counter for exposed agents vaccinated today

  ;; Calculate remaining vaccines after accounting for vaccines given during quarantine.
  let remaining_vaccines max list 0 (today-vaccinate-dose - quar-vax)

  ;; Define eligible persons for vaccination.
  let eligible-persons persons with [
    (epi-state = susceptible-code or epi-state = exposed-code or epi-state = immune-code)
      and not vaccinated? and not quarantine?
  ]

  let total-eligible count eligible-persons

  if total-eligible > 0 [
    ;; Separate eligible agents into susceptible/immune and exposed groups.
    let susceptible-eligible eligible-persons with [(epi-state = susceptible-code or epi-state = immune-code)]
    let exposed-eligible eligible-persons with [epi-state = exposed-code]

    let susceptible-count count susceptible-eligible
    let exposed-count count exposed-eligible

    ;; Distribute available vaccines proportionally between groups.
    let susceptible-vaccine round ((susceptible-count / total-eligible) * remaining_vaccines)
    let exposed-vaccine round ((exposed-count / total-eligible) * remaining_vaccines)

    ;; Vaccinate susceptible agents.
    if susceptible-count > 0 and susceptible-vaccine > 0 [
      ask n-of min (list susceptible-vaccine susceptible-count) susceptible-eligible [
        set epi-state immune-code
        set vaccinated? true
        set color violet
        set time-vaccinated 0
        set daily-vax-sus daily-vax-sus + 1
      ]
    ]
    ;; Vaccinate exposed agents.
    if exposed-count > 0 and exposed-vaccine > 0 [
      ask n-of min (list exposed-vaccine exposed-count) exposed-eligible [
        set vaccinated? true
        set color violet
        set time-vaccinated 0
        set daily-vax-exp daily-vax-exp + 1
        ;; Debugging: Uncomment to print which exposed agent was vaccinated.
        ;; print (word "Vaccinating exposed agent: " self " | Time exposed: " time-exposed)
      ]
    ]
    ;; Update daily vaccine dose counter.
    set daily-vax-dose daily-vax-exp + daily-vax-sus + quar-vax
    ;; Debug: Print total vaccinated agents today if needed.
    ;; print (word "Total vaccinated today: " daily-vax-dose)
  ]
end

;;; update-vaccine-immunity: Increases the counter for vaccinated agents and checks if
;;; the vaccine immune onset period (vax-immune-onset) has passed. If so, it assigns full
;;; immunity based on vaccine efficacy.
to update-vaccine-immunity
  ask persons with [vaccinated?] [
    set time-vaccinated time-vaccinated + 1  ;; Increment counter (in days)
    if time-vaccinated = vax-immune-onset [
      ifelse random-float 1 <= vax-dose-1-eff [
        set epi-state immune-code
        set immunity 1  ;; Full immunity achieved
      ][
        ;; Vaccine failed to confer immunity; immunity remains 0.
      ]
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                        QUARANTINE PROCEDURES                             ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; check-quarantine: Called each day if quarantine is enabled. It checks if today's
;;; date matches any of the quarantine dates provided via CSV or continuous mode. Based
;;; on the date, it determines the number of agents to quarantine.
to check-quarantine
  if not is-scenario? [
    foreach range length quarantine-dates [
      i ->
      let q_date item i quarantine-dates
      if time:is-equal? date q_date [
        let num-to-quarantine item i quarantine-data
        ;; Debug: Print day and number to quarantine.
        ; print (word "Day: " day " Number to quarantine: " num-to-quarantine)
        let infected-agents persons with [ epi-state = infectious-code ]
        ifelse any? infected-agents [
          ;; If there are infectious agents, attempt to quarantine their contacts.
          ask infected-agents [
            apply-quarantine-to-contacts self num-to-quarantine
          ]
        ][
          ;; Otherwise, apply quarantine to a random selection of agents.
          apply-quarantine num-to-quarantine
        ]
      ]
    ]
  ]
end

;;; apply-quarantine-to-contacts: For an infectious agent, this procedure quarantines a
;;; number of contacts based on the provided number (n). It distinguishes between exposed/infectious
;;; and susceptible/immune contacts and attempts to quarantine a proportion based on quar-eff.
to apply-quarantine-to-contacts [infected-person n]
  ask infected-person [
    let contact-list remove-duplicates contacts  ;; Remove duplicate contacts from the daily list
    let quarantined-count 0

    ;; Filter valid contacts (non-nobody and not already quarantined)
    let valid-contacts filter [agent -> agent != nobody and not quarantine?] contact-list
    let E_I_contacts filter [agent -> [epi-state] of agent = exposed-code or [epi-state] of agent = infectious-code] valid-contacts
    let S_IM_contacts filter [agent -> [epi-state] of agent = susceptible-code or [epi-state] of agent = immune-code] valid-contacts

    ;; Determine the number of exposed contacts to quarantine using a binomial random process.
    let quar_from_exp random-binomial (length E_I_contacts) quar-eff
    ;; Calculate how many susceptible contacts to quarantine, ensuring the total reaches n.
    let quar_from_sus max (list 0 (n - quar_from_exp))
    let selected_susceptible n-of min (list quar_from_sus length S_IM_contacts) S_IM_contacts

    ;; If additional susceptible agents are needed, select them from the overall population.
    let additional_needed max list 0 (n - (quar_from_exp + length selected_susceptible))
    if additional_needed > 0 [
      let additional_susceptible n-of min (list additional_needed count persons with [epi-state = susceptible-code or epi-state = immune-code])
                                          (persons with [epi-state = susceptible-code or epi-state = immune-code])
      set selected_susceptible sentence selected_susceptible additional_susceptible
    ]

    ;; Choose the actual contacts to quarantine.
    let quarantine-exposed n-of quar_from_exp E_I_contacts
    let quarantine-susceptible selected_susceptible

    ;; Apply quarantine to exposed contacts.
    foreach quarantine-exposed [contact ->
      ask contact [
        if not quarantine? [
          set epi-state quarantine-code
          set quarantine? true
          set color orange
          set quarantine-time 0
          set quarantined-count quarantined-count + 1
        ]
      ]
    ]

    ;; Apply quarantine to susceptible contacts.
    foreach quarantine-susceptible [contact ->
      ask contact [
        if not quarantine? [
          set epi-state quarantine-code
          set quarantine? true
          set color orange
          set quarantine-time 0
          set quarantined-count quarantined-count + 1
        ]
      ]
    ]

    ;; Update the daily quarantine counter.
    set daily-quarantine-count daily-quarantine-count + quarantined-count
    ; Debug: Print total quarantined count if needed.
    ; print (word "Total quarantined for this infection: " daily-quarantine-count)
  ]
end

;;; apply-quarantine: When there are no specific contacts (or if infected contacts are not present),
;;; this procedure selects n agents at random from the pool of valid contacts and quarantines them.
to apply-quarantine [n]
  if n > 0 [
    let valid-contacts persons with [
      (epi-state = susceptible-code or epi-state = exposed-code) and not quarantine?
    ]
    let valid-contacts-list sort valid-contacts
    if length valid-contacts-list >= n [
      let quarantined-count 0
      let E_I_contacts filter [agent -> [epi-state] of agent = exposed-code or [epi-state] of agent = infectious-code] valid-contacts-list
      let S_IM_contacts filter [agent -> [epi-state] of agent = susceptible-code or [epi-state] of agent = immune-code] valid-contacts-list

      let quar_from_exp random-binomial (length E_I_contacts) quar-eff
      let quar_from_sus max (list 0 (n - quar_from_exp))
      let selected_susceptible n-of min (list quar_from_sus length S_IM_contacts) S_IM_contacts

      let quarantine-exposed n-of quar_from_exp E_I_contacts
      let quarantine-susceptible selected_susceptible

      foreach quarantine-exposed [contact ->
        ask contact [
          if not quarantine? [
            set epi-state quarantine-code
            set quarantine? true
            set color orange
            set quarantine-time 0
            set quarantined-count quarantined-count + 1
          ]
        ]
      ]
      foreach quarantine-susceptible [contact ->
        ask contact [
          if not quarantine? [
            set epi-state quarantine-code
            set quarantine? true
            set color orange
            set quarantine-time 0
            set quarantined-count quarantined-count + 1
          ]
        ]
      ]
      set daily-quarantine-count daily-quarantine-count + quarantined-count
    ]
  ]
end

;;; update-quarantine-status: Called each day to update agents in quarantine.
;;; It increments their quarantine time and, based on the duration, either
;;; applies vaccination (if not already immune) or releases them (making them immune).
to update-quarantine-status
  ask persons with [quarantine?] [
    set quarantine-time quarantine-time + 1

    if quarantine-time = 2 [
      ifelse serology? or vaccinated? [
        ;; If the agent already has serology positive or is vaccinated, release them.
        set epi-state immune-code
        set quarantine? false
        set color violet
      ][
        ;; Otherwise, if vaccination is enabled, vaccinate the agent while in quarantine.
        if is-vaccination? [
          set vaccinated? true
          set immunity 1
          set color violet
          set quar-vax quar-vax + 1
          set time-vaccinated 100  ;; Large value to simulate immediate immunity onset after quarantine
        ]
      ]
    ]

    if quarantine-time = quarantine-duration [
      ;; After quarantine duration is reached, the agent is marked as recovered/immune.
      set epi-state immune-code
      set quarantine? false
      set color violet
      set immunity 1
    ]
  ]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                     DISEASE TRANSMISSION PROCEDURES                      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; assign-transmission-rate: Sets the base transmission probability for an agent.
to assign-transmission-rate
  set base-transmission-rate transmission-prob
end

;;; assign-latent-period: Assigns a latent period (in days) using a log-normal distribution
;;; based on the behavioral latent period parameter. This period represents the time
;;; between exposure and becoming infectious.
to assign-latent-period
  let mean-log-latent ln(behav-latent-period)
  let sd-latent 0.25
  let latent-period-day round exp (random-normal mean-log-latent sd-latent)
  set latent-period latent-period-day
end

;;; assign-inf-period: Assigns the infectious period (in days) for an agent.
to assign-inf-period
  set infectious-period inf-dur
end

;;; interact: Moves agents that are not in quarantine. Each agent turns randomly
;;; and moves forward by 2 units. This simulates natural movement in the shelter.
to interact
  ask persons [
    if not quarantine? [
      rt random 360
      lt random 360
      forward 2
      ; Alternative movement: setxy random-xcor random-ycor (if random relocation is desired)
    ]
  ]
end

;;; expose: Infectious agents interact with nearby agents. Each infectious agent
;;; performs a number of contacts (up to its daily contact limit) and, if the contacted
;;; agent is susceptible, exposes them (i.e., changes their state to exposed) with a
;;; probability based on the transmission rate.
to expose
  ask persons with [epi-state = infectious-code] [
    if not quarantine? [
      ;; Identify nearby agents (within interaction-radius) for possible contact.
      let nearby-agents sort [self] of other persons in-radius interaction-radius
      repeat num-daily-contacts [
        ;; Ensure each contact happens only once per day.
        let available-persons filter [target -> not member? target daily-contacts] nearby-agents
        if not empty? available-persons [
          let chosen-person one-of available-persons
          ;; Record this interaction in both agents' contact lists.
          set daily-contacts lput chosen-person daily-contacts
          set contacts lput chosen-person contacts
          ask chosen-person [
            if not member? myself contacts [
              set contacts lput myself contacts
              set daily-contacts lput myself daily-contacts
            ]
          ]
          ;; If the contacted agent is susceptible, attempt disease exposure.
          if [epi-state] of chosen-person = susceptible-code [
            exposure self chosen-person
          ]
        ]
      ]
    ]
  ]
end

;;; exposure: Called when an infectious agent interacts with a susceptible agent.
;;; With a probability equal to the transmission rate of the infectious agent,
;;; the susceptible agent becomes exposed (and is colored yellow).
to exposure [infector target]
  ask target [
    if epi-state = susceptible-code [
      let transmission-rate [base-transmission-rate] of infector
      if random-float 1 < transmission-rate [
        set epi-state exposed-code
        set color yellow
        assign-latent-period  ;; Assign a latent period for the new exposed agent
        ;; Track exposed cases by age group.
        if age-group = "1"   [ set exposed-cases-pop-1 exposed-cases-pop-1 + 1 ]
        if age-group = "2"   [ set exposed-cases-pop-2 exposed-cases-pop-2 + 1 ]
        if age-group = "3"   [ set exposed-cases-pop-3 exposed-cases-pop-3 + 1 ]
        if age-group = "18+" [ set exposed-cases-pop-4 exposed-cases-pop-4 + 1 ]
        if age-group = "5"   [ set exposed-cases-pop-5 exposed-cases-pop-5 + 1 ]
      ]
    ]
  ]
end

;;; make-infectious: Checks all exposed agents to see if their latent period has elapsed.
;;; If yes, changes their state to infectious (colored red), assigns them an infectious period,
;;; and updates the new case counter for their age group.
to make-infectious
  ask persons with [epi-state = exposed-code] [
    set time-exposed time-exposed + 1

    ;; Post-Exposure Prophylaxis (PEP): If the agent was vaccinated shortly after exposure,
    ;; there is a chance (vax-dose-1-eff) that the agent recovers immediately.
    if is-PEP? and vaccinated? and ((time-exposed - time-vaccinated) <= 3) [
      if random-float 1 < vax-dose-1-eff [
        set epi-state recovered-code
        set immunity 1
        set color green
      ]
    ]

    if time-exposed >= latent-period [
      set epi-state infectious-code
      set color red
      assign-inf-period
      assign-transmission-rate
      ;; Update daily new case counters based on age group.
      if age-group = "1"   [ set new-cases-pop-1 new-cases-pop-1 + 1 ]
      if age-group = "2"   [ set new-cases-pop-2 new-cases-pop-2 + 1 ]
      if age-group = "3"   [ set new-cases-pop-3 new-cases-pop-3 + 1 ]
      if age-group = "18+" [ set new-cases-pop-4 new-cases-pop-4 + 1 ]
      if age-group = "5"   [ set new-cases-pop-5 new-cases-pop-5 + 1 ]
    ]
  ]
end

;;; recover: Processes infectious agents to determine recovery. If the agent has been
;;; infectious for at least its infectious period, the procedure either places it into
;;; quarantine (if enabled) or recovers it immediately.
to recover
  ask persons with [epi-state = infectious-code] [
    set time-infected time-infected + 1
    ;; If the infectious period is complete, process recovery.
    if time-infected >= infectious-period [
      ifelse is-quarantine? [
        if not isolation? [
          ;; If quarantine is enabled, check historical quarantine data (if available) to
          ;; decide whether to quarantine contacts.
          ifelse not is-scenario? [
            let last-quar-date last quarantine-dates
            if time:is-after? date last-quar-date and continue-quarantine? [
              apply-quarantine-to-contacts-conti self
            ]
          ][
            apply-quarantine-to-contacts-conti self
          ]
          ;; Mark the agent as in isolation/quarantine.
          set epi-state quarantine-code
          set isolation? true
          set isolation-time 0
          set color orange
          set inf-quar inf-quar + 1
          set daily-quarantine-count daily-quarantine-count + inf-quar
          set total-quar total-quar + inf-quar
          set vax-started? true  ;; Start vaccination once the first case is detected.
        ]
        if isolation? [
          set isolation-time isolation-time + 1
        ]
        if isolation-time >= 4 [
          ;; Once isolation time reaches threshold, agent recovers.
          set epi-state recovered-code
          set color green
          set isolation? false
          set immunity 1
          update_recovery_stats_by_age_group
        ]
      ][
        ;; If quarantine is disabled, recover immediately after the infectious period.
        set epi-state recovered-code
        set color green
        set immunity 1
        set vax-started? true
        update_recovery_stats_by_age_group
      ]
    ]
  ]
end

;;; update_recovery_stats_by_age_group: Updates the recovered count for the corresponding
;;; age group when an agent recovers.
to update_recovery_stats_by_age_group
  if age-group = "1"   [ set recovered-pop-1 recovered-pop-1 + 1 ]
  if age-group = "2"   [ set recovered-pop-2 recovered-pop-2 + 1 ]
  if age-group = "3"   [ set recovered-pop-3 recovered-pop-3 + 1 ]
  if age-group = "18+" [ set recovered-pop-4 recovered-pop-4 + 1 ]
  if age-group = "5"   [ set recovered-pop-5 recovered-pop-5 + 1 ]
end

;;; apply-quarantine-to-contacts-conti: Similar to apply-quarantine-to-contacts, this procedure
;;; is used when an infectious agent recovers and its contacts need to be quarantined as part of
;;; a continuous intervention strategy.
to apply-quarantine-to-contacts-conti [infected-person]
  ask infected-person [
    let contact-list remove-duplicates contacts
    let quarantined-count 0
    let max_quarantine 0
    let effective-count 0

    ifelse not is-scenario? [
      set max_quarantine last quarantine-data  ;; Get the maximum number to quarantine from data
    ][
      set max_quarantine quar_val
    ]

    let valid-contacts filter [agent -> agent != nobody] contact-list
    let E_I_contacts filter [agent -> [epi-state] of agent = exposed-code or [epi-state] of agent = infectious-code] valid-contacts
    let S_IM_contacts filter [agent -> [epi-state] of agent = susceptible-code or [epi-state] of agent = immune-code] valid-contacts

    let quar_from_exp random-binomial (length E_I_contacts) quar-eff
    let quar_from_sus max (list 0 (max_quarantine - quar_from_exp))
    let selected_susceptible n-of min (list quar_from_sus length S_IM_contacts) S_IM_contacts
    let additional_needed max list 0 (max_quarantine - (quar_from_exp + length selected_susceptible))
    if additional_needed > 0 [
      let additional_susceptible n-of min (list additional_needed count persons with [epi-state = susceptible-code or epi-state = immune-code])
                                          (persons with [epi-state = susceptible-code or epi-state = immune-code])
      set selected_susceptible sentence selected_susceptible additional_susceptible
    ]
    let quarantine-exposed n-of quar_from_exp E_I_contacts
    let quarantine-susceptible selected_susceptible

    foreach quarantine-exposed [contact ->
      ask contact [
        if not quarantine? [
          set epi-state quarantine-code
          set quarantine? true
          set color orange
          set quarantine-time 0
          set quarantined-count quarantined-count + 1
          set effective-count effective-count + 1
        ]
      ]
    ]
    foreach quarantine-susceptible [contact ->
      ask contact [
        if not quarantine? [
          set epi-state quarantine-code
          set quarantine? true
          set color orange
          set quarantine-time 0
          set quarantined-count quarantined-count + 1
        ]
      ]
    ]
    set daily-quarantine-count daily-quarantine-count + quarantined-count
    set total-quar total-quar + quarantined-count + inf-quar
    set effective-capture-count effective-capture-count + effective-count
  ]
  ; Debug: Uncomment to print a separator for quarantine logging.
  ; print(word "____________________________________________________________")
end

;;; Plotting Procedures: These procedures update the NetLogo plots to display epidemic data.
;;; Each plot is updated with relevant daily values.

to plot-effective-capture
  set-current-plot "Quarantine Effective Capture"
  set-current-plot-pen "capture"
  plot effective-capture-count
  set-current-plot-pen "quar"
  plot daily-quarantine-count
end

to plot-effective-capture-per
  set-current-plot "Quarantine Effective Capture Per"
  set-current-plot-pen "capture-per"
  let capture-per 0
  if daily-quarantine-count > 0 [
    set capture-per ((effective-capture-count / daily-quarantine-count) * 100)
  ]
  plot capture-per
end

;;; plot-epidata: (Optional) Plots overall epidemic data (exposed, infectious, recovered).
to plot-epidata
  let exposed count persons with [epi-state = exposed-code]
  let inf count persons with [epi-state = infectious-code]
  let rec count persons with [epi-state = recovered-code]
  set-current-plot "Epidemic Data"
  set-current-plot-pen "Expose"
  plot exposed
  set-current-plot-pen "Infected"
  plot inf
  set-current-plot-pen "Recover"
  plot rec
end

;;; plot-cases: Plots the number of daily new cases for each population group.
to plot-cases
  set-current-plot "Daily Cases"
  set-current-plot-pen "1"
  plot daily-cases-pop-1
  set-current-plot-pen "2"
  plot daily-cases-pop-2
  set-current-plot-pen "3"
  plot daily-cases-pop-3
  set-current-plot-pen "4"
  plot daily-cases-pop-4
  set-current-plot-pen "5"
  plot daily-cases-pop-5
end

;;; plot-vaccination: Updates the plot with the daily number of vaccinations given.
to plot-vaccination
  set-current-plot "Daily Vaccinations"
  set-current-plot-pen "vax dose"
  plot daily-vax-count
end

;;; plot-quarantine: Updates the plot with the cumulative number of quarantined agents.
to plot-quarantine
  set-current-plot "Quarantine People"
  set-current-plot-pen "Quarantine"
  plot total-quar
end

;;; plot-contacts: Plots both the expected and actual daily contacts for infectious agents.
to plot-contacts
  set-current-plot "Daily contacts"
  set-current-plot-pen "expected daily contacts"
  let total-conc sum [num-daily-contacts] of persons with [not quarantine? and epi-state = infectious-code]
  let total-pop count persons with [not quarantine? and epi-state = infectious-code]
  if total-pop > 0 [
    let exp-conc total-conc / total-pop
    plot exp-conc
  ]
  set-current-plot-pen "daily contacts"
  let daily-con sum [length daily-contacts] of persons with [not quarantine? and epi-state = infectious-code]
  let combine-pop count persons with [not quarantine? and epi-state = infectious-code]
  if combine-pop > 0 [
    let daily-cont daily-con / combine-pop
    plot daily-cont
  ]
end

;;; plot-immune: Plots the number of persons with partial immunity (if any).
to plot-immune
  let immune-people count persons with [immunity != 1]
  set-current-plot "Persons with Immunity"
  set-current-plot-pen "immunity"
  plot immune-people
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                              UTILITY REPORTERS                           ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;; random-binomial: Custom reporter that simulates a binomial random process.
;;; Given n trials and probability p, it returns the number of successes.
to-report random-binomial [n p]
  let successes 0
  repeat n [
    if random-float 1 < p [
      set successes successes + 1
    ]
  ]
  report successes
end

;;; current-quarantine-count: Returns the number of agents currently quarantined or in isolation.
to-report current-quarantine-count
  report count persons with ([quarantine? or isolation?])
end

;;; daily-vax-count: Returns the total number of vaccines administered today.
to-report daily-vax-count
  report daily-vax-dose
end

;;; daily-quar-people: Returns the daily quarantine count.
to-report daily-quar-people
  report daily-quarantine-count
end

;;; Exposed and New Cases Reporters: Return daily exposed and new case counts by age group.
to-report daily-exposed-pop-1
  report exposed-cases-pop-1
end
to-report daily-exposed-pop-2
  report exposed-cases-pop-2
end
to-report daily-exposed-pop-3
  report exposed-cases-pop-3
end
to-report daily-exposed-pop-4
  report exposed-cases-pop-4
end
to-report daily-exposed-pop-5
  report exposed-cases-pop-5
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                    DAILY NEW CASES REPORTERS (BY AGE GROUP)                ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report daily-cases-pop-1
  report new-cases-pop-1
end
to-report daily-cases-pop-2
  report new-cases-pop-2
end
to-report daily-cases-pop-3
  report new-cases-pop-3
end
to-report daily-cases-pop-4
  report new-cases-pop-4
end
to-report daily-cases-pop-5
  report new-cases-pop-5
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                         DAILY RECOVERED REPORTERS                          ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report daily-recovered-pop-1
  report recovered-pop-1
end
to-report daily-recovered-pop-2
  report recovered-pop-2
end
to-report daily-recovered-pop-3
  report recovered-pop-3
end
to-report daily-recovered-pop-4
  report recovered-pop-4
end
to-report daily-recovered-pop-5
  report recovered-pop-5
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                         TOTAL CASES REPORTERS                              ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report total-cases-pop-1
  report pop-1-total-cases
end
to-report total-cases-pop-2
  report pop-2-total-cases
end
to-report total-cases-pop-3
  report pop-3-total-cases
end
to-report total-cases-pop-4
  report pop-4-total-cases
end
to-report total-cases-pop-5
  report pop-5-total-cases
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                       SUSCEPTIBLE COUNT REPORTERS                          ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report susc-pop-1
  report susceptible-pop-1
end
to-report susc-pop-2
  report susceptible-pop-2
end
to-report susc-pop-3
  report susceptible-pop-3
end
to-report susc-pop-4
  report susceptible-pop-4
end
to-report susc-pop-5
  report susceptible-pop-5
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                      DAILY QUARANTINE COUNT REPORTER                       ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
to-report daily-quar-report
  report daily-quarantine-count
end
@#$#@#$#@
GRAPHICS-WINDOW
85
10
499
425
-1
-1
7.961
1
10
1
1
1
0
0
0
1
-25
25
-25
25
0
0
1
ticks
30.0

BUTTON
14
194
77
227
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
14
238
77
271
NIL
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
13
274
76
307
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
518
10
575
63
Days
day
17
1
13

PLOT
78
428
508
632
Daily Cases
NIL
NIL
0.0
10.0
0.0
2.0
true
false
"" ""
PENS
"1" 1.0 0 -5298144 true "" ""
"2" 1.0 0 -7500403 true "" ""
"3" 1.0 0 -12087248 true "" ""
"4" 1.0 0 -955883 true "" ""
"5" 1.0 0 -6459832 true "" ""

MONITOR
964
10
1050
63
Vaccinated
count persons with [vaccinated?]
17
1
13

MONITOR
580
10
655
63
Total Pop
count persons
17
1
13

PLOT
517
634
955
813
Daily Vaccinations
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"vax dose" 1.0 0 -11221820 true "" ""

MONITOR
874
10
960
63
Quarantine
count persons with [quarantine?]
17
1
13

MONITOR
1054
10
1249
63
NIL
date
17
1
13

SWITCH
1280
301
1418
334
continue-vax?
continue-vax?
0
1
-1000

SWITCH
1271
502
1436
535
continue-quarantine?
continue-quarantine?
0
1
-1000

TEXTBOX
1257
281
1447
300
Continue Vaccination After Data?
13
0.0
1

TEXTBOX
1264
481
1459
499
Continue Qurantine After Data?
13
0.0
1

PLOT
81
636
507
812
Quarantine People
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"Quarantine" 1.0 0 -16777216 true "" ""

PLOT
962
631
1386
811
Daily contacts
NIL
NIL
0.0
10.0
0.0
0.0
true
true
"" ""
PENS
"daily contacts" 1.0 0 -3844592 true "" ""
"expected daily contacts" 1.0 0 -7500403 true "" ""

SWITCH
1017
202
1144
235
is-quarantine?
is-quarantine?
0
1
-1000

TEXTBOX
1024
180
1147
198
Enable Quarantine?
13
0.0
1

SWITCH
1018
268
1155
301
is-vaccination?
is-vaccination?
0
1
-1000

TEXTBOX
1026
251
1149
269
Enable Vaccination?
13
0.0
1

INPUTBOX
753
222
867
282
behav-latent-period
8.0
1
0
Number

INPUTBOX
896
217
978
277
quar-eff
0.5
1
0
Number

INPUTBOX
534
215
609
275
pop-1
0.0
1
0
Number

INPUTBOX
534
300
611
360
pop-2
0.0
1
0
Number

TEXTBOX
538
173
611
207
Population
14
0.0
1

TEXTBOX
536
460
640
478
Age-group 18+
13
0.0
1

TEXTBOX
768
202
867
220
Latent Period\t\t
14
0.0
1

TEXTBOX
881
198
996
216
Quarantine Efficacy
13
0.0
1

TEXTBOX
855
65
1075
93
List of Parameters
25
0.0
1

INPUTBOX
535
389
610
449
pop-3
0.0
1
0
Number

INPUTBOX
537
478
612
538
pop-4
1544.0
1
0
Number

INPUTBOX
538
560
613
620
pop-5
0.0
1
0
Number

INPUTBOX
901
299
981
359
inf-dur
5.0
1
0
Number

INPUTBOX
757
389
863
449
vax-immune-onset
10.0
1
0
Number

INPUTBOX
643
303
732
364
initial-inf-pop-2
0.0
1
0
Number

INPUTBOX
644
220
733
280
initial-inf-pop-1
0.0
1
0
Number

INPUTBOX
900
388
987
448
vax-dose-1-eff
0.93
1
0
Number

INPUTBOX
756
481
866
541
per-seroprevalence
0.65
1
0
Number

INPUTBOX
893
134
975
194
sero-immunity
1.0
1
0
Number

INPUTBOX
760
559
870
619
quarantine-duration
14.0
1
0
Number

INPUTBOX
906
479
990
539
inf-dur-index
9.0
1
0
Number

MONITOR
752
10
871
63
Infected  People
count persons with [epi-state = infectious-code]
17
1
13

INPUTBOX
646
389
731
449
initial-inf-pop-3
0.0
1
0
Number

INPUTBOX
647
478
734
538
initial-inf-pop-4
1.0
1
0
Number

INPUTBOX
647
562
736
622
initial-inf-pop-5
0.0
1
0
Number

TEXTBOX
649
179
736
197
Initial Infected
14
0.0
1

TEXTBOX
537
283
687
301
Age-group: 
13
0.0
1

TEXTBOX
536
369
603
387
Age-group:
13
0.0
1

TEXTBOX
540
543
690
561
Age-group:
13
0.0
1

TEXTBOX
536
197
603
215
Age-group:
13
0.0
1

TEXTBOX
646
202
729
220
Age-group:
13
0.0
1

TEXTBOX
647
284
720
302
Age-group:
13
0.0
1

TEXTBOX
649
369
720
387
Age-group:
13
0.0
1

TEXTBOX
646
459
741
477
Age-group:18+
13
0.0
1

TEXTBOX
652
545
724
563
Age-group:
13
0.0
1

TEXTBOX
885
280
987
298
Infection Duration
13
0.0
1

TEXTBOX
888
462
1038
480
Index case Inf duration
13
0.0
1

TEXTBOX
895
112
983
130
Sero-Immunity
13
0.0
1

TEXTBOX
900
368
975
386
Vax Efficacy
13
0.0
1

TEXTBOX
743
370
888
388
Vax Immunity Onset dur
13
0.0
1

TEXTBOX
760
463
868
481
percent with sero
13
0.0
1

MONITOR
657
10
749
63
Total Cases
total-cases-pop-4
17
1
13

INPUTBOX
1186
177
1342
272
vax-raw-dates
2024-07-13,2024-07-14,2024-07-15,2024-07-16,2024-07-17,2024-07-18,2024-07-19,2024-07-20,2024-07-21,2024-07-22,2024-07-23,2024-07-24,2024-07-25,2024-07-26,2024-07-27,2024-07-28,2024-07-29,2024-07-30
1
1
String

INPUTBOX
1196
389
1350
478
quar-raw-dates
2024-07-12,2024-07-13
1
1
String

TEXTBOX
1282
138
1404
160
Vaccination Data\t
15
0.0
1

TEXTBOX
1292
343
1407
362
Quarantine Data\t
15
0.0
1

INPUTBOX
641
109
747
169
start-date
2024-07-05
1
0
String

TEXTBOX
519
131
669
149
Model Start Date:
14
0.0
1

INPUTBOX
1344
180
1503
272
vax-raw-values
89,90,53,53,54,74,74,74,15,15,17,19,19,19,19,18,18,32
1
1
String

TEXTBOX
1199
158
1305
176
Vaccination Dates
13
0.0
1

TEXTBOX
1367
157
1517
175
Vaccination Values
13
0.0
1

INPUTBOX
1357
385
1508
477
quar-raw-values
9,10
1
1
String

TEXTBOX
1218
366
1324
384
Quarantine Dates
13
0.0
1

TEXTBOX
1380
364
1488
382
Quarantine Values
13
0.0
1

INPUTBOX
904
563
1001
623
R0
17.0
1
0
Number

TEXTBOX
941
547
967
565
R0
12
0.0
1

INPUTBOX
754
302
862
362
transmission-prob
0.9
1
0
Number

SWITCH
1029
129
1132
162
is-PEP?
is-PEP?
0
1
-1000

INPUTBOX
1094
393
1150
453
quar_val
10.0
1
0
Number

INPUTBOX
1029
393
1084
453
vax_val
50.0
1
0
Number

SWITCH
1029
357
1149
390
is-scenario?
is-scenario?
1
1
-1000

TEXTBOX
1029
319
1179
353
Vax and quar without historical data\t
14
0.0
1

MONITOR
1392
11
1571
68
Daily Quarantine People
daily-quarantine-count
17
1
14

MONITOR
1259
10
1388
67
Total Quarantine
current-quarantine-count
17
1
14

PLOT
1968
255
2338
445
Quarantine Effective Capture Per
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"capture-per" 1.0 0 -16777216 true "" ""

TEXTBOX
770
108
859
126
Seed for Model
12
0.0
1

INPUTBOX
785
130
845
190
seed
10.0
1
0
Number

TEXTBOX
1006
110
1156
128
Post Exposure Prophylaxis
13
0.0
1

TEXTBOX
1279
86
1429
106
Historical Data
16
0.0
0

MONITOR
1578
13
1666
66
Daily Cases
daily-cases-pop-4
17
1
13

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.4.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="measles_sim" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100"/>
    <metric>day</metric>
    <metric>daily-cases-pop-4</metric>
    <metric>total-cases-pop-4</metric>
    <metric>daily-quarantine-count</metric>
    <metric>total-quar</metric>
    <metric>current-quarantine-count</metric>
    <enumeratedValueSet variable="per-seroprevalence">
      <value value="0.65"/>
      <value value="0.78"/>
      <value value="0.88"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="is-quarantine?">
      <value value="true"/>
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vax_val">
      <value value="0"/>
      <value value="50"/>
      <value value="200"/>
    </enumeratedValueSet>
    <steppedValueSet variable="quar-eff" first="0.1" step="0.1" last="1"/>
    <steppedValueSet variable="seed" first="1" step="1" last="1000"/>
  </experiment>
  <experiment name="real_data_exp" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100"/>
    <metric>day</metric>
    <metric>daily-cases-pop-4</metric>
    <metric>total-cases-pop-4</metric>
    <metric>daily-quarantine-count</metric>
    <metric>total-quar</metric>
    <metric>current-quarantine-count</metric>
    <enumeratedValueSet variable="per-seroprevalence">
      <value value="0.896"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vax-immune-onset">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quar_val">
      <value value="10"/>
    </enumeratedValueSet>
    <steppedValueSet variable="seed" first="1" step="1" last="10000"/>
    <enumeratedValueSet variable="continue-quarantine?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="is-scenario?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vax-raw-dates">
      <value value="&quot;2024-07-13,2024-07-14,2024-07-15,2024-07-16,2024-07-17,2024-07-18,2024-07-19,2024-07-20,2024-07-21,2024-07-22,2024-07-23,2024-07-24,2024-07-25,2024-07-26,2024-07-27,2024-07-28,2024-07-29,2024-07-30&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="transmission-prob">
      <value value="0.9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-inf-pop-4">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sero-immunity">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="continue-vax?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quar-raw-dates">
      <value value="&quot;2024-07-12,2024-07-13&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quarantine-duration">
      <value value="14"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mean-contacts">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="is-PEP?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vax-raw-values">
      <value value="&quot;89,90,53,53,54,74,74,74,15,15,17,19,19,19,19,18,18,32&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="is-quarantine?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vax_val">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="is-vaccination?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-dur">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vax-dose-1-eff">
      <value value="0.93"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="pop-4">
      <value value="1544"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quar-raw-values">
      <value value="&quot;9,10&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behav-latent-period">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="start-date">
      <value value="&quot;2024-07-05&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="quar-eff">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="inf-dur-index">
      <value value="9"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
