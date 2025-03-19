# Preparation
## Add a new process
Replace `YOUR_LIBRARY_NAME` with your library path to add a new process in the files below:
```
detail_delay_and_power_simulation.m
delay_and_power_partial_simulation_leakage_calculation_fast.m
detail_delay_and_power_simulation_ng_bndry_fast.m
detail_delay_and_power_simulation_ng_bndry_fast_IRdrop.m
initial_measure_simulation_fast.m
detail_delay_and_power_simulation_starting_point_fast.m
leakage_calculation_fast.m
xk_initial_delay_estimation_ng_fast.m
```

## Add a new design
In the file `po1vt_lo_k10_0p1_t3_ws1_01_1vd_bndry_aw2.m`, please set the following variables:
```
LDH = dlmread('<ldh_of_your_design');    # Reading the logic path histogram of your design
total_gates = XXX;                       # Total number of gates in your design
total_stage = XX;                        # Number of bins divided in LDH
L = XXX;                                 # Channel length of the transistor (unit: m)
W = XXX;                                 # Channel width of the transistor (unit: m)
wire_len_per_net = XXX;                  # Average wire length per net (unit: um)
C_total = XXX;                           # C total from spef (pF)
total_length = XXX;                      # Total net length of your design
fanout = XXX;                            # Fan-out (raw from spef)
```


# Evaluation Mode
In the file `po1vt_lo_k10_0p1_t3_ws1_01_1vd_bndry_aw2.m`, please set the following variables:
- Room Temperature (85℃): `temp_set = 'hi'`
- Low Temperature (-196℃): `temp_set = 'lo'`
- Aging: `aging = 1`
- IR drop: `flag_IR_drop = 1`
- Interconnect model selection: `accurate_wire = 0/1/2`
- Multi-Vdd: `vdSelect = '2vd'`
- Multi-Vth: `vtSelect = '2vt'`
- Activity factor: `acti = 0.1`
- Vth variation: `vt_var_mode = 1`


# Execution

Run the following command in MATLAB:

    po1vt_lo_k10_0p1_t3_ws1_01_1vd_bndry_aw2

# Paper Link
This project is associated with the paper: _A Comparative Analysis of Low Temperature and Room Temperature Circuit Operation_. The paper is published at https://ieeexplore.ieee.org/document/10787912. This paper utilizes PROCEED-LT to evaluate and optimize VLSI benchmarks to show the benefits of low-temperature circuit operations compared to room temperature.

