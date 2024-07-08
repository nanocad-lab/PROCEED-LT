# Preparation

Replace `YOUR_LIBRARY_NAME` with your library path in the files below:

    detail_delay_and_power_simulation.m
    delay_and_power_partial_simulation_leakage_calculation_fast.m
    detail_delay_and_power_simulation_ng_bndry_fast.m
    detail_delay_and_power_simulation_ng_bndry_fast_IRdrop.m
    initial_measure_simulation_fast.m
    detail_delay_and_power_simulation_starting_point_fast.m
    leakage_calculation_fast.m
    xk_initial_delay_estimation_ng_fast.m


# Execution

Run the following command in MATLAB:

    po1vt_lo_k10_0p1_t3_ws1_01_1vd_bndry_aw2


# Evaluation Mode
- Room Temperature (85℃): `temp_set = 'hi'`
- Low Temperature (-196℃): `temp_set = 'lo'`
- Aging: `aging = 1`
- IR drop: `flag_IR_drop = 1`
- Interconnect model selection: `accurate_wire = 0/1/2`
- Multi-Vdd: `vdSelect = '2vd'`
- Multi-Vth: `vtSelect = '2vt'`
- Activity factor: `acti = 0.1`
- Vth variation: `vt_var_mode = 1`

