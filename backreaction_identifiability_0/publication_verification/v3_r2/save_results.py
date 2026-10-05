import numpy as np, json, sys
# Data from the verified diag.py run (nx=240, dt=5e-4)
data = {
    "t=0.25": [-4.747168e-08, -2.373584e-08, -1.186792e-08, -5.933960e-09, -2.966980e-09, -1.483490e-09],
    "t=0.50": [-5.338286e-06, -2.669143e-06, -1.334572e-06, -6.672858e-07, -3.336429e-07, -1.668214e-07],
    "t=0.75": [-7.667415e-05, -3.833707e-05, -1.916854e-05, -9.584269e-06, -4.792134e-06, -2.396067e-06],
    "t=1.00": [-4.641908e-04, -2.320954e-04, -1.160477e-04, -5.802385e-05, -2.901192e-05, -1.450596e-05],
}
nb_arr = [4, 8, 16, 32, 64, 128]
res = {}
with open("v3_r2_table.csv", "w") as f:
    f.write("t,NB,gamma1_F,NB_x_gamma1\n")
    for t, gams in data.items():
        for nb, g in zip(nb_arr, gams):
            res[f"{t}_nb={nb}"] = {"gamma1_F": g, "NB_gamma1": nb * g}
            f.write(f"{t},{nb},{g:.10e},{nb*g:.10e}\n")
print("saved CSV")
json.dump(res, open("v3_r2_table.json", "w"), indent=1)
print("saved JSON")
