import tbmodels
import numpy as np
import math
import matplotlib.pyplot as plt

# use all files
model = tbmodels.Model.from_wannier_files(
    hr_file='silicon_hr.dat',
    wsvec_file='silicon_wsvec.dat',
    win_file='silicon.win'
)

k_space = np.linspace(-0.5, 0.5, 1000)
eigvals = []
for k in k_space:
    eigvals += model.eigenval([k, k, 0]),

plt.plot(k_space, eigvals)
plt.xlabel("k")
plt.ylabel("E/eV")
plt.show()