"""Created on Mon Apr 4 10:07:22 2022"""

import os

from make_lightcurve import make_lightcurves
from matplotlib import pyplot as plt

source_name = "190114C"
t_05 = 0.704
t_95 = 117.056
trigger_time = 569192227.62590

path_ = "GRB190114873"
cwd = os.getcwd()

dat = [f for f in os.listdir(f"{cwd}/{path_}") if f.endswith(".dat")]
dat = [i.split(".")[0] for i in dat]

dat_NaI = [i.split(".")[0] for i in dat if "n" in i]
dat_BGO = list(set(dat) ^ set(dat_NaI))

os.chdir(f"{cwd}/{path_}")

fit = [f for f in os.listdir(f"{cwd}/{path_}") if f.endswith(".fits")][0]
f, ax = make_lightcurves(
    src_name=f"GRB{source_name}",
    met_time=trigger_time,
    start=t_05,
    stop=t_95,
    nai_detector_list=dat_NaI,
    bgo_detector_list=dat_BGO,
    lat_gtrspgen=f"{cwd}/{path_}/{fit}",
)

times = [0.704, 2.304, 3.648, 4.480, 5.120, 15.296, 18.432, 80, 117.056]

# ax[0].axhline(0, color='r', ls='--', zorder=100)
colors = ["r", "g", "gold", "b", "maroon", "cyan", "tab:orange", "tab:blue", "tab:green", "tab:red", "tab:purple"]

for ax_i in ax:
    for index, ((start, end), color) in enumerate(zip(zip(times[:-1], times[1:]), colors)):
        if end != times[-1]:
            ax_i.axvline(end, color=colors[index + 1], ls="--")
        ax_i.axvspan(start, end, color=color, alpha=0.1)

# ax[2].set_ylim(bottom=-1150, top=2350)
f.subplots_adjust(hspace=0.030)

# plt.show()
[plt.savefig(f"GRB{source_name}_lightcurve.{i}", dpi=600) for i in ["png", "pdf"]]
plt.close()

os.chdir(cwd)
