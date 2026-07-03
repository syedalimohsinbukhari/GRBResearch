"""Created on Fri May 1 19:08:22 2026"""

import os

from matplotlib import pyplot as plt

from grb_research import update_style
from make_lightcurve import make_lightcurves

update_style()

source_name = "131014A"
t_05 = 0.960
t_95 = 4.160
trigger_time = 403420143.2

path_ = "GRB131014215"
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

# the light curve division roughly follows
# https://ui.adsabs.harvard.edu/abs/2009Sci...323.1688A/abstract
# -0.192, 6.976
times = [t_05, 2.432, t_95]
colors = ["r", "g", "gold", "b", "orange", "cyan"]

for ax_i in ax:
    for index, ((start, end), color) in enumerate(zip(zip(times[:-1], times[1:]), colors)):
        if end != times[-1]:
            ax_i.axvline(end, color=colors[index + 1], ls="--")
        ax_i.axvspan(start, end, color=color, alpha=0.1)

f.subplots_adjust(hspace=0.030)

ax[-1].set_xlim(-0.512, 7.168)

# plt.show()
[plt.savefig(f"GRB{source_name}_lightcurve.{i}", dpi=600) for i in ["png", "pdf"]]
plt.close()

os.chdir(cwd)
