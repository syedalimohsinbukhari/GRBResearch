"""Created on Mar 29 12:44:19 2026"""

import os
import pathlib
import shutil

HOME = os.environ['HOME']

codes_for_paper = pathlib.Path(f"{HOME}/PycharmProjects/GRBResearchWork/codes-for-paper")
GRBResearch = pathlib.Path(f"{HOME}/PycharmProjects/GRBResearchWork/GRBResearch")

#######################################################################################################################
# SECTION 5
#######################################################################################################################

all_safe_unsafe = codes_for_paper / "all-safe-unsafe/all-safe-unsafe.png"
# all_safe_unsafe__ALT = codes_for_paper / "all-safe-unsafe/all-safe-unsafe__ALT.png"
amati_plot = codes_for_paper / "amati_relationship/amati_relationship.png"
butterfly_all = codes_for_paper / "butterfly_plots/butterfly_all.png"
peak_energy_best__all = codes_for_paper / "model_parameters/peak_energy_best__all.png"

combined = [all_safe_unsafe, butterfly_all, peak_energy_best__all, amati_plot]

section_5 = GRBResearch / "images/section5"

if not os.path.exists(section_5):
    os.makedirs(section_5, exist_ok=True)

for c in combined:
    shutil.copy2(c, section_5)
