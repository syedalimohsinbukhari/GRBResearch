"""Created on Apr 23 11:29:24 2026"""

import json

from src.grb_research import GRBCatalog, find_project_root, ModelSet
from src.grb_research.grb_constants import short_to_long

SOURCE_ROOT = find_project_root()
result_file = SOURCE_ROOT / "results.json"
grb_name = ["080916C", "131014A", "140206B", "231129C"]

with open(result_file, "r") as f:
    data = json.load(f)

grb = GRBCatalog.from_iterable(grb_list=grb_name, data=data, name_mapping=short_to_long)

with open("best_names.txt", "w") as f:
    for idx, _ in enumerate(grb):
        f.write(f" GRB{grb_name[idx]} ".center(51, "=") + "\n")
        p: ModelSet = _.get_all_best_models()
        for model_ in p:
            interval_str = model_.interval.to_string()
            params_str = ", ".join(f"{p.name} {p.value:.4g}({p.error:.4g})" for p in model_.parameters)
            f.write(f"{interval_str}: {model_.name}\n{params_str}\n\n")

            if not model_.name.endswith("_BB"):
                bb_name = model_.name + "_BB"
                bb_model = model_.interval.models.get(bb_name)
                if bb_model is not None:
                    bb_params_str = ", ".join(f"{p.name} {p.value:.4g}({p.error:.4g})" for p in bb_model.parameters)
                    unconstrained = [p.name for p in bb_model.parameters if p.is_unconstrained]
                    f.write(f"  +BB: {bb_model.name}\n  {bb_params_str}\n")
                    if unconstrained:
                        f.write(f"  Unconstrained: {', '.join(unconstrained)}\n")
                        unc_error = [f'{p.relative_error:.3%}' for p in bb_model.parameters if p.is_unconstrained]
                        f.write(f"  Error: {', '.join(unc_error)}\n\n")
                    else:
                        f.write(f"  cstat={bb_model.cstat:.4g}, dof={bb_model.dof}\n")
                        f.write(f"  cstat difference={model_.cstat - bb_model.cstat:.4g}\n\n")
