"""Created on Sep 20 12:49:29 2025"""

import json
import os
import sys
from datetime import datetime

import src.grb_research.grb_utils as utils
import src.grb_research.safe_good_best as sgb

timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
log_filename = f"cstat_run_{timestamp}.log"

req_dir = ['GRB080916009', 'GRB131014215', 'GRB140206275', 'GRB231129779']
# req_dir.append('GRB131014215GBM')
# req_dir.append('GRB131014215GBM40keV')

orig_stdout, orig_stderr = sys.stdout, sys.stderr
with open(log_filename, "w", buffering=1) as log_file:
    sys.stdout = log_file
    sys.stderr = log_file
    ep_ext = "T90"
    tr_count, tex_count, brk_count = 0, 0, 0
    try:
        res_safe, res_unsafe, res_total = {}, {}, {}
        res_best, res_marginal = {}, {}
        cwd_ = os.getcwd()
        outer_dirs = utils.get_directories_in_current_folder()
        outer_dirs = [x for x in outer_dirs if x in req_dir]
        for out_ in outer_dirs:
            inner_dirs = utils.get_directories_in_current_folder(f"{cwd_}/{out_}")
            for in_ in inner_dirs:
                cw_test = in_.split("__")[0].split("/")[-1]
                if "0" in cw_test:
                    ep_ext = "T90"
                elif "A" in cw_test or "B" in cw_test:
                    ep_ext = f"EX{tex_count}"
                elif "X" in cw_test or "Y" in cw_test or "Z" in cw_test:
                    ep_ext = f"BR{brk_count + 1}"
                else:
                    ep_ext = f"TR{tr_count + 1}"
                print(f"\n[RUN] Started at {timestamp} on directory {cwd_}/{out_}/{in_}\n")
                cwd = f"{cwd_}/{out_}/{in_}"
                candidates = [m + ".fit" for m in sorted(sgb.ALLOWED_MODELS)]
                existing = [f for f in candidates if os.path.exists(os.path.join(cwd, f))]
                mapping = sgb.collect_model_cstat([os.path.join(cwd, f) for f in existing])
                for k, v in mapping.items():
                    print(f"{k}: {v[0]:.4f}/{v[1]}")
                mapping = {k: v[0] for k, v in mapping.items()}

                # Step 1 — find the best single model
                simple_models = ["PL", "CPL", "BAND", "SBPL"]
                try:
                    base_filtered = sgb.filter_models_by_error(
                        c_stats=mapping, folder_path=cwd, candidates=simple_models
                    )
                    if base_filtered:
                        best_simple, best_c = sgb.pick_best_single_model(base_filtered)
                        print(f"Best single model: {best_simple} (cstat={best_c})")
                    else:
                        print("Best single model: none")
                        best_simple = None
                except Exception as e:
                    print(f"Best single model unavailable ({e})")
                    best_simple = None

                # -- use exact match, not substring ------------------------------
                bb_map = {"PL": "PL_BB", "CPL": "CPL_BB", "BAND": "BAND_BB", "SBPL": "SBPL_BB"}

                pl_bb_map = {
                    "PL": None,  # PL has no direct PL+BB, the nearest is CPL_PL_BB
                    "CPL": "CPL_PL_BB",
                    "BAND": "BAND_PL_BB",
                    "SBPL": "SBPL_PL_BB",
                }

                # Step 2 — BB comparison
                best_current = best_simple
                bb_candidate = bb_map.get(best_simple)

                if bb_candidate:
                    bb_candidate_filtered = sgb.filter_models_by_error(
                        c_stats=mapping, folder_path=cwd, candidates=[bb_candidate]
                    )

                    new_candidates = list(bb_candidate_filtered.keys()) + [best_current]

                    if bb_candidate_filtered:
                        try:
                            best_bb, best_bb_c = sgb.pick_best_model(
                                c_stats=mapping,
                                candidates=new_candidates,
                                group_name="+BB",
                                folder_path=cwd,
                                is_separate_group=0,
                            )
                            print(f"Best +BB model: {best_bb} (cstat={best_bb_c})")
                            best_current = best_bb
                        except Exception as e:
                            best_bb = None
                            print(f"Best +BB model unavailable ({e})")
                    else:
                        best_bb = None
                        print(f"Best +BB model: {bb_candidate} rejected by error criterion")
                else:
                    best_bb = None
                    print(f"No +BB candidate defined for {best_simple}")

                # Step 3 — PL+BB comparison
                # compares against a BB model if BB was accepted, else against the simpler model
                is_separate_group = 0 if best_current in bb_map.values() else 1
                pl_bb_candidate = pl_bb_map.get(best_simple)  # always derived from the simpler model base

                if pl_bb_candidate:
                    pl_bb_candidate_filtered = sgb.filter_models_by_error(
                        c_stats=mapping, folder_path=cwd, candidates=[pl_bb_candidate]
                    )

                    new_candidates = list(pl_bb_candidate_filtered.keys()) + [best_current]

                    if pl_bb_candidate_filtered:
                        try:
                            best_pl_bb, best_pl_bb_c = sgb.pick_best_model(
                                c_stats=mapping,
                                candidates=new_candidates,
                                group_name="+PL+BB",
                                folder_path=cwd,
                                is_separate_group=is_separate_group,
                            )
                            print(f"Best +PL+BB model: {best_pl_bb} (cstat={best_pl_bb_c})")
                            best_current = best_pl_bb
                        except Exception as e:
                            print(f"Best +PL+BB model unavailable ({e})")
                            best_pl_bb = None
                    else:
                        best_pl_bb = None
                        print(f"Best +PL+BB model: {pl_bb_candidate} rejected by error criterion")
                else:
                    best_pl_bb = None
                    print(f"No +PL+BB candidate defined for {best_simple}")

                safe_dict = sgb.filter_models_by_error(c_stats=mapping, folder_path=cwd, candidates=sgb.ALLOWED_MODELS)
                safe = list(safe_dict.keys())

                marginal_dict = sgb.filter_models_by_error(
                    c_stats=mapping, folder_path=cwd, candidates=sgb.ALLOWED_MODELS, par_constraint=0.5
                )
                _marginal = list(marginal_dict.keys())
                marginal = list(set(safe) ^ set(_marginal))

                if marginal:
                    new_list = safe + marginal
                    unsafe = [m for m in mapping if m not in new_list]
                else:
                    unsafe = [m for m in mapping if m not in safe]

                good_n = []
                if best_simple:
                    good_n.append(best_simple)
                if best_bb:
                    good_n.append(best_bb)
                if best_pl_bb:
                    good_n.append(best_pl_bb)

                safe.sort()
                unsafe.sort()

                # EXTRACT BEST_CURRENT MODEL FROM SAFE_DICT to BEST_DICT
                if best_current:
                    best_dict = {best_current: safe_dict[best_current]}
                    safe_dict.pop(best_current)

                print(f"BEST models: {sorted(best_dict.keys())}")
                sgb.list_par_err(
                    cwd_=cwd, fit_type=list(best_dict.keys()), string="BEST", result_dict=res_best, ep_ext=ep_ext
                )

                print(f"SAFE models: {sorted(safe_dict.keys())}")
                sgb.list_par_err(
                    cwd_=cwd, fit_type=list(safe_dict.keys()), string="SAFE", result_dict=res_safe, ep_ext=ep_ext
                )
                if marginal:
                    marginal_dict = {k: v for k, v in marginal_dict.items() if k in marginal}
                    print(f"MARGINAL models: {sorted(marginal_dict.keys())}")
                    sgb.list_par_err(
                        cwd_=cwd,
                        fit_type=list(marginal_dict.keys()),
                        string="MARGINAL",
                        result_dict=res_marginal,
                        ep_ext=ep_ext,
                    )
                print(f"UNSAFE models: {sorted(unsafe)}")
                sgb.list_par_err(cwd_=cwd, fit_type=unsafe, string="UNSAFE", result_dict=res_unsafe, ep_ext=ep_ext)
                print(f"[RUN] Finished at {datetime.now().strftime('%Y%m%d_%H%M%S')}")
                for i in [res_best, res_safe, res_marginal, res_unsafe]:
                    res_total = utils.deep_merge(d=res_total, u=i)
                pp = utils.flatten_results(res_total)
                with open("results.json", "w") as f:
                    json.dump(obj=utils.make_json_safe(res_total), fp=f, indent=4)
                if "0" in cw_test:
                    pass
                elif "A" in cw_test or "B" in cw_test:
                    tex_count += 1
                elif "X" in cw_test or "Y" in cw_test or "Z" in cw_test:
                    brk_count += 1
                else:
                    tr_count += 1

            tr_count = 0
            tex_count = 0
            brk_count = 0

    finally:
        sys.stdout = orig_stdout
        sys.stderr = orig_stderr
        print(f"Logs written to {log_filename}")
