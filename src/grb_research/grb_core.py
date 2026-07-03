"""Created on Dec 26 14:40:22 2025"""

import json
from dataclasses import dataclass

from .grb_constants import short_to_long
from .grb_enums import GoodnessOfFit
from .grb_model import Model, ModelSet
from .grb_time import EpisodeTypes, TimeInterval, TimeIntervalSet


@dataclass
class GRB:
    name: str
    intervals: TimeIntervalSet | None = None

    @classmethod
    def from_dictionary(cls, name: str, grb_data: dict) -> "GRB":
        """
        Create a GRB instance from a nested dictionary representation.

        Parameters
        ----------
        cls :
            The GRB class.
        name :
            The name for the created GRB.
        grb_data :
            Mapping where keys are time-interval strings and values are dictionaries mapping model names to model dictionaries.

        Returns
        -------
        GRB
            A `GRB` instance whose `intervals` attribute is a `TimeIntervalSet` populated from the provided keys and whose `models` attribute on
            each `TimeInterval` is set to a `ModelSet` constructed from the nested model dictionaries.

        Raises
        ------
        KeyError
            If an expected interval key referenced during model construction is missing in `data`.
        """
        grb = cls(name=name)
        grb.intervals = grb.__interval_container(grb_data.keys())
        for interval_ in grb.intervals:
            gm = []
            temp_data = grb_data[interval_.to_string()]
            for m_name in temp_data.keys():
                gm.append(Model.from_dictionary(m_name, temp_data[m_name], interval_))
            interval_.models = ModelSet(gm)

        return grb

    def __repr__(self):
        return f"GRB(name={self.name}, intervals={self.intervals})"

    def __str__(self):
        tr_eps = f"TR({self.intervals.n_trs})"
        ex1 = self.intervals.ex0
        ex2 = self.intervals.ex1
        if ex2 is not None:
            ex_eps = f"{ex1.kind}/{tr_eps}/{ex2.kind}"
        else:
            ex_eps = f"{ex1.kind}/{tr_eps}"
        return f"  GRBName: {self.name},\n" f"Intervals: {self.intervals.t90.kind}/{ex_eps}"

    @staticmethod
    def __interval_container(episode_keys):
        time_intervals = []
        for interval_str in episode_keys:
            interval = TimeInterval.from_string(interval_str)
            time_intervals.append(interval)
        return TimeIntervalSet(time_intervals)

    def get_all_best_models(self) -> ModelSet:
        """Get the best model for each interval."""
        m_total = []
        for i in self.intervals:
            for m in i.models:
                if m.status is GoodnessOfFit.BEST:
                    m_total.append(m)
        return ModelSet(m_total)

    def get_model(self, model_name, interval=None, tr_index=None) -> Model | ModelSet:
        all_models = []
        for interval_ in self.intervals:
            if interval and interval_.kind != interval:
                continue

            if tr_index is not None and interval_.kind in [EpisodeTypes.TR, EpisodeTypes.SP]:
                if interval_.index != tr_index:
                    continue

            model = interval_.models.get(model_name)
            if model:
                all_models.append(model)

        return ModelSet(all_models) if len(all_models) > 1 else all_models[0]

    def get_model_count(self, separate=False, interval_type=None):
        int_mask = [i for i in self.intervals if i.kind is interval_type] if interval_type else self.intervals

        if separate:
            m_count_safe, m_count_unsafe = {}, {}
            m_count_marginal, m_count_best = {}, {}
            for interval_ in int_mask:
                for m in interval_.models:
                    if m.status is GoodnessOfFit.SAFE:
                        if m.name not in m_count_safe:
                            m_count_safe[m.name] = 1
                        else:
                            m_count_safe[m.name] += 1
                    if m.status is GoodnessOfFit.UNSAFE:
                        if m.name not in m_count_unsafe:
                            m_count_unsafe[m.name] = 1
                        else:
                            m_count_unsafe[m.name] += 1
                    if m.status is GoodnessOfFit.MARGINAL:
                        if m.name not in m_count_marginal:
                            m_count_marginal[m.name] = 1
                        else:
                            m_count_marginal[m.name] += 1
                    if m.status is GoodnessOfFit.BEST:
                        if m.name not in m_count_best:
                            m_count_best[m.name] = 1
                        else:
                            m_count_best[m.name] += 1
            return {
                "SAFE": dict(m_count_safe),
                "UNSAFE": dict(m_count_unsafe),
                "MARGINAL": dict(m_count_marginal),
                "BEST": dict(m_count_best),
            }
        else:
            m_count = {}
            for interval_ in int_mask:
                for m in interval_.models:
                    if m.name not in m_count:
                        m_count[m.name] = 1
                    else:
                        m_count[m.name] += 1

            return m_count


@dataclass
class GRBCatalog:
    grb_list: list[GRB]

    def __post_init__(self):
        self._grb_list: dict[str, GRB] = {grb.name: grb for grb in self.grb_list}

    def __getitem__(self, key: str) -> GRB:
        return self._grb_list[key]

    def __iter__(self):
        return iter(self._grb_list.values())

    def __len__(self):
        return len(self._grb_list)

    def __repr__(self) -> str:
        if not self._grb_list:
            return "GRBCatalog(empty)"

        lines = ["GRBCatalog("]
        for grb in self._grb_list.values():
            lines.append(f"  {grb.name}")
        lines.append(")")
        return "\n".join(lines)

    def __str__(self):
        return self.__repr__()

    @classmethod
    def from_iterable(cls, grb_list: str | list[str], data: dict, name_mapping: dict) -> "GRBCatalog":
        """Construct a GRBCatalog object from iterable"""
        if isinstance(grb_list, str):
            grb_list = [grb_list]
        grb_ = [name_mapping[i] for i in grb_list]
        eps_ = [data[i] for i in grb_]
        grb_data_ = [GRB.from_dictionary(name=i, grb_data=j) for i, j in zip(grb_, eps_)]
        return GRBCatalog(grb_data_)

    def get_grb(self, name: str) -> GRB | None:
        """Get the GRB from the catalog by name."""
        return self._grb_list.get(name)


def prepare_grbs(grb_list, result_file, name_mapping=short_to_long, get_best=False):
    """Load results, build GRB catalog, and return objects and best-model lists."""
    with open(result_file, "r") as f:
        data = json.load(f)

    grb_list_long = [name_mapping[i] for i in grb_list]
    gc = GRBCatalog.from_iterable(grb_list=grb_list, data=data, name_mapping=name_mapping)
    grb_objs: list[GRB] = [gc.get_grb(name) for name in grb_list_long]
    grb_best = [g.get_all_best_models() for g in grb_objs]
    if get_best:
        return gc, grb_list_long, grb_objs, grb_best
    else:
        return gc, grb_list_long, grb_objs
