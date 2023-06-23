from typing import Any, Dict, Union, Generator, Tuple, TextIO, Optional
from pathlib import Path
from functools import cached_property

import networkx as nx
import numpy as np
import numpy.typing as npt

from optlis.shared import set_product
from optlis.static.problem_data import Instance as StaticInstance


class Instance(StaticInstance):

    EPSILON = 0.01
    NEUTRALIZING_SPEED = 0.3
    CLEANING_SPEED = 0.075

    def __init__(
        self, nodes, risk, degradation_rate, metabolization_rate, initial_concentration
    ):
        super().__init__()
        self.add_nodes_from(nodes)
        self._risk = risk
        self._degradation_rate = degradation_rate
        self._metabolization_rate = metabolization_rate
        self._initial_concentration = initial_concentration

    @property
    def node_resources(self):
        raise DeprecationWarning

    @property
    def resources(self):
        return dict(
            Qn=sum(nx.get_node_attributes(self, "Qn").values()),
            Qc=sum(nx.get_node_attributes(self, "Qc").values()),
        )

    # TODO: rename this
    @cached_property
    def nodes_duration(self):
        duration = np.zeros((len(self.nodes), len(self.time_periods)), dtype=np.int32)
        for i in self.nodes:
            for t in self.time_periods:
                duration[i][t] = self.cleaning_duration(i, t)
        return duration

    @cached_property
    def cleaning_start_times(self):
        """Returns a 2d vector with the latest start times."""
        nnodes, ntime_units = (len(self.nodes), len(self.time_units))
        nodes_duration = np.zeros(shape=(nnodes, ntime_units), dtype=np.int32)
        for i, t in set_product(self.nodes, self.time_units):
            nodes_duration[i][t] = self._cleaning_start_times(i, t)

        return nodes_duration

    def _cleaning_start_times(self, site, time):
        """Returns the latest start time for an op. if it finishes exactly at time t.
        If it is not possible for the task to finish exactly at time t, returns 0.
        """
        for s in self.time_units:
            v = max(self.initial_concentration(site, p) for p in self.products)
            tt = s
            while v > self.EPSILON:
                v -= self.CLEANING_SPEED
                tt += 1
            else:
                if tt == time:
                    return s
        return 0

    @cached_property
    def neutralizing_start_times(self):
        """Returns a 3d vector with the latest start times."""
        nnodes, nproducts, ntime_units = (
            len(self.nodes),
            len(self.products),
            len(self.time_units),
        )
        duration = np.zeros(shape=(nnodes, nproducts, ntime_units), dtype=np.int32)
        for i, p, t in set_product(self.nodes, self.products, self.time_units):
            duration[i][p][t] = self._neutralizing_duration(i, p, t)

        return duration

    def _neutralizing_duration(self, site, product, time):
        """Returns the latest start time for an op. if it finishes exactly at time t."""
        for s in self.time_units:
            v = self.initial_concentration(site, product)
            tt = s
            while v > self.EPSILON:
                v -= v * self.NEUTRALIZING_SPEED
                tt += 1
            else:
                if tt == time:
                    return s
        return 0

    @cached_property
    def risk(self):
        raise DeprecationWarning

    @cached_property
    def products_risk(self):
        return np.array(self._risk, dtype=np.float64)

    @property
    def products(self):
        nproducts = len(self._risk)
        return list(range(nproducts))

    def initial_concentration(self, i, p):
        return self._initial_concentration[i][p]

    @cached_property
    def degradation_rates(self):
        return np.array(self._degradation_rate, dtype=np.float64)

    def degradation_rate(self, p):
        raise DeprecationWarning

    def metabolization_rate(self, p, q):
        raise DeprecationWarning

    @cached_property
    def metabolizing_rates(self):
        nproducts = len(self.products)
        rates = np.zeros(shape=(nproducts, nproducts), dtype=np.float64)

        for (p, q) in set_product(self.products, self.products):
            rates[p][q] = self._metabolization_rate[p][q]

        return rates

    @cached_property
    def time_units(self):
        return np.array(range(101), dtype=np.int32)

    @property
    def time_periods(self):
        raise DeprecationWarning

    def c_struct(self) -> "c_instance":
        nnodes = len(self.nodes)
        ntasks = len(self.tasks)
        nproducts = len(self.products)
        ntime_units = len(self.time_units)

        return c_instance(
            c_size_t(nnodes),
            c_size_t(ntasks),
            c_size_t(nproducts),
            c_size_t(ntime_units),
            np.array(
                [self.resources["Qn"], self.resources["Qc"]], dtype=np.int32
            ).ctypes.data_as(POINTER(c_int32)),
            self.tasks.ctypes.data_as(POINTER(c_int32)),
            self.cleaning_start_times.ctypes.data_as(POINTER(c_int32)),
            self.neutralizing_start_times.ctypes.data_as(POINTER(c_int32)),
            self.products_risk.ctypes.data_as(POINTER(c_double)),
            self.degradation_rates.ctypes.data_as(POINTER(c_double)),
            self.metabolizing_rates.ctypes.data_as(POINTER(c_double)),
        )


def load_instance(path):
    """Loads an instance from a file."""
    nodes = []
    risk = []
    degradation_rate = []
    metabolization_rate = {}
    initial_concentration = {}

    with open(path, "r") as f:
        lines = f.readlines()

    assert lines[0].startswith("# format: dynamic")
    instance_data = (l for l in lines if not l.startswith("#"))

    nproducts = int(next(instance_data))

    # Parses products' risk
    for _ in range(nproducts):
        line = next(instance_data)
        id_, risk_ = line.split()
        risk.append(float(risk_))

    # Parses products' degradation rate
    for _ in range(nproducts):
        line = next(instance_data)
        id_, degradation_rate_ = line.split()
        degradation_rate.append(float(degradation_rate_))

    # Parses products'
    for _ in range(nproducts):
        line = next(instance_data)
        id_, *metabolization_rate_ = line.split()
        metabolization_rate[int(id_)] = tuple(float(r) for r in metabolization_rate_)

    nnodes = int(next(instance_data))
    for _ in range(nnodes):
        line = next(instance_data)
        nid, ntype, Qn, Qc, D = line.split()
        nodes.append(
            (
                int(nid),
                {
                    "type": int(ntype),
                    "Qn": int(Qn),
                    "Qc": int(Qc),
                    "D": int(D),
                },
            )
        )

    nconcentration = int(next(instance_data))
    for _ in range(nconcentration):
        line = next(instance_data)
        id_, *initial_concentration_ = line.split()
        initial_concentration[int(id_)] = tuple(
            float(c) for c in initial_concentration_
        )

    instance = Instance(
        nodes, risk, degradation_rate, metabolization_rate, initial_concentration
    )
    instance.time_horizon = int(next(instance_data))

    return instance


def export_instance(instance: Instance, outfile_path: Union[str, Path]) -> None:
    """Exports a problem instance to a file."""
    with open(outfile_path, "w") as outfile:
        _write_instance(instance, outfile)


def _write_instance(instance: Instance, outfile: TextIO) -> None:
    """Writes a problem instance to a file."""
    outfile.write("# format: dynamic\n")

    # Write product risk
    outfile.write(f"{len(instance.products)}\n")
    for pid in instance.products:
        outfile.write(f"{pid} {instance.products_risk[pid]:.2f}\n")

    # Write product degradation rate
    for pid in instance.products:
        outfile.write(f"{pid} {instance.degradation_rates[pid]:.2f}\n")

    # Write product metabolization matrix
    for pid in instance.products:
        outfile.write(f"{pid}")
        for sid in instance.products:
            outfile.write(f" {instance.metabolizing_rates[pid][sid]:.2f}")
        outfile.write("\n")

    # Write nodes information
    outfile.write(f"{len(instance.nodes)}\n")
    for id, data in instance.nodes(data=True):
        type, Qn, Qc, D = (data["type"], data["Qn"], data["Qc"], data["D"])
        outfile.write(f"{id} {type} {Qn} {Qc} {D}\n")

    # Write initial concentration
    outfile.write(f"{len(instance.nodes)}\n")
    for id in instance.nodes:
        outfile.write(f"{id}")
        for pid in instance.products:
            outfile.write(f" {instance.initial_concentration(id, pid):.2f}")
        outfile.write("\n")

    T = instance.time_units[-1]
    outfile.write(f"{T}\n")


from optlis.dynamic.models.ctypes import (
    c_instance,
    c_int32,
    c_size_t,
    c_double,
    POINTER,
)
