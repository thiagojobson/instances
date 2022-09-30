from pathlib import Path
from io import StringIO

import pytest
import numpy as np

from optlis import load_instance, import_solution
from optlis.utils import _write_solution, _write_instance


def test_Graph():
    """Tests the Graph class."""
    inst = load_instance(Path("data/instances/example.dat"))
    assert np.array_equal(inst.depots, np.array([0]))
    assert np.array_equal(inst.tasks, np.array([1, 2, 3, 4, 5, 6, 7, 8]))
    # Same time horizon as defined in the instance file
    assert inst.time_horizon == 56
    assert list(inst.time_periods)[-1] == 55
    # We set G.time_horizon = None, now the time horizon is estimated
    # by a formula (see G.time_periods)
    inst.time_horizon = None
    assert list(inst.time_periods)[-1] == 75


@pytest.mark.parametrize(["d", "precedence"], [
  (0, [(8, 7), (8, 6), (8, 5), (8, 4), (8, 3), (8, 2), (8, 1),
       (7, 6), (7, 5), (7, 4), (7, 3), (7, 2), (7, 1),
       (6, 4), (6, 3), (6, 2), (6, 1),
       (5, 4), (5, 3), (5, 2), (5, 1),
       (4, 3), (4, 2), (4, 1),
       (3, 2), (3, 1),
       (2, 1)]),
  (0.1, [(8, 6), (8, 5), (8, 4), (8, 3), (8, 2), (8, 1),
         (7, 4), (7, 3), (7, 2), (7, 1),
         (6, 3), (6, 2), (6, 1),
         (5, 3), (5, 2), (5, 1),
         (4, 2), (4, 1),
         (3, 1)]),
  (0.3, [(8, 3), (8, 2), (8, 1),
         (7, 2), (7, 1),
         (6, 1),
         (5, 1)]),
  (1, [])])
def test_Graph__precedence(d, precedence):
    """Tests the precedence generation with the relaxation threshold."""
    inst = load_instance(Path("data/instances/example.dat"))
    assert set(inst.precedence(d=d)) == set(precedence)


def test_load_instance(instance_grid3x3_data):
    """Tests the function to load an instance from a file."""
    inst = load_instance(Path("data/instances/example.dat"))
    assert list(inst.nodes(data=True)) == instance_grid3x3_data[0]
    assert set(inst.edges) == set(instance_grid3x3_data[1])


def test_export_instance():
    """Tests the function to export a problem instance to a text file."""
    inst = load_instance(Path("data/instances/example.dat"))
    outfile = StringIO()
    _write_instance(inst, outfile)
    outfile.seek(0)
    assert outfile.read() == "\n".join(("9",
                                        "0 0 0 1 0.0",
                                        "1 1 5 0 0.1",
                                        "2 1 5 0 0.2",
                                        "3 1 5 0 0.3",
                                        "4 1 5 0 0.4",
                                        "5 1 5 0 0.5",
                                        "6 1 5 0 0.5",
                                        "7 1 5 0 0.6",
                                        "8 1 5 0 0.7",
                                        "12",
                                        "0 3",
                                        "0 1",
                                        "1 4",
                                        "1 2",
                                        "2 5",
                                        "3 6",
                                        "3 4",
                                        "4 7",
                                        "4 5",
                                        "5 8",
                                        "6 7",
                                        "7 8",
                                        "55\n"))


def test_import_solution():
    """"Tests the function to import a solution from a text file."""
    assert import_solution("data/solutions/example.sol") == dict(
        # Completion date
        cd_1=56,
        cd_2=50,
        cd_3=42,
        cd_4=36,
        cd_5=30,
        cd_6=22,
        cd_7=16,
        cd_8=10,
        makespan=56,
        overall_risk=85.2,
        sd_1=50,
        sd_2=42,
        sd_3=36,
        sd_4=30,
        sd_5=22,
        sd_6=16,
        sd_7=10,
        sd_8=1,
        y_0_8_1=1,
        y_1_0_56=1,
        y_2_1_50=1,
        y_3_2_42=1,
        y_4_3_36=1,
        y_5_4_30=1,
        y_6_5_22=1,
        y_7_6_16=1,
        y_8_7_10=1,
    )


def test_export_solution():
    """"Tests the function to export a solution to a text file."""
    instance_path = "instance.dat"
    solution = {
        "makespan": 2,
        "overall_risk": 0.5,
        "cd_1": 1,
        "sd_1": 2,
        "y_0_1_1": 1,
        "y_1_0_2": 1,
    }
    outfile = StringIO()
    _write_solution(solution, instance_path, outfile)
    outfile.seek(0)
    # Ignores first line
    outfile.readline()
    # Checks the written variables
    assert outfile.read() == "\n".join(("makespan = 2",
                                        "overall_risk = 0.5",
                                        "cd_1 = 1",
                                        "sd_1 = 2",
                                        "y_0_1_1 = 1",
                                        "y_1_0_2 = 1\n"))
