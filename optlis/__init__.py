import argparse
from pathlib import Path

from .utils import (
    Instance,
    load_instance,
    export_instance,
    import_solution,
    export_solution,
)
from .generate import generate_instance

# Top level parser
parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(dest="subcommand", help="choose the subcommand")

# Generic solver command line parser (for cplex and ils)
solver_parser = argparse.ArgumentParser(add_help=False)
solver_parser.add_argument("instance-path", type=Path, help="problem instance path")
solver_parser.add_argument(
    "--relaxation",
    type=float,
    default=0.0,
    help="relaxation threshold (in range [0, 1], default 0.0)",
)
solver_parser.add_argument(
    "--perturbation",
    type=float,
    default=0.5,
    help="perturbation strength (in range [0, 1], default 0.5)",
)
solver_parser.add_argument(
    "--tt-off", dest="travel_times", action="store_false", help="disable travel times"
)
solver_parser.add_argument(
    "--log-path", type=Path, help="path to write the execution log"
)
solver_parser.add_argument("--sol-path", type=Path, help="path to write the solution")
solver_parser.add_argument(
    "--dynamic",
    dest="dynamic",
    action="store_true",
    help="set the problem type to dynamic",
)

# CPLEX command parser
cplex_parser = subparsers.add_parser(
    "cplex", parents=[solver_parser], help="solve a problem instance with cplex"
)
cplex_parser.add_argument(
    "--time-limit", type=int, help="maximum time limit for the execution (in seconds)"
)

# ILS command parser
ils_parser = subparsers.add_parser(
    "ils", parents=[solver_parser], help="solve a problem instance with ils"
)
ils_parser.add_argument(
    "--evaluations",
    type=int,
    default=0,
    help="max number of evaluation calls (default |v_d| * 10,000)",
)
ils_parser.add_argument(
    "--runs", type=int, default=35, help="number of repetitions to perform (default 35)"
)
ils_parser.add_argument(
    "--parallel",
    type=int,
    default=4,
    help="number of parallel processes to spawn (default 4)",
)
ils_parser.add_argument(
    "--seed",
    type=int,
    default=0,
    help="seed for the random number generator (default 0)",
)
ils_parser.add_argument(
    "--tuning",
    dest="tuning",
    action="store_true",
    help="activate the tuning mode (disable multiple runs)",
)

# Generate command parser
generate_parser = subparsers.add_parser(
    "generate", help="generate the benchmark instances"
)
generate_parser.add_argument(
    "export-dir", type=Path, help="directory to export instances (must exist)"
)
generate_parser.add_argument(
    "--seed",
    type=int,
    default=0,
    help="seed for the random number generator (default 0)",
)


def main() -> None:
    import optlis.generate
    from optlis.solvers import cplex, ils
    from optlis.dynamic.models import milp

    args = vars(parser.parse_args())

    if args["dynamic"]:
        if args["subcommand"] == "cplex":
            milp.from_command_line(args)
        else:
            raise NotImplementedError
        return

    if args["subcommand"] == "generate":
        optlis.generate.from_command_line(args)

    elif args["subcommand"] == "ils":
        ils.from_command_line(args)

    elif args["subcommand"] == "cplex":
        cplex.from_command_line(args)


if __name__ == "__main__":
    main()
