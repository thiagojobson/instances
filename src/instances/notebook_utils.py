import networkx as nx
import matplotlib.pyplot as plt

from instances import load_instance
from instances.utils import import_solution


def y_bar(G, sol={}):
    """Generates the y values (the sum of the risk of not cleaned sites)
       for plotting a solution."""
    dest = G.destinations
    r = nx.get_node_attributes(G, "r")
    cd = {i: sol.get(f"cd_{i}") for i in dest}
    assert None not in cd.values(), "some of the completion dates are None"
    makespan = sol.get("makespan", 0)
    for t in range(1, makespan + 1):
        yield sum(r[i] for i in dest if t < cd[i])


def plot_solutions(instance_path, sol_paths=[],
                   alpha=0.8, labels=[]):
    """Plot solutions sobreposed over a graph of risk per time period."""
    G = load_instance(instance_path)
    print(f"Instance {instance_path}")
    fig, ax = plt.subplots()
    for i, sol_path in enumerate(sol_paths):
        try:
            sol = import_solution(sol_path)
        except FileNotFoundError:
            continue
        x = list(range(sol.get("makespan", 0)))
        y = list(y_bar(G, sol))
        try:
            label = labels[i]
        except IndexError:
            label = f"sol {i}"
        print(f"{label} makespan: {sol.get('makespan')}, area: {sum(y)}")
        ax.fill_between(x, y, alpha=alpha, label=label)
    ax.set(xlabel="time", ylabel="accumulated risk")
    ax.legend(loc='upper right')
    plt.show()


def makespan(sol={}):
    """Returns the makespan of a solution (model 1's objective function)."""
    return sol["makespan"]


def weighted_sum_completion_dates(G, sol={}):
    """Returns the weighted sum of risks * completion times
       (model 2's objective function)."""
    dest = G.destinations
    r = nx.get_node_attributes(G, "r")
    cd = {i: sol.get(f"cd_{i}") for i in dest}
    assert None not in cd.values(), "some of the completion dates are None"
    return sum(r[i] * cd[i] for i in dest)


def overall_risk(G, sol={}):
    """Returns the overall risk generated by a solution."""
    return sum(y_bar(G, sol))


if __name__ == "__main__":
    # Prints some stats for two diff instances
    instance_name = "g_n16_p7_q2_r05"
    G = load_instance(f"data/instances/{instance_name}.dat")
    sol1 = import_solution(f"data/solutions/m1/{instance_name}.sol")
    sol2 = import_solution(f"data/solutions/m2/{instance_name}.sol")

    # The makespan
    mks1 = makespan(sol1)
    mks2 = makespan(sol2)
    print(f"Makespan 1 = {mks1}\nMakespan 2 = {mks2}")

    # The sum of completion dates
    sum_cd1 = weighted_sum_completion_dates(G, sol1)
    sum_cd2 = weighted_sum_completion_dates(G, sol2)
    print(f"\nSum of cd 1 = {sum_cd1}\nSum of cd 2 = {sum_cd2}")

    # The overall risk
    overall_risk1 = sum(y_bar(G, sol1))
    overall_risk2 = sum(y_bar(G, sol2))
    print(f"\nOverall risk 1 = {overall_risk1}\nOverall risk 2 = {overall_risk2}")
