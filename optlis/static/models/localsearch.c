#include<stdint.h>
#include<limits.h>
#include<stdbool.h>
#include<stdio.h>

struct instance {
  size_t nnodes;
  size_t ntasks;
  size_t nresources;
  int32_t *tasks;
  int32_t *node_durations;
  double *node_risks;
  int32_t *travel_times;
};

struct solution {
  size_t ntasks;
  int32_t *task_list;
  double *task_risks;
  double objective;
  int32_t *start_times;
  int32_t *finish_times;
  int32_t found_at;
  double relaxation_threshold;
};

struct budget {
  int32_t max;
  int32_t consumed;
};

bool can_relocate(size_t i, size_t j, double relaxation_threshold,
              const double *task_risks);
void relocate(size_t i, size_t j, int32_t *task_list, double *task_risks);
double overall_risk(size_t ntasks, const double *task_risks,
                    const int32_t *finish_times);
double earliest_finish_time(const struct instance *inst, struct solution *sol);
bool try_improve_solution(const struct instance *inst, struct solution *sol,
                          struct budget *budget);
void local_search(const struct instance *inst, struct solution *sol,
                  struct budget *budget);

bool try_improve_solution(const struct instance *inst, struct solution *sol,
                          struct budget *budget) {
  double current_objective = sol->objective;

  for (size_t i = 0; i < sol->ntasks; ++i) {
    for (size_t j = i + 1; j < sol->ntasks; ++j) {

      if (budget->consumed >= budget->max)
        return false;

      if (can_relocate(i, j, sol->relaxation_threshold, sol->task_risks)) {
        relocate(i, j, sol->task_list, sol->task_risks);

        earliest_finish_time(inst, sol);
        budget->consumed++;

        if (sol->objective < current_objective) {
          sol->found_at = budget->consumed;
          return true;
        }

        relocate(j, i, sol->task_list, sol->task_risks);
        sol->objective = current_objective;
      }
    }
  }

  return false;
}

void local_search(const struct instance *inst, struct solution *sol,
                  struct budget *budget) {
  // FIXME: The perturbation procedure doesnt update the obj function
  earliest_finish_time(inst, sol);
  while(try_improve_solution(inst, sol, budget));
}

void print_info(const struct instance *inst, const struct solution *sol,
                const struct budget *budget) {
  printf("Instance:\n");
  printf("%ld nodes %ld tasks %ld resources\n", inst->nnodes,
         inst->ntasks, inst->nresources);
  printf("\nSolution:\n");
  printf("%lf objective\n %lf threshold", sol->objective, sol->relaxation_threshold);

  printf("\nBudget:\n");
  printf("%d consumed %d max\n", budget->consumed, budget->max);
}

void relocate(size_t i, size_t j, int32_t *task_list, double *task_risks) {

  if (i == j) {
      return;
  }

  // Relocate for `task_list`
  int32_t task_to_relocate = task_list[i];

  // Shift elements to the right if moving to a higher index
  if (i < j) {
    for (int k = i; k < j; k++) {
      task_list[k] = task_list[k + 1];
    }
  } 
  // Shift elements to the left if moving to a lower index
  else {
    for (int k = i; k > j; k--) {
      task_list[k] = task_list[k - 1];
    }
  }

  task_list[j] = task_to_relocate;

  // Relocate for `task_risks`
  double risk_to_relocate = task_risks[i];

  // Shift elements to the right if moving to a higher index
  if (i < j) {
    for (int k = i; k < j; ++k) {
      task_risks[k] = task_risks[k + 1];
    }
  } 
  // Shift elements to the left if moving to a lower index
  else {
    for (int k = i; k > j; --k) {
      task_risks[k] = task_risks[k - 1];
    }
  }

  task_risks[j] = risk_to_relocate;
}

// FIXME: if the solution is unfeasible, it will return false even for valid swaps
bool can_relocate(size_t i, size_t j, double relaxation_threshold,
                  const double *task_risks) {

  if (relaxation_threshold >= 1)
    return true;

  double lowest_risk = INT_MAX;
  double highest_risk = 0;

  for (size_t l = i; l <= j; ++l) {

    if (task_risks[l] < lowest_risk)
      lowest_risk = task_risks[l];

    if (task_risks[l] > highest_risk)
      highest_risk = task_risks[l];
  }

  return (highest_risk <= lowest_risk + relaxation_threshold);
}

double overall_risk(size_t nnodes, const double *node_risks,
                    const int32_t *finish_times) {
  double value = 0;

  for (size_t i = 0; i < nnodes; ++i)
    value += node_risks[i] * finish_times[i];

  return value;
}

double earliest_finish_time(const struct instance *inst, struct solution *sol) {
  // NOTE: assumes the only depot is node 0 and the max number of resources is 64
  size_t nnodes = inst->nnodes;
  int resources[64] = {0};

  int cycle = 0;
  for (size_t i = 0; i < sol->ntasks; ++i) {

    // Loop through resources to find the one that provides the EFT for task i
    int est = 0; // earliest start time
    int eft = INT_MAX; // earliest finish time
    size_t efr = 0;  // earliest finishing resource
    size_t task_id = sol->task_list[i];

    for (size_t k = 0; k < inst->nresources; ++k) {
      // Gets the current task of this resource
      size_t res_task_id = resources[k];
      // Gets the finish time of the current task
      int finish_at = sol->finish_times[res_task_id];

      // Avoids time travel back to past lol
      if (finish_at < cycle)
        finish_at = cycle;

      int finish_time = finish_at + inst->travel_times[res_task_id*nnodes + task_id]
                        + inst->node_durations[task_id];

      if (finish_time < eft) {
          est = finish_at;
          eft = finish_time;
          efr = k;
      }
    }

    // Allocates resource `k` to task `i`
    sol->start_times[task_id] = est;
    sol->finish_times[task_id] = eft;
    resources[efr] = task_id;

    if (sol->start_times[task_id] > cycle)
      cycle = sol->start_times[task_id];
  }

  sol->objective = overall_risk(nnodes, inst->node_risks, sol->finish_times);
  return sol->objective;
}
