import os
import psutil
import sys
from time import time
from pysat.solvers import Glucose4
from pysat.pb import PBEnc

def get_file_names(dataset_folder):
    base = os.path.basename(dataset_folder)
    if base.lower().startswith("graph"):
        return {
            "domain": os.path.join(dataset_folder, "dom.txt"),
            "var": os.path.join(dataset_folder, "var.txt"),
            "ctr": os.path.join(dataset_folder, "ctr.txt")
        }
    elif base.lower().startswith("scen"):
        return {
            "domain": os.path.join(dataset_folder, "DOM.TXT"),
            "var": os.path.join(dataset_folder, "VAR.TXT"),
            "ctr": os.path.join(dataset_folder, "CTR.TXT")
        }
    else:
        raise ValueError("Not a valid dataset: " + dataset_folder)

def read_domain(file):
    domain = []
    with open(file) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            values = list(map(int, parts[2:]))
            domain.append(values)
    return domain 

def read_var(file, domain):
    var = {}
    with open(file) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            idx = int(parts[0])
            subset_idx = int(parts[1])
            var[idx] = domain[subset_idx]
    return var # domain subset for each variable

def create_var_map(var):
    var_map = {}
    counter = 1
    for i, vals in var.items():
        for v in vals:
            var_map[(i, v)] = counter
            counter += 1
    return var_map # dict mapping (i, v) to variable number

def build_constraints(solver, var, var_map, ctr_file):
    # Exactly One for each variable
    for i, vals in var.items():
        solver.add_clause([var_map[(i, v)] for v in vals])
        for j in range(len(vals)):
            for k in range(j+1, len(vals)):
                solver.add_clause([-var_map[(i, vals[j])], -var_map[(i, vals[k])]])

    # Distance constraints
    with open(ctr_file) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            i, j = int(parts[0]), int(parts[1])
            vals_i = var.get(i, [])
            vals_j = var.get(j, [])
            if '>' in parts:
                # format like: i j > d  (but original used parts[4]; keep same logic)
                # we try to be robust: find last token as the number
                distance = int(parts[-1])
                for vi in vals_i:
                    for vj in vals_j:
                        if abs(vi - vj) <= distance:
                            solver.add_clause([-var_map[(i, vi)], -var_map[(j, vj)]])
            elif '=' in parts:
                target = int(parts[-1])
                for vi in vals_i:
                    for vj in vals_j:
                        if abs(vi - vj) == target:
                            # note: original code used implication (-var_i v var_j)
                            # keep the same meaning: if i takes vi then j must take vj (one of them)
                            solver.add_clause([-var_map[(i, vi)], var_map[(j, vj)]])

def create_label_var_map(labels, start_index):
    label_var_map = {}
    current = start_index
    for lb in labels:
        label_var_map[lb] = current
        current += 1
    return label_var_map
    
def build_label_constraints(solver, var_map, label_var_map):
    # If (i,v) true then label v is used
    for (i, v), varnum in var_map.items():
        lb_varnum = label_var_map[v]
        solver.add_clause([-varnum, lb_varnum])

def add_conditional_atmost_k(solver, label_var_map, K):
    """
    Build PBEnc.leq for label_vars <= K and add clauses to solver as:
      (act -> clause)
    Return act_var (positive int) used as activation literal.
    Implementation detail:
      - We call PBEnc.leq to get encoding clauses (which may introduce aux vars),
      - Compute current max variable number from those clauses,
      - Choose act = max_var + 1,
      - Add clauses [-act] + clause  so that when act is assumed true the clauses must hold.
    """
    label_vars = list(label_var_map.values())
    # Build encoding (this creates clauses possibly using new aux variables)
    enc = PBEnc.leq(lits=label_vars, weights=[1]*len(label_vars), bound=K, encoding=1)
    # Determine max var used in this encoding
    max_var = 0
    for clause in enc.clauses:
        for lit in clause:
            if abs(lit) > max_var:
                max_var = abs(lit)
    # Choose activation var number
    act = max_var + 1
    # Add each clause guarded by -act (act -> clause  <=>  -act v clause)
    for clause in enc.clauses:
        guarded = [-act] + list(clause)
        solver.add_clause(guarded)
    return act

def solve_and_print(solver, var_map):
    # Use solver.get_model() as list of assigned literals; check membership
    if solver.solve():
        model = solver.get_model()
        model_set = set(model)
        assignment = {}
        for (i, v), varnum in var_map.items():
            if varnum in model_set:
                assignment[i] = v
        print("Solution:")
        print("{" + ", ".join(f"{assignment[i]}" for i in sorted(assignment.keys())) + "}")
        return assignment
    else:
        print("Cannot find solution.")
        return None

def verify_solution_simple(assignment, var, ctr_file):
    if assignment is None:
        return False
    with open(ctr_file) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            i, j = int(parts[0]), int(parts[1])
            if i not in assignment or j not in assignment:
                return False
            vi = assignment[i]
            vj = assignment[j]
            if (vi not in var[i]) or (vj not in var[j]):
                return False
            if '>' in parts:
                distance = int(parts[-1])
                if abs(vi - vj) <= distance:
                    return False
            elif '=' in parts:
                value = int(parts[-1])
                if abs(vi - vj) != value:
                    return False
    return True

def main():
    start_time = time()
    if len(sys.argv) < 2:
        print("Use: python main.py <dataset_folder>")
        return

    dataset_folder = os.path.join("dataset", sys.argv[1])

    try:
        files = get_file_names(dataset_folder)
    except ValueError as e:
        print(e)
        return

    domain = read_domain(files["domain"])
    var = read_var(files["var"], domain)
    var_map = create_var_map(var)

    solver = Glucose4()
    build_constraints(solver, var, var_map, files["ctr"])

    # Solve initial problem (no label limit)
    assignment = solve_and_print(solver, var_map)
    if assignment is None:
        solver.delete()
        return
    if verify_solution_simple(assignment, var, files["ctr"]):
        print("Correct solution!")
        num_labels = len(set(assignment.values()))
        print("Number of labels used: ", num_labels)
    else:
        print("Incorrect solution!")
        solver.delete()
        return

    end_time = time()
    print(f"Time for initial solve: {end_time - start_time:.2f} seconds")
    process = psutil.Process(os.getpid())
    print(f"Memory used: {process.memory_info().rss / 1024**2:.2f} MB")

    # Build label vars starting right after var_map max var
    start_index = max(var_map.values()) + 1
    label_var_map = create_label_var_map(domain[0], start_index)
    build_label_constraints(solver, var_map, label_var_map)

    # Now incremental optimization: try to reduce number of labels
    print("\n=== Start incremental optimization of labels ===")
    UB = num_labels
    best_assignment = assignment

    # We'll create for each k (UB-1 down to 1) one conditional encoding (act -> atmost_k)
    # and try to solve with act assumed true. If SAT, update UB and record solution.
    # Note: we add each encoding once. Encoded clauses are inert until their act is assumed.

    acts = {}  # map k -> activation var

    for k in range(UB - 1, 0, -1):
        print(f"\nPreparing <= {k} labels encoding and trying it...")
        act = add_conditional_atmost_k(solver, label_var_map, k)
        acts[k] = act

        # Try solving with the activation assumption (act = True)
        sat = solver.solve(assumptions=[act])
        if sat:
            # fetch and record new assignment
            model = solver.get_model()
            model_set = set(model)
            new_assignment = {}
            for (i, v), varnum in var_map.items():
                if varnum in model_set:
                    new_assignment[i] = v

            if verify_solution_simple(new_assignment, var, files["ctr"]):
                new_num_labels = len(set(new_assignment.values()))
                print(f"SAT: found solution with {new_num_labels} labels (asked <= {k})")
                UB = new_num_labels
                best_assignment = new_assignment
                # We can continue trying smaller k (loop continues)
            else:
                print("SAT but solution invalid on verification (shouldn't happen).")
        else:
            print(f"UNSAT for <= {k} labels (as expected if optimal > {k}).")

        # print resource usage
        process = psutil.Process(os.getpid())
        print(f"Memory used: {process.memory_info().rss / 1024**2:.2f} MB")

    print("\n=== Optimization finished ===")
    print("Best number of labels =", UB)
    print("Best solution (variable->label):")
    print(best_assignment)

    solver.delete()

if __name__ == "__main__":
    main()
