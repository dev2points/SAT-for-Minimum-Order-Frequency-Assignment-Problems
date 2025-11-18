import os
import psutil
import sys
from time import time
# from pysat.solvers import Cadical195
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
    # Exactly One
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
                distance = int(parts[4])
                for vi in vals_i:
                    for vj in vals_j:
                        if abs(vi - vj) <= distance:
                            solver.add_clause([-var_map[(i, vi)], -var_map[(j, vj)]])
            elif '=' in parts:
                target = int(parts[4])
                for vi in vals_i:
                    for vj in vals_j:
                        if abs(vi - vj) == target:
                            solver.add_clause([-var_map[(i, vi)], var_map[(j, vj)]])
                            
def create_label_var_map(labels, start_index):
    label_var_map = {}
    current = start_index
    for lb in labels:
        label_var_map[lb] = current
        current += 1
    return label_var_map
    
def build_label_constraints(solver, var_map, label_var_map):
    for (i, v), varnum in var_map.items():
        lb_varnum = label_var_map[v]
        solver.add_clause([-varnum, lb_varnum])

def add_limit_label_constraints(solver, label_var_map, max_labels):
    label_vars = list(label_var_map.values())
    
    atmost_k = PBEnc.leq(lits=label_vars, weights=[1]*len(label_vars),
                         bound=max_labels, encoding = 1) 

    for clause in atmost_k.clauses:
        solver.add_clause(clause)


def solve_and_print(solver, var_map):
    if solver.solve():
        model = solver.get_model()
        assignment = {}
        for (i, v), varnum in var_map.items():
            if model[varnum-1] > 0:
                assignment[i] = v
        print("Solution:")
        print("{" + ", ".join(f"{v}" for i, v in sorted(assignment.items())) + "}")
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
            if(vi not in var[i]) or (vj not in var[j]):
                return False
            if '>' in parts:
                distance = int(parts[4])
                if abs(vi - vj) <= distance:
                    return False
            elif '=' in parts:
                value = int(parts[4])
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
    # solver = Cadical195()
    build_constraints(solver, var, var_map, files["ctr"])

    assignment = solve_and_print(solver, var_map)
    if assignment is None:
        return
    if verify_solution_simple(assignment, var, files["ctr"]):
        print("Correct solution!")
        num_lables = len(set(assignment.values()))
        print("Number of lables used: ", num_lables)
    else:
        print("Incorrect solution!")
        return
    end_time = time()
    print(f"Time taken: {end_time - start_time:.2f} seconds")
    process = psutil.Process(os.getpid())
    print(f"Memory used: {process.memory_info().rss / 1024**2:.2f} MB")
    lable_var_map = create_label_var_map(domain[0], len(var_map) + 1)
    build_label_constraints(solver, var_map, lable_var_map)

    while True:
        start_time = time()
        solver.delete()
        solver = Glucose4()
        build_constraints(solver, var, var_map, files["ctr"])
        build_label_constraints(solver, var_map, lable_var_map)
        add_limit_label_constraints(solver, lable_var_map,num_lables - 1)
        print(f"\nSolving with at most {num_lables - 1} labels...")
        build_time = time()
        print(f"Build solver complete with {build_time - start_time:.2f} seconds")
        assignment = solve_and_print(solver, var_map)
        if assignment is None:
            break
        if verify_solution_simple(assignment, var, files["ctr"]):
            print("Correct solution!")
            new_num_lables = len(set(assignment.values()))
            print("Number of lables used: ", new_num_lables)
            num_lables = new_num_lables 
            
        else:
            print("Incorrect solution!")
            break

        end_time = time()
        print(f"Time taken: {end_time - start_time:.2f} seconds")
        process = psutil.Process(os.getpid())
        print(f"Memory used: {process.memory_info().rss / 1024**2:.2f} MB")

    solver.delete()

if __name__ == "__main__":
    main()
