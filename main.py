import os
import sys
from pysat.solvers import Glucose3

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
        raise ValueError("Tên folder phải bắt đầu bằng 'graph' hoặc 'scen'")

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
    return var

def create_var_map(var):
    var_map = {}
    counter = 1
    for i, vals in var.items():
        for v in vals:
            var_map[(i, v)] = counter
            counter += 1
    return var_map

def build_constraints(solver, var, var_map, ctr_file):
    # One-Hot: mỗi đỉnh chọn đúng 1 giá trị
    for i, vals in var.items():
        solver.add_clause([var_map[(i, v)] for v in vals])
        for j in range(len(vals)):
            for k in range(j+1, len(vals)):
                solver.add_clause([-var_map[(i, vals[j])], -var_map[(i, vals[k])]])

    # Ràng buộc khoảng cách
    with open(ctr_file) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            i, j = int(parts[0]), int(parts[1])
            vals_i = var.get(i, [])
            vals_j = var.get(j, [])
            if '>' in parts:
                distance = int(parts[-1])
                for vi in vals_i:
                    for vj in vals_j:
                        if abs(vi - vj) <= distance:
                            solver.add_clause([-var_map[(i, vi)], -var_map[(j, vj)]])
            elif '=' in parts:
                target = int(parts[-1])
                for vi in vals_i:
                    for vj in vals_j:
                        if abs(vi - vj) != target:
                            solver.add_clause([-var_map[(i, vi)], -var_map[(j, vj)]])

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

def verify_solution_simple(assignment, ctr_file):
    if assignment is None:
        return False
    with open(ctr_file) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            i, j = int(parts[0]), int(parts[1])
            if i not in assignment or j not in assignment:
                continue
            vi = assignment[i]
            vj = assignment[j]
            if '>' in parts:
                distance = int(parts[-1])
                if abs(vi - vj) <= distance:
                    return False
            elif '=' in parts:
                value = int(parts[-1])
                if abs(vi - vj) != value:
                    return False
    return True

# -------------------------------
# Main
# -------------------------------
def main():
    if len(sys.argv) < 2:
        print("Sử dụng: python main.py <dataset_folder>")
        return

    # tự động thêm dataset/ vào đầu
    dataset_folder = os.path.join("dataset", sys.argv[1])

    try:
        files = get_file_names(dataset_folder)
    except ValueError as e:
        print(e)
        return

    domain = read_domain(files["domain"])
    var = read_var(files["var"], domain)
    var_map = create_var_map(var)

    solver = Glucose3()
    build_constraints(solver, var, var_map, files["ctr"])
    assignment = solve_and_print(solver, var_map)

    # Kiểm tra kết quả
    if verify_solution_simple(assignment, files["ctr"]):
        print("\nCorrect solution!")
        print("Number of lables used: ", len(set(assignment.values())))
    else:
        print("\nIncorrect solution!")

    solver.delete()

if __name__ == "__main__":
    main()
