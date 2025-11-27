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

def create_order_var_map(var,var_map, solver):
    counter = solver.nof_vars() + 1
    order_var_map = {}

    for u, labels in var.items():
        for i in labels:
            order_var_map[(u,i)] = counter
            counter += 1

    # Monotonicity constraints 
    for u, labels in var.items():
        #(1)
        first_i = labels[0]
        solver.add_clause([-var_map[(u, first_i)], order_var_map[(u, first_i)]])   # x -> y
        solver.add_clause([-order_var_map[(u, first_i)], var_map[(u, first_i)]])   # y -> x
        for idx in range(len(labels)-1):
            solver.add_clause([-order_var_map[(u, labels[idx])], order_var_map[(u, labels[idx+1])]]) #(4)
        solver.add_clause([order_var_map[(u, labels[-1])]]) #(3)
        #(2)
        for idx, i in enumerate(labels):
            if idx > 0:
                solver.add_clause([-var_map[(u, i)], order_var_map[(u, labels[idx])]])
                solver.add_clause([-var_map[(u, i)], -order_var_map[(u, labels[idx - 1])]])  
                solver.add_clause([-order_var_map[(u, labels[idx])], order_var_map[(u, labels[idx - 1])], var_map[(u, i)]])
                

    

    return order_var_map # dict mapping (u,i) to order variable number

def build_constraints(solver, var, var_map, ctr_file):
    # Exactly One
    for i, vals in var.items():
        lit = [var_map[(i, v)] for v in vals]
        solver.add_clause(lit)
        for j in range(len(vals)):
            for k in range(j+1, len(vals)):
                solver.add_clause([-var_map[(i, vals[j])], -var_map[(i, vals[k])]])

    # Order encoding of distance constraints
    order_var_map = create_order_var_map(var,var_map, solver)

    with open(ctr_file) as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue

            u, v = int(parts[0]), int(parts[1])
            vals_u = var.get(u, [])
            vals_v = var.get(v, [])
            distance = int(parts[4])
            if '=' in parts:
                for iu in vals_u:
                    for jv in vals_v:
                        if abs(iu - jv) == distance:
                            solver.add_clause([-var_map[(u, iu)],  var_map[(v, jv)]])
                            solver.add_clause([-var_map[(v, jv)],  var_map[(u, iu)]])
                
            elif '>' in parts:
                # (5)
                for iu in vals_u:
                    if (iu - distance < vals_v[0] and iu + distance > vals_v[-1]):
                        solver.add_clause([-var_map[(u, iu)]]) #(5)
                    elif (iu - distance < vals_v[0]):
                        for jv in vals_v:
                            if jv - iu >= distance:
                               solver.add_clause([-var_map[(u, iu)], order_var_map[(v, jv)]]) #(6)
                               break
                    elif iu + distance > vals_v[-1]:
                        T = iu + distance - 1
                        # tìm nhãn gần nhất <= T
                        candidates = [t for t in vals_v if t <= T]
                        if candidates:
                            t_high = candidates[-1]
                            solver.add_clause([-var_map[(u, iu)], -order_var_map[(v, t_high)]]) #(7)
                    else : # (8)
                        limit_low  = iu - distance - 1
                        limit_high = iu + distance - 1

                        clause = [-var_map[(u, iu)]]
                        low_candidates = [t for t in vals_v if t <= limit_low]
                        if low_candidates:
                            t_low = low_candidates[-1]
                            clause.append(order_var_map[(v, t_low)])
                        high_candidates = [t for t in vals_v if t <= limit_high]
                        if high_candidates:
                            t_high = high_candidates[-1]
                            clause.append(-order_var_map[(v, t_high)])
                        if len(clause) > 1:
                            solver.add_clause(clause)       

                
                
                #     if (iu - distance < vals_v[0]):
                #         for jv in vals_v:
                #             if jv - iu >= distance:
                #                solver.add_clause([-var_map[(u, iu)], order_var_map[(v, jv)]]) #(6)
                #                break
                #     if iu + distance > vals_v[-1]:
                #         T = iu + distance - 1
                #         # tìm nhãn gần nhất <= T
                #         candidates = [t for t in vals_v if t <= T]
                #         if candidates:
                #             t_high = candidates[-1]
                #             solver.add_clause([-var_map[(u, iu)], -order_var_map[(v, t_high)]]) #(7)
                    # else : # (8)
                    #     limit_low  = iu - distance - 1
                    #     limit_high = iu + distance - 1

                    #     clause = [-var_map[(u, iu)]]

                    #         # ---- Vế trái: y_{v, j-d-1} ----
                    #     low_candidates = [t for t in vals_v if t <= limit_low]
                    #     if low_candidates:
                    #         t_low = low_candidates[-1]
                    #         clause.append(order_var_map[(v, t_low)])

                    #     # ---- Vế phải: ¬y_{v, j+d-1} ----
                    #     high_candidates = [t for t in vals_v if t <= limit_high]
                    #     if high_candidates:
                    #         t_high = high_candidates[-1]
                    #         clause.append(-order_var_map[(v, t_high)])

                    #     # Nếu có ít nhất 1 vế, ta mới thêm mệnh đề
                    #     if len(clause) > 1:
                    #         solver.add_clause(clause)        
                    

                    
                        




            # elif '>' in parts:
                # low_index = 0
                # high_index = 1
                # for i in vals_u:
                #     #print(f"vertice {u} label {i}")
                #     low = i - distance            # min label v cannot take           
                #     high = i + distance          # high label v cannot take
                #     #print(f" vertice..: {v} distance: {distance} low: {low}, high: {high}")
                #     while(low_index < len(vals_v)):
                #         if(vals_v[low_index] >= low and low_index > 0):
                #             # R(u,i) = 1 -> R(v, low) = 1
                #             solver.add_clause([-var_map[(u, i)], order_var_map[(v, vals_v[low_index - 1])]])
                #             #print(f"add clause: ({u},{i}) -> ({v},{vals_v[low_index - 1]})")
                #             break
                #         low_index += 1

                #     if(vals_v[-1] <= high):
                #         solver.add_clause([-var_map[(u, i)], -order_var_map[(v, vals_v[-1])]])
                #         #print(f"add clause: ({u},{i}) -> -({v},{vals_v[-1]})")
                #         continue
                #     while(high_index < len(vals_v)):
                #         if(vals_v[high_index] > high):
                #             # R(u,i) = 1 -> R(v, high - 1) = 0
                #             solver.add_clause([-var_map[(u, i)], -order_var_map[(v, vals_v[high_index - 1])]])
                #             #print(f"add clause: ({u},{i}) -> -({v},{vals_v[high_index - 1]})")
                #             break
                #         high_index += 1
                    
                    

                # low_index = 0
                # high_index = 1
                # for i in vals_v:
                #     low = i - distance       # min label u cannot take                  
                #     high = i + distance          # high label u cannot take
                #     while(low_index < len(vals_u)):
                #         if(vals_u[low_index] >= low and low_index > 0):
                #             solver.add_clause([-var_map[(v, i)], order_var_map[(u, vals_u[low_index - 1])]])
                #             break
                #         low_index += 1
                #     if(vals_u[-1] <= high):
                #         solver.add_clause([-var_map[(v, i)], -order_var_map[(u, vals_u[-1])]])
                #         continue
                #     while(high_index < len(vals_u)):
                #         if(vals_u[high_index] > high):
                #             # R(u,i) = 1 -> R(v, high - 1) = 0
                #             solver.add_clause([-var_map[(v, i)], -order_var_map[(u, vals_u[high_index - 1])]])
                #             break
                #         high_index += 1

                


    

    # for i, vals in var.items():
    #     lits = [var_map[(i, v)] for v in vals]
    #     eq1 = PBEnc.equals(lits=lits, bound=1, encoding=1)
    #     for clause in eq1.clauses:
    #         solver.add_clause(clause)

    # for i, vals in var.items():
    #     lits = [var_map[(i, v)] for v in vals]

    #     # at least 1
    #     solver.add_clause(lits)

    #     # at most 1
    #     am1 = PBEnc.leq(lits=lits, weights=[1]*len(lits), bound=1, encoding=1)
    #     for clause in am1.clauses:
    #         solver.add_clause(clause)


    # # Distance constraints
    # with open(ctr_file) as f:
    #     for line in f:
    #         parts = line.strip().split()
    #         if not parts:
    #             continue
    #         i, j = int(parts[0]), int(parts[1])
    #         vals_i = var.get(i, [])
    #         vals_j = var.get(j, [])
    #         if '>' in parts:
    #             distance = int(parts[4])
    #             # limit range
    #             for vi in vals_i:
    #                 allowed = [ var_map[(j,v)] for v in vals_j if v < vi-distance or v > vi+distance ]
    #                 if allowed:
    #                     solver.add_clause([-var_map[(i,vi)]] + allowed)
    #                 else:
    #                     # nếu không có nhãn hợp lệ thì loại bỏ giá trị này
    #                     solver.add_clause([-var_map[(i,vi)]])

    #             # # sequence counter
    #             # for vi in vals_i:
    #             #     allowed = [ var_map[(j,v)] for v in vals_j if v < vi-distance or v > vi+distance ]
    #             #     if allowed:
    #             #         atmost1 = PBEnc.atleast(lits = allowed, bound=1, encoding=1)
    #             #         solver.add_clause([-var_map[(i,vi)]])
    #             #         for clause in atmost1.clauses:
    #             #             solver.add_clause(clause)
    #             #     else:
    #             #         # nếu không có nhãn hợp lệ thì loại bỏ giá trị này
    #             #         solver.add_clause([-var_map[(i,vi)]])
                

    #             # pairwise
    #             for vi in vals_i:
    #                 for vj in vals_j:
    #                     if abs(vi - vj) <= distance:
    #                         solver.add_clause([-var_map[(i, vi)], -var_map[(j, vj)]])
            # elif '=' in parts:
            #     target = int(parts[4])
            #     for vi in vals_i:
            #         for vj in vals_j:
            #             if abs(vi - vj) == target:
            #                 solver.add_clause([-var_map[(i, vi)], var_map[(j, vj)]])
            #                 solver.add_clause([-var_map[(j, vj)], var_map[(i, vi)]])

    
    
                            
def create_label_var_map(labels, start_index):
    label_var_map = {}
    current = start_index
    for lb in labels:
        label_var_map[lb] = current
        current += 1
    return label_var_map
    
# ánh xạ biến active -> biến xác nhận label được sử dụng    
def build_label_constraints(solver, var_map, label_var_map):
    for (i, v), varnum in var_map.items():
        lb_varnum = label_var_map[v]
        solver.add_clause([-varnum, lb_varnum])

def add_limit_label_constraints(solver, label_var_map, max_labels):
    label_vars = list(label_var_map.values())
    
    atmost_k = PBEnc.leq(lits=label_vars, bound=max_labels, encoding = 4) 

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
    lable_var_map = create_label_var_map(domain[0], solver.nof_vars() + 1)
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
        print(f"Build solver complete in {build_time - start_time:.2f} seconds")
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
