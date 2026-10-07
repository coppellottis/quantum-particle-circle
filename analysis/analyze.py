import pandas as pd
import os
import numpy as np
from functions import blckbstr_error, error, safe_block_k
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

def analyze(args):
    sim_name, N, a, beta, position = args

    file_name = f"{sim_name}/N{N}_a{a:.4f}_beta{beta:.4f}.csv"

    try:
        data = pd.read_csv("data/" + file_name)
    except FileNotFoundError:
        raise FileNotFoundError(f"File non trovato: {file_name}")

    tag = f"N={N} a={a:.4f} beta={beta:.4f}"

    row_result = {
        "N": N,
        "a": a,
        "beta": beta
    }

    # 1. Energia Termodinamica (E)
    row_result["E"] = data["E"].mean()
    row_result["err_E"], tau_e = error(data["E"], label=f"{tag} obs=E")

    # 2. Winding Number (Q) - Dovrebbe mediare a zero per simmetria
    #row_result["Q"] = data["Q"].mean()
    #row_result["err_Q"], tau_q = error(data["Q"], label=f"{tag} obs=Q")

    # 3. Winding Number Quadrato (Q^2) - Misura la delocalizzazione topologica
    #k_q = safe_block_k(tau_q, len(data["Q"]), label=f"{tag} obs=Q")
    #row_result["Q2"] = (data["Q"]**2).mean()
    #row_result["err_Q2"] = blckbstr_error(data["Q"], lambda x: np.mean(x**2), k_q)

    # 4. Calore Specifico (Cv) derivato dalle fluttuazioni dell'energia
    #k_e = safe_block_k(tau_e, len(data["E"]), label=f"{tag} obs=E")
    #row_result["Cv"] = np.var(data["E"]) * (beta**2)
    #row_result["err_Cv"] = blckbstr_error(data["E"], lambda x: np.var(x) * (beta**2), k_e)

    return row_result

if __name__ == "__main__":
    sim_name = input("Simulation: ").strip()
    filename = f"data/{sim_name}/metadata.csv"

    metadata = pd.read_csv(filename)
    jobs = []

    # Estraiamo N, a, beta leggendo le fasi impostate nel main.c
    for i in range(metadata.shape[0]):
        N = int(metadata["N"][i])
        a = float(metadata["a"][i])
        beta = float(metadata["beta"][i])

        jobs.append((sim_name, N, a, beta))

    n_processes = 4
    
    # Round-robin per allocare dinamicamente le righe del terminale ai worker
    jobs = [(*job, 1 + i % n_processes) for i, job in enumerate(jobs)]

    os.makedirs(f"results/{sim_name}", exist_ok=True)
    results_list = []

    # Generazione concorrente con salvataggio sicuro dei crash
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [executor.submit(analyze, job) for job in jobs]
        
        with tqdm(total=len(jobs), desc="N completati", position=0, dynamic_ncols=True) as outer:
            for fut in as_completed(futures):
                results_list.append(fut.result())
                outer.update(1)

    # Assemblaggio finale e salvataggio
    df_results = pd.DataFrame(results_list)
    df_results = df_results.sort_values(by="N")
    df_results.to_csv(f"results/{sim_name}/results.csv", index=False)
    
    print(f"\nAnalisi completata! File salvato in results/{sim_name}/results.csv")