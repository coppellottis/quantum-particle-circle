
import pandas as pd
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm

from functions import error


def analyze(args):
    sim_name, N, a, beta = args

    file_name = f"data/{sim_name}/N{N}_a{a:.4f}_beta{beta:.4f}.csv"

    try:
        data = pd.read_csv(file_name)
    except FileNotFoundError:
        raise FileNotFoundError(f"File non trovato: {file_name}")

    tag = f"N={N} a={a:.4f} beta={beta:.4f}"

    row_result = {
        "N": N,
        "a": a,
        "beta": beta
    }

    # 1. Energia termodinamica
    row_result["E"] = data["E"].mean()
    row_result["err_E"], tau_e = error(
        data["E"], label=f"{tag} obs=E"
    )

    # 2. Winding number
    row_result["Q"] = data["Q"].mean()
    row_result["err_Q"], tau_q = error(
        data["Q"], label=f"{tag} obs=Q"
    )

    # 3. Winding number quadrato
    # Inserire qui il calcolo di Q^2 e relativo errore.

    # 4. Calore specifico
    # Inserire qui il calcolo di Cv e relativo errore.

    return row_result


def main():

    sim_name = input("Simulation: ").strip()
    filename = f"data/{sim_name}/metadata.csv"

    if not os.path.exists(filename):
        raise FileNotFoundError(f"Metadata non trovato: {filename}")

    metadata = pd.read_csv(filename)

    # Generazione dei job
    jobs = [
        (
            sim_name,
            int(row["N"]),
            float(row["a"]),
            float(row["beta"])
        )
        for _, row in metadata.iterrows()
    ]

    # Numero di processi
    n_processes = 2

    os.makedirs(f"results/{sim_name}", exist_ok=True)

    results_list = []
    failed_jobs = []

    # Esecuzione parallela
    with ProcessPoolExecutor(
        max_workers=n_processes
    ) as executor:

        futures = {
            executor.submit(analyze, job): job
            for job in jobs
        }

        with tqdm(
            total=len(jobs),
            desc="Analisi completate",
            dynamic_ncols=True
        ) as pbar:

            for future in as_completed(futures):

                job = futures[future]

                try:
                    result = future.result()
                    results_list.append(result)

                except Exception as e:
                    failed_jobs.append((job, str(e)))
                    print(f"\nErrore nel job {job}: {e}")

                finally:
                    pbar.update(1)

    # Salvataggio dei risultati
    if results_list:

        df_results = pd.DataFrame(results_list)

        df_results = df_results.sort_values(
            by=["N", "a", "beta"]
        )

        output_file = f"results/{sim_name}/results.csv"

        df_results.to_csv(output_file, index=False)

        print(f"\nAnalisi completata!")
        print(f"Risultati salvati in: {output_file}")

    else:
        print("\nNessun risultato valido.")

    # Riepilogo degli errori
    if failed_jobs:

        print(f"\nJob falliti: {len(failed_jobs)}/{len(jobs)}")

        for job, error_msg in failed_jobs:
            print(f"  {job}: {error_msg}")

    else:
        print("\nTutti i job completati correttamente.")


if __name__ == "__main__":
    main()
