###############################################################################
# Imports
###############################################################################
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ieda.ieda import IED
from ieda.MOL.PDB import PDB
from ieda.QM.molden import Molden
from ieda.QM.overlap_matrix import OverlapMatrix

###############################################################################
# Grupos de testes
###############################################################################
# Lista base de códigos
codes = ["1b69", "1f6u", "1g5j", "1hvn", "2kxh", "2l2k", "2luf", "2px9", "2wy2"]

# "1g5j", "2px9"

# --- Lista de arquivos molden ---
tests_molden_zip = [f"{code}/molden.zip" for code in codes]

# --- Lista de arquivos intmat ---
tests_s_zip = [f"./{code}/intmat.zip" for code in codes]

# --- Lista de arquivos PDB ---
tests_pdb = [f"{code}/{code}.pdb" for code in codes]

###############################################################################
# Função: Benchmark de cálculo (CPU vs GPU)
###############################################################################
def benchmark_compute(df_label, moldens, smats, pdbs, n_runs=1):
    """
    Executa o cálculo da matriz IED n_runs vezes para cada sistema
    e retorna média e desvio sem ler repetidamente os arquivos.
    """
    results = []

    for m, s, p in zip(moldens, smats, pdbs):
        system_name = p.split("/")[-1].split(".")[0]
        print(f"⚙️ Processando {system_name} ({df_label})")

        molden = Molden()
        molden.read_file(m)
        smatrix = OverlapMatrix()
        smatrix.read_from_multiwfn(s)
        pdb = PDB(p)
        n_atoms = pdb.data.inid.size

        for device, use_gpu in [("CPU", False), ("GPU", True)]:
            times = []
            for run in range(n_runs):
                print(f"   🔹 {device} | Run {run+1}/{n_runs} ...")
                ieda = IED()
                t0 = time.perf_counter()
                ieda.matrix(
                    pdb=pdb,
                    molden=molden,
                    Smatrix=smatrix,
                    path_out=f"./{system_name}_{df_label}_{device.lower()}_run{run+1}",
                    gpu=use_gpu,
                    write=False
                )
                t_elapsed = time.perf_counter() - t0
                times.append(t_elapsed)

            # Calcula média e desvio
            mean_time = sum(times) / n_runs
            std_time = (sum((t - mean_time)**2 for t in times) / (n_runs-1))**0.5 if n_runs > 1 else 0.0

            results.append({
                "system": system_name,
                "type": df_label,
                "device": device,
                "matrix_time": mean_time,
                "matrix_std": std_time,
                "n_atoms": n_atoms
            })

    return results

###############################################################################
# Função: Plotagem CPU vs GPU + Speedup
###############################################################################

def plot_results_cpu_gpu(df):
    # ==========================
    # Pré-processamento
    # ==========================
    # Filtra CPU e GPU e cria rótulo
    df_compute = df[df["device"].isin(["CPU", "GPU"])].copy()
    df_compute["label"] = df_compute["system"] + " (" + df_compute["n_atoms"].astype(str) + ")"

    # Ordena pelo número de átomos
    df_compute = df_compute.sort_values("n_atoms").reset_index(drop=True)

    # Define ordem dos rótulos
    order_labels = df_compute["label"].unique()

    sns.set(style="whitegrid")

    # ==========================
    # Gráfico: CPU vs GPU (linha)
    # ==========================
    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=df_compute,
        x="label",
        y="matrix_time",
        hue="device",
        marker="o",
        palette="Set1",
        errorbar="sd",
        linewidth=2,
        hue_order=["CPU", "GPU"]
    )
    plt.xticks(rotation=0, ha="center")
    plt.xlabel("Sistema")
    plt.ylabel("Tempo da etapa matrix (s)")
    plt.title("Benchmark: CPU vs GPU (média ± desvio padrão)")
    plt.legend(title="Dispositivo")
    plt.tight_layout()
    plt.savefig("cpu_vs_gpu_tempo_seaborn_error.png", dpi=300)
    plt.close()

    # ==========================
    # Gráfico: Speedup CPU/GPU (barra)
    # ==========================
    pivot = (
        df_compute.groupby(["label", "device"])["matrix_time"]
        .agg(["mean", "std"])
        .reset_index()
    )
    pivot_wide = pivot.pivot(index="label", columns="device", values=["mean", "std"]).reindex(order_labels)

    mean_cpu = pivot_wide["mean"]["CPU"]
    mean_gpu = pivot_wide["mean"]["GPU"]
    std_cpu = pivot_wide["std"]["CPU"]
    std_gpu = pivot_wide["std"]["GPU"]

    speedup_mean = mean_cpu / mean_gpu
    speedup_std = speedup_mean * np.sqrt((std_cpu / mean_cpu) ** 2 + (std_gpu / mean_gpu) ** 2)

    df_speedup = pd.DataFrame({
        "Sistema": order_labels,
        "Ganho": speedup_mean,
        "Desvio": speedup_std
    })

    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x="Sistema", y="Ganho", data=df_speedup, palette="Set2", alpha=0.9, order=order_labels)
    ax.errorbar(
        x=np.arange(len(df_speedup)),
        y=df_speedup["Ganho"],
        yerr=df_speedup["Desvio"],
        fmt='none',
        ecolor='black',
        capsize=4,
        linewidth=1.5
    )
    plt.xticks(rotation=0, ha="center")
    plt.xlabel("Sistema (ordenado por nº de átomos)")
    plt.ylabel("Ganho (CPU / GPU)")
    plt.title("Ganho de desempenho da GPU sobre a CPU (média ± desvio)")
    plt.axhline(1.0, color="gray", linestyle="--", linewidth=1)
    for i, (m, s) in enumerate(zip(df_speedup["Ganho"], df_speedup["Desvio"])):
        plt.text(i, m + s + 0.05, f"{m:.2f}×", ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    plt.savefig("gpu_speedup_bar.png", dpi=300)
    plt.close()

    # ==========================
    # Gráfico: Tempos CPU (barra)
    # ==========================
    df_cpu = df_compute[df_compute["device"] == "CPU"]
    cpu_stats = df_cpu.groupby("label")["matrix_time"].agg(["mean", "std"]).reindex(order_labels).reset_index()
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x="label", y="mean", data=cpu_stats, palette="Blues", alpha=0.9, order=order_labels)
    ax.errorbar(
        x=np.arange(len(cpu_stats)),
        y=cpu_stats["mean"],
        yerr=cpu_stats["std"],
        fmt='none',
        ecolor='black',
        capsize=4,
        linewidth=1.5
    )
    plt.xticks(rotation=0, ha="center")
    plt.xlabel("Sistema (ordenado por nº de átomos)")
    plt.ylabel("Tempo CPU (s)")
    plt.title("Tempos médios da CPU com desvio padrão")
    plt.tight_layout()
    plt.savefig("cpu_times_bar.png", dpi=300)
    plt.close()

    # ==========================
    # Gráfico: Tempos GPU (barra)
    # ==========================
    df_gpu = df_compute[df_compute["device"] == "GPU"]
    gpu_stats = df_gpu.groupby("label")["matrix_time"].agg(["mean", "std"]).reindex(order_labels).reset_index()
    plt.figure(figsize=(10, 6))
    ax = sns.barplot(x="label", y="mean", data=gpu_stats, palette="Greens", alpha=0.9, order=order_labels)
    ax.errorbar(
        x=np.arange(len(gpu_stats)),
        y=gpu_stats["mean"],
        yerr=gpu_stats["std"],
        fmt='none',
        ecolor='black',
        capsize=4,
        linewidth=1.5
    )
    plt.xticks(rotation=0, ha="center")
    plt.xlabel("Sistema (ordenado por nº de átomos)")
    plt.ylabel("Tempo GPU (s)")
    plt.title("Tempos médios da GPU com desvio padrão")
    plt.tight_layout()
    plt.savefig("gpu_times_bar.png", dpi=300)
    plt.close()

    print("✅ Gráficos salvos:")
    print("   • cpu_vs_gpu_tempo_seaborn_error.png")
    print("   • gpu_speedup_bar.png")
    print("   • cpu_times_bar.png")
    print("   • gpu_times_bar.png")

###############################################################################
# Execução principal
###############################################################################
if __name__ == "__main__":
    results = []

    # Benchmark de cálculo

    results += benchmark_compute("caatinga", tests_molden_zip, tests_s_zip, tests_pdb, n_runs=3)

    # Salvar resultados
    df = pd.DataFrame(results)
    df.to_csv("benchmark_results.csv", index=False)
    print("✅ Resultados salvos em benchmark_results.csv")
    df = pd.read_csv("benchmark_results.csv")
    # Gerar gráficos
    print(df)
    plot_results_cpu_gpu(df)
