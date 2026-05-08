import os
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))
import subprocess
import csv
from pathlib import Path
import matplotlib.pyplot as plt

BASE_DIR = Path(".")
DATA_DIR = BASE_DIR / "data" / "Datapoints"
OUTPUT_DIR = BASE_DIR / "outputs"
PLOT_DIR = BASE_DIR / "plots"

OUTPUT_DIR.mkdir(exist_ok=True)
PLOT_DIR.mkdir(exist_ok=True)


def read_output(file_path):
    lines = file_path.read_text().splitlines()
    cost = int(lines[0])
    time_ms = float(lines[3])
    memory_kb = float(lines[4])
    return cost, time_ms, memory_kb


def get_problem_size(input_file):
    # 👉 简单办法：直接调用你的 parser（推荐）
    from sequence_alignment.input_parser import read_problem
    s, t = read_problem(str(input_file))
    return len(s) + len(t)


def main():
    files = sorted(DATA_DIR.glob("*.txt"))

    results = []

    for f in files:
        print(f"Running {f.name}...")

        basic_out = OUTPUT_DIR / f"{f.stem}_basic.txt"
        eff_out = OUTPUT_DIR / f"{f.stem}_eff.txt"

        subprocess.run(["python", "basic.py", str(f), str(basic_out)])
        subprocess.run(["python", "efficient.py", str(f), str(eff_out)])

        _, bt, bm = read_output(basic_out)
        _, et, em = read_output(eff_out)

        size = get_problem_size(f)

        results.append({
            "file": f.name,
            "size": size,
            "basic_time": bt,
            "basic_mem": bm,
            "eff_time": et,
            "eff_mem": em
        })

    # 保存 CSV
    csv_file = OUTPUT_DIR / "results.csv"
    with open(csv_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=results[0].keys())
        writer.writeheader()
        writer.writerows(results)

    print("Saved CSV to", csv_file)

    # 排序
    results.sort(key=lambda x: x["size"])

    x = [r["size"] for r in results]

    # 时间图
    plt.figure()
    plt.plot(x, [r["basic_time"] for r in results], marker='o', label="Basic")
    plt.plot(x, [r["eff_time"] for r in results], marker='o', label="Efficient")
    plt.xlabel("Problem Size (m+n)")
    plt.ylabel("Time (ms)")
    plt.legend()
    plt.savefig(PLOT_DIR / "time.png")
    plt.close()

    # 内存图
    plt.figure()
    plt.plot(x, [r["basic_mem"] for r in results], marker='o', label="Basic")
    plt.plot(x, [r["eff_mem"] for r in results], marker='o', label="Efficient")
    plt.xlabel("Problem Size (m+n)")
    plt.ylabel("Memory (KB)")
    plt.legend()
    plt.savefig(PLOT_DIR / "memory.png")
    plt.close()

    print("Plots saved in plots/")


if __name__ == "__main__":
    main()