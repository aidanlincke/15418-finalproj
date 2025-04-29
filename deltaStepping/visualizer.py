import numpy as np
import matplotlib.pyplot as plt
import struct


def read_problem(filepath):
    with open(filepath, "rb") as f:
        header = f.read(6 * 4)
        rows, cols, start_row, start_col, end_row, end_col = struct.unpack("iiiiii", header)

        data = f.read(rows * cols * 4)
        map_flat = np.frombuffer(data, dtype=np.float32)
        map = map_flat.reshape((rows, cols))

    return map, (start_row, start_col), (end_row, end_col)


def read_plan(filepath):
    with open(filepath, "rb") as f:
        length_bytes = f.read(4)
        length = struct.unpack("i", length_bytes)[0]

        plan = [
            struct.unpack("ii", f.read(8))
            for _ in range(length)
        ]

    return plan

if __name__ == "__main__":
    map, (start_row, start_col), (end_row, end_col) = read_problem("problems/dc.bin")

    plan = read_plan('plans/dc.bin')
    plan_rows = [r for r, _ in plan]
    plan_cols = [c for _, c in plan]

    plt.figure(figsize=(10, 10))
    plt.imshow(map, cmap='gray', origin='upper')

    # plt.plot(plan_cols, plan_rows, color='blue', linewidth=1, label="Path")

    plt.scatter([start_col], [start_row], color='green', s=50, label="Start")
    plt.scatter([end_col], [end_row], color='red', s=50, label="Goal")

    plt.legend()
    plt.axis("off")
    plt.tight_layout()
    plt.savefig("path_plot.png", dpi=300) 
    plt.show()