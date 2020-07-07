
import matplotlib.pyplot as plt

one_core_x = []
one_core_y = []
all_cores_x = []
all_cores_y = []
with open("times1.txt") as f:
  for line in f:
    items = [int(item) for item in line.split(",")]
    one_core_x.append(items[1])
    one_core_y.append(items[0])

with open("times16.txt") as f:
  for line in f:
    items = [int(item) for item in line.split(",")]
    all_cores_x.append(items[1])
    all_cores_y.append(items[0])

plt.plot(one_core_x, one_core_y,"bo", label="one core")
plt.plot(all_cores_x, all_cores_y,"go", label="16 cores")
plt.legend(loc="upper left")
plt.savefig("plot.png")
