import matplotlib.pyplot as plt
import numpy as np

fig, axs = plt.subplots(2, figsize=(10, 10))
axs[0].set_xlim(10, 20)

axs[1].sharex(axs[0])
plt.show()
