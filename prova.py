import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure(figsize=(10, 10))
axs = fig.subplot_mosaic("""AB""")

fig.subplots_adjust(wspace=0.0)


axs["B"].sharex(axs["A"])
axs["B"].sharey(axs["A"])
axs["A"].tick_params(top=False, labeltop=False,
                                          bottom=True, labelbottom=True,
                                          left=True, labelleft=True,
                                          right=False, labelright=False)
axs["B"].tick_params(top=False, labeltop=False,
                                          bottom=True, labelbottom=True,
                                          left=False, labelleft=False,
                                          right=True, labelright=True)

axs["A"].set(xlim=(10, 20), ylim=(10, 20), aspect=1)
axs["B"].set(xlim=(10, 20), ylim=(10, 20), aspect=1)
for ax in axs.values():
    ax.plot(np.linspace(10, 20, 100), np.linspace(10, 20, 100), color='k')
print(axs["A"].get_xlim())
axs["A"].set_xlabel('pippo')
plt.show()
for ax in axs:
    print(ax)
    print(axs[ax].get_xlim())
    print(axs[ax].get_xlabel())
