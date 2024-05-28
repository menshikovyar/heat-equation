import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

file_path = 'U_etalon.csv'
U_etalon = np.loadtxt(file_path)

file_path = 'U_x.csv'
U_x = np.loadtxt(file_path)

fig, ax = plt.subplots(figsize=(3, 3))
ax.set_ylim([0, 20])
ax.plot(U_etalon[0])


def update(frame):

    plt.cla()
    ax.set_ylim([0, 20])
    ax.plot(U_etalon[frame], label="U etalon")
    ax.plot(U_x[frame], label="U found")
    ax.legend()

animation = FuncAnimation(
    fig=fig,
    func = update, 
    frames = 96,
    interval = 5
)

plt.show()

