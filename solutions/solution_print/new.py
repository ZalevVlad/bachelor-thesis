from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import numpy as np
from struct import unpack

folder = "./b_results/test_polinome/test4/"
folder_rz = folder


def solution_rz_read(filename: str):
    # Прочесть сетки по r и z и решение
    # Вернет одномерные массивы r и z и двумерный массив q
    f = open(folder_rz + filename, 'rb')
    sizes = []
    sizes.append(int.from_bytes(f.read(4), 'little'))
    sizes.append(int.from_bytes(f.read(4), 'little'))
    print(sizes)

    r_mesh = unpack(str(sizes[0]) + "d", f.read(8*sizes[0]))
    z_mesh = unpack(str(sizes[1]) + "d", f.read(8*sizes[1]))
    q_size = sizes[0]*sizes[1]
    q = unpack(str(q_size)+"d", f.read(8*q_size))
    print(r_mesh)
    print(z_mesh)
    q = np.array(q)
    q = q.reshape(sizes[1], sizes[0])
    print(q.shape)
    f.close()
    return r_mesh, z_mesh, q


def paint_layer_uniform(x, y, q, lims: tuple, lvls):
    x0, x1, y0, y1 = lims

    fig, ax = plt.subplots(1, 1)  # Создадим фигуру
    fig.set_size_inches(5, 5)
    # Отношение размеров координатных осей
    ax.set_aspect(1)
    # Повернем метки по x на 90 градусов, чтобы они не перекрывали друг друга
    plt.xticks(rotation=90)

    if (lvls == 0):  # Если не указано сколько уровней, то не указываем
        c = ax.contourf(x, y, q)
    else:
        c = ax.contourf(x, y, q, levels=lvls)

    # sf = ScalarFormatter()
    # sf.set_powerlimits((0, 0))
    # ax.xaxis.set_major_formatter(sf)

    plt.colorbar(c, extendrect=True, drawedges=True)
    plt.show()


def paint_layer_uniform_3d(x, y, q, lims: tuple, lvls):
    x0, x1, y0, y1 = lims

    fig = plt.figure(figsize=(6, 6))
    ax_3d = fig.add_subplot(projection='3d')
    ax_3d.set_xlabel('x')
    ax_3d.set_ylabel('y')
    ax_3d.set_zlabel('z')

    xgrid, ygrid = np.meshgrid(x, y)
    ax_3d.plot_wireframe(xgrid, ygrid, q)

    plt.show()


# Отрисовка решения двумерной задачи
folder_rz = "./b_results/test_polinome/test4/"
folder_rz = "./b_results/meshs/mesh1/"
r, z, q = solution_rz_read("solution_rz_q")
rad = 0.02
print(q)
paint_layer_uniform(r, z, q, (0, 0, 0, 0), 100)
paint_layer_uniform_3d(r, z, q, (0, 0, 0, 0), 0)
