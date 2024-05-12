import matplotlib.pyplot as plt
import numpy as np
from struct import unpack

def paint_layer_uniform(x,y,q,lims: tuple,lvls):
    x0,x1,y0,y1 = lims

    fig,ax=plt.subplots(1,1) # Создадим фигуру
    fig.set_size_inches(5,5)
    ax.set_aspect(1) # Отношение размеров координатных осей
    plt.xticks(rotation = 90) # Повернем метки по x на 90 градусов, чтобы они не перекрывали друг друга

    if(lvls==0): # Если не указано сколько уровней, то не указываем
        c = ax.contour(x,y,q)
    else:
        c = ax.contour(x,y,q,levels = lvls)
    ax.clabel(c, inline=True, fontsize=8)
    # sf = ScalarFormatter()
    # sf.set_powerlimits((0, 0))
    # ax.xaxis.set_major_formatter(sf)
    ax.set_xlabel('r, м')
    ax.set_ylabel('z, м')
    # plt.colorbar(c,extendrect = True, drawedges = True, label = "V, Дж/Кл")
    plt.show()

def solution_rz_read(filename: str):
    # Прочесть сетки по r и z и решение
    # Вернет одномерные массивы r и z и двумерный массив q
    f = open(filename, 'rb')
    sizes = []
    sizes.append(int.from_bytes(f.read(4), 'little'))
    sizes.append(int.from_bytes(f.read(4), 'little'))
    print(sizes)

    r_mesh = unpack(str(sizes[0])+ "d", f.read(8*sizes[0]))
    z_mesh = unpack(str(sizes[1])+ "d", f.read(8*sizes[1]))
    q_size = sizes[0]*sizes[1]
    q = unpack(str(q_size)+"d",f.read(8*q_size))
    # print(r_mesh)
    # print(z_mesh)
    r_mesh = np.array(r_mesh)
    z_mesh = np.array(z_mesh)
    q = np.array(q)
    q = q.reshape(sizes[1],sizes[0])
    print("Q SHAPE: ",q.shape)
    f.close()
    return r_mesh, z_mesh, q

def solution_xyz_read(filename: str) -> tuple: 
    f = open(filename, 'rb')
    sizes = []
    sizes.append(int.from_bytes(f.read(4), 'little'))
    sizes.append(int.from_bytes(f.read(4), 'little'))
    sizes.append(int.from_bytes(f.read(4), 'little'))
    print(sizes)

    x_mesh = unpack(str(sizes[0])+ "d", f.read(8*sizes[0]))
    y_mesh = unpack(str(sizes[1])+ "d", f.read(8*sizes[1]))
    z_mesh = unpack(str(sizes[2])+ "d", f.read(8*sizes[2]))
    q_size = sizes[0]*sizes[1]*sizes[2]
    q = unpack(str(q_size)+"d",f.read(8*q_size))
    # print(x_mesh)
    # print(y_mesh)
    # print(z_mesh)
    x_mesh = np.array(x_mesh)
    y_mesh = np.array(y_mesh)
    z_mesh = np.array(z_mesh)
    q = np.array(q)
    q = q.reshape(sizes[2],sizes[1],sizes[0])
    print("Q SHAPE: ",q.shape)
    f.close()
    return x_mesh,y_mesh,z_mesh,q

# Вернет индекс начала отрезка, который содержит точку num
# x - сортированный по неубыванию массив точек
def low_border(x : np.ndarray, num):
    lb = x.searchsorted(num, side='left')
    rb = x.searchsorted(num, side='right')
    if lb==rb:
        return(lb-1)
    else:
        return(max(0,rb-2))