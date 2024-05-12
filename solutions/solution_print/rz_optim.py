import help_functions as hf
import matplotlib.pyplot as plt
import numpy as np

def low_border(x,n):
    return x.searchsorted(n)
r,z,q = hf.solution_rz_read("./rz_optim/optim/q_rz")
r=r[1:]
q = q[:,1:]
# hf.paint_layer_uniform(r,z,q,(0,0,0,0),0)
q_a = 1/(2*np.pi*r*1e-2)

fig,ax = plt.subplots(figsize=(6,6))
ax.scatter(r,q[0],color='black', marker = 'o', s=2)
ax.scatter(r,q_a, color='red', marker = 'p', s=2)
ax.set_xlim(0,10)
ax.set_ylim(0,10)
plt.show()  

# fig,ax = plt.subplots(figsize=(6,6))
# ax.plot(r,q[0],color='green')
# ax.plot(r,r_a, color='red')
# ax.set_xlim(0,10)
# ax.set_ylim(0,10)
# plt.show()

fig,ax = plt.subplots(figsize=(6,6))
ax.plot(r,abs(q[0]-q_a)/abs(q[0])*100,color='blue')
ax.set_xlim(0,2000)
ax.set_ylim(0,100)
plt.show()

point = low_border(r,1000)
print(r[point], q[0][point],q_a[point])

#Тест rz optim : h0
fig,ax = plt.subplots(figsize=(6,6))
ax.set_xlim(0,10)
ax.set_ylim(0,10)
for i in range(1,6):
    r,z,q = hf.solution_rz_read("./rz_optim/h0/1e-"+str(i))
    r=r[1:]
    q = q[:,1:]
    q_a = 1/(2*np.pi*r*1e-2)
    ax.scatter(r,q[0], marker = 'o', s=2)
    point = low_border(r,1)
    print(r[point], q[0][point],q_a[point], abs(q_a[point]-q[0][point])/abs(q[0][point]))
plt.show()