#MKO Left Neuman Up Neuman Down Neuman Right Dirichle
import math
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot  as  plt
Dx = 1
Dy = 1

M = 20
h = (1.0 - 0)/(M)
pi =math.pi

def Iij(x,y):
  return (-1/pi)*(Dy*pi*pi-Dx)*(np.exp(-(x+h))- np.exp(-x))*(np.sin(pi*(y+h)) - np.sin(pi*y))

A = np.zeros(((M)*(M), (M)*(M)) , 'double' )
f = np.zeros((M)*(M), 'double')

for k in range(0, (M)*(M)):
    xi = ( k % (M) )*h
    yj = ( k // (M) )*h
    f[k] = Iij(xi,yj)
    if (k == 0): #Левый нижний
      A[k][k] = Dx+Dy
      A[k][k+1] = -Dx
      A[k][k+M] = -Dy
      f[k] += (Dx/pi)*np.sin(pi*h)
    elif (k == M-1): #Правый нижний
      A[k][k] = 3*Dx+Dy
      A[k][k-1] = -Dx
      A[k][k+M] = -Dy
      f[k] += 2*Dx*np.exp(-1)*(np.sin(pi*(yj+h)) - np.sin(pi*yj))/(h*pi)
    elif (k < M-1): #Нижняя
      A[k][k] = 2*Dx+Dy
      A[k][k-1] = -Dx
      A[k][k+1] = -Dx
      A[k][k+M] = -Dy
    elif (k == (M)*(M) - 1): #Правый верхний
      A[k][k] = 3*Dx+Dy
      A[k][k-1] = -Dx
      A[k][k-M] = -Dy
      f[k] += 2*Dx*np.exp(-1)*(np.sin(pi*(yj+h)) - np.sin(pi*yj))/(h*pi)
    elif (k == (M-1)*M): #Левый верхний
      A[k][k] = Dy+Dx
      A[k][k+1] = -Dx
      A[k][k-M] = -Dy
      f[k] += (Dx/pi)*(np.sin(pi*(yj+h)) - np.sin(pi*(yj)))
    elif (k > (M-1)*M ): #Верхняя
      A[k][k] = 2*Dx+Dy
      A[k][k-1] = -Dx
      A[k][k+1] = -Dx
      A[k][k-M] = -Dy
    elif (k % (M) == 0): #Левая
      A[k][k] = 2*Dy+Dx
      A[k][k+1] = -Dx
      A[k][k+M] = -Dy
      A[k][k-M] = -Dy
      f[k] += (Dx/pi)*(np.sin(pi*(yj+h)) - np.sin(pi*(yj)))
    elif ( (k+1) % (M) == 0): #Правая
      A[k][k] = 2*Dy+3*Dx
      A[k][k-1] = -Dx
      A[k][k+M] = -Dy
      A[k][k-M] = -Dy
      f[k] += 2*Dx*np.exp(-1)*(np.sin(pi*(yj+h)) - np.sin(pi*yj))/(h*pi)
    else: #Внутренние
      A[k][k] = 2*Dx+2*Dy
      A[k][k-1] = -Dx
      A[k][k+1] = -Dx
      A[k][k+M] = -Dy
      A[k][k-M] = -Dy

C = np.linalg.solve(A, f)

C_exact = np.zeros((M)*(M), 'double')
for k in range(0, (M)*(M)):
    xi = (k % (M)) * h + 0.5*h
    yj = (k // (M)) * h + 0.5*h
    C_exact[k] = np.exp(-xi)*np.cos(pi*yj)

print(np.linalg.norm(C - C_exact, np.inf) )

x = np.zeros((M),'float')
for i in range(0, (M)):
  x[i] = i*h + 0.5*h
y = np.zeros((M),'float')
for j in range(0, (M)):
  y[j] = j*h + 0.5*h
C1= np.reshape(C_exact, (M, M))
x1,y1 = np.meshgrid(x,y)
my_map = plt.get_cmap('rainbow')
plt.contourf(x1,y1,C1,cmap=my_map)
plt.colorbar()
plt.show()

C= np.reshape(C, (M, M))
X=np.linspace(0, 1, M)
Y=np.linspace(0, 1, M)
X1,Y1=np.meshgrid(X, Y)
fig = plt.figure()
ax = plt.axes(projection = "3d")
ax.plot_surface(X1,Y1,C, cmap= "viridis")
ax.set_xlabel("x")
ax.set_ylabel("y")
plt.show()

X1,Y1=np.meshgrid(X, Y)
fig = plt.figure()
ax = plt.axes(projection = "3d")
ax.plot_surface(X1,Y1,C1, cmap= "viridis")
ax.set_xlabel("x")
ax.set_ylabel("y")