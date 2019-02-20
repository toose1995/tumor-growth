import numpy as np
import math as math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

rhow = 0.01 #day^-1
rhog = 0.005 #day^-1
D = 0.01 #cm^2/day
chat = 4.00e4
R1 = 2.0 #cm
R2 = 4.0 #cm
R3 = 6.0 #cm
R4 = 8.0 #cm
R5 = 10.0 #cm
phi = math.pi/8.0

n = 100 # number r nodes
m = 100 # number theta nodes

deltaR = R5/n #cm
deltaT = (math.pi /2.0)/m #rad
#m = m -1 # eliminate north border ( equals south border )

#we use horizontal numbering
def rowChecker(alpha, n, m):
    rowNumber = m-1-(alpha//n)
    return rowNumber;

def columnChecker(alpha, n, m):
    columnNumber = alpha%n
    return columnNumber;

def colorChecker(rowNumber, columnNumber, deltaR, deltaT, R1, R2, R3, R4, phi):
    r = deltaR*columnNumber;
    theta = deltaT*rowNumber;
    if ((r >= R4) and ((theta <= phi) or (theta >= math.pi/2 - phi))):
        if ((r == R4) or (theta == phi) or (theta == math.pi/2 - phi)):
            rho = (rhow + rhog)/2.0
        else:
            rho = rhog
    
    elif ((r >= R2) and (r <= R3) and (theta >= phi) and (theta <= math.pi/2 - phi)):
        if ((r==R2) or (r==R3) or (theta==phi) or (theta == math.pi/2 - phi)):
            rho = (rhow + rhog)/2.0
        else:
            rho = rhog
    
    elif (r <= R1):
        if (r==R1):
            rho = (rhow + rhog)/2.0
        else:
            rho = rhog
            
    else:
        rho = rhow
    
    return rho;

rc = np.zeros(n)
re = np.zeros(n)
rw = np.zeros(n)
M = np.zeros((n*m,n*m))
S = np.zeros((n*m,n*m))

for i in range(0,n):
    rc[i] = deltaR*i + deltaR/2.0
    re[i] = rc[i] + deltaR/2.0
    rw[i] = rc[i] - deltaR/2.0
    
for alpha in range(0,n*m):
    row = rowChecker(alpha, n, m)
    column = columnChecker(alpha, n, m)
    rho = colorChecker(row, column, deltaR, deltaT, R1, R2, R3, R4, phi)
    
    # M-matrix
    M[alpha,alpha] = rc[column]*deltaR*deltaT
    # S-matrix
    
    #internal nodes
    if ((row > 0) and (row < m - 1) and (column > 0) and (column < n - 1)):
        S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2.0*D*deltaR/(rc[column]*deltaT) - D*(re[column]+rw[column])*deltaT/deltaR #center
        S[alpha, alpha - 1] = D*rw[column]*deltaT/deltaR #west
        S[alpha, alpha + 1] = D*re[column]*deltaT/deltaR #east
        S[alpha, alpha + n] = D*deltaR/(rc[column]*deltaT) #north
        S[alpha, alpha - n] = D*deltaR/(rc[column]*deltaT) #south
    #east boundary 
    elif (column == n-1):
        if (row > 0) and (row < m-1): #not a corner
            S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2*D*deltaR/(rc[column]*deltaT) - D*rw[column]*deltaT/deltaR
            S[alpha, alpha - 1] = D*rw[column]*deltaT/deltaR
            S[alpha, alpha + n] = D*deltaR/(rc[column]*deltaT)
            S[alpha, alpha - n] = D*deltaR/(rc[column]*deltaT)
        elif (row == 0): # top right corner
            S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2*D*deltaR/(rc[column]*deltaT) - D*rw[column]*deltaT/deltaR
            S[alpha, alpha - 1] = D*rw[column]*deltaT/deltaR
            S[alpha, n - 1] = D*deltaR/(rc[column]*deltaT)
            S[alpha, alpha - n] = D*deltaR/(rc[column]*deltaT)
        else: #bottom right corner            
            S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2*D*deltaR/(rc[column]*deltaT) - D*rw[column]*deltaT/deltaR
            S[alpha, alpha - 1] = D*rw[column]*deltaT/deltaR
            S[alpha, alpha + n] = D*deltaR/(rc[column]*deltaT)
            S[alpha, n*m-1] = D*deltaR/(rc[column]*deltaT)
    #west boundary 
    elif (column == 0):
        if (row > 0) and (row < m-1): #not a corner
            S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2.0*D*deltaR/(rc[column]*deltaT) -D*re[column]*deltaT/deltaR
            S[alpha, alpha + 1] = D*re[column]*deltaT/deltaR
            S[alpha, alpha + n] = D*deltaR/(rc[column]*deltaT)
            S[alpha, alpha - n] = D*deltaR/(rc[column]*deltaT)
        elif (row == 0): #top left corner
             S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2*D*deltaR/(rc[column]*deltaT) - D*re[column]*deltaT/deltaR
             S[alpha, alpha + 1] = D*re[column]*deltaT/deltaR
             S[alpha, 0] = D*deltaR/(rc[column]*deltaT)
             S[alpha, alpha - n] = D*deltaR/(rc[column]*deltaT)
        else: #bottom left corner
            S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2*D*deltaR/(rc[column]*deltaT) - D*re[column]*deltaT/deltaR
            S[alpha, alpha + 1] = D*re[column]*deltaT/deltaR
            S[alpha, alpha + n] = D*deltaR/(rc[column]*deltaT)
            S[alpha, (m-1)*n] = D*deltaR/(rc[column]*deltaT)
    #bottom boundary without corners
    elif ((row == m - 1) and (column > 0) and (column < n - 1)):
        S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2*D*deltaR/(rc[column]*deltaT) - D*(re[column]+rw[column])*deltaT/deltaR
        S[alpha, alpha - 1] = D*rw[column]*deltaT/deltaR
        S[alpha, alpha + 1] = D*re[column]*deltaT/deltaR
        S[alpha, alpha + n] = D*deltaR/(rc[column]*deltaT)
        S[alpha, alpha + (m-1)*n] = D*deltaR/(rc[column]*deltaT)
    #top boundary without corners
    elif ((row == 0) and (column > 0) and (column < n - 1)):
        S[alpha, alpha] = rho*rc[column]*deltaR*deltaT - 2*D*deltaR/(rc[column]*deltaT) - D*(re[column]+rw[column])*deltaT/deltaR
        S[alpha, alpha - 1] = D*rw[column]*deltaT/deltaR
        S[alpha, alpha + 1] = D*re[column]*deltaT/deltaR
        S[alpha, alpha - (m-1)*n] = D*deltaR/(rc[column]*deltaT)
        S[alpha, alpha - n] = D*deltaR/(rc[column]*deltaT)
        
# Time integration via Backward Euler

inverseM = np.zeros((n*m,n*m))
I = np.identity(n*m)
A = np.zeros((n*m,n*m))



#Initial condition
c0 = np.zeros((n*m,1))
c1 = np.zeros((n*m,1))
dt = 0.1;

for alpha in range(0,n*m):
    column = columnChecker(alpha, n, m)
    c0[alpha] = chat*math.exp(-rc[column]*rc[column])

for i in range(0,n*m):
    inverseM[i,i] = 1.0/M[i,i]

A = I - dt*(inverseM@S)

eigvals = np.linalg.eigvals(inverseM@S)

area = 0
t = 0.0

while (area < 0.25*0.25*math.pi*R5*R5):
    area = 0
    t = t + dt
    c1 = np.linalg.solve(A,c0)
    c0 = c1
    for i in range(0,n*m):
        if c1[i] > chat:
            column = columnChecker(i,n,m)
            area = area + rc[column]*deltaR*deltaT;
    
print(t)

#Plot 
           
fig = plt.figure()
ax = fig.gca(projection='3d')

ax1 = np.linspace(deltaR/2.0,R5-deltaR/2.0,n)
ax2 = np.linspace(deltaT/2,4*phi-deltaT/2.0,m)
ax1, ax2 = np.meshgrid(ax1, ax2)
xx = np.zeros((n,m))
yy = np.zeros((n,m))

for i in range(0,n):
    for j in range(0,m):
        xx[i,j] = ax1[i,j]*math.cos(ax2[i,j])
        yy[i,j] = ax1[i,j]*math.sin(ax2[i,j])
    
ax1 = xx
ax2 = yy

Z = np.zeros((m,n))

for i in range(0,n):
    for j in range(0,m):
        Z[j][i] = c0[i + j*n,0]


# Plot the surface.
surf = ax.plot_surface(ax1, ax2, Z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)

# Customize the z axis.
#ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()
