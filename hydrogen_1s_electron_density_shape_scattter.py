import matplotlib.pyplot as plt
import numpy as np
#Probability of 1s
def prob_1s(x,y,z):
    r=np.sqrt(np.square(x)+np.square(y)+np.square(z))
    #Remember.. probability is psi squared!
    return np.square(np.exp(-r)/np.sqrt(np.pi))
#Random coordinates
x=np.linspace(0,1,30)
y=np.linspace(0,1,30)
z=np.linspace(0,1,30)
elements = []
probability = []
for ix in x:
    for iy in y:
        for iz in z:
            #Serialize into 1D object
            elements.append(str((ix,iy,iz)))
            probability.append(prob_1s(ix,iy,iz))
            
#Ensure sum of probability is 1
probability = probability/sum(probability)
#Getting electron coordinates based on probabiliy
coord = np.random.choice(elements, size=100000, replace=True, p=probability)
elem_mat = [i.split(',') for i in coord]
elem_mat = np.matrix(elem_mat)
x_coords = [float(i.item()[1:]) for i in elem_mat[:,0]] 
y_coords = [float(i.item()) for i in elem_mat[:,1]] 
z_coords = [float(i.item()[0:-1]) for i in elem_mat[:,2]]
#Plotting
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x_coords, y_coords, z_coords, alpha=0.05, s=2)
ax.set_title("Hydrogen 1s density")
plt.show()