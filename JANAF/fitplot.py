import matplotlib.pylab as plt
import numpy as np

# Open the file and read everything as a string into 'f'
f = open('fitout.dat','r')
#Separate the lines
lines = f.readlines()[0:]
#Close the file
f.close()

# Initialise lists to add the data to
temp = []
rawpt = []

for line in lines:
    # Split our line into columns
    data = line.split()
    # Add the columns to the lists
    temp.append(float(data[0]))
    rawpt.append(float(data[1]))

f = open('fitout2.dat','r')
lines = f.readlines()
f.close()
temp2 = []
fitpt = []

for line in lines:
    # Split our line into columns
    data = line.split()
    # Add the columns to the lists
    temp2.append(float(data[0]))
    fitpt.append(float(data[1]))

temp  = np.array(temp)
temp2 = np.array(temp2)
rawpt = np.array(rawpt)
fitpt = np.array(fitpt)

# plot raw data
plt.scatter(temp,rawpt,color='red',marker='o',label='data')

# plot model
plt.plot(temp2,fitpt,color='blue',label='model')

# Add x and y labels to the graph
plt.xlabel('T [K]')
plt.ylabel('dG, Kp or Pvap')

# Add a legend without the frame
plt.legend(frameon=False,loc='lower right')

# Show plot on screen
plt.show()
