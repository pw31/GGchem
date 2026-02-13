import matplotlib.pylab as plt
import numpy as np

#fits = [['fitout.dat','fitout2.dat']]
fits = [['liquid.out','liquid2.out'],['solid.out','solid2.out']]
for fit in fits:

  # Open the file and read everything as a string into 'f'
  f = open(fit[0],'r')
  lines = f.readlines()[0:]
  f.close()
  temp = []
  rawpt = []
  for line in lines:
    data = line.split()
    temp.append(float(data[0]))
    rawpt.append(float(data[1]))

  f = open(fit[1],'r')
  lines = f.readlines()
  f.close()
  temp2 = []
  fitpt = []
  for line in lines:
    data = line.split()
    temp2.append(float(data[0]))
    fitpt.append(float(data[1]))

  temp  = np.array(temp)
  temp2 = np.array(temp2)
  rawpt = np.array(rawpt)
  fitpt = np.array(fitpt)

  # plot raw data
  plt.scatter(temp,rawpt,marker='o')

  # plot model
  plt.plot(temp2,fitpt,label=fit[0])

  aa = -1.7811484E+07
  bb =  6.0492334E+03
  yy = (aa + bb*temp2)/1000.0
  plt.plot(temp2,yy)
  
# Add x and y labels to the graph
plt.xlabel('T [K]')
plt.ylabel('dG, Kp or Pvap')

# Add a legend without the frame
plt.legend(frameon=False,loc='lower right')

# Show plot on screen
plt.show()
