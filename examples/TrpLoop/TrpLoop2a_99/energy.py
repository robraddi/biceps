import numpy as np
a=np.loadtxt("population.txt")
b=[]
for i in a:
	b.append(i/99791.0)
for j in b:	
	f=-np.log(j) # units in kT
	print f
