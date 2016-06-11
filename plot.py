# Peter Han, A10187618, Physics 239, Murphy
# Python 3.5.1 :: Anaconda 4.0.0 (64-bit)
import matplotlib.pyplot as plt
import numpy as np
from operator import itemgetter

# A = Z + N
# A : massNumber, Z : atomicNumber, N : neutronNumber
# https://en.wikipedia.org/wiki/Semi-empirical_mass_formula#The_formula
def bindingEnergy(A):
	Z = zMinimum(A)
	N = A - Z
	volume = 15.75 * A
	surface = -17.8 * A ** (2 / 3)
	coulomb = -0.711 * Z ** 2 / A ** (1 / 3)
	asymmetric = -23.7 * (A - 2 * Z) ** 2 / A
	pairing = 11.18 * pair(Z, N) / A ** 0.5
	return volume + surface + coulomb + asymmetric + pairing

# helper function to bindingEnergy
def zMinimum(A):
	return round(A / (0.5 * 0.711 / 23.7 * A ** (2 / 3) + 2))

# helper function to bindingEnergy
def pair(Z, N):
	zParity = Z % 2
	nParity = N % 2
	if zParity == 0 and nParity == 0:
		return 1
	elif zParity != 0 and nParity != 0:
		return -1
	else:
		return 0

def getData(filepath):
	try:
		f = open(filepath, "r")
	except IOError as e:
		print("Handling run-time error:", e)
		raise SystemExit(1)
	data = []
	for line in f:
		# remove weird lines from data file
		if line[0] != "#" and line[len(line) - 5] != "-":
			data.append(line.split())
	f.close()
	return data

def getBindingEnergyPerNucleon(data):
	subset = [] # temporary holding list of entries to find minimum
	collect = [] # list of entries with minimum binding energy per nucleon
	count = 1 # compare with mass number

	# append entries that have the same mass number to subset list
	# when mass number increments, find the minimum from the subset list
	# and append that to collect
	# clear the subset list and loop until end of data
	for i in range(len(data)):
		if int(data[i][0]) != count:
			# get the minimum mass from the subset list
			collect.append(min(subset, key = itemgetter(1)))
			subset = []
			count += 1
		subset.append([int(data[i][0]), float(data[i][3]), float(data[i][5])])
	a = np.array(collect)
	x = a[:, 0].tolist()
	y = a[:, 2].tolist()
	return (x, y)

def getError(model, data):
	x = []
	y = []
	for i in range(258):
		if data[i] != 0:
			x.append(i + 1)
			y.append(abs(model[i] - data[i]) / data[i])
	return (x, y)

# data
file = "C:/Users/Peter/Dropbox/Programming/Python/binding-energy/binding.txt"
data = getData(file)
xData, yData = getBindingEnergyPerNucleon(data)

# model
xModel = []
yModel = []
for i in range(280):
	k = i + 1
	xModel.append(k)
	yModel.append(bindingEnergy(k) / k)

# error
xError, yError = getError(yModel, yData)
logyError = np.log10(yError)
print("Average error from 56Fe to 258Md:", np.average(yError[54:]))

# plotting data on top of model
plt1 = plt
plt1.plot(xData, yData, ".r")
plt1.plot(xModel, yModel, "-k")
plt1.xlim([0, 300])
plt1.ylim([0, 9])
plt1.xlabel("Nuclear Mass Number, A")
plt1.ylabel("Binding Energy Per Nucleon (MeV)")
plt1.title("Nuclear Liquid Drop Model Predictions")
plt1.savefig("plot1.png")
plt1.show()

# plotting error
plt2 = plt
plt2.plot(xError, logyError, "-ko")
plt2.xlim([0, 300])
plt2.ylim([-6, 1])
plt2.xlabel("Nuclear Mass Number, A")
plt2.ylabel("log10 %Error Between Model And Data")
plt2.title("Nuclear Liquid Drop Model Error")
plt2.savefig("plot2.png")
plt2.show()
