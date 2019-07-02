import csv
import math
from math import log
import os

# Parameter initialization / FIXED PARAMETERS
############################  General parameters  ############################ 
tau = 40						# cell cycle [min]
time_effect = log(2)/tau		# doubling time effect [1/min]
nA = 6.023e23					# Avogadro's number: # particles/mol
Vcell =  1.1e-9					# typical volume of E. coli (microliters/cell). Source: Bionumbers
Vext = 1e-3						# culture medium volume [ul] in a microfluidics device
Vc = Vcell/Vext					# Volume Ratio, for nM to numbr of molecules conversion
Ccells_OD600_1=8e5				# cells/microliter for OD=1.
									# Source: Agilent, E. coli Cell Culture Concentration from OD600 Calculator
factor_units=1e15/(Vcell*nA)	# Conversion factor from number of particles to
									# concentration (1 nMolar=nanomols/liter)
									# 1 particle/cell x 1/Vcell cell/microliter
									# x 1/nA mol/particle x
									# 1e6 microliter/liter x 1e9 nanomol/l

# Burst size and plasmid number
b = 20							# burst size, from Weber
pN_luxI = 17					# plasmid number of luxI in pBR322 (15-20 copies/cell)
pN_luxR = 10					# plasmid number of luxI in pACYC184 (10 copies/cell)

# LuxI
alphaI = 0.01					# leakage of repressor
#pN = pNluxI*factor_units		# copy number [nM/min]
dmI = time_effect + 0.23		# degradation rate of LuxI mRNA [1/min]. Typical mRNA degradation rate for E coli
dI = time_effect + 0.01			# degradation rate of LuxI [1/min]
pI = 3.09						# translation rate of LuxI #mRNA [1/min]. b*dmI from Weber.  [3.0928 - 6.1856]
kI = 1.03						# average transcription rate of LuxI [1/min] from calculator [1.0309 - 10.3093] from our rates calculator based on Alberts, and
									# Ref 2: 50/b from Weber

# LuxR
alphaR = 0.01					# ratio between unactivated vs activated of expresion LuxR
dmR = time_effect + 0.23		# degradation rate of LuxR mRNA [1/min]
dR = time_effect + log(2)/5		# degradation rate of LuxR [1/min]
pR =  2.38 #b*dmR;				# translation rate of LuxR [1/min] [2.381 - 4.7619] from our rates calculator, Alberts
kR = 0.79 #200/b;				# transcription rate of LuxR [1/min] [0.79365 - 7.9365] from our rates calculator
cR = kR*pN_luxR					# constitutive expression of LuxR (always ConstR = initCond of plasmids in complete model)

# Monomer
kd1 = 100						# dissociation constant of R to A [nM], Urbanowski etal. 2004
k_1 = 10						# unbinding rate LuxR to AHL [1/min]
k1 = k_1/kd1					# binding rate LuxR to AHL [1/min]

dRA = time_effect + log(2)/5	# degradation rate of (LuxR.A) [1/min]. Buchler et al. 2004 Monomer half-life is just few minutes.

# Dimer
kd2 = 20						# or 10 dissociation cte (LuxR.A) [nM], Buchler et al. 2003
									# Koren, R. & Hammes, G. G. (1976) Biochemistry 15, 1165�1171.
									# Northrup, S. H. & Erickson, H. P. (1992) Proc. Natl. Acad. Sci. USA 89, 3338�3342.
k_2 = 1							# dissociation rate dimer (LuxR.A)2 [1/min]
k2 = k_2/kd2					# binding rate LuxR to AHL [1/min]

kdLux = 100						# dissociation cte (LuxR.A)2 to promoter [nM], Bucler et al [1 1000]nM

dRA2 =time_effect				# Ron Weiss et al. A synthetic multicellular system for programmed pattern formation
								# Buchler et al 2004 Degradation rate for (LuxR.AHL)_2. Corresponding to dilution if cell half-life = tau minutes
# AHL and AHLe
kA = 0.04						# synthesis rate of AHL by LuxI [1/min] from Bionumbers

# with QS
D = 2
								# kinetic rate of AHL external transport [1/min] across the cell membrane, calculated

dA = time_effect + 0.04			# [0.05 0.03 min^-1]Degradation from Bionumbers online
dAe = 0.00000164				# antes 0.04. 0.0164, Degradation rate for AHL. From Kaufmann etal. 2005. Similar to You etal. Nature, 2004
								# [0.05 0.03 min^-1] Degradation from Bionumbers online
								# 0.000282 Degradation rate for external AHL. Fitted using half-life of 180 minutes, from Englmann etal. 2007
								# In Kauffmann & Sartorio, 2005 they use 0.0018

# Arguments of the langevin program.
Ncells = 240
ruido = 0.15
ahl_e_0 = 0/Vc
STO = 1							# STO = 0 is deterministic simulation, and STO = 1 is stochastic.
HISTO = 1						# HISTO = 1 for obtaining the steady state histogram
TEMPO =  1						# TEMPO = 1 For obtaining the temporal response of the means and std
TEMPOT =  1						# TEMPOT = 1 For obtaining the temporal response of all cells

print('Saving params...')
# Params list
param_out = [dI, pI, kI, pN_luxI, dmI, kdLux, alphaI, dR, pR, cR, dmR, k_1, kd1, k_2, kd2, dRA2, kA, dA, D, Vcell, Vext, dRA, dAe]

# Write a file named param.sT, with the struct param_out
filename = 'param.dat'
with open(filename, 'w') as myfile:
	wr = csv.writer(myfile, quoting=csv.QUOTE_NONE)
	wr.writerow(param_out)

print('Params saved in ' + filename)
# Excecuting the external C++ program langevin with 4 cores, and with the files param.dat as input
command = 'mpirun -np 1 ./langevin param.dat ' + str(Ncells) + ' ' + str(ruido) + ' ' + str(ahl_e_0) + ' ' + str(STO) + ' ' + str(HISTO) + ' ' + str(TEMPO) + ' ' + str(TEMPOT)
print(command)
os.system(command)