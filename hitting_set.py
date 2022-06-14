#!/usr/bin/python3
import sys
from math import ceil
import scipy.optimize

def compute(V, init = 10):
	if len(V) <= 1:
		return 1
	def h(x):
		return sum([x**(-v) for v in V])-1
	return scipy.optimize.brenth(h,1,init)

def compute_conditional(V, bn1, bn2):
	if sum([bn2**(-v) for v in V]) > 1:
		return bn+1000
	elif bn1 > 0 and sum([bn1**(-v) for v in V]) < 1:
		return 0
	else:
		return compute(V)

def compute_conditional2(V, bn1, bn2):
	return compute(V)

def Ceil(x): return ceil(x*100000)/100000

def str_vec(V):
	X = [f"{x:.3f}" for x in V]
	return "("+", ".join(X)+")"

def calc(Psi, deg, upper_bn, verbose = False):
	bn = 0

	if verbose:	print("Rule B1")
	for d2 in range(1,deg+1):
		for m2 in range(d2, 8+1): # range is 8 since Psi[x]=Psi[8] for every x>8
			Delta1 = 1-(Psi[m2]-Psi[m2-d2])
			d3 = max(0,3-d2)
			Delta2 = d2-(Psi[m2]-Psi[max(m2-d2**2,0)+d3])
			bn2 = compute_conditional([Delta1,Delta2], bn, upper_bn)
			if verbose:
				print(f"({Delta1:.3f}, {Delta2:.3f}) {Ceil(bn2):.5f} d2 = {d2} m2 = {m2} ")
			bn = max(bn, bn2)
			if bn > upper_bn:
				return bn
	if verbose:
		print(f"max bn: {Ceil(bn):.5f}")

	if verbose: print("\nRule B2")
	# Since Psi is monotone, we can skip the cases when m2 < 3
	v = [1,2-max(Psi[m2]-Psi[m2-3] for m2 in range(3,len(Psi)) )]
	bn2 = compute_conditional(v, bn, upper_bn)
	if verbose:	print(str_vec(v), Ceil(bn2))
	bn = max(bn, bn2)
	if bn > upper_bn: return bn

	for m2 in range(4,len(Psi)):
		# Since Psi is monotone, we can skip the cases when m2 < 4
		v = [ 2-(Psi[m2]-Psi[max(m2-5,0)]), 2-(Psi[m2]-Psi[m2-4]) ]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(m2, str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)
		if bn > upper_bn: return bn

	# We don't need to check the vector (1,1-Delta(-2)) since its branching number is less than 2

	if deg in [3,4]:
		if verbose: print("\nRule B3")
		v = [min(2,1+Psi[4]),Psi[2]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)
		if bn > upper_bn: return bn

	if verbose and deg > 3:	print("\nRule B6")

	if deg in [4,5]:
		if verbose: print(" +Rule R3 or R4")
		v = [1,Psi[deg+1]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

	if deg == 5:
		if verbose: print(" +Rule B1, d2(x)=1")
		v = [1,1+Psi[deg-1],1+Psi[deg+1]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

		if verbose: print(" +Rule B1, d2(x)=2")
		v = [1,1+Psi[deg-2],2+Psi[deg-3]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

	if deg == 4:
		if verbose: print(" +Rule B1, d2(x)=1")
		v = [1,1+Psi[deg-1],1+Psi[deg+2]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

		v = [1+Psi[2],1+Psi[deg-1],1+Psi[deg+1]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

		if verbose: print(" +Rule B1, d2(x)=2")
		v = [1,1+Psi[deg-2],2+Psi[2]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

	if deg in [4,5]:
		if verbose: print(" +Rule B1, d2(x)=3")
		v = [1,1+Psi[deg-3],3+Psi[1]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

	if deg == 5:
		if verbose: print(" +Rule B1, d2(x)=4")
		v = [1,1+Psi[deg-4],4+Psi[1]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

	if deg == 4:
		if verbose: print(" +Rule B2")
		v = [1,1+Psi[deg],2+Psi[deg-3]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

	if deg == 6:
		v = [1,Psi[deg]]
		bn2 = compute_conditional(v, bn, upper_bn)
		if verbose:	print(str_vec(v), Ceil(bn2))
		bn = max(bn, bn2)

	if verbose:
		print(f"\nMax branching number: {Ceil(bn):.5f}", bn)

	return bn

def W_to_Psi(W):
	Psi = [0]
	for x in W:
		Psi.append(Psi[-1]+x)
	return Psi

def frange(a, b, step):
	x = (a+b)/2
	yield x
	while x < b:
		x += step
		yield x
	x = (a+b)/2
	while x > a:
		x -= step
		yield x

def frange2(a, step):
	return frange(max(0,a-step*4),a+step*4,step)

mode = 1
deg = 4
params = 8

for x in sys.argv[1:]:
	if x == "-x":
		mode = 2
	elif x[0] == "-":
		params = int(x[1:])
	else:
		deg = int(x)

if mode == 1:
	Psi3 = [0, 0.2337409973144531, 0.4515090942382813, 0.6475447654724121, 0.8149245262145997, 0.9563242912292481, 1.0559596061706544, 1.1075743675231935, 1.1075743675231935, 1.1075743675231935]
	Psi4 = [0, 0.25637664794921877, 0.4892762184143067, 0.6895000457763673, 0.8463504791259766, 1.0463504791259766, 1.2151004791259765, 1.3338504791259764, 1.3963504791259764, 1.3963504791259764]
	Psi5 = [0, 0.23786649703979496, 0.4564464569091797, 0.6487490653991699, 0.806581974029541, 0.9213740348815919, 1.0286317825317384, 1.0843505859375, 1.0843505859375, 1.0843505859375]
	Psi6 = [0, 0.20819206237792973, 0.4082907676696778, 0.5888676643371583, 0.7442749023437502, 0.8683744430541993, 0.9552198410034181, 1.0, 1.0, 1.0]

	Psi = [Psi3, Psi4, Psi5, Psi6][deg-3]

	compute_conditional = compute_conditional2
	bn = calc(Psi, deg, 999, True)
	sys.exit(0)

bn = 999
W = None
for a in frange(0, 0.4, 0.1):
	for b in frange(0, 0.4, 0.1):
		for c in frange(0, 0.4, 0.1):
			for d in frange(0, 0.4, 0.1):
				for e in frange(0, 0.4, 0.1):
					for f in frange(0, 0.4, 0.1):
						for g in frange(0, 0.4, 0.1) if params >= 7 else [0]:
							for h in frange(0, 0.4, 0.1) if params >= 8 else [0]:
								W2 = [a,b,c,d,e,f,g,h]+[0]
								Psi = W_to_Psi(W2)
								bn2 = calc(Psi, deg, bn)
								if bn2 < bn:
									bn = bn2
									W = W2
									print(f"{a:.1f} {b:.1f} {c:.1f} {d:.1f} {e:.1f} {f:.1f} {g:.1f} {h:.1f} {Ceil(bn)}")

delta = 0.1
for i in range(10):
	print(W_to_Psi(W))
	print()
	X = W
	delta = delta/4
	for a in frange2(X[0], delta):
		for b in frange2(X[1], delta):
			for c in frange2(X[2], delta):
				for d in frange2(X[3], delta):
					for e in frange2(X[4], delta):
						for f in frange2(X[5], delta):
							for g in frange2(X[6], delta) if params >= 7 else [0]:
								for h in frange2(X[7], delta) if params >= 8 else [0]:
									W2 = [a,b,c,d,e,f,g,h]+[0]
									Psi = W_to_Psi(W2)
									bn2 = calc(Psi, deg, bn)
									if bn2 < bn:
										bn = bn2
										W = W2
										print(f"{a:.4f} {b:.4f} {c:.4f} {d:.4f} {e:.4f} {f:.4f} {g:.4f} {h:.4f} {Ceil(bn)}")


print()
Psi = W_to_Psi(W)
compute_conditional = compute_conditional2
bn = calc(Psi, deg, 999, True)
print(Psi)
