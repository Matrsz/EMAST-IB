import numpy as np
import matplotlib.pyplot as plt
import control as ct
import warnings
warnings.filterwarnings('ignore')
from scipy.signal import ellip, bessel, butter, cheby1

FILTRO = False

Is = 1069

def satelite_TF(rm2=0.34, f2=0.3):
    ns = 1
    ds = np.array([Is, 0, 0])

    tfS = ct.tf(ns,ds)

    sys = tfS

    rm1 = 0.04
    I0 = Is*(1-rm1)

    f1 = 2.0
    psi1 = 0.001
    w1 = 2*np.pi*f1

    psi1m = psi1*np.sqrt(Is/I0)
    w1m = w1*np.sqrt(Is/I0)

    n1 = Is/I0*np.array([1, 2*psi1*w1, w1*w1])
    d1 = np.array([1, 2*psi1m*w1m, w1m*w1m])

    tf1 = ct.tf(n1,d1) 
    sys = sys * tf1

    I0 = Is*(1-rm2)

    psi2 = 0.001

    w2 = 2*np.pi*f2
    psi2m = psi2*np.sqrt(Is/I0)
    w2m = w2*np.sqrt(Is/I0)

    n2 = Is/I0*np.array([1, 2*psi2*w2, w2*w2])
    d2 = np.array([1, 2*psi2m*w2m, w2m*w2m])

    tf2 = ct.tf(n2,d2)
    sys = sys * tf2
    return sys


fs =  [2*np.pi*1e-3, 2*np.pi*1e2]
mag,phase,w = ct.bode(satelite_TF(0.5,1.0),fs,plot=True,Hz=True,dB=True,deg=True)
mag,phase,w = ct.bode(satelite_TF(0.5,0.1),fs,plot=True,Hz=True,dB=True,deg=True)
mag,phase,w = ct.bode(satelite_TF(0.1,1.0),fs,plot=True,Hz=True,dB=True,deg=True)
plt.show()

ct.nichols(satelite_TF(0.5,1.0),[0.01,10])
ct.nichols(satelite_TF(0.5,0.1),[0.01,10])
ct.nichols(satelite_TF(0.1,1.0),[0.01,10])
plt.show()


def lead_lag_TF(k=3.17):
    psi = 1.
    Wn = np.sqrt ( k/Is )
    kd = 2*psi*Wn*Is
    Td = kd/k

    nc = np.array([k*Td, k])
    dc = np.array([Td/10, 1])

    tfC = ct.tf(nc,dc)
    sys = tfC

    Tdelay = 1
    ndelay = 1
    ddelay = np.array([Tdelay, 1])

    tfDelay = ct.tf(ndelay,ddelay)

    sys = sys*tfDelay
    return sys

fs =  [2*np.pi*1e-3, 2*np.pi*1e2]

mag,phase,w = ct.bode(lead_lag_TF(1), fs, plot=True,Hz=True,dB=True,deg=True)
mag,phase,w = ct.bode(lead_lag_TF(10),fs,plot=True,Hz=True,dB=True,deg=True)
mag,phase,w = ct.bode(lead_lag_TF(100),fs,plot=True,Hz=True,dB=True,deg=True)
plt.show()


def get_TF(rm2=0.34, f2=0.3, k=3.17, filtro=False, fc=0.1):
    sys_sat = satelite_TF(rm2, f2)
    sys_leadlag = lead_lag_TF(k)
    sys = sys_sat*sys_leadlag
    if filtro:
      n, d = cheby1(4, 1, fc, analog=True)
      sys_filtro = ct.tf(n,d)
      sys = sys*sys_filtro
    return sys


sys = get_TF(0.34, 0.3, 3.17)

mag,phase,w = ct.bode(sys,plot=True,Hz=True,dB=True,deg=True)

def findfirst(predicate, lst):
    for i, val in enumerate(lst):
        if predicate(val):
            return i
    return None  # Return None if no match is found

def first_local_maximum(lst):
    for i in range(1, len(lst) - 1):
        if lst[i] > lst[i - 1] and lst[i] > lst[i + 1]:
            return i # Return the index and value of the local maximum
    return None  # Return None if no local maximum is found

def cumple_requisitos(rm2, f2, k, filtro=False):
    sys = get_TF(rm2, f2, k, filtro)
    mag,phase,w = ct.bode(sys,plot=False)
    ico = findfirst(lambda x: x<1, mag)
    mag_dB = 10*np.log10(mag)
    phase_deg = 180*phase/np.pi
    phase_margin = phase_deg[ico]+180
    phase_margin = phase_margin % 360
    if phase_margin > 180:
        phase_margin -= 360
    imax = first_local_maximum(mag)
    flex_gain = mag_dB[imax]
    return flex_gain<=-6 and abs(phase_margin)>=30  , w[ico]*2*np.pi, phase_margin, flex_gain

print(cumple_requisitos(0.34, 0.3, 3.17))


def binary_search_predicate(predicate, low, high, tolerance=1e-3):
    while high - low > tolerance:
        mid = np.sqrt(low * high)  # Geometric mean
        if predicate(mid):
            low = mid    # Search the upper half
        else:
            high = mid   # Search the lower half
    return mid if predicate(mid) else low  # Return the midpoint as the approximate result


## F = 0.1
F = 0.1
R = 0.1
predicate = lambda k: cumple_requisitos(R, F, k)[0]
k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)

result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.3
predicate = lambda k: cumple_requisitos(R, F, k)[0]

k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.5
predicate = lambda k: cumple_requisitos(R, F, k)[0]
k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

## F = 0.5
F = 0.5
R = 0.1
predicate = lambda k: cumple_requisitos(R, F, k)[0]

k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.3
predicate = lambda k: cumple_requisitos(R, F, k)[0]

k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.5
predicate = lambda k: cumple_requisitos(R, F, k)[0]
k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")


## F = 1
F = 1.0
R = 0.1
predicate = lambda k: cumple_requisitos(R, F, k)[0]

k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.3
predicate = lambda k: cumple_requisitos(R, F, k)[0]

k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.5
predicate = lambda k: cumple_requisitos(R, F, k)[0]
k = binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

def find_best_k(R, F):
    predicate = lambda k: cumple_requisitos(R, F, k)[0]
    return binary_search_predicate(predicate, low=0.1, high=100, tolerance=1e-3)

#fig, [a1, a2] = plt.subplots(2,1)
#
#Rs = np.linspace(0.1, 0.5, 100)
#F = 0.5
#
#a1.plot(Rs, [find_best_k(R, F) for R in Rs])
#
#Fs = np.linspace(0.1, 1, 100)
#R = 0.3
#
#a2.plot(Fs, [find_best_k(R, F) for F in Fs])
#plt.show()

print("Con filtro!!!!!")

## F = 0.1
F = 0.1
R = 0.1
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]
k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)

result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.3
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]

k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.5
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]
k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

## F = 0.5
F = 0.5
R = 0.1
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]

k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.3
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]

k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.5
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]
k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")


## F = 1
F = 1.0
R = 0.1
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]

k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.3
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]

k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")

R = 0.5
predicate = lambda k: cumple_requisitos(R, F, k, True)[0]
k = binary_search_predicate(predicate, low=0.1, high=100000, tolerance=1e-3)
result, fco, phase_margin, flex_gain = cumple_requisitos(R, F, k, True)
print(f"Con f = {F} Hz y r = {R}\n k = {k}, fco = {fco} Hz, margen fase = {phase_margin} deg y g modo flex = {flex_gain} dB")
