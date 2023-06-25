import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

N, beta, q, x, dq, dx = np.loadtxt('plot_data/datiplot.dat', unpack=True)

BETA = np.linspace(np.min(beta), np.max(beta), 1000)

#====================================================================================
# Plot
#====================================================================================

plt.figure(1)
plt.title("topological susceptibility, gauge U(1)", fontsize=1)
plt.xlabel(r'$\beta$', fontsize=15)
plt.ylabel(r'$\chi \, a^2$', fontsize=15)
plt.grid()

#plt.errorbar(beta, q, dq, fmt='.', color='r', label='charge')
plt.errorbar(beta, x, dx, fmt='.', color='k', label='susc')
plt.plot(BETA, 1/(4 * np.pi**2 * BETA), 'b', label=r'$\frac{1}{4 \pi \beta}$')
plt.legend(loc='best')



BETA = np.linspace(1/np.min(beta), 1/np.max(beta), 1000)

#====================================================================================
# Fit
#====================================================================================

def f(x, m , q):
    return m*x + q

#x = x[1:]
#dx = dx[1:]
#beta = beta[1:]


pars, covm = curve_fit(f, 1/beta, x, sigma=dx)

print('m  = %.5f +- %.5f ' % (pars[0], np.sqrt(covm.diagonal()[0])))
print('q  = %.5f +- %.5f ' % (pars[1], np.sqrt(covm.diagonal()[1])))


chisq = sum(((x - f(1/beta, *pars))/dx)**2.)
ndof  = len(x) - len(pars)
print(f'chi quadro = {chisq:.3f} ({ndof:d} dof)')


c = np.zeros((len(pars),len(pars)))

for i in range(0, len(pars)):
    for j in range(0, len(pars)):
       c[i][j] = (covm[i][j])/(np.sqrt(covm.diagonal()[i])*np.sqrt(covm.diagonal()[j]))
print(c) #correlation matrix

#====================================================================================
# Fit's plot
#====================================================================================

fig1 = plt.figure(2)
frame1=fig1.add_axes((.1,.35,.8,.6))

frame1.set_title("Topological susceptibility, gauge U(1)", fontsize=15)
plt.ylabel(r'$\chi \, a^2$', fontsize=15)
plt.grid()


plt.errorbar(1/beta, x, dx, fmt='.', color='k', label='data')
plt.plot(BETA, f(BETA, *pars), color='blue', label='best fit') 
plt.legend(loc='best')


#residuals' graphic
frame2=fig1.add_axes((.1,.1,.8,.2))

# residuals
ff = (x - f(1/beta, *pars))/dx
frame2.set_ylabel('Residui Normalizzati')
plt.xlabel(r'$1/\beta$', fontsize=15)

plt.plot(BETA, 0*BETA, color='red', linestyle='--', alpha=0.5)
plt.plot(1/beta, ff, '.', color='black')
plt.grid()


plt.show()
