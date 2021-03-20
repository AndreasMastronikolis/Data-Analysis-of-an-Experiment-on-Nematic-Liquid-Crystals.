'''
Abstract:

This project constitutes an illustration of how I analyze data I collect from a 
laboratory experiment; which in this case aims to calculate the response times, 
threshold voltage and viscosity (at different temperatures) of the Nematic Liquid 
Crystal 5CB. The data analyzed were not measured by me, they were generated with 
pseudo-random algorithm (to simulate errors). The script that performs the 
analysis is written in python and produces analytic plots of my findings

Author: Andreas Mastronikolis, Date: March 2019
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.optimize import curve_fit
from scipy import constants as cs
from scipy import integrate

Un_Angle = 1 # Uncertainty in the Analyzer Angle (Precision of the Analyzer) in Degrees.
Un_Voltage = 0.02 # Uncertainty in the Voltage reading in the photodetector (Precision of the Oscilloscope reading) in Volts.
Un_Temperature = 0.3 # Uncertainty in the Temperature.
LC_Thickness = 5 * 10**(-6) # SI
Lambda_0 = 632.812 * 10**(-9) # SI
epsilon_0 = cs.epsilon_0

## ------------------- Input Command ------------------- ##

# The variable Command shall take a keyword that will make the script execute a specific task.

Command = input("Command Prompt: ")

def Cos_Func(X, Delta, I):
    return I*(np.cos(np.radians(X - Delta)))**2

def Sin_Func(X, Delta, I):
    return I*(np.sin(np.radians(2*X - 2*Delta)))**2

def Line_Fit(X, Slope, Offset):
    return Slope * X + Offset

def S_Line(x, eta):
    return eta * x

def Local_Function_Fit(T, Beta, Critical_Temp, A):
    return A*(np.abs(Critical_Temp - T))**Beta

def Polynomial_Fit(x, A, B, C, D, E):
        return A * x**3 + B * x**2 + C * x + D + E * x**4
    

if Command == 'Photodiode Callibration' or Command == 'All':
    
    # Uncertainties and other numerical information
    
    Un_Angle = 1.5 # Uncertainty in the Analyzer Angle (Precision of the Analyzer) in Degrees.
    Un_Voltage = 0.05 # Uncertainty in the Voltage reading in the photodetector (Precision of the Oscilloscope reading) in Volts.
        
    ## File A1 contains the data acquired from the task (A1) in the Lab Script.

    A1 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\A1.csv', delimiter = ',', dtype=np.dtype(float))
    
    A1_Curve_Fit, A1_Cov_Fit = curve_fit(Cos_Func, A1[:,0], A1[:,1], p0 = [55, 2], sigma = np.full(np.shape(A1[:,1]), Un_Voltage), absolute_sigma = True)
    A1_Delta_Error = np.sqrt(np.diag(A1_Cov_Fit))
    
    A1_Lin_Fit, A1_Cov_Lin_Fit = curve_fit(Line_Fit, Cos_Func(A1[:,0], A1_Curve_Fit[0], A1_Curve_Fit[1]), A1[:,1], p0 = [1, 0], sigma = np.full(np.shape(A1[:,1]), Un_Voltage), absolute_sigma = True)
    A1_Lin_Error = np.sqrt(np.diag(A1_Cov_Lin_Fit))
    
    ## Analysis for (A1) Data and Visualization
    
    Figure = plt.figure(figsize = (10,9))
    Sizes = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    
    Data_And_Fit = Figure.add_subplot(Sizes[0], label = 'Data_And_Fit')
    Data_And_Fit.grid(alpha = 1)
    Data_And_Fit.set_ylabel(r'Photodetector Voltage [Volts]', fontsize = 13)
    Data_And_Fit.set_xlabel(r'$\Theta$ [Degrees$^{\circ}$]', fontsize = 13)
    Data_And_Fit.set_title(r"Voltage Response vs. Analyzer's Angular Position ($\Theta$)")
    Data_And_Fit.errorbar(A1[:,0], A1[:,1], xerr = Un_Angle, yerr = Un_Voltage, fmt = 's r', markersize = 5.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Measured Data', alpha = 0.6)
    Domain = np.linspace(0,175,num=200,endpoint=True)
    Data_And_Fit.plot(Domain, Cos_Func(Domain, A1_Curve_Fit[0], A1_Curve_Fit[1]), '-.k', label = 'Best Fit', lw = 2, alpha = 1)
    Data_And_Fit.text(0.01, 0.050, r'$\Delta = ($' + ('%.2f' % A1_Curve_Fit[0]) + r' $\pm$ ' + ('%.2f' % A1_Delta_Error[0]) + r')$^{\circ}$', fontsize=11, transform=plt.gcf().transFigure)
    Data_And_Fit.text(0.01, 0.025, r'V$_0 = ($' + ('%.2f' % A1_Curve_Fit[1]) + r' $\pm$ ' + ('%.2f' % A1_Delta_Error[1]) + r') Volts', fontsize=11, transform=plt.gcf().transFigure)
    Data_And_Fit.legend()
    
    Residuals = Figure.add_subplot(Sizes[1], label = 'Residuals')
    Residuals.set_ylabel('Residuals', fontsize = 13)
    Residuals.set_xlabel(r'$\Theta$ [Degrees$^{\circ}$]', fontsize = 13)
    Residuals.errorbar(A1[:,0], Cos_Func(A1[:,0], A1_Curve_Fit[0], A1_Curve_Fit[1]) - A1[:,1], yerr = Un_Voltage, fmt = 's r', markersize = 5.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Measured Data', alpha = 0.6)
    Residuals.axhline(0, color = 'black')
    Residuals.grid(alpha = 1)
    
    ## Linearity Check of the Photodetector.
    
    Figure = plt.figure(figsize=(10,9))
    Photodetector_Data_And_Fit = Figure.add_subplot(Sizes[0], label = 'Data_And_Fit')
    Photodetector_Data_And_Fit.grid(alpha = 1)
    Photodetector_Data_And_Fit.set_ylabel(r"Photodetector's Voltage Response [Volts]", fontsize = 13)
    Photodetector_Data_And_Fit.set_title(r"Photodetector's Voltage Response vs. Input Intensity", fontsize = 13)
    Photodetector_Data_And_Fit.set_xlabel(r'Input Intensity [Volts]', fontsize = 13)
    
    Photodetector_Data_And_Fit.errorbar(Cos_Func(A1[:,0], A1_Curve_Fit[0], A1_Curve_Fit[1]), A1[:,1], yerr = Un_Voltage, fmt = 's r', markersize = 5.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Observed Data Points', alpha = 0.6)
    Range_Cos_Func = np.linspace(0, 1.75, num = 100, endpoint = True)
    Photodetector_Data_And_Fit.plot(Range_Cos_Func, Line_Fit(Range_Cos_Func, A1_Lin_Fit[0], A1_Lin_Fit[1]), color = 'black', ls = '-.', label = 'Best Fit (Linear)', alpha = 1)
    Photodetector_Data_And_Fit.text(0.01, 0.050, r'$G = ($' + ('%.2f' % A1_Lin_Fit[0]) + r' $\pm$ ' + ('%.2f' % A1_Lin_Error[0]) + r')$^{\circ}$', fontsize=11, transform=plt.gcf().transFigure)
    Photodetector_Data_And_Fit.text(0.01, 0.025, r'$z = ($' + ('%.2f' % A1_Lin_Fit[1]) + r' $\pm$ ' + ('%.2f' % A1_Lin_Error[1]) + r') Volts', fontsize=11, transform=plt.gcf().transFigure)
    Photodetector_Data_And_Fit.legend()
    
    # Residual Plot for the Linearity Test.
    
    Residuals = Figure.add_subplot(Sizes[1], label = 'Residuals')
    Residuals.set_ylabel('Residuals', fontsize = 13)
    Residuals.grid(alpha = 1)
    Residuals.axhline(0, color = 'black') # Here, we set a flat horizontal line y = 0; as it wish to present positive or negative deviations of the data points from the fit.

    
    X_Values = Cos_Func(A1[:,0], A1_Curve_Fit[0], A1_Curve_Fit[1]) # This is an array of the totality of the X-Values we wish to plot.
    Y_Values = Line_Fit(Cos_Func(A1[:,0], A1_Curve_Fit[0], A1_Curve_Fit[1]), A1_Lin_Fit[0], A1_Lin_Fit[1]) - A1[:,1] # This is an array of the totality of the Y-Values we wish to plot.
    Residuals.errorbar(X_Values, Y_Values, yerr = Un_Voltage, fmt = 's r', markersize = 6.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Measured Data', alpha = 0.6)
    
    
    
if Command == 'A2' or Command == 'All':
    
    # File A2 contains the data acquired from the task (A2) in the Lab Script.
    
    A2 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\A2.csv', delimiter = ',', dtype=np.dtype(float))
    A1 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\A1.csv', delimiter = ',', dtype=np.dtype(float))
    A2_Curve_Fit, A2_Cov_Fit = curve_fit(Sin_Func, A2[:,0], A2[:,1], p0 = [0, 2], sigma = np.full(np.shape(A1[:,1]), Un_Voltage), absolute_sigma = True)
    A2_Delta_Error = np.sqrt(np.diag(A2_Cov_Fit))
    
    ## Analysis for (A2) Data and Visualization.
    
    Figure = plt.figure(figsize = (10,9))
    Sizes = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    
    Data_And_Fit = Figure.add_subplot(Sizes[0])
    Data_And_Fit.grid(alpha = 1)
    Data_And_Fit.set_ylabel(r'Photodetector Voltage [Volts]', fontsize = 13)
    Data_And_Fit.set_title(r"Voltage Response vs. LC Angular Position ($\Theta'$)")
    Data_And_Fit.set_xlabel(r"$\Theta'$ [Degrees$^{\circ}$]", fontsize = 13)
    Data_And_Fit.errorbar(A2[:,0], A2[:,1], xerr = Un_Angle, yerr = Un_Voltage, fmt = 's r', markersize = 5.0, mec = 'black', mfc = 'white', capsize = 2, elinewidth = 1, label = 'Measured Data', alpha = 0.5)
    Domain = np.linspace(0,175,num=200,endpoint=True)
    Data_And_Fit.plot(Domain, Sin_Func(Domain, A2_Curve_Fit[0], A2_Curve_Fit[1]), '--k', label = 'Best Fit', lw = 2, alpha = 1)
    Data_And_Fit.text(0.01, 0.050, r'$\Delta = ($' + ('%.2f' % A2_Curve_Fit[0]) + r' $\pm$ ' + ('%.2f' % A2_Delta_Error[0]) + r')$^{\circ}$', fontsize = 11, transform=plt.gcf().transFigure)
    Data_And_Fit.text(0.01, 0.025, r'V$_0 = ($' + ('%.2f' % A2_Curve_Fit[1]) + r' $\pm$ ' + ('%.2f' % A2_Delta_Error[1]) + r') Volts', fontsize = 11, transform=plt.gcf().transFigure)
    Data_And_Fit.legend(fontsize = 13)
    
    Residuals = Figure.add_subplot(Sizes[1], label = 'Residuals')
    Residuals.set_ylabel('Residuals', fontsize = 13)
    Residuals.errorbar(A2[:,0], Sin_Func(A2[:,0], A2_Curve_Fit[0], A2_Curve_Fit[1]) - A2[:,1], yerr = Un_Voltage, fmt = 's r', markersize = 7.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, alpha = 0.6)
    Residuals.axhline(0, color = 'black', ls = '--')
    Residuals.grid(alpha = 1)
    

if Command == 'Viscosity of 5CB' or Command == 'All':
    
    A33 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\A3_3.csv', delimiter = ',', dtype=np.dtype(float))
    A33[0,0] = 32.5
    B3 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\B3.csv', delimiter = ',', dtype = np.dtype(float), encoding=None)
    B3[0,0] = 20
    B2 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\B2.csv', delimiter = ',', dtype = np.dtype(float), encoding=None)
    B2[0,0] = 20
    
    A3_1 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\A3_1.csv', delimiter = ',', dtype=np.dtype(float))
    A3_2 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\A3_2.csv', delimiter = ',', dtype=np.dtype(float))

    B1 = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\B1.csv', delimiter = ',', dtype=np.dtype(float))
    B1[0,0] = 20
    
    A3_3_Curve_Fit, A3_3_Cov_Fit = curve_fit(Polynomial_Fit, A33[:,0], A33[:,1], p0 = [1, 1, 1, 1, 1], sigma = np.full(np.shape(A33[:,1]), 0.1), absolute_sigma=True) # Polynomail Fit, 4th degree.
    A3_3_Error = np.sqrt(np.diag(A3_3_Cov_Fit)) # Errors on the parameters.
    
    A3_3_Curve_Fit_Sqrt, A3_3_Cov_Fit_Sqrt = curve_fit(Local_Function_Fit, A33[:,0], A33[:,1], p0 = [0.5, 32.5, 5], sigma = np.full(np.shape(A33[:,1]), 0.1), absolute_sigma=True)
    A3_3_Error_Sqrt = np.sqrt(np.diag(A3_3_Cov_Fit_Sqrt)) # Errors on the parameters.
    
    ## Finding Beta
    
    def Chi_Squared_Pol(X, Y):
        return np.sum( (Y - Polynomial_Fit(X, A3_3_Curve_Fit[0], A3_3_Curve_Fit[1], A3_3_Curve_Fit[2], A3_3_Curve_Fit[3], A3_3_Curve_Fit[4]))**2 / 0.1**2 ) / 3
    
    def Chi_Squared_Sqrt(X, Y):
        return np.sum( (Y - Local_Function_Fit(X, A3_3_Curve_Fit_Sqrt[0], A3_3_Curve_Fit_Sqrt[1], A3_3_Curve_Fit_Sqrt[2]))**2 / 0.1**2 ) / 3
    
    print('Polynomial Fit {}'.format(Chi_Squared_Pol(A33[:,0], A33[:,1])))
    print('Sqrt Fit {}'.format(Chi_Squared_Sqrt(A33[:,0], A33[:,1])))
    
    Figure = plt.figure(figsize = (10,9))
    Sizes = gridspec.GridSpec(2, 1, height_ratios=[3, 2])
    
    Data_And_Fit = Figure.add_subplot(Sizes[0])
    Data_And_Fit.grid(alpha = 1)
    Data_And_Fit.set_ylabel(r'Dielectric Anisotropy [Dimensionless]', fontsize = 13)
    Data_And_Fit.set_title(r'Dielectric Anisotropy vs. Temperature')
    Data_And_Fit.set_xlabel(r'Temperature T [$^{\circ}$C]', fontsize = 13)
    Data_And_Fit.errorbar(A33[:,0], A33[:,1], yerr = 0.1, xerr = Un_Temperature, fmt = 's k', markersize = 3.0, mec = 'black', mfc = 'white', capsize = 2, elinewidth = 1, label = 'Measured Anisotropies')
    Domain = np.linspace(18.5, 32.5, num = 100, endpoint=True)
    Data_And_Fit.plot(Domain, Polynomial_Fit(Domain, A3_3_Curve_Fit[0], A3_3_Curve_Fit[1], A3_3_Curve_Fit[2], A3_3_Curve_Fit[3], A3_3_Curve_Fit[4]), ls = '-.', label = 'Best Fit (Polynomial)', color = 'orange', alpha = 1)
    Data_And_Fit.plot(Domain, Local_Function_Fit(Domain, A3_3_Curve_Fit_Sqrt[0], A3_3_Curve_Fit_Sqrt[1], A3_3_Curve_Fit_Sqrt[2]), ls = '--', label = 'Best Fit (Landau Critical Exponent Fit)', lw = 3, alpha = 0.4, color = 'green')
    Data_And_Fit.legend()
    
    Residuals = Figure.add_subplot(Sizes[1], label = 'Residuals')
    Residuals.set_ylabel('Residuals', fontsize = 13)
    Residuals.errorbar(A33[:,0], Polynomial_Fit(A33[:,0], A3_3_Curve_Fit[0], A3_3_Curve_Fit[1], A3_3_Curve_Fit[2], A3_3_Curve_Fit[3], A3_3_Curve_Fit[4]) - A33[:,1], yerr = 0.1, fmt = 's k', markersize = 7.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, alpha = 0.5, label = 'Deviation from the Polynomial Fit')
    Residuals.errorbar(A33[:,0], Local_Function_Fit(A33[:,0], A3_3_Curve_Fit_Sqrt[0], A3_3_Curve_Fit_Sqrt[1], A3_3_Curve_Fit_Sqrt[2]) - A33[:,1], yerr = 0.1, fmt = '^ k', markersize = 7.0, mec = 'black', mfc = 'green', capsize = 2, elinewidth = 1, alpha = 0.5, label = 'Deviation from the Landau Fit')
    Residuals.axhline(0, color = 'black')
    Residuals.legend()
    Residuals.grid(alpha = 1)
    
    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (10,9))
    X_Dictionary = {'X_20': [], 'X_22': [], 'X_24': [], 'X_26': [], 'X_28': [], 'X_29': [], 'X_30': [], 'X_31': [], 'X_32': [], 'X_33': []}
    Y_Dictionary = {'Y_20': [], 'Y_22': [], 'Y_24': [], 'Y_26': [], 'Y_28': [], 'Y_29': [], 'Y_30': [], 'Y_31': [], 'Y_32': [], 'Y_33': []}
    Fit_Dictionary = []
    Error_Dictionary = []
    V = [2, 3, 4, 5]
    
    
    for T in B3[:,0]:
        for v in range(len(V)):
            i, = np.where(B3[:,0] == T)
            i = i[0]
            X_Dictionary[list(X_Dictionary.keys())[i]].append(LC_Thickness**2 * (1000 / B2[i,1] + 1000 / B3[i,v + 1]) ) # Conversion to SI here.
            Y_Dictionary[list(Y_Dictionary.keys())[i]].append(epsilon_0 * V[v]**2 * Polynomial_Fit(T, A3_3_Curve_Fit[0], A3_3_Curve_Fit[1], A3_3_Curve_Fit[2], A3_3_Curve_Fit[3], A3_3_Curve_Fit[4]))
    
    for i in range(len(B3[:,0])):
        Curve_Fit, Cov_Fit = curve_fit(S_Line, np.array(X_Dictionary[list(X_Dictionary.keys())[i]]), np.array(Y_Dictionary[list(Y_Dictionary.keys())[i]]), p0 = [1], sigma = epsilon_0 * 1 * np.array(V)**2, absolute_sigma=True)
        Error = np.sqrt(np.diag(Cov_Fit))
        Fit_Dictionary.append(Curve_Fit[0])
        Error_Dictionary.append(Error[0])
    
    Visc_Method_2 = []
    Visc_Errors_Meth_2 = []
    
    for i in range(len(B3[:,0])):
        
        def Viscosity_C(Anisotropy, V_Thresh, Time_Off, Thickness):
            return epsilon_0 * Anisotropy * V_Thresh**2 * Time_Off / Thickness**2
        
        Dielectric_Anisotropy = Local_Function_Fit(B3[i,0], A3_3_Curve_Fit_Sqrt[0], A3_3_Curve_Fit_Sqrt[1], A3_3_Curve_Fit_Sqrt[2])
        Viscosity = Viscosity_C(Dielectric_Anisotropy, B1[i,1] * 10**(-3), B2[i,1] * 10**(-3), LC_Thickness)
        # Viscosity = epsilon_0 *  * B1[i,1]**2 * 10**(-6) * B2[i,1] * 10**(-3) / LC_Thickness**2
        Visc_Method_2.append(Viscosity)
        Error = np.abs(Viscosity_C(Dielectric_Anisotropy + 0.1, (B1[i,1] + 5) * 10**(-3), (B2[i,1] + 0.1) * 10**(-3), LC_Thickness + 0.1*LC_Thickness) - Viscosity)
        Visc_Errors_Meth_2.append(Error)
        
    
    ax[0].errorbar(B3[:,0],  np.array(Fit_Dictionary), yerr = Error_Dictionary, xerr = 0.3, fmt = 's k', markersize = 7.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Method I', alpha = 0.7)
    ax[0].errorbar(B3[:,0],  np.array(Visc_Method_2), yerr = np.array(Visc_Errors_Meth_2), xerr = 0.3, fmt = 's k', markersize = 7.0, mec = 'black', mfc = 'green', capsize = 2, elinewidth = 1, label = 'Method II', alpha = 0.7)
    ax[0].grid(alpha = 1)
    ax[0].set_ylabel(r'Viscosity [Pa s]', fontsize = 13)
    ax[0].set_title(r'Viscosity vs. Temperature')
    ax[0].set_xlabel(r'Temperature T [$^{\circ}$C]', fontsize = 13)
    ax[0].legend(fontsize = 13)
    
    Final_Visc_Values = (np.array(Fit_Dictionary) / np.array(Error_Dictionary)**2 + np.array(Visc_Method_2) / np.array(Visc_Errors_Meth_2)**2) / (1 / np.array(Visc_Errors_Meth_2)**2 + 1/ np.array(Error_Dictionary)**2)
    Final_Errors = np.sqrt(1 / (np.array(Error_Dictionary)**(-2) + np.array(Visc_Errors_Meth_2)**(-2)))
    
    ax[1].errorbar(B3[:,0], Final_Visc_Values, yerr = Final_Errors, xerr = 0.3, fmt = 'h k', markersize = 8.0, mec = 'black', mfc = 'lightskyblue', capsize = 4, elinewidth = 2, alpha = 0.6, label = 'Final Values of the Viscosity')
    ax[1].grid(alpha = 1)
    ax[1].legend(fontsize = 13)
    ax[1].set_ylabel(r'Viscosity [Pa s]', fontsize = 13)
    ax[1].set_title(r'Viscosity vs. Temperature')
    ax[1].set_xlabel(r'Temperature T [$^{\circ}$C]', fontsize = 13)       

    
if Command == 'Threshold Voltage' or Command == 'All':
    
    V_Thresh = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\Voltage_Threshold_Fake_Data.csv', delimiter = ',', dtype=np.dtype(float))
    Uncertainty_on_Data = 0.03
    
    Figure = plt.figure(figsize = (10,9))
    Sizes = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    
    Data_And_Fit = Figure.add_subplot(Sizes[0])
    Data_And_Fit.grid(True)
    Data_And_Fit.errorbar(V_Thresh[:,0],  V_Thresh[:,1], yerr = Uncertainty_on_Data, fmt = 's k', markersize = 7.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Measured Transmissions', alpha = 0.6)
    
    Curve_Fit, Curve_Cov_Fit = curve_fit(Line_Fit, V_Thresh[5:17,0], V_Thresh[5:17,1], sigma = np.full(np.shape(V_Thresh[5:17,1]), Uncertainty_on_Data), absolute_sigma = True, p0 = [-1,1])
    Error_E = np.sqrt(np.diag(Curve_Cov_Fit))
    
    Domain = np.linspace(min(V_Thresh[:,0]), max(V_Thresh[:,0]), num = 2, endpoint = True)
    Mean = np.average(V_Thresh[:5,1], weights = np.full(np.shape(V_Thresh[:5,1]), Uncertainty_on_Data))
    Standard_Dev = np.std(V_Thresh[:5,1])
    Data_And_Fit.scatter( (Mean - Curve_Fit[1]) / Curve_Fit[0], Mean, zorder = 2, label = 'Intercept Point', color = 'black', marker = 'x')
    Data_And_Fit.axhline(Mean, ls = '-.', color = 'steelblue', zorder = 1, label = r'Linear Fit for $V < V_{th}$', lw = 2)
    Data_And_Fit.plot(Domain, Line_Fit(Domain, Curve_Fit[0], Curve_Fit[1]), ls = '-.', color = 'steelblue', zorder = 1, label = r'Linear Fit for $V > V_{th}$', lw = 2)
    Data_And_Fit.axvline((Mean - Curve_Fit[1]) / Curve_Fit[0], ls = '-', ymin = 0, ymax = 1, label = 'Threshold Voltage', lw = 4, alpha = 0.5)
    Data_And_Fit.set_title('Transmission vs. Voltage', fontsize = 14)
    Data_And_Fit.set_xlabel('Voltage [Volts]', fontsize = 13)
    Data_And_Fit.set_ylabel('Transmission (a.u.)', fontsize = 13)
    Threshold_V = (Mean - Curve_Fit[1]) / Curve_Fit[0]
    Un = Mean * np.sqrt( (Standard_Dev**2 + Error_E[1]**2)/(Mean - Curve_Fit[1])**2 + Error_E[0]**2 / Curve_Fit[0]**2 )
    print('The threshold Voltage is: {:.2f} +- {:.2f}'.format(Threshold_V, Un))
    Data_And_Fit.legend(fontsize = 12)
    
    Residuals = Figure.add_subplot(Sizes[1], label = 'Residuals')
    Residuals.axhline(0, color = 'black', ls = '--')
    Residuals.grid(True)
    Residuals.set_ylabel('Residuals', fontsize = 13)
    Residuals.errorbar(V_Thresh[5:17,0], Line_Fit(V_Thresh[5:17,0], Curve_Fit[0], Curve_Fit[1]) - V_Thresh[5:17,1], yerr = Uncertainty_on_Data, fmt = 's k', markersize = 7.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 2, label = 'Measured Transmissions', alpha = 0.5)
    Residuals.errorbar(V_Thresh[:5,0], Mean - V_Thresh[:5,1], yerr = Uncertainty_on_Data, fmt = 's k', markersize = 7.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 2, label = 'Measured Transmissions', alpha = 0.5)
    
if Command == 'TOTO2' or Command == 'All':
    
    Time_On_Start = 0
    Time_On_End = 31
    
    TOTO = np.genfromtxt(r'D:\Physics\University of Manchester\Third Year (BSc)\Semester 6\Liquid Crystals Lab\Data\Fake_Time_On_Time_Off_2.csv', delimiter = ',', dtype=np.dtype(float))
    
    def Time_On_Fit_F(T, Time_On):
        return 1 - np.exp( - T / Time_On)
    
    def Time_Off_Fit_F(T, Time_Off):
        return np.exp(- (T- TOTO[Time_On_End,0]) / Time_Off)
    
    Time_On_Fit, Time_On_Cov_Fit = curve_fit(Time_On_Fit_F, TOTO[Time_On_Start:Time_On_End,0], TOTO[Time_On_Start:Time_On_End,1], sigma = np.full(np.shape(TOTO[Time_On_Start:Time_On_End,1]), 0.01), absolute_sigma = True, p0 = [1])
    Error_E_On = np.sqrt(np.diag(Time_On_Cov_Fit))
    
    Time_Off_Fit, Time_Off_Cov_Fit = curve_fit(Time_Off_Fit_F, TOTO[Time_On_End:,0], TOTO[Time_On_End:,1], sigma = np.full(np.shape(TOTO[Time_On_End:,1]), 0.01), absolute_sigma = True, p0 = [0.04])
    Error_E_Off = np.sqrt(np.diag(Time_On_Cov_Fit))
    
    Figure = plt.figure(figsize = (10,9))
    Sizes = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    
    Data_And_Fit = Figure.add_subplot(Sizes[0])
    Data_And_Fit.errorbar(TOTO[:,0], TOTO[:,1], yerr = 0.03, fmt = 's k', markersize = 5.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Measured Intensities', alpha = 0.5)
    Data_And_Fit.grid(True)
    Data_And_Fit.set_ylabel(r'Normalized Intensity [Volts]', fontsize = 13)
    Data_And_Fit.set_xlabel(r'Time [milliSeconds]', fontsize = 13)
    Data_And_Fit.set_title(r"Normalized Intensity vs. Time")
    
    Time_On_Domain = np.linspace(TOTO[Time_On_Start,0], TOTO[Time_On_End,0], num = 100, endpoint = True)
    Data_And_Fit.plot(Time_On_Domain, Time_On_Fit_F(Time_On_Domain, Time_On_Fit[0]), ls = '--', lw = 2, alpha = 0.8, color = 'darkslategrey')
    
    Time_Off_Domain = np.linspace(TOTO[Time_On_End,0], max(TOTO[:,0]), num = 100, endpoint=True)
    Data_And_Fit.plot(Time_Off_Domain, Time_Off_Fit_F(Time_Off_Domain, Time_Off_Fit[0]), ls = '--', lw = 2, alpha = 0.8, color = 'darkslategrey')

    Residuals = Figure.add_subplot(Sizes[1], label = 'Residuals')
    Residuals.set_ylabel('Residuals', fontsize = 13)
    Residuals.errorbar(TOTO[Time_On_Start:Time_On_End,0], Time_On_Fit_F(TOTO[Time_On_Start:Time_On_End,0], Time_On_Fit[0]) - TOTO[Time_On_Start:Time_On_End,1], yerr = 0.01, fmt = 's k', markersize = 5.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Measured Data', alpha = 0.5)
    Residuals.errorbar(TOTO[Time_On_End:,0], Time_Off_Fit_F(TOTO[Time_On_End:,0], Time_Off_Fit[0]) - TOTO[Time_On_End:,1], yerr = 0.01, fmt = 's k', markersize = 5.0, mec = 'black', mfc = 'orange', capsize = 2, elinewidth = 1, label = 'Measured Data', alpha = 0.5)
    Residuals.axhline(0, color = 'black', ls = '--')
    Residuals.grid(alpha = 1)