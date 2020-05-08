import numpy as np
# from numpy import array, zeros, round, 
from numpy.random import rand

from scipy.integrate import ode
from scipy.optimize import fsolve

def hovorka_parameters(BW):
    """
    PATIENT PARAMETERS
    BW - body weight in kilos
    """

    # Patient-dependent parameters:
    V_I = 0.12*BW              # Insulin volume [L]
    V_G = 0.16*BW              # Glucose volume [L]
    F_01 = 0.0097*BW           # Non-insulin-dependent glucose flux [mmol/min]
    EGP_0 = 0.0161*BW          # EGP extrapolated to zero insulin concentration [mmol/min]

    # Patient-independent(?) parameters:
    S_IT = 51.2e-4             # Insulin sensitivity of distribution/transport [L/min*mU]
    S_ID = 8.2e-4              # Insulin sensitivity of disposal [L/min*mU]
    S_IE = 520e-4              # Insluin sensitivity of EGP [L/mU]

    tau_G = 40                 # Time-to-maximum CHO absorption [min]
    tau_I = 55                 # Time-to-maximum of absorption of s.c. injected short-acting insulin [min]

    A_G = 0.8                  # CHO bioavailability [1]
    k_12 = 0.066               # Transfer rate [min]

    k_a1 = 0.006               # Deactivation rate of insulin on distribution/transport [1/min]
    k_b1 = S_IT*k_a1           # Activation rate of insulin on distribution/transport
    k_a2 = 0.06                # Deactivation rate of insulin on dsiposal [1/min]
    k_b2 = S_ID*k_a2           # Activation rate of insulin on disposal
    k_a3 = 0.03                # Deactivation rate of insulin on EGP [1/min]
    k_b3 = S_IE*k_a3           # Activation rate of insulin on EGP

    k_e = 0.138                # Insulin elimination from Plasma [1/min]

    # Summary of the patient's values:
    P = [tau_G, tau_I, A_G, k_12, k_a1, k_b1, k_a2, k_b2, k_a3, k_b3, k_e, V_I, V_G, F_01, EGP_0]

    return P


def hovorka_model(t, x, u, D, P): ## This is the ode version
    """HOVORKA DIFFERENTIAL EQUATIONS
    # t:    Time window for the simulation. Format: [t0 t1], or [t1 t2 t3 ... tn]. [min]
    # x:    Initial conditions
    # u:    Amount of insulin insulin injected [mU/min]
    # D:    CHO eating rate [mmol/min]
    # P:    Model fixed parameters
    #
    # Syntax :
    # [T, X] = ode15s(@Hovorka, [t0 t1], xInitial0, odeOptions, u, D, p);
    """
    # TODO: update syntax in docstring

    # u, D, P = args

    # Defining the various equation names
    D1 = x[ 0 ]               # Amount of glucose in compartment 1 [mmol]
    D2 = x[ 1 ]               # Amount of glucose in compartment 2 [mmol]
    S1 = x[ 2 ]               # Amount of insulin in compartment 1 [mU]
    S2 = x[ 3 ]               # Amount of insulin in compartment 2 [mU]
    Q1 = x[ 4 ]               # Amount of glucose in the main blood stream [mmol]
    Q2 = x[ 5 ]               # Amount of glucose in peripheral tissues [mmol]
    I =  x[ 6 ]                # Plasma insulin concentration [mU/L]
    x1 = x[ 7 ]               # Insluin in muscle tissues [1], x1*Q1 = Insulin dependent uptake of glucose in muscles
    x2 = x[ 8 ]               # [1], x2*Q2 = Insulin dependent disposal of glucose in the muscle cells
    x3 = x[ 9 ]              # Insulin in the liver [1], EGP_0*(1-x3) = Endogenous release of glucose by the liver
    C = x[10]

    # Unpack data
    tau_G = P[ 0 ]               # Time-to-glucose absorption [min]
    tau_I = P[ 1 ]               # Time-to-insulin absorption [min]
    A_G = P[ 2 ]                 # Factor describing utilization of CHO to glucose [1]
    k_12 = P[ 3 ]                # [1/min] k_12*Q2 = Transfer of glucose from peripheral tissues (ex. muscle to the blood)
    k_a1 = P[ 4 ]                # Deactivation rate [1/min]
    k_b1 = P[ 5 ]                # [L/(mU*min)]
    k_a2 = P[ 6 ]                # Deactivation rate [1/min]
    k_b2 = P[ 7 ]                # [L/(mU*min)]
    k_a3 = P[ 8 ]                # Deactivation rate [1/min]
    k_b3 = P[ 9 ]               # [L/(mU*min)]
    k_e = P[ 10 ]                # Insulin elimination rate [1/min]
    V_I = P[ 11 ]                # Insulin distribution volume [L]
    V_G = P[ 12 ]                # Glucose distribution volume [L]
    F_01 = P[ 13 ]               # Glucose consumption by the central nervous system [mmol/min]
    EGP_0 = P[ 14 ]              # Liver glucose production rate [mmol/min]

    # If some parameters are not defined
    if len(P) == 15:
        ka_int = 0.073
        R_cl = 0.003
        # R_thr = 9
        R_thr = 14
    elif len(P) == 18:
        R_cl = P[16]
        ka_int = P[15]
        R_thr = P[17]

    # Certain parameters are defined
    U_G = D2/tau_G             # Glucose absorption rate [mmol/min]
    U_I = S2/tau_I             # Insulin absorption rate [mU/min]

    # Constitutive equations
    G = Q1/V_G                 # Glucose concentration [mmol/L]

    # if (G>=4.5):
    #     F_01c = F_01           # Consumption of glucose by the central nervous system [mmol/min
    # else:
    #     F_01c = F_01*G/4.5     # Consumption of glucose by the central nervous system [mmol/min]

    F_01s = F_01/0.85
    F_01c = F_01s*G / (G + 1)

    # if (G>=9):
        # F_R = 0.003*(G-9)*V_G  # Renal excretion of glucose in the kidneys [mmol/min]
    # else:
        # F_R = 0                # Renal excretion of glucose in the kidneys [mmol/min]

    if (G >= R_thr):
        F_R = R_cl*(G - R_thr)*V_G  # Renal excretion of glucose in the kidneys [mmol/min]
    else:
        F_R = 0                # Renal excretion of glucose in the kidneys [mmol/min]

    # Mass balances/differential equations
    xdot = np.zeros (11);

    xdot[ 0 ] = A_G*D-D1/tau_G                                # dD1
    xdot[ 1 ] = D1/tau_G-U_G                                  # dD2
   
    xdot[ 2 ] = u-S1/tau_I                                    # dS1
    xdot[ 3 ] = S1/tau_I-U_I                                  # dS2
   
    xdot[ 4 ] = -(F_01c+F_R) - x1*Q1 + k_12*Q2 + U_G + max(EGP_0*(1-x3), 0)   # dQ1
    xdot[ 5 ] = x1*Q1-(k_12+x2)*Q2                            # dQ2

    xdot[ 6 ] = U_I/V_I-k_e*I                                 # dI
    
    xdot[ 7 ] = k_b1*I-k_a1*x1                                # dx1
    xdot[ 8 ] = k_b2*I-k_a2*x2                                # dx2
    xdot[ 9 ] = k_b3*I-k_a3*x3                                # dx3

    # ===============
    # CGM delay
    # ===============
    xdot[10] = ka_int*(G - C)


    return xdot


def hovorka_model_tuple(x, *pars):
    """HOVORKA DIFFERENTIAL EQUATIONS without time variable
    # t:    Time window for the simulation. Format: [t0 t1], or [t1 t2 t3 ... tn]. [min]
    # x:    Initial conditions
    # u:    Amount of insulin insulin injected [mU/min]
    # D:    CHO eating rate [mmol/min]
    # P:    Model fixed parameters
    #
    """
    # TODO: update syntax in docstring
    # import numpy as np

    u, D, P = pars

    t = 0

    xdot = hovorka_model(t, x, u, D, P)

    return xdot


def run_simulation():
    
    # Setting up meals (should be improved)
    meals = np.array([40, 80, 60, 30])
    meal_times = np.array([8*60, 12*60, 18*60, 22*60]) 
    meal_vector = np.zeros(1440)

    meal_vector[meal_times] = meals * 1000/180

    boluses = meal_vector

    # Initial values for parameters
    init_basal = 6.43
    P = hovorka_parameters(70)
    initial_pars = (init_basal, 0, P)

    # Bolus carb factor -- [g/U]
    carb_factor = 25

    # Initial value
    X0 = fsolve(hovorka_model_tuple, np.zeros(11), args=initial_pars)
    
    # Simulation setup
    integrator = ode(hovorka_model)
    # integrator.set_integrator('vode', method='bdf', order=4)
    integrator.set_integrator('dopri5')
    integrator.set_initial_value(X0, 0)

    action = init_basal

    blood_glucose_value = []

    for i in range(1440):

        if meal_vector[i] > 0:
            insulin_rate = action + np.round(max(meal_vector[i] * (180 / carb_factor), 0), 1)
        else:
            insulin_rate = action

        # Updating the carb and insulin parameters in the model
        integrator.set_f_params(insulin_rate, meal_vector[i], P)

        # Integration step
        integrator.integrate(integrator.t + 1)

        blood_glucose_value.append(int(integrator.y[-1]*18))

    # Returning blood glucose value
    return blood_glucose_value


if __name__ == "__main__":
        
    # app.run()
    bg = run_simulation()

    # import matplotlib.pyplot as plt

    # plt.figure()
    # plt.plot(bg)
    # plt.show()