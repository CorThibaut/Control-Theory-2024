import numpy as np

import matplotlib.pyplot as plt
from IPython.display import display, clear_output

from package_DBR import Bode, Process

#-----------------------------------        
def LL_RT(MV,Kp,Tlead,Tlag,Ts,PV,PVInit=0,method='EBD'):
    
    """
    The function "LL_T" needs to be included in a "for or while loop".
    
    :MV: input vector
    :Kp: process gain
    :Tlead: lead time constant [s]
    :Tlag: lag time constant [s]
    :Ts: sampling period [s]
    :PV: output vector
    :PVInit: (optional: default value is 0)
    :method: discretization method (optional: default value is 'EBD')
        EBD: Euler Backward difference
        EFD: Euler Forward difference
    
    The function "LEADLAG_RT" appends a value to the output vector "PV".
    The appended value is obtained from a recurrent equation that depends on the discretization method.
    """    
    
    if (Tlag != 0):
        K = Ts/Tlag
        if len(PV) == 0:
            PV.append(PVInit)
        else: # MV[k+1] is MV[-1] and MV[k] is MV[-2]
            if method == 'EBD':
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+(Tlead/Ts))*MV[-1])-((Tlead/Ts)*MV[-2])))
            elif method == 'EFD':
                PV.append((1-K)*PV[-1] + K*Kp*(((Tlead/Ts)*MV[-1]) + ((1-(Tlead/Ts))*MV[-2])))         
            else:
                PV.append((1/(1+K))*PV[-1] + (K*Kp/(1+K))*(((1+(Tlead/Ts))*MV[-1])-((Tlead/Ts)*MV[-2])))
    else:
        PV.append(Kp*MV[-1])
        
#----------------------------------------------------------------------------------------------------------
def PID_RT(SP, PV, Man , MVMan, MVFF, Kc, Ti, Td, alpha, Ts, MVMin,  MVMax, MV, MVP, MVI, MVD, E, ManFF=False, PVInit=0, method='EBD-EBD'):
   
    """
    The function "PID_RT" needs to be included in a "for or while loop".
    
    :SP: SP (or SetPoint) vector
    :PV: PV (or Process Value) vector
    :Man: Man (or Manual controller mode) vector (True or False)
    :MVMan: MVMan (or Manual value for MV) vector
    :MVFF: MVFF (or Feedforward) vector
    
    :Kc: controller gain
    :Ti: integral time constant [s]
    :Td: derivative time constant [s]
    :alpha: Tfd = alpha*Td where Tfd is the derivative filter time constant [s]
    :Ts: sampling period [s]
    
    :MVMin: minimum value for MV (used for saturation and anti wind-up)
    :MVMax: maximum value for MV (used for saturation and anti wind-up)
    
    :MV: MV (or Manipulated Value) vector
    :MVP: MVP (or Proportional part of MV) vector
    :MVI: MVI (or Integral part of MV) vector
    :MVD: MVD (or Derivative part of MV) vector
    :E: E (or control Error) vector
    
    :ManFF: Activated FF in manual mode (optional: default boolean value is False)
    :PVInit: Initial value for PV (optional: default value is 0): used if PID_RT is ran first in the sequence and no value of PV is available yet.
    
    :method: discretization method (optional: default value is 'EBD')
        EBD-EBD: EBD for integral action and EBD for derivative action
    
    The function "PID_RT" appends new values to the vectors "MV", "MVP", "MVI" and "MVD".
    The appended values are based on the PID algorithm, the controller mode, and feedforward. Note that saturation of "MV" within the limits [MVMin, MVMax] is implemented with anti wind-up.
    """

    #Initialisation of E
    if len(PV)==0:
        E.append(SP[-1]-PVInit)
    else:
        E.append(SP[-1]-PV[-1])
        
    #Proprotional action
    MVP.append(Kc*E[-1])
    
    #Initialisation of integral action
    if len(MVI)==0:
        MVI.append((Kc*Ts/Ti)*E[-1])
    else :
        MVI.append(MVI[-1] + (Kc*Ts/Ti)*E[-1])
                
    #Initialisation of derivative action
    if len(MVD)==0:
        MVD.append(0)
    else:
        MVD.append((Td/(Td+Ts)*MVD[-1]) + ((Kc*Td)/(Td+Ts))*(E[-1]-E[-2]))
          
    #Manual mode + anti wind-up
    if Man[-1] == True:
        if ManFF :
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1]
        else : 
            MVI[-1] = MVMan[-1] - MVP[-1] - MVD[-1] - MVFF[-1]
            
    #Automatic mode + anti wind-up
    else:
        if method == "EBD-EBD": # MV[k-1] is MV[-1] and MV[k] is MV[-2]
            if (Ti !=0 and Td+Ts !=0 and 1>(Td/(Td+Ts))>0):
                if (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1] > MVMax):
                    MVI[-1] = MVMax - MVP[-1] - MVD[-1] - MVFF[-1]
                elif (MVP[-1] + MVI[-1] + MVD[-1] + MVFF[-1] < MVMin):
                    MVI[-1] = MVMin - MVP[-1] - MVD[-1] - MVFF[-1]         
            else:
                MV.append(Kc*E[-1])
                
    #Addition of all actions
    MV.append(MVP[-1]+MVI[-1]+MVD[-1]+MVFF[-1])

#-----------------------------------------------------------------------------------------------------
def IMCTuning(k, Tlag1, Tlag2=0, theta=0, gamma=0.5, process='SOPDT'):
    
    """
    IMCTuning(K, Tlag1, Tlag2=0, theta=0, gamma=0.5, process='SOPDT')
    
    The function "IMCTuning" computes the IMC PID tuning parameters for SOPDT processes.
    For a SOPDT process we impleent the line I of the IMC tuning table (with T3=0)
    
    :K: process gain
    :Tlag1: first (or main) lag time constant [s]
    :Tlag2: second lag time constant [s]
    :theta: process delay
    :gamma: loop response time as a ratio of T1
    """
    
    Kc = ((Tlag1+Tlag2)/((gamma*Tlag1) + theta))/k
    Ti = Tlag1 + Tlag2
    Td = (Tlag1*Tlag2)/(Tlag1+Tlag2)
    
    return Kc, Ti, Td


#-----------------------------------------------------------------------------------------------------
# def Margins()


#-----------------------------------------------------------------------------------------------------
class Controller():
        
    def __init__(self, parameters):
        
        self.parameters = parameters
        self.parameters['Kc'] = parameters['Kc'] if 'Kc' in parameters else 2.0
        self.parameters['Ti'] = parameters['Ti'] if 'Ti' in parameters else 100.0
        self.parameters['Td'] = parameters['Td'] if 'Td' in parameters else 10.0
        self.parameters['alpha'] = parameters['alpha'] if 'alpha' in parameters else 1.0