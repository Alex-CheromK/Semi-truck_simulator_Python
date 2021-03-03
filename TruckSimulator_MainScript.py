# This main script uses the initial conditions specified
# by the user to simulate the motion of a truck-trailer system.
# Simulation, controller, and truck-trailer properties may be
# defined in the "Classes" module.

import numpy as np
from scipy.integrate import solve_ivp
import Classes
from Functions import TruckDynModel
from copy import deepcopy
import os
import time
import matplotlib.pyplot as plt

## Initialize states and results class objects
# Initialize state vector
InitStateVec = np.array([
    0, # Truck x position in global frame (m)
    0, # Truck y position in global frame (m)
    45*np.pi/180, # Truck angle (rad)
    0, # Truck forward velocity in truck frame (m/s)
    0, # Trailer angle relative to truck (rad)
    0, # Front wheel angle (rad)
    0]) # Front wheel angular velocity (rad/s)

# Initialize simulation class objects
Sim = Classes.SimulationProperties()
Control = Classes.ControllerProperties(Sim)
Global = Classes.GlobalProperties()
Truck = Classes.TruckProperties()
Trailer = Classes.TrailerProperties()
TimeStep = []
Sample = []

# Initialize simulation start time
StartTime = 0

## Run discrete-time simulation
for SampleNum in range(Control.NumSamples - 1):
    # Print simulation time for progress tracking
    print(StartTime)

    # Define control inputs for the current sampling period
    Control.ThrustForce = 100 # N
    Control.WheelTorque = 1 # N-m
    # Note: When testing controllers, the above terms may be calculated
    # from a control function.
    
    # Solve dynamics problem for the current smapling period
    Sol = solve_ivp(TruckDynModel,
                    [StartTime,StartTime + Control.SamplingPeriod],
                    InitStateVec,
                    args=(Control,Global,Truck,Trailer),
                    t_eval=np.linspace(StartTime,StartTime + Control.SamplingPeriod,Sim.NumSavesPerPeriod))
    
    # Save results at the end of the current controller sampling period
    if SampleNum == 0:
        TruckDynModel(
            Sol.t[0],
            Sol.y[:,0],
            Control,Global,Truck,Trailer)
        TimeStep.append(Classes.SampleResults(Sol.t[0],Sol.y[:,0],deepcopy(Truck),deepcopy(Trailer)))
        Sample.append(Classes.SampleResults(Sol.t[0],Sol.y[:,0],deepcopy(Truck),deepcopy(Trailer)))
        pass

    # Save simulation results
    for i in range(1,Sim.NumSavesPerPeriod):
        TruckDynModel(
            Sol.t[i],
            Sol.y[:,i],
            Control,Global,Truck,Trailer)
        TimeStep.append(Classes.SampleResults(Sol.t[i],Sol.y[:,i],deepcopy(Truck),deepcopy(Trailer)))
        pass

    # Save sampling period results
    TruckDynModel(
        Sol.t[-1],
        Sol.y[:,-1],
        Control,Global,Truck,Trailer)
    Sample.append(Classes.SampleResults(Sol.t[-1],Sol.y[:,-1],deepcopy(Truck),deepcopy(Trailer)))
    
    
    # Update ODE solver inputs for the next sampling period
    StartTime = StartTime + Control.SamplingPeriod
    InitStateVec = Sol.y[:,-1]
    pass

# Save both time step and sample results into a single structure
Results = Classes.SimulationResults(TimeStep,Sample)

## Generate output plots
# Extract plot arrays from results class objects
Time = np.zeros(Control.NumSamples)
TruckLongVel = np.zeros(Control.NumSamples)
TruckLatVel = np.zeros(Control.NumSamples)
TruckWheelAngle = np.zeros(Control.NumSamples)
TrailerAngleRel2Truck = np.zeros(Control.NumSamples)
for SampleNum in range(Control.NumSamples):
    Time[SampleNum] = Results.Sample[SampleNum].Time
    TruckLongVel[SampleNum] = Results.Sample[SampleNum].Truck.LongVel
    TruckLatVel[SampleNum] = Results.Sample[SampleNum].Truck.LatVel
    TruckWheelAngle[SampleNum] = Results.Sample[SampleNum].Truck.WheelAngle*180/np.pi
    TrailerAngleRel2Truck[SampleNum] = Results.Sample[SampleNum].Trailer.AngleRel2Truck*180/np.pi
    pass

# Generate time-series plots
fig, axs = plt.subplots(nrows=2,ncols=2,sharex=True)
axs[0,0].plot(Time,TruckLongVel)
axs[1,0].plot(Time,TruckLatVel)
axs[0,1].plot(Time,TruckWheelAngle)
axs[1,1].plot(Time,TrailerAngleRel2Truck)
axs[0,0].set_ylabel('Truck longitudinal speed (m/s)')
axs[1,0].set_ylabel('Truck lateral speed (m/s)')
axs[0,1].set_ylabel('Truck wheel angle (deg)')
axs[1,1].set_ylabel('Trailer angle relative to truck (deg)')
axs[1,0].set_xlabel('Time (sec)')
axs[1,1].set_xlabel('Time (sec)')
plt.show()