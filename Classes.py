# This module defines classes that store all the information
# that is necessary for simulating the truck-trailer system.

import numpy as np

class SimulationProperties:
    def __init__(self):
        # Simulation end time (sec)
        self.EndTime = 60
        # Time interval for saving results (sec).
        # Note: Sim.SavePeriod must be at most half of the controller sampling
        # period. It will automatically be reduced later on.
        self.SavePeriod = 0.1
        # Compute the number of saved samples
        self.NumSavePoints = int(self.EndTime/self.SavePeriod) + 1
        self.NumSavesPerPeriod = 0
        # Playback rate of the animation. Example: Enter "2" to play the animation
        # at twice the rate of the simulated physics.
        # Enter "0" for no animation.
        self.AnimPlayRate = 10
        pass
    pass

class ControllerProperties:
    def __init__(self,Sim):
        # Controller sampling period (sec)
        self.SamplingPeriod = 10
        # Compute number of controller sample points
        self.NumSamples = int(Sim.EndTime/self.SamplingPeriod) + 1
        # Initial control values
        self.ThrustForce = 0
        self.WheelTorque = 0
        # Set Sim.SavePeriod to at most half of Control.SamplingPeriod
        Sim.SavePeriod = min(Sim.SavePeriod,self.SamplingPeriod/2)
        Sim.NumSavesPerPeriod = int(self.SamplingPeriod/Sim.SavePeriod) + 1
        pass
    pass

class GlobalProperties:
    def __init__(self):
        # Gravitational acceleration (m/s^2)
        self.GravAccel = 9.81
        # Rolling friction coefficient between tires and road
        self.DynFricCoeff = 0.002
        # Density of air (kg/m^3)
        self.AirDensity = 1.2
        pass
    pass

class TruckProperties:
    def __init__(self):
        # Mass (kg)
        self.Mass = 1500
        # Moment of inertia about vertical axis (kg-m^2)
        self.MOI = 10000
        # Length from CG to front wheels (m)
        self.CG2FrontLength = 1
        # Length from CG to rear wheels (m)
        self.CG2RearLength = 1
        # Front wheel system moment of inertia about vertical axis (kg-m^2)
        self.WheelMOI = 1
        # Gear ratio from steering wheel to front wheels
        self.WheelGearRatio = 1
        # Wheel system linear rotational stiffness coefficient (N/rad)
        self.WheelStiffCoeff = 0
        # Wheel system pivot offset due to caster angle (m)
        self.WheelOffset = 0.25
        # Wheel system linear rotational damping coefficient (N/(rad/s))
        self.WheelDampCoeff = 10
        # Truck-trailer system aerodynamic drag coefficient
        self.DragCoeff = 5
        # Truck-trailer system reference drag area (m^2)
        self.DragArea = 3
        # Artifical maximum torque for limiting wheel rotation beyond maximum
        # allowable value (N-m)
        self.WheelStopStiff = 100
        # Maximum allowable wheel angle (deg)
        self.MaxWheelAngle = 10
        pass
    pass

class TrailerProperties:
    def __init__(self):
        # Mass (kg)
        self.Mass = 3000
        # Moment of inertia about vertical axis (kg-m^2)
        self.MOI = 20000
        # Length from CG to pivot point (m)
        self.Pivot2CG_Length = 6
        # Length from CG to rear wheels (m)
        self.Pivot2RearLength = 12
        pass
    pass

class TimeStepResults:
    def __init__(self,Time,StateVec,Truck,Trailer):
        self.Time = Time
        self.StateVec = StateVec
        self.Truck = Truck
        self.Trailer = Trailer
        pass
    pass

class SampleResults:
    def __init__(self,Time,StateVec,Truck,Trailer):
        self.Time = Time
        self.StateVec = StateVec
        self.Truck = Truck
        self.Trailer = Trailer
        pass
    pass

class SimulationResults:
    def __init__(self,TimeStep,Sample):
        self.TimeStep = TimeStep
        self.Sample = Sample
        pass
    pass