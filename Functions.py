import numpy as np

# This function uses the output of the "TruckSyms.m" script to compute
# numerical values for simulation.
# The function arguments are class objects containing global, truck, and
# trailer properties.
# The function outputs are...
# vAx_dot - Truck forward acceleration
# lambda1 - Truck front wheel lateral force
# lambda2 - Truck rear wheel lateral force
# lambda3 - Trailer rear wheel lateral force
def TruckKinetics(Global,Truck,Trailer):
    # Rename all incoming global, truck, and trailer properties to equivalent
    # symbolic variables (purely for visual purposes)
    IzzA = Truck.MOI
    IzzB = Trailer.MOI
    mA = Truck.Mass
    mB = Trailer.Mass
    g = Global.GravAccel
    mu_k = Global.DynFricCoeff
    L1 = Truck.CG2FrontLength
    L2 = Truck.CG2RearLength
    L3 = Trailer.Pivot2CG_Length
    L4 = Trailer.Pivot2RearLength
    phi = Truck.WheelAngle
    phi_dot = Truck.WheelAngVel
    T = Truck.ThrustForce
    vAx = Truck.LongVel
    vAy = Truck.LatVel
    omegaA = Truck.AngVel
    omegaB = Trailer.AngVel
    gamma = Trailer.AngleRel2Truck
    C_D = Truck.DragCoeff
    rho = Global.AirDensity
    A = Truck.DragArea

    # Define relavant position vectors
    r_A_1 = np.array([L1,0,0])
    r_A_2 = np.array([-L2,0,0])
    r_2_3 = np.array([-L4*np.cos(gamma),-L4*np.sin(gamma),0])

    # Define wheel normal forces
    N1 = L2/(L1 + L2)*mA*g
    N2 = L1/(L1 + L2)*mA*g + (L4 - L3)/L4*mB*g
    N3 = L3/L4*mB*g
    
    # Define relevant angular velocity and velocity vectors
    omegaA_vec = np.array([0,0,omegaA])
    omegaB_vec = np.array([0,0,omegaB])
    vA_vec = np.array([vAx,vAy,0])
    v1_vec = vA_vec + np.cross(omegaA_vec,r_A_1)
    v2_vec = vA_vec + np.cross(omegaA_vec,r_A_2)
    v3_vec = \
        vA_vec + \
        np.cross(omegaB_vec,r_2_3) + \
        np.cross(omegaA_vec,r_A_2 + r_2_3)

    # Define unit vectors that align with friction forces
    u_phi = np.array([np.cos(phi),np.sin(phi),0])
    u_gamma = np.array([np.cos(gamma),np.sin(gamma),0])
    i_hat = np.array([1,0,0])

    # Define the aerodynamic drag force
    F_D = -(1/2)*C_D*rho*A*np.sqrt(vAx**2 + vAy**2)*vA_vec
    F_Dx = F_D[0]
    # Define friction forces at the three wheel sets
    Ff1 = -np.tanh(np.dot(v1_vec,u_phi))*mu_k*N1*u_phi
    Ff1x,Ff1y = Ff1[0],Ff1[1]
    Ff2 = -np.tanh(np.dot(v2_vec,i_hat))*mu_k*N2*i_hat
    Ff2x = Ff2[0];
    Ff3 = -np.tanh(np.dot(v3_vec,u_gamma))*mu_k*N3*u_gamma
    Ff3x,Ff3y = Ff3[0],Ff3[1]

    # Compute function outputs unp.sing the analytical expressions derived from
    # "TruckSyms.m"
    vAx_dot = (F_Dx*L1**2*L4**3*np.cos(phi)**3 + F_Dx*L2**2*L4**3*np.cos(phi)**3 + Ff1x*L1**2*L4**3*np.cos(phi)**3 + Ff1x*L2**2*L4**3*np.cos(phi)**3 + Ff2x*L1**2*L4**3*np.cos(phi)**3 + Ff2x*L2**2*L4**3*np.cos(phi)**3 + L1**2*L4**3*T*np.cos(phi)**2 + L2**2*L4**3*T*np.cos(phi)**2 + Ff1y*L1**2*L4**3*np.cos(phi)**2*np.sin(phi) + Ff1y*L2**2*L4**3*np.cos(phi)**2*np.sin(phi) + 2*F_Dx*L1*L2*L4**3*np.cos(phi)**3 + 2*Ff1x*L1*L2*L4**3*np.cos(phi)**3 + 2*Ff2x*L1*L2*L4**3*np.cos(phi)**3 - IzzA*L4**3*phi_dot*vAx*np.sin(phi) - L2**3*L4**3*mB*omegaA**2*np.cos(phi)**3 + IzzB*L1**2*vAx**2*np.cos(gamma)*np.cos(phi)**3 + IzzB*L2**2*vAx**2*np.cos(gamma)*np.cos(phi)**3 + 2*L1*L2*L4**3*T*np.cos(phi)**2 + Ff3x*L1**2*L4**3*np.cos(gamma)**2*np.cos(phi)**3 + Ff3x*L2**2*L4**3*np.cos(gamma)**2*np.cos(phi)**3 - IzzB*L1**2*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 - IzzB*L2**2*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 + 2*Ff1y*L1*L2*L4**3*np.cos(phi)**2*np.sin(phi) + 2*IzzB*L1*L2*vAx**2*np.cos(gamma)*np.cos(phi)**3 + L1**2*L4**3*mA*omegaA*vAy*np.cos(phi)**3 + L2**2*L4**3*mA*omegaA*vAy*np.cos(phi)**3 + L1**2*L4**3*mB*omegaA*vAy*np.cos(phi)**3 + L2**2*L4**3*mB*omegaA*vAy*np.cos(phi)**3 + 2*Ff3x*L1*L2*L4**3*np.cos(gamma)**2*np.cos(phi)**3 - L1**2*L3**2*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 - L2**2*L3**2*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 - 2*L1*L2**2*L4**3*mB*omegaA**2*np.cos(phi)**3 - L1**2*L2*L4**3*mB*omegaA**2*np.cos(phi)**3 + L2**3*L3*L4**2*mB*omegaA**2*np.cos(phi)**3 - 2*IzzB*L1*L2*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 + Ff3y*L1**2*L4**3*np.cos(gamma)*np.cos(phi)**3*np.sin(gamma) + Ff3y*L2**2*L4**3*np.cos(gamma)*np.cos(phi)**3*np.sin(gamma) - L2**2*L4**3*mA*phi_dot*vAx*np.sin(phi) + L1**2*L3**2*mB*vAx**2*np.cos(gamma)*np.cos(phi)**3 + L2**2*L3**2*mB*vAx**2*np.cos(gamma)*np.cos(phi)**3 - L2**2*L4**3*mA*omegaA*vAx*np.cos(phi)**2*np.sin(phi) + 2*L1*L2*L4**3*mA*omegaA*vAy*np.cos(phi)**3 + 2*L1*L2*L4**3*mB*omegaA*vAy*np.cos(phi)**3 - IzzB*L1*L4**2*phi_dot*vAx*np.cos(phi)*np.sin(gamma) - IzzB*L2*L4**2*phi_dot*vAx*np.cos(phi)*np.sin(gamma) - L1**2*L3*L4**3*mB*omegaA**2*np.cos(gamma)*np.cos(phi)**3 - L1**2*L3*L4**3*mB*omegaB**2*np.cos(gamma)*np.cos(phi)**3 - L2**2*L3*L4**3*mB*omegaA**2*np.cos(gamma)*np.cos(phi)**3 - L2**2*L3*L4**3*mB*omegaB**2*np.cos(gamma)*np.cos(phi)**3 - 2*L1*L2*L3**2*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 + L1**2*L3*L4*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 + L2**2*L3*L4*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 - L1**2*L3*L4**2*mB*omegaA*vAy*np.cos(phi)**3 - L2**2*L3*L4**2*mB*omegaA*vAy*np.cos(phi)**3 - L2**3*L3*L4**2*mB*omegaA**2*np.cos(gamma)**2*np.cos(phi)**3 + 2*L1*L2**2*L3*L4**2*mB*omegaA**2*np.cos(phi)**3 + L1**2*L2*L3*L4**2*mB*omegaA**2*np.cos(phi)**3 + 2*Ff3y*L1*L2*L4**3*np.cos(gamma)*np.cos(phi)**3*np.sin(gamma) + 2*L1*L2*L3**2*mB*vAx**2*np.cos(gamma)*np.cos(phi)**3 - L1**2*L3*L4*mB*vAx**2*np.cos(gamma)*np.cos(phi)**3 - L2**2*L3*L4*mB*vAx**2*np.cos(gamma)*np.cos(phi)**3 - L1*L2*L4**3*mA*omegaA*vAx*np.cos(phi)**2*np.sin(phi) - 2*L1*L2*L3*L4**3*mB*omegaA**2*np.cos(gamma)*np.cos(phi)**3 - 2*L1*L2*L3*L4**3*mB*omegaB**2*np.cos(gamma)*np.cos(phi)**3 + 2*L1*L2*L3*L4*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)**3 - 2*L1**2*L3*L4**3*mB*omegaA*omegaB*np.cos(gamma)*np.cos(phi)**3 - 2*L2**2*L3*L4**3*mB*omegaA*omegaB*np.cos(gamma)*np.cos(phi)**3 - 2*L1*L2*L3*L4**2*mB*omegaA*vAy*np.cos(phi)**3 + L1**2*L3*L4**2*mB*omegaA*vAy*np.cos(gamma)**2*np.cos(phi)**3 + L2**2*L3*L4**2*mB*omegaA*vAy*np.cos(gamma)**2*np.cos(phi)**3 - 2*L1*L2*L3*L4*mB*vAx**2*np.cos(gamma)*np.cos(phi)**3 + IzzB*L1*L4*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) + IzzB*L2*L4*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) - 2*L1*L2**2*L3*L4**2*mB*omegaA**2*np.cos(gamma)**2*np.cos(phi)**3 - L1**2*L2*L3*L4**2*mB*omegaA**2*np.cos(gamma)**2*np.cos(phi)**3 - L1**2*L3*L4**2*mB*omegaA*vAx*np.cos(gamma)*np.cos(phi)**3*np.sin(gamma) - L2**2*L3*L4**2*mB*omegaA*vAx*np.cos(gamma)*np.cos(phi)**3*np.sin(gamma) - 4*L1*L2*L3*L4**3*mB*omegaA*omegaB*np.cos(gamma)*np.cos(phi)**3 + 2*L1*L2*L3*L4**2*mB*omegaA*vAy*np.cos(gamma)**2*np.cos(phi)**3 - L1*L3*L4**2*mB*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) + L1*L3**2*L4*mB*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) - L2*L3*L4**2*mB*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) + L2*L3**2*L4*mB*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) - 2*L1*L2*L3*L4**2*mB*omegaA*vAx*np.cos(gamma)*np.cos(phi)**3*np.sin(gamma))/(L4*(IzzA*L4**2*np.cos(phi) + IzzB*L1**2*np.cos(phi)**3 + IzzB*L2**2*np.cos(phi)**3 - IzzA*L4**2*np.cos(phi)**3 + L1**2*L4**2*mA*np.cos(phi)**3 + L1**2*L3**2*mB*np.cos(phi)**3 + L1**2*L4**2*mB*np.cos(phi)**3 + L2**2*L3**2*mB*np.cos(phi)**3 + L2**2*L4**2*mB*np.cos(phi)**3 + 2*IzzB*L1*L2*np.cos(phi)**3 - IzzB*L1**2*np.cos(gamma)**2*np.cos(phi)**3 - IzzB*L2**2*np.cos(gamma)**2*np.cos(phi)**3 + L2**2*L4**2*mA*np.cos(phi) + 2*L1*L2*L4**2*mA*np.cos(phi)**3 + 2*L1*L2*L3**2*mB*np.cos(phi)**3 + 2*L1*L2*L4**2*mB*np.cos(phi)**3 - 2*L1**2*L3*L4*mB*np.cos(phi)**3 - 2*L2**2*L3*L4*mB*np.cos(phi)**3 - L1**2*L3**2*mB*np.cos(gamma)**2*np.cos(phi)**3 - L2**2*L3**2*mB*np.cos(gamma)**2*np.cos(phi)**3 - 2*IzzB*L1*L2*np.cos(gamma)**2*np.cos(phi)**3 - 4*L1*L2*L3*L4*mB*np.cos(phi)**3 + IzzB*L1*L4*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) + IzzB*L2*L4*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) - 2*L1*L2*L3**2*mB*np.cos(gamma)**2*np.cos(phi)**3 + 2*L1**2*L3*L4*mB*np.cos(gamma)**2*np.cos(phi)**3 + 2*L2**2*L3*L4*mB*np.cos(gamma)**2*np.cos(phi)**3 + 4*L1*L2*L3*L4*mB*np.cos(gamma)**2*np.cos(phi)**3));
    lambda1 = ((F_Dx*IzzA*L1*L4**3*np.sin(2*phi))/2 - Ff1y*IzzB*L2**3*L4*np.cos(phi)**2 - Ff1y*IzzB*L1**3*L4*np.cos(phi)**2 + (F_Dx*IzzA*L2*L4**3*np.sin(2*phi))/2 + (Ff1x*IzzA*L1*L4**3*np.sin(2*phi))/2 + (Ff1x*IzzA*L2*L4**3*np.sin(2*phi))/2 + (Ff2x*IzzA*L1*L4**3*np.sin(2*phi))/2 + (Ff2x*IzzA*L2*L4**3*np.sin(2*phi))/2 + (IzzA*IzzB*L4*vAx**2*np.sin(2*gamma))/2 - L1**3*L4**3*T*mA*(np.sin(phi) - np.sin(phi)**3) - L1**3*L4**3*T*mB*(np.sin(phi) - np.sin(phi)**3) - L2**3*L4**3*T*mB*(np.sin(phi) - np.sin(phi)**3) + L2**3*L4**3*mA**2*phi_dot*vAx - Ff1y*L1**3*L4**3*mA*np.cos(phi)**2 - Ff1y*L2**3*L4**3*mA*np.cos(phi)**2 - Ff1y*L1**3*L4**3*mB*np.cos(phi)**2 - Ff1y*L2**3*L4**3*mB*np.cos(phi)**2 + (F_Dx*L2**3*L4**3*mA*np.sin(2*phi))/2 + (Ff1x*L2**3*L4**3*mA*np.sin(2*phi))/2 + (Ff2x*L2**3*L4**3*mA*np.sin(2*phi))/2 + IzzA*L1*L4**3*T*(np.sin(phi) - np.sin(phi)**3) + IzzA*L2*L4**3*T*(np.sin(phi) - np.sin(phi)**3) - IzzB*L1**3*L4*T*(np.sin(phi) - np.sin(phi)**3) - IzzB*L2**3*L4*T*(np.sin(phi) - np.sin(phi)**3) - IzzB*L1**2*L4**2*T*np.cos(phi)*np.sin(gamma) - IzzB*L2**2*L4**2*T*np.cos(phi)*np.sin(gamma) + (IzzB*L2**2*L4*mA*vAx**2*np.sin(2*gamma))/2 - (IzzA*L3*L4**2*mB*vAx**2*np.sin(2*gamma))/2 + (IzzA*L3**2*L4*mB*vAx**2*np.sin(2*gamma))/2 - 3*IzzB*L1*L2**2*L4*T*(np.sin(phi) - np.sin(phi)**3) - 3*IzzB*L1**2*L2*L4*T*(np.sin(phi) - np.sin(phi)**3) + Ff1y*IzzB*L1**3*L4*np.cos(gamma)**2*np.cos(phi)**2 + Ff1y*IzzB*L2**3*L4*np.cos(gamma)**2*np.cos(phi)**2 + IzzB*L1**2*L4**2*T*np.cos(phi)**3*np.sin(gamma) + IzzB*L2**2*L4**2*T*np.cos(phi)**3*np.sin(gamma) - 3*Ff1y*IzzB*L1*L2**2*L4*np.cos(phi)**2 - 3*Ff1y*IzzB*L1**2*L2*L4*np.cos(phi)**2 + IzzA*IzzB*L1*L4*phi_dot*vAx + IzzA*IzzB*L2*L4*phi_dot*vAx + L2**3*L4**3*mA*mB*phi_dot*vAx - (IzzA*L2**2*L4**3*mB*omegaA**2*np.sin(2*phi))/2 + L2**3*L4**3*mA**2*omegaA*vAx*np.cos(phi)**2 - (L2**4*L4**3*mA*mB*omegaA**2*np.sin(2*phi))/2 - 2*L1*L2**2*L4**3*T*mA*(np.sin(phi) - np.sin(phi)**3) - 3*L1**2*L2*L4**3*T*mA*(np.sin(phi) - np.sin(phi)**3) - 3*L1*L2**2*L4**3*T*mB*(np.sin(phi) - np.sin(phi)**3) - 3*L1**2*L2*L4**3*T*mB*(np.sin(phi) - np.sin(phi)**3) + 2*L1**3*L3*L4**2*T*mB*(np.sin(phi) - np.sin(phi)**3) - L1**3*L3**2*L4*T*mB*(np.sin(phi) - np.sin(phi)**3) + 2*L2**3*L3*L4**2*T*mB*(np.sin(phi) - np.sin(phi)**3) - L2**3*L3**2*L4*T*mB*(np.sin(phi) - np.sin(phi)**3) + (L2**3*L4**3*mA**2*omegaA*vAy*np.sin(2*phi))/2 + L1*L2**2*L4**3*mA**2*phi_dot*vAx - 3*Ff1y*L1*L2**2*L4**3*mA*np.cos(phi)**2 - 3*Ff1y*L1**2*L2*L4**3*mA*np.cos(phi)**2 - 3*Ff1y*L1*L2**2*L4**3*mB*np.cos(phi)**2 - 3*Ff1y*L1**2*L2*L4**3*mB*np.cos(phi)**2 + 2*Ff1y*L1**3*L3*L4**2*mB*np.cos(phi)**2 - Ff1y*L1**3*L3**2*L4*mB*np.cos(phi)**2 + 2*Ff1y*L2**3*L3*L4**2*mB*np.cos(phi)**2 - Ff1y*L2**3*L3**2*L4*mB*np.cos(phi)**2 + IzzA*L1*L4**3*mA*phi_dot*vAx + IzzA*L2*L4**3*mA*phi_dot*vAx + IzzB*L2**3*L4*mA*phi_dot*vAx + IzzA*L1*L4**3*mB*phi_dot*vAx + IzzA*L2*L4**3*mB*phi_dot*vAx + (F_Dx*L1*L2**2*L4**3*mA*np.sin(2*phi))/2 + (Ff1x*L1*L2**2*L4**3*mA*np.sin(2*phi))/2 + (Ff2x*L1*L2**2*L4**3*mA*np.sin(2*phi))/2 - 2*Ff1y*L1**3*L3*L4**2*mB*np.cos(gamma)**2*np.cos(phi)**2 + Ff1y*L1**3*L3**2*L4*mB*np.cos(gamma)**2*np.cos(phi)**2 - 2*Ff1y*L2**3*L3*L4**2*mB*np.cos(gamma)**2*np.cos(phi)**2 + Ff1y*L2**3*L3**2*L4*mB*np.cos(gamma)**2*np.cos(phi)**2 + IzzB*L2**3*mA*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzB*L2**3*L4*mA*phi_dot*vAx*np.cos(gamma)**2 + IzzB*L2**3*L4*mA*omegaA*vAx*np.cos(phi)**2 + 6*Ff1y*L1*L2**2*L3*L4**2*mB*np.cos(phi)**2 - 3*Ff1y*L1*L2**2*L3**2*L4*mB*np.cos(phi)**2 + 6*Ff1y*L1**2*L2*L3*L4**2*mB*np.cos(phi)**2 - 3*Ff1y*L1**2*L2*L3**2*L4*mB*np.cos(phi)**2 + IzzB*L1**3*L4*T*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) + IzzB*L2**3*L4*T*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) + (IzzA*L1*L4**3*mA*omegaA*vAy*np.sin(2*phi))/2 + (IzzA*L2*L4**3*mA*omegaA*vAy*np.sin(2*phi))/2 + (IzzA*L1*L4**3*mB*omegaA*vAy*np.sin(2*phi))/2 + (IzzA*L2*L4**3*mB*omegaA*vAy*np.sin(2*phi))/2 + IzzB*L1*L2**2*L4*mA*phi_dot*vAx - 2*IzzA*L1*L3*L4**2*mB*phi_dot*vAx + IzzA*L1*L3**2*L4*mB*phi_dot*vAx - 2*IzzA*L2*L3*L4**2*mB*phi_dot*vAx + IzzA*L2*L3**2*L4*mB*phi_dot*vAx + Ff3x*L2**3*L4**3*mA*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) + 2*IzzB*L1*L2*L4**2*T*np.cos(phi)**3*np.sin(gamma) - (IzzA*L1*L2*L4**3*mB*omegaA**2*np.sin(2*phi))/2 - IzzB*L2**3*mA*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) + 3*Ff1y*IzzB*L1*L2**2*L4*np.cos(gamma)**2*np.cos(phi)**2 + 3*Ff1y*IzzB*L1**2*L2*L4*np.cos(gamma)**2*np.cos(phi)**2 + IzzA*IzzB*L1*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + IzzA*IzzB*L2*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzA*IzzB*L1*L4*phi_dot*vAx*np.cos(gamma)**2 - IzzA*IzzB*L2*L4*phi_dot*vAx*np.cos(gamma)**2 + L2**3*L4**3*mA*mB*omegaA*vAx*np.cos(phi)**2 + (L2**3*L4**3*mA*mB*omegaA*vAy*np.sin(2*phi))/2 + L1*L2**2*L4**3*mA*mB*phi_dot*vAx - 2*L2**3*L3*L4**2*mA*mB*phi_dot*vAx + L2**3*L3**2*L4*mA*mB*phi_dot*vAx + Ff3x*IzzA*L1*L4**3*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) + Ff3x*IzzA*L2*L4**3*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) + (IzzA*L2**2*L3*L4**2*mB*omegaA**2*np.sin(2*phi))/2 - Ff1y*IzzB*L1**2*L4**2*np.cos(phi)*np.sin(gamma)*np.sin(phi) - Ff1y*IzzB*L2**2*L4**2*np.cos(phi)*np.sin(gamma)*np.sin(phi) + 2*L1*L2**2*L4**3*mA**2*omegaA*vAx*np.cos(phi)**2 + L1**2*L2*L4**3*mA**2*omegaA*vAx*np.cos(phi)**2 - (L2**2*L3*L4**2*mA*mB*vAx**2*np.sin(2*gamma))/2 + (L2**2*L3**2*L4*mA*mB*vAx**2*np.sin(2*gamma))/2 - IzzA*IzzB*L4*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma) - (L1*L2**3*L4**3*mA*mB*omegaA**2*np.sin(2*phi))/2 + (L2**4*L3*L4**2*mA*mB*omegaA**2*np.sin(2*phi))/2 - 2*IzzB*L1*L2*L4**2*T*np.cos(phi)*np.sin(gamma) + 6*L1*L2**2*L3*L4**2*T*mB*(np.sin(phi) - np.sin(phi)**3) - 3*L1*L2**2*L3**2*L4*T*mB*(np.sin(phi) - np.sin(phi)**3) + 6*L1**2*L2*L3*L4**2*T*mB*(np.sin(phi) - np.sin(phi)**3) - 3*L1**2*L2*L3**2*L4*T*mB*(np.sin(phi) - np.sin(phi)**3) - IzzA*IzzB*L1*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) - IzzA*IzzB*L2*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) + (L1*L2**2*L4**3*mA**2*omegaA*vAy*np.sin(2*phi))/2 - 2*L1**3*L3*L4**2*T*mB*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) + L1**3*L3**2*L4*T*mB*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) - 2*L2**3*L3*L4**2*T*mB*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) + L2**3*L3**2*L4*T*mB*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) + (L1*L2**2*L4**3*mA*mB*omegaA*vAy*np.sin(2*phi))/2 - (L2**3*L3*L4**2*mA*mB*omegaA*vAy*np.sin(2*phi))/2 - 2*L1*L2**2*L3*L4**2*mA*mB*phi_dot*vAx + L1*L2**2*L3**2*L4*mA*mB*phi_dot*vAx + (L1*L2**3*L3*L4**2*mA*mB*omegaA**2*np.sin(2*phi))/2 + Ff3y*IzzA*L1*L4**3*np.cos(gamma)*np.cos(phi)*np.sin(gamma)*np.sin(phi) + Ff3y*IzzA*L2*L4**3*np.cos(gamma)*np.cos(phi)*np.sin(gamma)*np.sin(phi) - L2**3*L3**2*mA*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) - IzzB*L2**3*L4*mA*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 - 6*Ff1y*L1*L2**2*L3*L4**2*mB*np.cos(gamma)**2*np.cos(phi)**2 + 3*Ff1y*L1*L2**2*L3**2*L4*mB*np.cos(gamma)**2*np.cos(phi)**2 - 6*Ff1y*L1**2*L2*L3*L4**2*mB*np.cos(gamma)**2*np.cos(phi)**2 + 3*Ff1y*L1**2*L2*L3**2*L4*mB*np.cos(gamma)**2*np.cos(phi)**2 + IzzB*L1*L2**2*mA*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + IzzA*L1*L3**2*mB*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + IzzA*L2*L3**2*mB*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzB*L1*L2**2*L4*mA*phi_dot*vAx*np.cos(gamma)**2 + 2*IzzA*L1*L3*L4**2*mB*phi_dot*vAx*np.cos(gamma)**2 - IzzA*L1*L3**2*L4*mB*phi_dot*vAx*np.cos(gamma)**2 + 2*IzzA*L2*L3*L4**2*mB*phi_dot*vAx*np.cos(gamma)**2 - IzzA*L2*L3**2*L4*mB*phi_dot*vAx*np.cos(gamma)**2 + 2*IzzB*L1*L2**2*L4*mA*omegaA*vAx*np.cos(phi)**2 + IzzB*L1**2*L2*L4*mA*omegaA*vAx*np.cos(phi)**2 + 3*IzzB*L1*L2**2*L4*T*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) + 3*IzzB*L1**2*L2*L4*T*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) - (IzzA*L1*L3*L4**2*mB*omegaA*vAy*np.sin(2*phi))/2 - (IzzA*L2*L3*L4**2*mB*omegaA*vAy*np.sin(2*phi))/2 + Ff3x*L1*L2**2*L4**3*mA*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) + (IzzA*L1*L2*L3*L4**2*mB*omegaA**2*np.sin(2*phi))/2 - IzzB*L2**2*L4*mA*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma) + IzzA*L3*L4**2*mB*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma) - IzzA*L3**2*L4*mB*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma) - 2*Ff1y*IzzB*L1*L2*L4**2*np.cos(phi)*np.sin(gamma)*np.sin(phi) - IzzB*L1*L2**2*mA*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) - IzzA*L1*L3**2*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) - IzzA*L2*L3**2*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) + Ff3y*L2**3*L4**3*mA*np.cos(gamma)*np.cos(phi)*np.sin(gamma)*np.sin(phi) + L2**3*L3**2*mA*mB*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + 2*L2**3*L3*L4**2*mA*mB*phi_dot*vAx*np.cos(gamma)**2 - L2**3*L3**2*L4*mA*mB*phi_dot*vAx*np.cos(gamma)**2 + 2*L1*L2**2*L4**3*mA*mB*omegaA*vAx*np.cos(phi)**2 + L1**2*L2*L4**3*mA*mB*omegaA*vAx*np.cos(phi)**2 - 2*L2**3*L3*L4**2*mA*mB*omegaA*vAx*np.cos(phi)**2 + L2**3*L3**2*L4*mA*mB*omegaA*vAx*np.cos(phi)**2 + Ff3y*L1*L2**2*L4**3*mA*np.cos(gamma)*np.cos(phi)*np.sin(gamma)*np.sin(phi) + 2*L2**3*L3*L4**2*mA*mB*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 - L2**3*L3**2*L4*mA*mB*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 - L2**3*L3*L4**3*mA*mB*omegaA**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - L2**3*L3*L4**3*mA*mB*omegaB**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + L1*L2**2*L3**2*mA*mB*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + L2**3*L3*L4*mA*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) + 2*L1*L2**2*L3*L4**2*mA*mB*phi_dot*vAx*np.cos(gamma)**2 - L1*L2**2*L3**2*L4*mA*mB*phi_dot*vAx*np.cos(gamma)**2 - 4*L1*L2**2*L3*L4**2*mA*mB*omegaA*vAx*np.cos(phi)**2 + 2*L1*L2**2*L3**2*L4*mA*mB*omegaA*vAx*np.cos(phi)**2 - 2*L1**2*L2*L3*L4**2*mA*mB*omegaA*vAx*np.cos(phi)**2 + L1**2*L2*L3**2*L4*mA*mB*omegaA*vAx*np.cos(phi)**2 - 6*L1*L2**2*L3*L4**2*T*mB*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) + 3*L1*L2**2*L3**2*L4*T*mB*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) - 6*L1**2*L2*L3*L4**2*T*mB*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) + 3*L1**2*L2*L3**2*L4*T*mB*np.cos(gamma)**2*np.cos(phi)**2*np.sin(phi) - (L1*L2**2*L3*L4**2*mA*mB*omegaA*vAy*np.sin(2*phi))/2 - IzzA*L1*L3*L4*mB*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzA*L2*L3*L4*mB*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzA*L2**2*L3*L4**2*mB*omegaA**2*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) + L2**2*L3*L4**2*mA*mB*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma) - L2**2*L3**2*L4*mA*mB*vAx**2*np.cos(gamma)*np.cos(phi)**2*np.sin(gamma) - L2**4*L3*L4**2*mA*mB*omegaA**2*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) - L1*L2**2*L3**2*mA*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) - 2*IzzB*L1*L2**2*L4*mA*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 - IzzB*L1**2*L2*L4*mA*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 - IzzA*L1*L3*L4**3*mB*omegaA**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzA*L1*L3*L4**3*mB*omegaB**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzA*L2*L3*L4**3*mB*omegaA**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzA*L2*L3*L4**3*mB*omegaB**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + IzzA*L1*L3*L4*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) + IzzA*L2*L3*L4*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) - L2**3*L3*L4*mA*mB*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + IzzB*L2**2*L4**2*mA*omegaA*vAx*np.cos(phi)*np.sin(gamma)*np.sin(phi) + IzzA*L1*L3*L4**2*mB*omegaA*vAy*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) + IzzA*L2*L3*L4**2*mB*omegaA*vAy*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) - L1*L2**2*L3*L4*mA*mB*vAx**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - 2*L2**3*L3*L4**3*mA*mB*omegaA*omegaB*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzA*L1*L2*L3*L4**2*mB*omegaA**2*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) + 4*L1*L2**2*L3*L4**2*mA*mB*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 - 2*L1*L2**2*L3**2*L4*mA*mB*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 + 2*L1**2*L2*L3*L4**2*mA*mB*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 - L1**2*L2*L3**2*L4*mA*mB*omegaA*vAx*np.cos(gamma)**2*np.cos(phi)**2 - L1*L2**2*L3*L4**3*mA*mB*omegaA**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) - L1*L2**2*L3*L4**3*mA*mB*omegaB**2*np.cos(gamma)*np.cos(phi)*np.sin(phi) + L1*L2**2*L3*L4*mA*mB*vAx**2*np.cos(gamma)**3*np.cos(phi)*np.sin(phi) + L2**3*L3*L4**2*mA*mB*omegaA*vAy*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) - 2*IzzA*L1*L3*L4**3*mB*omegaA*omegaB*np.cos(gamma)*np.cos(phi)*np.sin(phi) - 2*IzzA*L2*L3*L4**3*mB*omegaA*omegaB*np.cos(gamma)*np.cos(phi)*np.sin(phi) + IzzB*L1*L2*L4**2*mA*omegaA*vAx*np.cos(phi)*np.sin(gamma)*np.sin(phi) - L1*L2**3*L3*L4**2*mA*mB*omegaA**2*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) - L2**3*L3*L4**2*mA*mB*omegaA*vAx*np.cos(gamma)*np.cos(phi)*np.sin(gamma)*np.sin(phi) - 2*L1*L2**2*L3*L4**3*mA*mB*omegaA*omegaB*np.cos(gamma)*np.cos(phi)*np.sin(phi) - IzzA*L1*L3*L4**2*mB*omegaA*vAx*np.cos(gamma)*np.cos(phi)*np.sin(gamma)*np.sin(phi) - IzzA*L2*L3*L4**2*mB*omegaA*vAx*np.cos(gamma)*np.cos(phi)*np.sin(gamma)*np.sin(phi) + L1*L2**2*L3*L4**2*mA*mB*omegaA*vAy*np.cos(gamma)**2*np.cos(phi)*np.sin(phi) - L1*L2**2*L3*L4**2*mA*mB*omegaA*vAx*np.cos(gamma)*np.cos(phi)*np.sin(gamma)*np.sin(phi))/(L4*(L1 + L2)*(IzzA*L4**2*np.cos(phi) + IzzB*L1**2*np.cos(phi)**3 + IzzB*L2**2*np.cos(phi)**3 - IzzA*L4**2*np.cos(phi)**3 + L1**2*L4**2*mA*np.cos(phi)**3 + L1**2*L3**2*mB*np.cos(phi)**3 + L1**2*L4**2*mB*np.cos(phi)**3 + L2**2*L3**2*mB*np.cos(phi)**3 + L2**2*L4**2*mB*np.cos(phi)**3 + 2*IzzB*L1*L2*np.cos(phi)**3 - IzzB*L1**2*np.cos(gamma)**2*np.cos(phi)**3 - IzzB*L2**2*np.cos(gamma)**2*np.cos(phi)**3 + L2**2*L4**2*mA*np.cos(phi) + 2*L1*L2*L4**2*mA*np.cos(phi)**3 + 2*L1*L2*L3**2*mB*np.cos(phi)**3 + 2*L1*L2*L4**2*mB*np.cos(phi)**3 - 2*L1**2*L3*L4*mB*np.cos(phi)**3 - 2*L2**2*L3*L4*mB*np.cos(phi)**3 - L1**2*L3**2*mB*np.cos(gamma)**2*np.cos(phi)**3 - L2**2*L3**2*mB*np.cos(gamma)**2*np.cos(phi)**3 - 2*IzzB*L1*L2*np.cos(gamma)**2*np.cos(phi)**3 - 4*L1*L2*L3*L4*mB*np.cos(phi)**3 + IzzB*L1*L4*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) + IzzB*L2*L4*np.cos(phi)**2*np.sin(gamma)*np.sin(phi) - 2*L1*L2*L3**2*mB*np.cos(gamma)**2*np.cos(phi)**3 + 2*L1**2*L3*L4*mB*np.cos(gamma)**2*np.cos(phi)**3 + 2*L2**2*L3*L4*mB*np.cos(gamma)**2*np.cos(phi)**3 + 4*L1*L2*L3*L4*mB*np.cos(gamma)**2*np.cos(phi)**3));
    lambda2 = 0 # Not used in the simulation
    lambda3 = 0 # Not used in the simulation

    # Return outputs
    return vAx_dot,lambda1,lambda2,lambda3

# This function used incoming state and controller data to compute and
# output the rates of change of the state vector for integration by
# Python's ODE solver.
# The function arguments are...
# Time - A scalar containing the current simulation time.
# StateVec - A column vector containing the states of the truck-trailer
# system. These states include (in order from top-to-bottom): The truck's
# CG x and y positions in the global frame, the truck's angle, the truck's
# forward velocity measured in its local frame, the trailer angle relative
# to the truck, and the angle and angular velocity of the truck's front
# wheels.
# Control - A class object containing control inputs including the truck's
# thrust force and the torque applied to the truck's steering wheel.
# The remaining arguments are class objects containing global, truck, and
# trailer parameters.
# The function outputs are the rate of the change of the state vector and
# updated truck and trailer class objects for saving results.
def TruckDynModel(Time,StateVec,Control,Global,Truck,Trailer):
    # Extract input information from the control class object
    Truck.WheelTorque = Control.WheelTorque
    Truck.ThrustForce = Control.ThrustForce

    # Extract states information from the state vector
    Truck.CG_PosVec = StateVec[0:1]
    Truck.Angle = StateVec[2]
    Truck.LongVel = StateVec[3]
    Trailer.AngleRel2Truck = StateVec[4]
    Truck.WheelAngle = StateVec[5]
    Truck.WheelAngVel = StateVec[6]

    # Compute kinematic variables of the truck and trailer based on known
    # truck states.
    Truck.LatVel = \
        (Truck.CG2RearLength*Truck.LongVel*np.sin(Truck.WheelAngle))/ \
        (np.cos(Truck.WheelAngle)*(Truck.CG2FrontLength + Truck.CG2RearLength))
    Truck.AngVel = \
        (Truck.LongVel*np.sin(Truck.WheelAngle))/(np.cos(Truck.WheelAngle)* \
        (Truck.CG2FrontLength + Truck.CG2RearLength))
    Trailer.Angle = Truck.Angle + Trailer.AngleRel2Truck
    Trailer.AngVel = \
        -(Truck.LongVel*(Trailer.Pivot2RearLength*np.sin(Truck.WheelAngle) + \
        Truck.CG2FrontLength*np.cos(Truck.WheelAngle)* \
        np.sin(Trailer.AngleRel2Truck) + Truck.CG2RearLength* \
        np.cos(Truck.WheelAngle)*np.sin(Trailer.AngleRel2Truck)))/ \
        (Trailer.Pivot2RearLength*np.cos(Truck.WheelAngle)* \
        (Truck.CG2FrontLength + Truck.CG2RearLength))
    Trailer.LongVel = \
        Truck.LongVel + \
        Trailer.Pivot2CG_Length*Truck.AngVel*np.sin(Trailer.AngleRel2Truck) + \
        Trailer.Pivot2CG_Length*Trailer.AngVel*np.sin(Trailer.AngleRel2Truck)
    Trailer.LatVel = \
        Truck.LatVel - \
        Truck.AngVel*(Truck.CG2RearLength + \
        Trailer.Pivot2CG_Length*np.cos(Trailer.AngleRel2Truck)) - \
        Trailer.Pivot2CG_Length*Trailer.AngVel*np.cos(Trailer.AngleRel2Truck)

    # Compile truck and trailer velocities in to vectors.
    Truck.VelVec = np.array([
        Truck.LongVel,
        Truck.LatVel])
    Trailer.VelVec = np.array([
        Trailer.LongVel,
        Trailer.LatVel])

    # Solve truck-trailer kinetics problem by unp.sing the analytical expressions
    # derived symbolically to compute the wheel lateral forces and truck
    # forward acceleration.
    Truck.LongAcc,Truck.FrontLatForce,Truck.RearLatForce,Trailer.RearLatForce = TruckKinetics(Global,Truck,Trailer)

    # Solve steering wheel kinetics problem. The following components include:
    # Applied steering wheel torque...
    # Steering wheel linear stiffness...
    # Restoring torque resulting from caster angle...
    # Steering wheen linear damping...
    # Artificial force to prevent wheel rotation past max desired angle.
    Truck.WheelAngAcc = (1/Truck.WheelMOI)* \
        (Truck.WheelGearRatio*Truck.WheelTorque - \
        Truck.WheelStiffCoeff*Truck.WheelAngle - \
        (Truck.FrontLatForce*np.cos(Truck.WheelAngle))*Truck.WheelOffset - \
        Truck.WheelDampCoeff*Truck.WheelAngVel - \
        Truck.WheelStopStiff* \
        ((0.5*np.tanh(Truck.WheelAngle*180/np.pi - Truck.MaxWheelAngle) + 0.5) + \
        (-0.5*np.tanh(Truck.WheelAngle*180/np.pi + Truck.MaxWheelAngle) + 0.5)))

    # Define rotation and transformation matrices of the truck and trailer.
    Truck.RotMat = np.array([
        [np.cos(Truck.Angle),-np.sin(Truck.Angle)],
        [np.sin(Truck.Angle),np.cos(Truck.Angle)]])
    Trailer.RotMat = np.array([
        [np.cos(Trailer.Angle),-np.sin(Trailer.Angle)],
        [np.sin(Trailer.Angle),np.cos(Trailer.Angle)]])
    Truck.TransMat_FW = np.transpose(Truck.RotMat)
    Truck.TransMat_BW = np.transpose(Truck.TransMat_FW)

    # Compute position vectors relevant to the truck and trailer for plotting
    # purposes.
    Truck.RearPosVec = \
        Truck.CG_PosVec + Truck.RotMat@np.array([-Truck.CG2RearLength,0])
    Trailer.CG_PosVec = \
        Truck.CG_PosVec + \
        Truck.RotMat@np.array([-Truck.CG2RearLength,0]) + \
        Trailer.RotMat@np.array([-Trailer.Pivot2CG_Length,0])
    Trailer.RearPosVec = \
        Truck.CG_PosVec + \
        Truck.RotMat@np.array([-Truck.CG2RearLength,0]) + \
        Trailer.RotMat@np.array([-Trailer.Pivot2RearLength,0])

    # Compile the state dervative vector.
    StateDerivVec = np.hstack((
        Truck.TransMat_BW@Truck.VelVec,
        Truck.AngVel,
        Truck.LongAcc,
        Trailer.AngVel,
        Truck.WheelAngVel,
        Truck.WheelAngAcc))

    # Display simulation time for tracking progress
    # print(Time)

    # Return function output
    return StateDerivVec

