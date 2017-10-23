
#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/18/2017
# Latest update: 09/01/2017
#
# Comment: This class processes the input data and provides the initial information for the simulation. 
#
#####################################################################

class Initialization:

    def __init__(self):
        # -- Import libs/classes
        import numpy as np
        import Discretization_Class

        self.Disc = Discretization_Class.Discretization()

        self.Disc.Discretize()

        print(" ========== Initialization Class ==========")
        print(" Initialization ...")

        self.Q  = np.zeros(self.Disc.N_Cells,     dtype=np.float64 )
        self.V  = np.zeros(self.Disc.N_Cells,     dtype=np.float64 )
        self.L  = np.zeros(self.Disc.N_Cells,     dtype=np.float64 )
        self.Z  = np.zeros(self.Disc.N_Cells,     dtype=np.float64 )
        self.Z_F= np.zeros(self.Disc.N_Cells*2+1, dtype=np.float64 )
        self.M  = np.zeros(self.Disc.N_Cells,     dtype=np.float64 )
        self.B  = np.zeros(self.Disc.N_Cells,     dtype=np.float64 )
        self.X  = np.zeros(self.Disc.N_Cells,     dtype=np.float64 )
        self.X_F= np.zeros(self.Disc.N_Cells*2+1, dtype=np.float64 )
        self.S_Cell= np.zeros(self.Disc.N_Cells,     dtype=np.float64 )

        self.Total_Time = self.Disc.Total_Time
        self.Time_Step  = self.Disc.Time_Step
        self.h_dw = self.Disc.h_dw
        self.Q_Up = self.Disc.Q_Up
        self.N_Cells    = self.Disc.N_Cells
        self.Output_Dir = self.Disc.Output_Dir

        self.L[:]   = self.Disc.Length_Cell[:]
        self.Z[:]   = self.Disc.Z_Cell[:]
        self.Z_F[:] = self.Disc.Z_Full[:]
        self.M[:]   = self.Disc.Manning_Cell[:]
        self.B[:]   = self.Disc.Width_Cell[:]
        self.X[:]   = self.Disc.X_Disc[:]
        self.X_F[:] = self.Disc.X_Full[:]
        self.S_Cell[:] = self.Disc.S_Cell[:]

        for ii in range( self.Disc.N_Cells ):
            #self.V[ii] = (self.Disc.V_in)* (1+  ( float(ii)/(self.Disc.N_Cells))* (self.Disc.V_ratio)    )
            self.V[ii] = ( (0.33 -self.Z[ii])* self.B[ii] )*self.L[ii]

        for ii in range( self.Disc.N_Cells ):
            self.Q[ii] = 0
            #self.Q[ii] = self.Disc.Q_Up * (self.Disc.N_Cells - ii)/self.Disc.N_Cells
        #self.Q[0] = self.Disc.Q_Up

        del self.Disc

        print(" ========== Initialization Class Ends. ==========")
        print()

    # <modify> Use geometry function and elevation functino, instead of direct definition in the solver function, to make the code a more general one. 
    #  def Geometry_func(self, V, L): # Returns the cross-section area: A
    #    return V/L

    #  def Elevation_func(self, A, B , Z):  # Returns the free surface elevation: eta
    #    return A/B + Z
