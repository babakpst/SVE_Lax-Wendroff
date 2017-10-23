

#####################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    07/31/2017
# Latest update: 09/20/2017
#
# Comment: This class conducts the numerical solution of the Saint-Venant equation
#
#####################################################################

class Solver:

    def __init__(self):
        # -- Import libs/classes
        import numpy as np
        import sys
        import matplotlib.pyplot as plt
        import os

        import Initialization_Class as Initial
        import Visualization_Class  as Visual

        Gravity = 9.81
        Draw = Visual.Visualization()
        Ex = Initial.Initialization()

        print(" ========== Solver Class ==========")
        print(" Solving the DE ...")

        print(" Allocating memory ...")
        print(" Initialization ... ")

        Plot = 200
        h_downstream = Ex.h_dw
        Q_upstream = Ex.Q_Up
        g = Gravity
        DT = Ex.Time_Step
        N_Cells = Ex.N_Cells
        nx = N_Cells
        N_Steps = int(Ex.Total_Time/Ex.Time_Step)
        nt = N_Steps

        h = np.zeros ( nx )
        uh = np.zeros ( nx )

        hm = np.zeros ( nx - 1 )
        uhm = np.zeros ( nx - 1 )

        x = np.zeros ( nx )
        z = np.zeros ( nx )

        L = np.zeros ( nx )
        M = np.zeros ( nx )
        B = np.zeros ( nx )

        S_f = np.zeros ( nx )
        S_0 = np.zeros ( nx )
        S   = np.zeros ( nx )
        S_f_m = np.zeros ( nx-1 )
        S_m   = np.zeros ( nx-1 )

        t = np.zeros ( nt + 1 )

        #h_array = np.zeros ( [ nx, nt + 1 ] )
        #uh_array = np.zeros ( [ nx, nt + 1 ] )

        x[:] = Ex.X[:]
        t = np.linspace ( 0, Ex.Total_Time, nt + 1 )

        z[:]   = Ex.Z[:]
        S_0[:] = -Ex.S_Cell[:]
        L[:]   = Ex.L[:]
        M[:]   = Ex.M[:]
        B[:]   = Ex.B[:]
        n_manning  = M[0]
        print(" n_manning: ", n_manning)
        
        dt = DT
        dx = L[0]  #x_length / float ( nx - 1 )
        B = B[0]

        print(" B: ", B)
        print(" DX: ", dx)
        print(" dt: ", dt)
        #for i in range(nx):
        #    print("  L:", L[i] )
        #    print("S_0:", S_0[i])
        #   print("  x:", x[i])

        print(" h downstream: ", h_downstream )

        h, uh = initial_conditions (nx, nt, h, uh, x, z )
        h, uh = boundary_conditions (nx, nt, h, uh, t[0], h_downstream, Q_upstream, B )

        #h_array[0:nx,0] = h[0:nx]
        #uh_array[0:nx,0] = uh[0:nx]

        print(" Time marching ... ")
        for it in range ( 1, nt + 1 ):
            print(it)

            if (it%Plot) == 1:

                fig = plt.figure()

                ax1 = fig.add_subplot(211)
                ax1.grid(True, color='k')   
                #ax1.plot(X_Arr, Q_Arr, label ="Water flow" , color = "c", linewidth = 2.0)
                ax1.fill_between (x, z[:], z[:] +h[:])
                #plt.fill_between ( x, z[:], h[:] )
                
                title_string = ( 'H(T) - Time = %8.2f' % ( t[it] ) )
                plt.title(title_string, fontsize = 16)

                plt.xlabel ( 'X',  fontsize=12 )
                plt.ylabel ( 'H(X,T)',  fontsize=12 )


                #plt.axis ( [ 0.0, 2000, 0, 10 ] )
                #plt.fill_between ( x, z[:], z[:]+h[:] )

                ax2 = fig.add_subplot(212)
                ax2.grid(True, color='k')   
                #ax1.plot(X_Arr, Q_Arr, label ="Water flow" , color = "c", linewidth = 2.0)
                ax2.plot (x, uh[:], label ="Water flow" , color = "c", linewidth = 2.0)
                
                title_string = ( 'UH(T) - Time = %8.2f' % ( t[it] ) )
                plt.title(title_string, fontsize = 16)

                plt.xlabel ( 'X',  fontsize=12 )
                plt.ylabel ( 'UH(X,T)',  fontsize=12 )

                mng = plt.get_current_fig_manager()
                mng.resize(*mng.window.maxsize())

                #plt.show ( )
                #plt.show(block=False) # <modify> See why the execution stops when the the command gets here. 

                FileName = os.path.join(Ex.Output_Dir, 'Time_' +str(t[it])+"_s" +'.jpg')
                #print(FileName)
                plt.savefig(FileName)
                #savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format=None, transparent=False, bbox_inches=None, pad_inches=0.1, frameon=None)
                plt.close(fig)  



            S_f[0:nx] = (n_manning**2.0) * ( (uh[0:nx]/h[0:nx]) * abs(uh[0:nx]/h[0:nx]) ) / (((B * h[0:nx]) /(B + 2.0 * h[0:nx]) )**(4.0/3.0) )
            S  [0:nx] = - g * h[0:nx] * ( S_0[0:nx] - S_f[0:nx] )

            #for i in range(nx):
            #    print("{:5d},  S_f: {:20.15f},  S: {:20.15f},  h: {:20.15f}".format(i,S_f[i], S[i], h[i]) )
            #print(" This is the end ----------------")

            hm[0:nx-1] = ( h[0:nx-1] + h[1:nx] ) / 2.0 - ( dt / 2.0 ) * ( uh[1:nx] - uh[0:nx-1] ) / dx
            uhm[0:nx-1] = ( uh[0:nx-1] + uh[1:nx] ) / 2.0 - ( dt / 2.0 ) * (  (uh[1:nx] ** 2.0) / h[1:nx]   + 0.5 * g * (h[1:nx] ** 2.0) - (uh[0:nx-1] ** 2.0 )/ h[0:nx-1] - 0.5 * g * (h[0:nx-1] ** 2.0 )) / dx - ( dt / 4.0 ) * ( S[1:nx] + S[0:nx-1] )
            #print(uhm)
            #print(hm)


            S_f_m[0:nx-1] = (n_manning**2.0) * ( (uhm[0:nx-1]/hm[0:nx-1]) * abs(uhm[0:nx-1]/hm[0:nx-1]) ) / (  ( (B * hm[0:nx-1]) /(B + 2 * hm[0:nx-1]) )**(4.0/3.0) )
            S_m  [0:nx-1] = - g * hm[0:nx-1] * ( S_0[0:nx-1] - S_f_m[0:nx-1] )

            h[1:nx-1] = h[1:nx-1] - dt * ( uhm[1:nx-1] - uhm[0:nx-2] ) / dx
            uh[1:nx-1] = uh[1:nx-1] - dt * ( (uhm[1:nx-1] ** 2.0)  / hm[1:nx-1] + 0.5 * g * (hm[1:nx-1] ** 2.0) - (uhm[0:nx-2] ** 2.0)  / hm[0:nx-2] - 0.5 * g * (hm[0:nx-2] ** 2.0) ) / dx - ( dt / 2.0 ) * ( S_m[1:nx-1] + S_m[0:nx-2] )


            h, uh = boundary_conditions ( nx, nt, h, uh, t[it], h_downstream, Q_upstream, B )

            #h_array[0:nx,it] = h[0:nx]
            #uh_array[0:nx,it] = uh[0:nx]


    print(" ========== Solver Class ==========")
    print()


def initial_conditions ( nx, nt, h, uh, x, z ):

#*****************************************************************************80
#
## INITIAL_CONDITIONS sets the initial conditions.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 December 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer NX, the number of spatial nodes.
#
#    Input, integer NT, the number of times steps.
#
#    Input, real H(NX,1), an array to hold the height.
#
#    Input, real UH(NX,1), an array to hold the mass velocity.
#
#    Input, real X(NX,1), the coordinates of the nodes.
#
#    Output, real H(NX,1), the initial height for all space.
#
#    Output, real UH(NX,1), the initial mass velocity for all space.
#
  import numpy as np

#  h[:] = 4.0 #+ np.sin ( 2.0 * np.pi * x )
  h[:] = 0.33 -z[:]
  uh = np.zeros ( nx )

  return h, uh


def boundary_conditions ( nx, nt, h, uh, t, h_downstream, Q_upstream, B ):

#*****************************************************************************80
#
## INITIAL_CONDITIONS sets the initial conditions.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    09 December 2016
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer NX, the number of spatial nodes.
#
#    Input, integer NT, the number of times steps.
#
#    Input, real H(NX), the height for all space.
#
#    Input, real UH(NX), the mass velocity for all space.
#
#    Input, real T, the current time.
#
#    Output, real H(NX), the height, with H(1) and H(NX) adjusted for
#    boundary conditions.
#
#    Output, real UH(NX), the mass velocity, with UH(1) and UH(NX)
#    adjusted for boundary conditions.
#
  bc = 4
#
#  Periodic boundary conditions on H and UH.
#
  if ( bc == 1 ):
    h[0] = h[nx-2]
    h[nx-1] = h[1]
    uh[0] = uh[nx-2]
    uh[nx-1] = uh[1]
#
#  Free boundary conditions on H and UH.
#
  elif ( bc == 2 ):
    h[0] = h[1]
    h[nx-1] = h[nx-2]
    uh[0] = uh[1]
    uh[nx-1] = uh[nx-2]
#
#  Reflective boundary conditions on UH, free boundary conditions on H.
#
  elif ( bc == 3 ):
    h[0] = h[1]
    h[nx-1] = h[nx-2]
    uh[0] = - uh[1]
    uh[nx-1] = - uh[nx-2]

  elif ( bc == 4 ):
    h[0] = h[1]
    h[nx-1] = h_downstream

    uh[0] = Q_upstream/B
    uh[nx-1] = uh[nx-2]

  return h, uh


