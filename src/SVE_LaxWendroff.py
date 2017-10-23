
##############################################################################
#
# Code developed by: Dr. Babak Poursartip
# Supervised by:     Dr. Ben R. Hodges
# 
# Start date:    08/18/2017
# Latest update: 10/10/2017
#
# Comment: This code is developed based on the following open source code:
#
#  Reference:
#
#    Cleve Moler, "The Shallow Water Equations", Experiments with MATLAB.
#    Developed by: John Burkardt
#
#    Machalinska et al. Lax-Wendroff and McCormack Scheme for Numerical 
#    Simulation of Unsteady Gradually and Rapidly Varied Open Channel Flow
#
##############################################################################

def main(arg):


    # Import built-in libraries =======================================================================
    # import numpy as np --delete

    # Import classes ==================================================================================
    import sys
    import Solver_Class
    import math  # 
    import os    #  You can create/del dir using this module
    import time  # time.sleep(2)
    from datetime import datetime 

    # Code begins =====================================================================================
    print()
    print("{:^80}".format("-------- 1D Saint-Venant Equation based on Lax-Wendroff scheme -------"))
    print("{:^80}".format("---------- Developers:: Babak Poursartip/Ben R. Hodges ----------"))
    print()
    print("{:^80}".format(" Simulation starts ..."))
    print()

    Results = Solver_Class.Solver()

    print("{:80}".format("---------- Simulation was conducted successfully ----------"))
    print()

if __name__ == '__main__':
    import sys    
    main(sys.argv)





