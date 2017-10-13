import math
from sympy import symbols, Eq, init_printing, sqrt
from IPython.display import display, Latex
init_printing(use_latex='mathjax')


class TimberBeam(object):
    """A class which contains the properties of a timber beam and methods to calculate its load bearing capacity according to Eurocode 5"""
    
        
    def __init__(self, ID, L, A, Iy, Iz, Ip, Wy, Wz, Wt, alfaY, alfaZ, alfaLTB):
        """Input units: L [mm], A [mm2], Iy[mm4], Iz[mm4], Ip[mm4], Wy [mm3], Wz [mm3], Wt [mm3], alfaY [-], alfaZ [-], alfaLTB [-]"""
        
        #Cross section properties
        self.ID = ID
        self.L = L
        self.Lsy = L*alfaY
        self.Lsz = L*alfaZ
        self.LLTB = L*alfaLTB
        self.A = A
        self.Iy = Iy
        self.Iz = Iz
        self.Ip = Ip
        self.Wy = Wy
        self.Wz = Wz
        self.Wt = Wt
        
        #Material properties (default C30)
        #Characteristic values
        self.E0k = 8000     # [MPa]
        self.G = 800        # [MPa]
        self.gamma = 460    # [kg/m3]
        self.fmk = 30       # [MPa]
        self.ft0k = 18      # [MPa]
        self.ft90k = 0.6    # [MPa]
        self.fc0k = 23      # [MPa]
        self.fc90k = 2.7    # [MPa]
        self.fvk = 3.0      # [MPa]
        
        self.kmod = 0.6
        self.partialCoeff = 1.35
        
        #Design values
        self.E0d = (self.kmod*self.E0k)/self.partialCoeff
        self.fmd = (self.kmod*self.fmk)/self.partialCoeff
        self.ft0d = (self.kmod*self.ft0k)/self.partialCoeff
        self.ft90d = (self.kmod*self.ft90k)/self.partialCoeff
        self.fc0d = (self.kmod*self.fc0k)/self.partialCoeff
        self.fc90d = (self.kmod*self.fc90k)/self.partialCoeff
        self.fvd = (self.kmod*self.fvk)/self.partialCoeff
        
        #Other
        self.isGlulam = False
        self.isRect = True
        
        
    def __str__(self):
        return "TimberBeam: {0}".format(self.ID)
        
        
        
    # ------------------------------------------------ MATERIAL PROPERTIES --------------------------------------------------
    def setMatProps(self, E0k, G, gamma, fmk, ft0k, ft90k, fc0k, fc90k, fvk, kmod, partialCoeff, isGlulam, isRect):
        """Define the material properties of the timber beam.
    
        Input units: E0k [MPa], G [MPa], gamma [kg/m3], fmk/ft0k/ft90k/fc0k/fc90k/fvk [MPa], kmod [-], partialCoeff [-], isGlulam [T/F], isRect [T/F]"""
        
        self.E0k = E0k
        self.G = G
        self.gamma = gamma
        self.fmk = fmk
        self.ft0k = ft0k
        self.ft90k = ft90k
        self.fc0k = fc0k
        self.fc90k = fc90k
        self.fvk = fvk
        
        self.E0d = (kmod*E0k)/partialCoeff
        self.fmd = (kmod*fmk)/partialCoeff
        self.ft0d = (kmod*ft0k)/partialCoeff
        self.ft90d = (kmod*ft90k)/partialCoeff
        self.fc0d = (kmod*fc0k)/partialCoeff
        self.fc90d = (kmod*fc90k)/partialCoeff
        self.fvd = (kmod*fvk)/partialCoeff 
        
        self.isGlulam = isGlulam
        self.isRect = isRect   
        
        
        
    # -------------------------------------------------- BEAM/COLUMN AXIAL ---------------------------------------------------
    def calcUtil_Axial(self, NormalForce, printParam, printCalc, outputResult):
    
        """Calculates the utilisation of a timber beam/column in tension/compression without column effect.
    
        Input units: N [N], printParam [T/F], printCalc [T/F], outputResult [T/F]"""
        
        
        #Containers and variables
        doc_list = []
        equation_list = []
        result_list = []
        unit_list = []
        
        if(NormalForce >= 0.0):
            suffix = "t"
            _fd = self.ft0d
        else:
            suffix = "c"
            _fd = self.fc0d
            
            
        #Print title and parameters
        if(printCalc):
            display(Latex("**AXIALLY LOADED BEAM/COLUMN**"))
            
            if(printParam):
                print("Input parameters:")
                print("N = {0} {1}".format(NormalForce, "N"))
                print("A = {0} {1}".format(self.A, "mm2"))
                print("f{0}0d = {1} {2}".format(suffix, _fd, "MPa"))
                  
    
        #Calculate utilisation
        u, N, A, fd = symbols('u_{0} N A f_{0}0d'.format(suffix))
        data_u = [(N,math.fabs(NormalForce)), (A,self.A), (fd,_fd)]
        eq_u = Eq( u, (N/A)/fd )
        result_u = eq_u.subs(data_u).rhs.evalf()
        unit_u = "-"
    
        doc_list.append("\nCalculate axial utilisation: ")
        equation_list.append(eq_u)
        result_list.append(result_u)
        unit_list.append(unit_u)
        

        #Print output
        if(printCalc):           
            self.printOutput(doc_list, equation_list, result_list, unit_list)
    
        if(outputResult):  
            return result_list[-1]
           
         
    
    # -------------------------------------- COLUMN IN COMPRESSION WITH COLUMN EFFECT ----------------------------------------
    def calcUtil_Axial_ColumnEffect(self, NormalForce, axis, printParam, printCalc, outputResult):
        """Calculates the utilisation of a timber column in compression with column effect.
    
        Input units: N [N], axis [0:y, 1:z], printParam [T/F], printCalc [T/F], outputResult [T/F]"""
        
        
        #Containers and variables
        doc_list = []
        equation_list = []
        result_list = []
        unit_list = []
        
        if(axis == 0):
            suffix = 'y'
            _ls = self.Lsy
            _I = self.Iy              
        else:
            suffix = 'z'
            _ls = self.Lsz
            _I = self.Iz 
            
        
        #Print title and parameters
        if(printCalc):
            display(Latex("**COLUMN IN COMPRESSION WITH COLUMN EFFECT ({0}-AXIS)**".format(suffix.upper())))
            
            if(printParam):
                print("Input parameters:")
                print("N = {0} {1}".format(NormalForce, "N"))
                print("E0d = {0} {1}".format(self.E0d, "MPa"))
                print("I{2} = {0} {1}".format(_I, "mm4", suffix))
                print("Ls{2} = {0} {1}".format(_ls, "mm", suffix))
                print("A = {0} {1}".format(self.A, "mm2"))
                print("fc0d = {0} {1}".format(self.fc0d, "MPa"))
                
    
        #Calculate relative slenderness
        lambda_rel, ls, I, A, fcd, pi, Ed  = symbols('lambda_rel{0} l_s{0} I_{0} A f_0cd pi E_0d'.format(suffix))
        eq_lambda_rel = Eq( lambda_rel, (ls/(sqrt(I/A))) * sqrt(fcd/(pi**2 * Ed)) )
    
        data_lambda_rel = [(Ed,self.E0d), (I,_I), (ls,_ls), (A,self.A), (pi,math.pi), (fcd, self.fc0d)]
        result_lambda_rel = eq_lambda_rel.subs(data_lambda_rel).rhs.evalf()
        unit_lambda_rel = "-"
    
        doc_list.append("\nCalculate relative slenderness: ")
        equation_list.append(eq_lambda_rel)
        result_list.append(result_lambda_rel)
        unit_list.append(unit_lambda_rel)
    
    
        #Calculate kc factor based on the relative slenderness
        kc, k, beta  = symbols('k_c{0} k beta'.format(suffix))
    
        if(result_lambda_rel <= 0.3):
            eq_kc = Eq( kc, 1.0 )
            result_kc = eq_kc.rhs
            unit_kc = "-"
        
            doc_list.append("\nCalculate strength reduction factor kc based on the relative slenderness: ")
            equation_list.append(eq_kc)
            result_list.append(result_kc)
            unit_list.append(unit_kc)
        
        else:
            eq_k = Eq( k, 0.5 * (1 + beta * (lambda_rel - 0.3) + lambda_rel**2) )
    
            _beta = 0.2
            if(self.isGlulam):
                _beta = 0.1
        
            data_k = [(lambda_rel, result_lambda_rel), (beta,_beta)]
            result_k = eq_k.subs(data_k).rhs.evalf()
            unit_k = "-"
    
            eq_kc = Eq( kc, 1 / (k + sqrt(k**2 - lambda_rel**2)) )
            data_kc = [(k, result_k), (lambda_rel,result_lambda_rel)]
            result_kc = eq_kc.subs(data_kc).rhs.evalf()
            unit_kc = "-"
        
            doc_list.append("\nCalculate strength reduction factor kc based on the relative slenderness (beta = {0}): ".format(_beta))
            equation_list.append(eq_k)
            result_list.append(result_k)
            unit_list.append(unit_k)
        
            doc_list.append(" ")
            equation_list.append(eq_kc)
            result_list.append(result_kc)
            unit_list.append(unit_kc)
    
    
        #Calculate utilisation
        uc, N  = symbols('u_c{0} N'.format(suffix))
        eq_uc = Eq( uc, (N/A) / (kc * fcd) )
        data_uc = [(N, math.fabs(NormalForce)), (A,self.A), (kc,result_kc), (fcd, self.fc0d)]
        result_uc = eq_uc.subs(data_uc).rhs.evalf()
        unit_uc = "-"
    
        doc_list.append("\nCalculate utilisation: ")
        equation_list.append(eq_uc)
        result_list.append(result_uc)
        unit_list.append(unit_uc)
    
    
        #Print output
        if(printCalc):        
            self.printOutput(doc_list, equation_list, result_list, unit_list)
    
        if(outputResult):  
            return result_list[-1]
        
    
    
    # ---------------------------------------------- BEAM/COLUMN BENDING -----------------------------------------------------
    def calcUtil_Bending(self, My, Mz, axis, printParam, printCalc, outputResult):
    
        """Calculates the utilisation of a timber beam/column in bending.
    
        Input units: My [Nmm], Mz [Nmm], axis [0:y, 1:z], printParam [T/F], printCalc [T/F], outputResult [T/F]"""
    
        
        #Containers and variables
        doc_list = []
        equation_list = []
        result_list = []
        unit_list = []
    
        km = 1.0
        if(self.isRect):
            km = 0.7
            
        suffix = "y"
        if(axis == 1):
            suffix = "z"
            
    
        #Print title and parameters
        if(printCalc):
            display(Latex("**BEAM/COLUMN IN BENDING ({0}-AXIS AS PRIMARY)**".format(suffix.upper())))
        
            if(printParam):
                print("Input parameters:")
                print("My = {0} {1}".format(My, "Nmm"))
                print("Mz = {0} {1}".format(Mz, "Nmm"))
                print("Wy = {0} {1}".format(self.Wy, "mm3"))
                print("Wz = {0} {1}".format(self.Wz, "mm3"))
                print("fmd = {0} {1}".format(self.fmd, "MPa"))
                

        #Define general symbols and specify data
        u_my, u_mz, M_y, W_y, M_z, W_z, f_md, k_m  = symbols('u_myz u_mzy M_y W_y M_z W_z f_md k_m')
        data_um = [(M_y, math.fabs(My)), (M_z, math.fabs(Mz)), (W_y,self.Wy), (W_z,self.Wz), (f_md,self.fmd), (k_m,km)]
    
        #Calculate utilisation about y-axis
        if(axis==0): 
            eq_umy = Eq( u_my, (M_y/W_y)/f_md + k_m * (M_z/W_z)/f_md )
    
            result_umy = eq_umy.subs(data_um).rhs.evalf()
            unit_umy = "-"
    
            doc_list.append("\nCalculate utilisation (km = {0}): ".format(km))
            equation_list.append(eq_umy)
            result_list.append(result_umy)
            unit_list.append(unit_umy)
        
        #Calculate utilisation about z-axis
        else:
            eq_umz = Eq( u_mz, k_m * (M_y/W_y)/f_md + (M_z/W_z)/f_md )
    
            result_umz = eq_umz.subs(data_um).rhs.evalf()
            unit_umz = "-"
    
            doc_list.append("\nCalculate utilisation (km = {0}): ".format(km))
            equation_list.append(eq_umz)
            result_list.append(result_umz)
            unit_list.append(unit_umz)
        
    
        #Print output
        if(printCalc):                 
            self.printOutput(doc_list, equation_list, result_list, unit_list)
    
        if(outputResult):  
            return result_list[-1]
     
    
    
    # ---------------------------------- BEAM/COLUMN COMBINED AXIAL FORCE AND BENDING ---------------------------------------
    def calcUtil_CombinedAxialBending(self, NormalForce, My, Mz, axis, columnEffect, printParam, printCalc, outputResult):
    
        """Calculates the utilisation of a timber beam/column under combined axial force and bending with/without column effect.
    
        Input units: N [N], My [Nmm], Mz [Nmm], axis [0:y, 1:z], columnEffect [T/F], printParam [T/F], printCalc [T/F], outputResult [T/F]"""
        
        
        #Containers and variables
        doc_list = []
        equation_list = []
        result_list = []
        unit_list = []
    
        if(NormalForce >= 0.0):
            suffix0 = "t"
            _fd = self.ft0d
        else:
            suffix0 = "c"
            _fd = self.fc0d
        
        if(columnEffect and axis==0):
            suffix1 = "y"
        elif(columnEffect and axis==1):
            suffix1 = "z"
        else:
            suffix1 = ""
               
        if(axis==0):
            suffix2 = "yz"
        else:
            suffix2 = "zy"
            
        
        #Print title and parameters
        if(printCalc):
            display(Latex("**BEAM/COLUMN UNDER COMBINED AXIAL FORCE AND BENDING**"))
            
            if (printParam):
                print("Input parameters:")
                print("N = {0} {1}".format(NormalForce, "N"))
                print("A = {0} {1}".format(self.A, "mm2"))
                print("f{0}0d = {1} {2}".format(suffix0, _fd, "MPa"))
                print("My = {0} {1}".format(My, "Nmm"))
                print("Mz = {0} {1}".format(Mz, "Nmm"))
                print("Wy = {0} {1}".format(self.Wy, "mm3"))
                print("Wz = {0} {1}".format(self.Wz, "mm3"))
                print("fmd = {0} {1}".format(self.fmd, "MPa"))
                
        
        #Calculate individual utilisation results
        #Axial
        if(columnEffect and NormalForce < 0.0):
            util_axial = self.calcUtil_Axial_ColumnEffect(NormalForce, axis, False, printCalc, True)
        else:
            util_axial = self.calcUtil_Axial(NormalForce, False, printCalc, True)
            
        #Bending
        util_bending = self.calcUtil_Bending(My, Mz, axis, False, printCalc, True)
        
        
        #Calculate utilisation
        u_tot, u_axial, u_bending = symbols('u_tot{0} u_{1} u_m{2}'.format(suffix1, suffix0 + suffix1, suffix2))
        
        data_u_tot = [(u_axial, util_axial), (u_bending, util_bending)]
        eq_u_tot = Eq( u_tot, u_axial + u_bending )
        result_u_tot = eq_u_tot.subs(data_u_tot).rhs.evalf()
        unit_u_tot = "-"
        
        doc_list.append("\nCalculate combined utilisation: ")
        equation_list.append(eq_u_tot)
        result_list.append(result_u_tot)
        unit_list.append(unit_u_tot)
        

        #Print output
        if(printCalc):                
            self.printOutput(doc_list, equation_list, result_list, unit_list)
    
        if(outputResult):  
            return result_list[-1]



    # ------------------------------------------ PRESTRESS FROM INITIAL BENDING ---------------------------------------------
    def calcUtil_BendingPrestress(self, curvatureRadius, printParam, printCalc, outputResult):
    
        """Calculates the utilisation of a timber beam/column that is initially bent into its position.
    
        Input units: curvatureRadius [m], printParam [T/F], printCalc [T/F], outputResult [T/F]"""
        
        
        #Containers and variables
        doc_list = []
        equation_list = []
        result_list = []
        unit_list = []
        

        #Calculate thickness
        _t = (self.Iy / self.Wy) * 2
            
        #Print title and parameters
        if(printCalc):
            display(Latex("**INITIALLY BENT BEAM/COLUMN**"))
            
            if(printParam):
                print("Input parameters:")
                print("r = {0} {1}".format(curvatureRadius, "m"))
                print("t = {0} {1}".format(_t, "mm"))
                print("E0k = {0} {1}".format(self.E0k, "MPa"))
                print("fmd = {0} {1}".format(self.fmd, "MPa"))
                  
    
        #Calculate utilisation
        u_b, r, t, E, fmd = symbols('u_bent r t E_0k f_md')
        data_u_b = [(r,curvatureRadius*1e3), (t,_t), (E, self.E0k), (fmd,self.fmd)]
        eq_u_b = Eq( u_b, ((E*t)/(2*r))/fmd )
        result_u_b = eq_u_b.subs(data_u_b).rhs.evalf()
        unit_u_b = "-"
    
        doc_list.append("\nCalculate utilisation from initial bending: ")
        equation_list.append(eq_u_b)
        result_list.append(result_u_b)
        unit_list.append(unit_u_b)
        

        #Print output
        if(printCalc):           
            self.printOutput(doc_list, equation_list, result_list, unit_list)
    
        if(outputResult):  
            return result_list[-1]
           
        
        
    # -------------------------------------------------- PRINT FUNCTION -----------------------------------------------------
    def printOutput(self, doc_list, equation_list, result_list, unit_list):
        for i in range(0,len(equation_list)):
        
            print(doc_list[i])
        
            display(equation_list[i])
        
            unit = unit_list[i]
            if(unit == "-"):
                unit = ""
            
            outputText = "= {0:.2f} {1}".format(result_list[i], unit)
            print("{0:>100}".format(outputText))
    
    