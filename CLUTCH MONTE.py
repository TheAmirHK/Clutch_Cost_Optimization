# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 09:57:40 2021

@author: amirh
"""

import numpy as np
import math
from pyomo.environ import * 
from pyomo.opt import SolverFactory
from scipy.stats import norm
from termcolor import colored
import matplotlib.pyplot as plt
import time
import timeit
from scipy import stats
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from smt.sampling_methods import LHS
import sympy as sp

# In[] "-----Dynamic Plotting Prepration -----"
xdata = []
ydata = []
plt.close()
plt.show() 
axes = plt.gca()
axes.set_xlim(0.005 , 1.005)
axes.set_ylim(1e1, 3000)
line, = axes.plot(xdata, ydata, 'r-')

# In[] 
##########################   Partial Derivation Calculator   #####################################################

number_of_varibles = 3
derivation_coefficient = list()
nominal = [55.2973, 22.86, 101.6] # Dimenssions' Nominal Values
x = [sp.symbols('x%d' % i) for i in range(number_of_varibles)] # Number of dimenssions
f = sp.acos((x[0]+x[1])/(x[2]-x[1])) # Functional requirement 

for i in x:        
    Derivation_coefficient_equation = sp.diff(f,i)  #  Derivation coefficient
    Derivation_coefficient_value = Derivation_coefficient_equation.subs((x[j], nominal[j]) for j in range (number_of_varibles)) 
    derivation_coefficient.append(sp.Abs(Derivation_coefficient_value))
derivation_coefficient = [float(i) for i in derivation_coefficient ]    

# In[]  "-----Monte Carlo -----"
lty = np.linspace (0.1 , 1 , 10) # Quality requirement in degree 
nmct = len(lty)
Result_list = list ()
for iter in range(nmct):
    start_time = time.time()
            
    ##########################   Optimisation Algorithm   #####################################################
    model = ConcreteModel(name="(Allocation)")
    
    #Rework = int(input("Would you do rework (Zero or One): "))
    
    #model.partitions = int(input("Number of partitions for storage: "))
    
    Rework = 1
    model.partitions = 3
    
    ty = lty[iter]
    #Scalars
    model.it = 1000 # iteration
    model.NParts = 3
    model.op = 3
    #model.partitions = 2
    NRR=np.zeros(model.NParts)
    al=0.0027 #Error type I
    be=0.00005 #Error type II
    Q = 100     # Number of demands 
    
    #Sets
    model.m = [0,1,2]
    model.n = [0,1,2]
    size = (model.op,) * model.NParts
    
    N = np.zeros((model.NParts,model.op,model.it))
    NR = np.zeros((model.NParts,model.op,model.it))
    NRR2 = np.zeros((model.NParts,model.op,model.it))
    
    #Model costs
#    Components_Manufacturing_cost= ((5, 8, 10), (3, 2.5, 2.95), (2.95, 3.15, 4))
    s = ((0.566, 0.133, 0.100),(0.166, 0.300, 0.208),(0.208, 0.133, 0.09))
#    s = ((0.0566, 0.0133, 0.0100),(0.0166, 0.0300, 0.0208),(0.0208, 0.0133, 0.009))

    Components_Manufacturing_cost= ((2, 3.15, 3.5), (3, 2.5, 2.95), (2.95, 3.15, 4)) 
    Components_Scraping_cost = (0.5,0.5,0.5) 
    Components_Inspection_cost = (1,1.5,1)
    Components_Reworking_cost = (1,1,1)
    Components_Inventory_cost = (2,2,2)
    
    Product_Scraping_cost = 10
    Product_Inspection_cost = 0.5
    Product_Assembly_cost = 3
    
    #s = ((0.166, 0.0633, 0.0600),(0.0666, 0.0800, 0.0708),(0.0708, 0.0633, 0.059))
    
    
    # Model initialization
    A = [(i,j,k) for i in range(model.NParts) for j in range(model.op) for k in range(model.it)]
    B = [(i,j) for i in range(model.NParts) for j in range(model.op)]
    C = [(i) for i in range(model.NParts)]
    D = [(i) for i in range(model.op**model.NParts)]
    
    
    g = np.zeros((model.NParts,model.op,model.it))
    cp = np.zeros((model.NParts,model.op,model.it))
    gp = np.zeros((model.NParts,model.op,model.it))
    grework = np.zeros((model.NParts,model.op,model.it))
    gprework = np.zeros((model.NParts,model.op,model.it))
    ps = np.zeros((model.NParts,model.op,model.it))
    z = np.zeros((model.NParts,model.op,model.it))
    AU = np.zeros((model.NParts,model.op,model.it))
    r = np.zeros((model.NParts,model.op,model.it))
    gr = np.zeros((model.NParts,model.op,model.it,model.partitions))
    
    Lower_bound = 0.01
    Upper_bound = 2
    u= (Lower_bound + Upper_bound)/2
    sigma = ((Lower_bound-Upper_bound)**2)/12
    t = np.random.uniform(Lower_bound, Upper_bound, (model.NParts, model.op, model.it))
    #t = np.random.normal(u, sigma, (model.NParts, model.op, model.it))
    
    np.random.seed(0)
    
    t = np.sort(t)[ : : +1]
    l = np.zeros((size))
    model.laa = 0
    "-----Define components ratio pass and regarding prcess capability -----"
    
    for i in range (model.NParts):
        for j in range(model.op):
            for k in range (model.it):
                g[i,j,k]  = math.erf(3*(t[i,j,k])/(3*np.array(s)[i,j]*math.sqrt(2)))  # Confirmity Ratio   
                grework[i,j,k] = Rework*(1 - g[i,j,k])/2 # Reworking Ratio
                ps[i,j,k] = norm.cdf((-t[i,j,k])/(3*np.array(s)[i,j])) # Scrap Ratio
                gprework[i,j,k] = g[i,j,k] + grework[i,j,k]*g[i,j,k] # New Confirmity Ratio after reworking
                gp[i,j,k] = gprework[i,j,k]*(1-al)+(1-gprework[i,j,k])*be # Confirmity ratio considering failures
                cp[i,j,k] = t[i,j,k]/(3*np.array(s)[i,j]) # process capability
                N[i,j,k]  = (Q/gp[i,j,k]) # Number of components to be processed without residual
                for o in range (model.partitions):
                    gr[i,j,k,o] = norm.cdf((6*cp[i,j,k]*(o-(model.partitions)/2+1))/(model.partitions)) - norm.cdf((6*cp[i,j,k]*(o-(model.partitions)/2))/(model.partitions))
                z[i,j,k] = sum(np.sqrt(gr[i,j,k,l]*(1-gr[i,j,k,l])) for l in range(model.partitions))         
                AU[i,j,k]= np.sqrt(1/(2*np.pi*Q))*z[i,j,k]
                r[i,j,k]=AU[i,j,k]*np.sqrt((AU[i,j,k]**2)+2)-AU[i,j,k]**2   
                NR[i,j,k] = int(round(Q/((gp[i,j,k])*(1-r[i,j,k]))))
                NRR2[i,j,k] = int(round(Q*(r[i,j,k])/((gp[i,j,k])*(1-r[i,j,k]))))                
                       
    
    "-----Define assembled product conformity regarding operations devation -----"             
    for i in range (model.op):
        for j in range (model.op):
            for k in range (model.op): 
                l[i,j,k] = math.erf(ty/(math.sqrt((derivation_coefficient[0]*np.array(s)[0,i])**2+(derivation_coefficient[1]*np.array(s)[1,j])**2+(derivation_coefficient[2]*np.array(s)[2,k])**2)*math.sqrt(2)))
    l = l.flatten('F')
    
    "-----Define Variables -----"            
    model.x = Var(A, within = Binary)
    model.y = Var(C, within = NonNegativeReals)
    model.indexx = Var(initialize=0)
    model.la = Var(initialize=0)
    model.la2 = Var(initialize=0)
    model.Assembly_Conformity = Var(initialize=0)
    model.count = Var(D, within = Binary)
    model.NRR = Var(C, within = NonNegativeReals)
    ##
    "-----Define Objective Functions -----"   
    def objective_rule(model):
        CM = sum((NR[i,j,k]*np.array(Components_Manufacturing_cost)[i,j] /gp[i,j,k]* model.x[i,j,k]) for i in range(model.NParts) for j in range(model.op) for k in range (model.it))
        CS = sum ((NR[i,j,k]*np.array(Components_Scraping_cost)[i]*(ps[i,j,k])/(gp[i,j,k])* model.x[i,j,k] ) for i in range(model.NParts) for j in range(model.op) for k in range (model.it)) 
        CI = sum((NR[i,j,k]*np.array(Components_Inspection_cost)[i] /(gp[i,j,k])* model.x[i,j,k]) for i in range(model.NParts) for j in range(model.op) for k in range (model.it))
        CR = sum((NR[i,j,k]*grework[i,j,k]*Components_Reworking_cost[i]* model.x[i,j,k]/(gp[i,j,k])) for i in range(model.NParts) for j in range(model.op) for k in range (model.it))
        CIP = Product_Inspection_cost *Q *(model.la)
        CCI = sum(Components_Inventory_cost[i]*model.NRR[i] for i in range (model.NParts))    
        CA = (Product_Assembly_cost*Q *(model.la))
        CSP = Product_Scraping_cost*Q*(1-model.la2)
    
#            expr = CM + CS + CI + CR + CIP + CCI + CA + CSP
        expr = CM + CS + CI + CIP + CA + CSP + Rework*CR
        return  expr 
    model.objective = Objective(rule=objective_rule, sense=minimize)
    #
    #
    "---------Define constriants---------"
    def Each_component_one_operation(model, i):
        return (sum(model.x[i,j,k] for j in range(model.op) for k in range (model.it))==1)
    model.constraint = Constraint(model.m, rule=Each_component_one_operation)
    
    def Take_the_operation_index_for_each_component(model, i):
        return (model.y[i] == sum(model.x[i,j,k]*(j) for j in range (model.op) for k in range(model.it)))
    model.constraint2 = Constraint(model.m, rule=Take_the_operation_index_for_each_component)
    
    def Calculate_the_index_for_vectorized_assembly_conformity(model):
        return (model.indexx == sum(model.op*(i+1)*model.y[i+1] for i in range (model.NParts-1)) + model.y[0])
    model.constraint3 = Constraint(rule=Calculate_the_index_for_vectorized_assembly_conformity)
    
    def Take_the_proper_position_of_the_vectorized_assembly_conformity(model):
        return (sum(i*model.count[i] for i in range (model.op**model.NParts))  == model.indexx)
    model.constraint4 = Constraint(rule=Take_the_proper_position_of_the_vectorized_assembly_conformity)
    
    def Just_one_position_of_the_vectorized_assembly_conformity_should_be_selected(model):
        return (sum(model.count[i] for i in range (model.op**model.NParts))  == 1  )
    model.constraint5 = Constraint(rule=Just_one_position_of_the_vectorized_assembly_conformity_should_be_selected)
    
    def Calculated_Assembly_conformity_with_failure(model):
        return (model.la  == sum (model.count[i] /(l[i]*(1-al)+(1-l[i])*be)  for i in range (model.op**model.NParts)))
    model.constraint6 = Constraint( rule=Calculated_Assembly_conformity_with_failure)
    
    def Calculated_Assembly_conformity_with_failure2(model):
        return (model.la2  == sum (model.count[i] * (l[i]*(1-al)+(1-l[i])*be)  for i in range (model.op**model.NParts)))
    model.constraint9 = Constraint( rule=Calculated_Assembly_conformity_with_failure2)
    
    #def constraint_rule49(model, i):
    #    return (model.NRR[i]  == sum ((NR[i,j,k]*model.x[i,j,k] - N[i,j,k]*model.x[i,j,k])  for j in range(model.op) for k in range (model.it)))
    #model.constraint7 = Constraint(model.m,  rule=constraint_rule49)
    
    def constraint_rule49(model, i):
        return (model.NRR[i]  == sum (NRR2[i,j,k]*model.x[i,j,k]  for j in range(model.op) for k in range (model.it)))
    model.constraint7 = Constraint(model.m,  rule=constraint_rule49)
    
    "-----process capability -----"
#        def proces_capability(model, i ):
#            return (sum(cp[i,j,k]*model.x[i,j,k] for j in range (model.op) for k in range (model.it))<=2)
#        model.constraint10 = Constraint(model.m, rule=proces_capability)
    #    
    
    
    "-----Functional Requirement -----"
    def Functional_requirement(model ):
        return (sum(((derivation_coefficient[i]*t[i,j,k])**2)*model.x[i,j,k] for i in range (model.NParts) for j in range (model.op) for k in range (model.it))<=(ty)**2)
    model.constraint8 = Constraint( rule=Functional_requirement)
    
    "-----Solve-----"
    opt = SolverFactory('cplex')
    results = opt.solve(model, tee=False)
#    Cost = np.zeros((2,2))
    Resources = np.zeros(3)
    Component_conformity = np.zeros(3)
    Component_tolerance = np.zeros(3)
#    Number_of_Scraps = np.zeros(3)
#    Assembly_comformity = np.zeros((2,2))
    print (colored('Functional requirement (mm)=','red'),ty)
    
#    a= np.zeros(model.NParts)
    for i in range (model.NParts):
        for j in range(model.op):
            for k in range (model.it):
                if value(model.x[i,j,k])==1:
                    Resources[i] = j+1 
                    Component_conformity[i] = "%.2f"  % value(gp[i,j,k]*100)
                    Component_tolerance[i] = "%.3f"  % t[i,j,k]
#                    Number_of_Scraps[i] = int(round(N[i,j,k])-Q)
                    print(colored('component','red'), i+1, colored('will be processed on','red'), colored(j+1,'green'))
                    print(colored('component','red'), i+1, colored('allocated tolerance T=','red'),"%.3f"  % t[i,j,k])
                    print(colored('component','red'), i+1, colored('allocated process capability CP=','red'),"%.3f"  % cp[i,j,k])
                    print(colored('component','red'), i+1, colored('conforimity ratio G= %','red'), "%.2f"  % value(gp[i,j,k]*100)) 
                    print(colored('component','red'), i+1, colored('number of components to be manufactured NR=','red'),int(round(NR[i,j,k])))
                    print(colored('component','red'), i+1, colored('number of scraps =','red'), int(round(N[i,j,k])-Q))                                                   
                    print(colored('component','red'), i+1, colored('number of residuals =','red'),int(round(NR[i,j,k])) - int(round(N[i,j,k])))
                    print(colored('component','red'), i+1, colored('number of reworkeds =','red'), int(int(round(NR[i,j,k]))*grework[i,j,k]))               

    #print(value(model.NRR[0]))
#    Assembly_comformity = [ty,value(model.la2 *100)]
    Cost = [ty,value(model.objective )]               
    print(colored('Production cost =','blue'),"%.1f" % value(model.objective ))
    print(colored('Assembly conformity rate = %','blue'), "%.2f"  % value(model.la2 *100))
     
    print("--- RunTime = %s seconds ---" % int(round((time.time() - start_time))))
    Result_list.extend(['ty', ty ,Resources,Component_tolerance, Component_conformity, "%.2f"  %  value(model.la2 *100),"%.1f" % value(model.objective )])
#    Cost_list.extend((Cost,Resources))
    
# In[] "-----Dynamic Plotting  -----"
    xdata.append(lty[iter])
    ydata.append(value(model.objective))
    line.set_xdata(xdata) # set x cordinates
    line.set_ydata(ydata) # set y cordinates
    plt.xlabel("Tolerance")
    plt.ylabel("Production Cost")
    plt.draw()
    plt.pause(1e-2)
Tolerances_list = np.zeros((nmct,model.NParts))
Resources_list = np.zeros((nmct,model.NParts))
Conformity_list = np.zeros((nmct,model.NParts))
Manufacturing_cost_list = [Result_list[7*(i+1)-1] for i in range(nmct)]
Assembly_conformity_list = [Result_list[7*(i+1)-2] for i in range(nmct)]
for i in range(nmct):
    for j in range(model.NParts):
        Resources_list[i][j] = Result_list[7*(i+1)-5][j]
        Tolerances_list[i][j] = Result_list[7*(i+1)-4][j]
        Conformity_list[i][j] = Result_list[7*(i+1)-3][j]
