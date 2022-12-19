


import pyomo
import pyomo.environ as pe
import itertools 
import numpy as np
import scipy.linalg as spL
from scipy.optimize import lsq_linear


import pyutilib.services
from pyomo.opt import TerminationCondition







def PCO_opt(obj_coefs,vcount_upper_bound,f_log_counts,S,K):

    '''  
    Computes a maximum flux solution with respect to obj_coefs

     n_react:   Number of reactions
     n_M_v:     Number of variable metabolites
     n_M_f:     Number of fixed metabolites



    INPUTS
    obj_coefs:              array of shape (n_react,) with objective coefficients for each reaction 
    vcount_upper_bound:     array of shape (n_M_v,) with upper bound on the log count for each variable metabolite
    f_log_counts:           array of shape (n_M_f,) with log count for each fixed metabolite
    S:                      array of shape (n_react, n_M_v+n_M_f) the stoichiometric matrix
    K:                      array of shape (n_react,) with the equilibrium constant for each reaction
    
    
    OUTPUTS
    y_sol:                  array of shape (n_react,) with steady state flux for each reaction at the solution
    alpha_sol:              array of shape (n_react,) with activity coefficient for each reaction at the solution
    n_sol:                  array of shape (n_M_v,) with the log count of each variable metabolite at the solution

    '''







   
    '''
    ##################################
    Set Optimization Parameters

    These parameters have been chosen heuristically to achieve good convergence for a specific implementation 
    and may need to be adjusted to achieve convergence in other implementations
    ##################################
    '''
    
    VarM_lbnd = -150  # lower bound on the log count of the metabolites
    Mb = 50          # Big M coefficient to relax flux sign constraints
    lmbda_a = 1e-50   # parameter for the relaxed signum function
    lmbda_b = 1e50    # parameter for the relaxed signum function
    


    #Solver Parameters

    max_iter = 10000  # maximum interations for the solver
    max_cpu_time = 800000

    tolerance = 1e-7
    acc_tolerance = 1e-8    




   
    '''
    ##################################
    
    Pre-Process Inputs for Optimization

    ##################################
    '''




    # Flip Stoichiometric Matrix
    S_T = S    # S_T is the Stoich matrix with rows as reactions, columns as metabolites
    S = np.transpose(S_T) # this now is the Stoich matrix with rows metabolites, and columns reactions
    n_react = np.shape(S)[1]

    #Set the System parameters 
    VarM = vcount_upper_bound
    FxdM = np.reshape(f_log_counts,(len(f_log_counts),1) )
    K_eq = np.reshape(K,(n_react,1))
    n_M = len(VarM) + len(FxdM)

    #Metabolite parameters
    n_M_f = len(f_log_counts)   #total number of fixed metabolites
    n_M_v = len(vcount_upper_bound)   #total number of variable metabolites
    n_M = n_M_f + n_M_v   # total number of metabolites

    ### construct parameter indices
    react_idx = np.arange(0,n_react)

    TotM_idx = np.arange(0,n_M)
    VarM_idx = np.arange(0,n_M_v) #variable metabolite indices
    FxdM_idx = np.arange(n_M_v,n_M) # fixed metabolite indices

    # Split S into the component corresponding to the variable metabolites S_v 
    # and the component corresponding to the fixed metabolites S_f
    S_v_T = np.delete(S_T,FxdM_idx,1)
    S_v = np.transpose(S_v_T)

    S_f_T = np.delete(S_T,VarM_idx,1)
    S_f = np.transpose(S_f_T)

    #find a basis for the nullspace of S_v
    Sv_N = spL.null_space(S_v)
    dSv_N = np.shape(Sv_N)[1] # the dimension of the nullspace


    beta_idx = np.arange(0,dSv_N)


    #Find the gradient direction
    beta_grad = np.matmul(obj_coefs,Sv_N)
    beta_grad = np.transpose(beta_grad)
    y_grad = np.matmul(Sv_N,beta_grad)
    y_grad = np.reshape(y_grad,(n_react,1))

    

    '''
    ##################################
     Compute the initial condition
    ##################################
    '''

    K_v = np.reshape(K,(len(K),1))

    #### adjust signs to match - y_grad
    S_v_T_sgn = S_v_T*np.sign(-y_grad)
    S_v_sgn = np.transpose(S_v_T_sgn)
    S_f_T_sgn = S_f_T*np.sign(-y_grad)



    ### append the identity for the slack variables
    S_id = np.eye(n_react)
    A = np.concatenate([S_v_T_sgn,-S_id],axis = 1)

    ### compute the right hand side
    v = np.matmul(S_f_T_sgn,FxdM) - np.log(K_v)*np.sign(-y_grad)

    #construct the bounds
    x_upper = np.concatenate([vcount_upper_bound,1000*np.ones(n_react)])
    x_lower = np.concatenate([-300*np.ones(n_M_v), np.zeros(n_react)])


    #compute the solution
    opt_out = lsq_linear(A,-np.ravel(v),bounds = (x_lower,x_upper))


    n_out = opt_out['x']
    n_out = n_out[0:n_M_v]
    n_out = np.reshape(n_out,(len(n_out),1))


    v_true = np.matmul(S_f_T,FxdM) - np.log(K_v)
    flux_out = np.matmul(S_v_T,n_out) + v_true


    n_ini = np.ravel(n_out)
    y_ini = -1e1*np.ravel(flux_out)
    beta_ini = np.ravel(np.matmul(np.transpose(Sv_N),y_ini) )
    y_ini = np.matmul(Sv_N,beta_ini)
    
    
    #Set the initial condition
    b_ini = np.matmul(S_v_T,np.reshape(n_ini,(len(n_ini),1))) + np.matmul(S_f_T,np.reshape(FxdM,(n_M_f,1) )    )
    b_ini = np.reshape(b_ini,len(b_ini))
    
   
    h_ini = np.sign(y_ini)*( np.log(2) - np.log( np.abs(y_ini) + np.sqrt( np.power(y_ini,2) + 4 ) )  )
    

    
    '''
    ##################################

    Construct the Optimization Problem

    ##################################
    '''
    

    #Pyomo Model Definition
    ######################
    m =  pe.ConcreteModel()

    #Input the model parameters

    #set the indices
    #####################
    m.react_idx = pe.Set(initialize = react_idx) 

    m.TotM_idx = pe.Set(initialize = TotM_idx)
    m.VarM_idx = pe.Set(initialize = VarM_idx)
    m.FxdM_idx = pe.Set(initialize = FxdM_idx)

    m.beta_idx = pe.Set(initialize = beta_idx)

    # Stochiometric matrix
    ###########################
    S_idx = list(itertools.product(np.arange(0,n_M),np.arange(0,n_react)))
    S_vals = list(np.reshape(S,[1,n_M*n_react])[0])
    S_dict = dict(list(zip(S_idx,S_vals)))

    m.S = pe.Param(S_idx ,initialize = S_dict,mutable = True)

    ## Nullspace basis Matrix
    ####################
    SvN_idx = list(itertools.product(np.arange(0,n_react),np.arange(0,dSv_N)))
    SvN_vals = list(np.reshape(Sv_N,[1,n_react*dSv_N])[0])
    SvN_dict = dict(list(zip(SvN_idx,SvN_vals)))

    m.SvN=pe.Param(SvN_idx, initialize = SvN_dict)

    # Reaction Equilibrium constants
    ##################
    K_dict = dict(list(zip(react_idx,K) ))
    m.K=pe.Param(m.react_idx, initialize = K_dict)

    # Fixed metabolite log counts
    FxdM_dict = dict( list( zip(FxdM_idx,f_log_counts) ) )
    m.FxdM = pe.Param(m.FxdM_idx, initialize = FxdM_dict)

    # Bounds on the log of the metabolites
    #########
    M_ubnd_dict = dict(list(zip(VarM_idx,vcount_upper_bound)))
    m.VarM_ubnd = pe.Param(m.VarM_idx, initialize = M_ubnd_dict)


    # Objective Coefficients
    ##################
    obj_c_dict = dict(list(zip(react_idx,obj_coefs) ))
    m.obj_coefs=pe.Param(m.react_idx, initialize = obj_c_dict)




    #SET the Variables
    #############################


    ## Variable metabolites (log)
    ######################

    Mini_dict = dict(list(zip(VarM_idx,n_ini)))
    m.VarM = pe.Var(VarM_idx,initialize = Mini_dict)

    # steady state fluxes
    yini_dict = dict(list(zip(react_idx,y_ini)))
    m.y = pe.Var(react_idx, initialize = yini_dict)

    # flux null space representation
    betaini_dict = dict(list(zip(beta_idx,beta_ini)))
    m.beta = pe.Var(beta_idx, initialize = betaini_dict)

    # Steady state condition RHS
    bini_dict = dict(list(zip(react_idx,b_ini)))
    m.b = pe.Var(react_idx, initialize = bini_dict)

    hini_dict = dict(list(zip(react_idx,h_ini)))
    m.h = pe.Var(react_idx, initialize = hini_dict)



    # y sign variable
    y_sign_ini_dict = dict(list(zip(react_idx,.5+.5*np.sign(y_ini) )))
    m.u = pe.Var(react_idx,bounds=(0,1),initialize = y_sign_ini_dict)


    
    
    # Set the Constraints
    #############################


    #flux null space constraint
    def flux_null_rep(m,i):
        return m.y[i]   ==  sum( m.SvN[(i,j)]*m.beta[j]  for j in m.beta_idx )
    m.fxnrep_cns = pe.Constraint(m.react_idx, rule = flux_null_rep)

    

    def steady_state_metab(m,j):
        return m.b[j]  == sum( m.S[(k,j)]*m.VarM[k] for k in m.VarM_idx ) + sum( m.S[(k,j)]*m.FxdM[k] for k in m.FxdM_idx ) 
    m.ssM_cns = pe.Constraint(m.react_idx, rule = steady_state_metab)

    
    
    
    def h_cns(m,i):
        return m.h[i] ==(m.y[i]*lmbda_b/(abs(m.y[i])*lmbda_b + lmbda_a))*(pe.log(2) -  pe.log(abs(m.y[i]) + pe.sqrt(m.y[i]**2 + 4 ) )  ) 
    m.h_cns = pe.Constraint(m.react_idx, rule = h_cns)
    
    
        
    
    def relaxed_reg_cns_upper(m,i):
        return ( m.b[i] - pe.log(m.K[i]) ) >= m.h[i] - Mb*(m.u[i])  
    m.rxr_cns_up = pe.Constraint(m.react_idx,rule = relaxed_reg_cns_upper)
    
    
    def relaxed_reg_cns_lower(m,i):
        return ( m.b[i] - pe.log(m.K[i]) ) <= m.h[i] + Mb*(1 - m.u[i])
    m.rxr_cns_low = pe.Constraint(m.react_idx,rule = relaxed_reg_cns_lower)
    
    
    
    def sign_constraint(m,i):
        return (pe.log(m.K[i]) - m.b[i])*m.y[i] >= 0
    m.sign_y_cns = pe.Constraint(m.react_idx, rule = sign_constraint)
    
    
    def y_sign_relax(m,i):
        return 2*m.u[i] - 1 == (m.y[i]/(abs(m.y[i]) + lmbda_a) ) 
    m.y_sign_relax_cns = pe.Constraint(m.react_idx,rule = y_sign_relax)
    
   
    
    # Variable metabolite upper and lower bounds
    def M_upper_cnstrnts(m,i):
        return  m.VarM[i] <= m.VarM_ubnd[i]
    m.VarM_ub_cns = pe.Constraint(m.VarM_idx,rule = M_upper_cnstrnts)


    def M_lower_cnstrnts(m,i):
        return  m.VarM[i] >= VarM_lbnd
    m.VarM_lb_cns = pe.Constraint(m.VarM_idx,rule = M_lower_cnstrnts)

    
    
    
    # Set the Objective function

    def _Obj(m):
        return sum( m.y[j]*m.obj_coefs[j]  for j in m.react_idx )
    m.Obj_fn = pe.Objective(rule = _Obj, sense = pe.maximize) 
    


    
    '''
    ##################################

    Compute a Solution

    ##################################
    '''



    #Set the solver to use
    opt=pe.SolverFactory('ipopt', solver_io='nl')


    #Set solver otpions
    opts = {'max_iter': max_iter,
          'max_cpu_time': max_cpu_time,
          'tol':tolerance,
          'acceptable_tol':acc_tolerance}



    ## Solve the Model
    status_obj = opt.solve(m, options=opts, tee=True)




    ## Extract the Solution
    n_sol=np.zeros(n_M_v)
    y_sol = np.zeros(n_react)

    for i in react_idx:
        y_sol[i] = pe.value(m.y[i])

    for i in VarM_idx:
        n_sol[i] = pe.value(m.VarM[i])
        
    
    E_regulation = np.ones(len(y_sol))
    unreg_rxn_flux = np.ravel( rxn_flux(n_sol, f_log_counts,S_T, K, E_regulation) )
    alpha_sol = y_sol/unreg_rxn_flux
    alpha_sol = np.ravel(alpha_sol)

    return(y_sol, alpha_sol, n_sol)








def rxn_flux(v_log_counts, f_log_counts, S, K, E_regulation):
    ''' 
    Computes reaction fluxes as proportional to the thermodynamic odds of each reaction represented in the stoichiometeric matrix S with equilibrium
    constants K, given log metabolite counts v_log_counts and f_log_counts.
    '''


    # Flip Stoichiometric Matrix
    S_T = S    # S_T is the Stoich matrix with rows as reactions, columns as metabolites
    S = np.transpose(S) # this now is the Stoich matrix with rows metabolites, and columns reactions
    v_log_counts = np.reshape(v_log_counts,(len(v_log_counts),1))
    f_log_counts = np.reshape(f_log_counts,(len(f_log_counts),1))
    tot_log_counts = np.concatenate((v_log_counts,f_log_counts))
    K = np.reshape(K,(len(K),1))
    E_regulation = np.reshape(E_regulation,(len(E_regulation),1))
    
    forward_odds = K*np.exp(-.25*np.matmul(S_T,tot_log_counts) )**4
    reverse_odds = np.power(K,-1)*np.exp(.25*np.matmul(S_T,tot_log_counts) )**4
    
    return E_regulation*(forward_odds - reverse_odds )
    




