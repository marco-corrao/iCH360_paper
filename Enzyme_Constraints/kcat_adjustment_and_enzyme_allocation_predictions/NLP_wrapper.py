"""
A simple wrapper around the package CadADi, useful to construct parametric NLPs in a user-friendly way, while mainitaining access to lower-level functionality.
"""
import numpy as np
import casadi as ca

class NLP():
    def __init__(self):
        self.x=[]
        self.g=[]
        self.lbg=[]
        self.ubg=[]
        self.lbx=[]
        self.ubx=[]
        self.f=0.
        self.x0=[] #initial guesses
        self.solver_fcn=None
        self.variables={}
        self.parameters={}

        self.default_x0=2
    def add_parameter(self,var,name):
        self.parameters[name]=var

    def add_variable(self,var,name,lb=None,ub=None,x0=None):
        self.x.append(var)
        self.variables[name]=var
        if lb is None:
            self.lbx.append([-ca.inf]*var.shape[0])
        elif isinstance(lb,float):
            self.lbx.append([lb]*var.shape[0])
        elif isinstance(lb,list):
            assert len(lb)==var.shape[0]
            self.lbx.append(lb)
        elif isinstance(lb,np.ndarray):
            lb_list=lb.flatten().tolist()
            assert len(lb_list)==var.shape[0]
            self.lbx.append(lb_list)
        else:
            print(f"Unable to parse lower bounds for {name}")

            #Now upper bounds
        if ub is None:
            self.ubx.append([ca.inf]*var.shape[0])
        elif isinstance(ub,float):
            self.ubx.append([ub]*var.shape[0])
        elif isinstance(ub,list):
            assert len(ub)==var.shape[0]
            self.ubx.append(ub)
        elif isinstance(ub,np.ndarray):
            ub_list=ub.flatten().tolist()
            assert len(ub_list)==var.shape[0]
            self.ubx.append(ub_list)
        else:
            print(f"Unable to parse upper  bounds for {name}")
        #initial guess
        if x0 is None:
            self.x0.append([self.default_x0]*var.shape[0])
        elif isinstance(x0,float):
            self.x0.append([x0]*var.shape[0])
        elif isinstance(x0,list):
            assert(len(x0))==var.shape[0]
            self.x0.append(x0)
        else:
            assert x0.shape[0]==var.shape[0]
            self.x0.append(x0)
            

    def add_constraint(self,expression,sense,rhs):
        self.g.append(expression)
        if isinstance(rhs,float):
            rhs=rhs*np.ones(expression.shape)
        else:
            assert rhs.shape==expression.shape
        if sense=='>':
            self.lbg.append(rhs)
            self.ubg.append([ca.inf]*rhs.shape[0])
        elif sense=='<':
            self.ubg.append(rhs)
            self.lbg.append([-ca.inf]*rhs.shape[0])
        elif sense=='=':
            self.ubg.append(rhs)
            self.lbg.append(rhs)
        else:
            raise NotImplemented
    def set_objective(self,expression,sense):
        """
        By default, IPOPT minimize the objective so we flip the sign in case of maximisation
        """
        if sense=='max':
            self.f=-expression
        elif sense=='min':
            self.f=expression

    def construct_nlp(self,opts={}):
        self.w=ca.vertcat(*self.x)
        self.par_names=list(self.parameters.keys())
        self.p=ca.vertcat(*[self.parameters[par_name] for par_name in self.par_names])
        ipopt_nlp = {'x':self.w, 
                    'f':self.f,
                    'g':ca.vertcat(*self.g),
                    'p':self.p
                    }
        self.solver_fcn=ca.nlpsol("solver", 'ipopt', ipopt_nlp, opts)

        #Now construct a retrival function for the optimal values of the optimisation variables
        self.retrieve_optimal_values=ca.Function('optimum_retrieval',
                                                 [self.w,self.p],
                                                 [var for var in self.variables.values()],
                                                 ['w','p'],
                                                 [name for name in self.variables.keys()]
                                                 )
    def add_evaluation_expression(self,expression,name):
        self.variables[name]=expression
        
    def solve(self,parameters={},warm_start=None,verbose=False):
        p_vector=[]
        for par_name in self.par_names:
            if par_name not in parameters.keys():
                print(f"Missing parameter {par_name}. Unable to run the solver")
                raise ValueError
            if isinstance(parameters[par_name],float):
                p_vector.append([parameters[par_name]])
            elif isinstance(parameters[par_name],list):
                p_vector.append(parameters[par_name])
        if len(p_vector)>0:
            p_vector=np.concatenate(p_vector)
        if warm_start is None:
            sol = self.solver_fcn(lbg=ca.vertcat(*self.lbg),
                            ubg=ca.vertcat(*self.ubg),
                            lbx=np.concatenate(self.lbx),
                            ubx=np.concatenate(self.ubx),
                            x0=np.concatenate(self.x0),
                            p=p_vector
                            )
        else:
            # In this case, we reconstruct the solver with the warm start options
            warm_start_opts={'ipopt':{'warm_start_init_point':'yes',
                                      'warm_start_bound_push':1e-8,
                                      'mu_init':1e-8,
                                      'warm_start_mult_bound_push':1e-8,
                                      'warm_start_slack_bound_push':1e-8,
                                      'warm_start_bound_push':1e-8,
                                      'max_iter':10000,
                                        'print_level':5 if verbose else 0,
                                      },
                                        'print_time':0
                                      }
            self.construct_nlp(opts=warm_start_opts)

            sol = self.solver_fcn(lbg=ca.vertcat(*self.lbg),
                            ubg=ca.vertcat(*self.ubg),
                            lbx=np.concatenate(self.lbx),
                            ubx=np.concatenate(self.ubx),
                            x0=warm_start['x'],
                            lam_g0=warm_start['lam_g'],
                            lam_x0=warm_start['lam_x'],
                            p=p_vector
                            )
    
        
        optimal_values=self.retrieve_optimal_values(w=sol['x'],p=p_vector)
        optimal_values_numpy={}
        for key,value in optimal_values.items():
            optimal_values_numpy[key]=np.array(value)
        return optimal_values_numpy,sol,self.solver_fcn.stats()
