import numpy as np
from scipy import optimize

#function usage as follows:
#
# conf, m = getM(counts=[0,1,1,2,8,... etc colony counts], 
#                cellN=#ofcells_in_culture, 16, r_max=170)
#
#           16 ~ number of independent cultures used in the experiment
#           r_max ~ max number of mutants / can leave at 170 if all mutant #s < 170
#
#     returns: 
#                m ~ mutation rate
#                conf[0], conf[1] ~ 95% confidence intervals



### Probability function of observing 'r' mutants in a given culture Given a rate m  $$p_0 = e^{-m}$$   $$p_r = \frac{m}{r} \sum_{i=0}^{r-1} \frac{p_i}{(r-i+1)} $$   Maximum likelihood function (this is what needs to be maximized.   $$f(r|m) = \prod_{i=1}^{c} f(r_i|m) = (p_0)^{c_0} (p_1)^{c_1} (p_2)^{c_2} (p_{r_{max}})^{c_{r_{max}}} $$    We also have the following MMS-MLE confidence intervals for the MSS-MLE method    $$CL_{+95\%} = ln(m) + 1.96 \sigma (e^{1.96 \sigma})^{-0.315}$$  $$CL_{-95\%} = ln(m) - 1.96 \sigma (e^{1.96 \sigma})^{+0.315}$$  with $$\sigma_{ln(m)} \approx \frac{1.225m^{-0.315}}{\sqrt{C}} $$

def Pr_ARRAY(m, r=0, P_R=[], r_max=50): 
    #where m is mutations
    #where 'r' is observed mutations?
    #r_max is the maximum value of p_r calculated (most number of cultures we would count)
    if r==0:
        P_R=np.zeros((r_max,np.shape(m)[0]))
        
        #print(P_R.shape)
        
        P_R[0]=np.e**-m
        return Pr_ARRAY(m, 1, P_R, r_max)
    
    elif r==r_max:
        return P_R.T
    
    else:
        quot   = np.arange(2,r+2,dtype=np.float)[::-1] **-1
        
        #print(P_R[:r])
        #print(quot)
        #print(np.sum(P_R[:r]*(quot),axis=1))
        #print(m.shape)
        #print(np.sum(P_R[:r].T*(quot),axis=0).shape)
        
        P_R[r] = m/r*np.sum(P_R[:r].T*quot,axis=1)
    
    return Pr_ARRAY(m, r+1, P_R, r_max)


def Pr(m, r=0, P_R=[], r_max=50): 
    #where m is mutations
    #where 'r' is observed mutations?
    #r_max is the maximum value of p_r calculated (most number of cultures we would count)
    if r==0:
        P_R=np.zeros(r_max)
        
        #print(P_R.shape)
        
        P_R[0]=np.e**-m
        return Pr(m, 1, P_R, r_max)
    
    elif r==r_max:
        return P_R
    
    else:
        quot   = np.arange(2,r+2,dtype=np.float)[::-1] **-1
        P_R[r] = m/r*P_R[:r].dot(quot)
    
    return Pr(m, r+1, P_R, r_max)



def mut_to_array(x=[0,0,0,0,0,0,0,0],r_max=50):
    #takes an array x to a histogram of x
    #e.g. x=[1,4,4,0] is transformed to y=[1,1,0,0,2,0]
    
    temp = np.zeros(r_max)
    for i in x:
        temp[i]+=1
        
    return temp
    
    
def f_MLF(P_R, C):
    return np.log(P_R).dot(C)

def MLF(m,Y,r_max=50):
    #m is the mutation rate probe
    #Y is the number of mutations matrix derived from mut-to-array
    
    #basically we are returning the 'cost' function
    #here we have log maximum likelihood estimator for data 'Y' and free parameter 'm'
    if isinstance(m, float): return f_MLF(Pr(m, r_max=r_max), Y)
    
    else: return f_MLF(Pr_ARRAY(m, r_max=r_max), Y)
    
    
def optimize_test(Y, x_1=.0001, x_2=3, r_max=50, mult=2):
		print("WARNING, no longer functions correctly missing LOG/EXP corrections")    
		m = np.logspace(x_1,x_2,num=10)
		m_arg = np.nanargmax(MLF(m,Y,r_max=r_max)) 
    #print(MLF(m,Y))
    #print(m_arg)
    #print(m)
    #print(x_1, x_2)
		x_1_new = np.log10(m[m_arg]-mult)
		x_2_new = np.log10(m[m_arg]+mult)    
		#print(x_1_new, x_2_new)    
		return x_1_new, x_2_new, mult
    

def fun(m,y,r_max):
		return -MLF(m,y,r_max)

def confidence(m,cells,c=8):
    sigma=1.96 * 1.225*m**-.315/np.sqrt(c)
    cl_upp=np.log(m) + sigma*(np.e**sigma)**-.315
    cl_low=np.log(m) - sigma*(np.e**sigma)**+.315
    return np.e**cl_upp/cells*1E7, np.e**cl_low/cells*1E7


def getM(x, cells, c=8, r_max=50, decay=.3, draws=60, mult=2):
		y=np.array(mut_to_array(x,r_max))
				
		#quick python module optimizer for 1 paramter optimization (easy anyway since it is convex)
		z = optimize.minimize_scalar(fun,bounds=[0,100],args=(y,r_max))
		m = z.x

		#print(m)
		
		#limit, m = optimizer(y,decay,r_max=r_max,decay=decay,draws=draws,mult=mult)
		#print(m)    
    #returns mutation rate (1E7 scaled) 
    #returns the upper and lower bound intervals
    
    #get format to 1 decimal place.
		return confidence(m,cells,c), m/cells*1E7
