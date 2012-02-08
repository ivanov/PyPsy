import numpy as np
import scipy.special as special

def ptoz(p):
    z = np.sqrt(2)*special.erfinv(p*2-1);
    return z

def afc2_p_to_dprime(ph, pf=None):
    #function [dP C lnB pC]=PAL_SDT_2AFC_PHFtoDP(pHF)

    if pf is None:
        pf = ph[:,1]
        ph = ph[:,0]

    zH=ptoz( ph )
    zF=ptoz( pf )

    dP=(zH-zF)/np.sqrt(2)
    C=-(zH+zF)/np.sqrt(2)
    lnB=(zF**2-zH**2)/2
    pC=(ph+(1.0-pf))/2

    return dP, C, lnB, pC

# p163
demo_prins_h = numpy.array( [0.6, 0.7, 0.8] )
demo_prins_f = numpy.array( [0.2, 0.2, 0.2] )
demores = afc2_p_to_dprime( demo_prins_h, demo_prins_f )

# p189
prins_tab62_h = numpy.array([0.61,0.69,0.79,0.88,0.97,0.99]) 
prins_tab62_f = numpy.array([0.53,0.42,0.33,0.18,0.06,0.03]) 
hw2res= afc2_p_to_dprime( prins_tab62_h, prins_tab62_f )

"""
Dear class,
My software is now working.
The problem with the 2AFC problem is that many psychophysicists are confused about it and factors of sqrt(2). So our discussion on Thursday will hopefully clarify almost everything you wanted to know about signal detection theory and 2AFC.  I'll try to undo all your confusion.  I had thought that the Exercise 2 had asked you to calculate the 2AFC dprime the normal way assuming there is no bias. Thus the assignment I had wanted you to do for  is the data on p. 189 is do calculate d' using PAL_SDT_2AFC_PHFtoDP and also using PAL_SDT_MAFC_PCtoDP.
Hint:   The answer to the latter problem is   dprimeNoBias= erfinv(2*pC-1)*2
 
 I much prefer you to calculate dprimeNoBias rather than the requested PcMax that is complicated.  Kingdom/Prins made the problem more confusing than it had to be (as seen by how simple my Hint is).
 For the plot one would plot the two dprimes (with and without bias). 
  
  Given that this is an extra thing to do (though simpler than the original assignment) you can have extra time say until noon Wednesday. Send me what you have at that point and then if you have further questions you can join me at 4:30 in my office at 420 Minor Hall (not the 520 I said before).  
   
   Stan
"""
