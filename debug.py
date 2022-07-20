import pyscal as pc
from calCD import calCD
from calCD import A
from multiprocessing import Pool
import cChargeAnalysis 
from multiprocessing import get_context
from copy import deepcopy
import dill
#file=files[40]

if __name__ == '__main__':
    a=cChargeAnalysis.ChargeAnalysis()
    dill.dumps(a)
    #calCD([3],a,[])