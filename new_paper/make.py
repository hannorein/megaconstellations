import os
import numpy as np
import rebound
import sys
sys.path.insert(0,'../model')
import mega
for name in ["SXODC", "Sunrise", "Stampede"]:
#for name in ["Sunrise", "Stampede"]:
    constellations = {name: mega.constellations_all[name]}
    os.system("rm mega_"+name+".bin")
    sims = mega.get_simulations(constellations)
    print("Done",name,": N=", sims[name].N)
