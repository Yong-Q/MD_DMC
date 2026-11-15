import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from transport_v3 import select_full

def test_full_curve_selects_linear():
 t=np.arange(15000.);r=select_full(t,2*t+1)
 assert r['status']=='candidate' and abs(r['fit']['alpha']-1)<.01

def test_full_curve_rejects_sqrt_and_ballistic():
 t=np.arange(15000.)
 assert select_full(t,np.sqrt(t+1))['status']=='no_candidate'
 assert select_full(t,t*t)['status']=='no_candidate'

def test_no_window_when_trajectory_short():
 t=np.arange(1000.);assert select_full(t,t+1)['status']=='no_candidate'
