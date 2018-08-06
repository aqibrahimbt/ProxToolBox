import sys
sys.path.append('../proxtoolbox/Problems/CT')
sys.path.append('..')
from ART import ART
import ART_alternating_proj_in

CT = ART(ART_alternating_proj_in.new_config)
CT.solve()
CT.show()
