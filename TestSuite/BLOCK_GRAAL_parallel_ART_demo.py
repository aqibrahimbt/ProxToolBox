import sys
sys.path.append('../proxtoolbox/Problems/CT')
sys.path.append('..')
from ART import ART
import BLOCK_GRAAL_parallel_ART_in

CT = ART(BLOCK_GRAAL_parallel_ART_in.new_config)
CT.solve()
CT.show()
