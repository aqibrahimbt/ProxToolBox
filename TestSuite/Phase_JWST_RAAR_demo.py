import sys
sys.path.append('../proxtoolbox/Problems/Phase')
sys.path.append('..')
import JWST_data_processor
import JWST_RAAR_in
from phase import Phase

JWST = Phase(JWST_RAAR_in.new_config)
JWST.solve()
JWST.show()
