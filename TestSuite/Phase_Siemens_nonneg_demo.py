import sys
sys.path.append('../proxtoolbox/Problems/Phase')
sys.path.append('..')
import Siemens_processor
import Siemens_nonneg_in
from phase import Phase

Siemens = Phase(Siemens_nonneg_in.new_config)
Siemens.solve()
Siemens.show()
