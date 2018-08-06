import sys
sys.path.append('../proxtoolbox/Problems/Phase')
sys.path.append('..')
import Siemens_processor
import Siemens_real_in
from phase import Phase

Siemens = Phase(Siemens_real_in.new_config)
Siemens.solve()
Siemens.show()
