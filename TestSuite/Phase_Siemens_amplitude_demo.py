import sys
sys.path.append('../proxtoolbox/Problems/Phase')
sys.path.append('..')
import Siemens_processor
import Siemens_amplitude_in
from phase import Phase

Siemens = Phase(Siemens_amplitude_in.new_config)
Siemens.solve()
Siemens.show()
