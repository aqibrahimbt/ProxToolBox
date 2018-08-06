import sys
sys.path.append('../proxtoolbox/Problems/Phase')
sys.path.append('..')
import Siemens_processor
import Near_field_living_worm_in
from phase import Phase

Siemens = Phase(Near_field_living_worm_in.new_config)
Siemens.solve()
Siemens.show()