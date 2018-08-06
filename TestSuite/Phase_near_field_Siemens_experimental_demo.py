import sys
sys.path.append('../proxtoolbox/Problems/Phase')
sys.path.append('..')
import Siemens_processor
import Near_field_Siemens_experimental_in
from phase import Phase

Siemens = Phase(Near_field_Siemens_experimental_in.new_config)
Siemens.solve()
Siemens.show()
#Siemens.compare_to_matlab()
#Siemens.compare_data_to_matlab()
