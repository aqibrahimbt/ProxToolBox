import sys
sys.path.append('../proxtoolbox/Problems/Phase')
sys.path.append('..')
import tasse_raar_in
from phase import Phase


tasse = Phase(tasse_raar_in.new_config)
tasse.solve()
tasse.show()