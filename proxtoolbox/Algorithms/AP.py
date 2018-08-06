from .SimpleAlgortihm import SimpleAlgorithm

class AP(SimpleAlgorithm):

    def evaluate(self, u):
        tmp_u = self.prox2.work(u)
        unew = self.prox1.work(tmp_u)
        return unew
