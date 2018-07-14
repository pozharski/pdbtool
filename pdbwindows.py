from matplotlib.figure import Figure
from matplotlib.pyplot import grid
class BWindow(Figure):
    def __init__(self, *args, **kwds):
        Figure.__init__(self, *args, **kwds)
    def set_model(self, model):
        self.model = model
    def plot(self):
        if 'resids' not in dir(self):
            self.resids, self.b0, self.b1, self.b2 = self.model.GetResidueBvectorByChain()
            self.resnums = dict([(k,[int(x) for x in v]) for k,v in self.resids.items()])
        self.add_subplot(111)
        self.axes[0].clear()
        self.pt = {}
        for chid in self.resids:
            self.pt[chid] = self.axes[0].plot(self.resnums[chid],self.b0[chid])
        grid(True)
        self.canvas.draw()
