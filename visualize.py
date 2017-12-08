"""
Animated visualization suite

Written with reference to
https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot
"""
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import Molec

class MolecularViewer:
    """
    No attempts are made to deallocate memory used by the simulation.
    This should be handled in a try finally block surrounding the
    instantiation and use of this class.
    """
    def __init__(self, lennard_jones_simulation, updates_per_display = 10,
                 num_cycles = None):
        self.sim = lennard_jones_simulation
        self.fig, self.ax = plt.subplots()
        self.stream = self.data_stream()
        self.xlim = [0, self.sim.get_box_width()]
        self.ylim = [0, self.sim.get_box_width()]

        self.anim = animation.FuncAnimation(self.fig, self.update,
                                            interval=1,
                                            init_func=self.setup_plot,
                                            blit=True)

    def setup_plot(self):
        data = self.sim.get_center_data(plottable = True)
        self.scat = self.ax.scatter(*data, animated=True)
        self.ax.axis([*self.xlim, *self.ylim])
        return self.scat,
        
    def data_stream(self):
        while True:
            self.sim.step()
            #import pdb; pdb.set_trace()
            data = self.sim.get_center_data(plottable = False)
            # print(data)
            yield data

    def update(self, i):
        """
        i: just a parameter required by matplotlib animation
        """
        data = next(self.stream)
        #print(data)
        self.scat.set_offsets(data)

        # TODO: add settings for different colors for a subset

        return self.scat,

    def show(self):
        plt.show()

if __name__ == "__main__":
    print("Getting initial data (solid phase)")
    initial_data = Molec.get_initial_data(100, 30, 8, 14, 1000)
    print(len(initial_data))
    print("Done")
    sim = Molec.LennardJones(initial_data, m = 1., rmin = 1., epsilon = 0.8,
                             T = 0.01, alpha = 0.2, h = 1e-4, box_width = 30)
    #import pdb; pdb.set_trace()
    try:
        viewer = MolecularViewer(sim)
        viewer.show()
    finally:
        sim.deconstruct()
