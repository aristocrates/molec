"""
Animated visualization suite

Written with reference to
https://stackoverflow.com/questions/9401658/how-to-animate-a-scatter-plot
"""
import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import Molec
import numpy as np
import os

class MolecularViewer:
    """
    No attempts are made to deallocate memory used by the simulation.
    This should be handled in a try finally block surrounding the
    instantiation and use of this class.
    """
    def __init__(self, lennard_jones_simulation, num_cycles = None):
        self.num_cycles = num_cycles

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
        counter = 0
        target = self.num_cycles
        if self.num_cycles is None:
            target = -1
        while counter != target:
            counter += 1
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
        #print(data.shape)
        self.scat.set_offsets(data)

        # TODO: add settings for different colors for a subset

        return self.scat,

    def show(self):
        plt.show()

def zero_pad(n, length):
    ans = str(n)
    ans = "0"*(length - len(ans)) + ans
    return ans
        
def write_plot(x, y, n, digits, fileformat, directory, xlim, ylim):
    plt.scatter(x, y)
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.savefig(os.path.join(directory, "figure" + zero_pad(n, digits)
                             + "." + fileformat))
    plt.close()

def molecule_gif(lennard_jones_simulation, image_dir, nframes = 1000,
                 steps_per_frame = 10, verbose = True, fileformat = "png"):
    assert(nframes >= 0)
    assert(int(nframes) == nframes)
    digits_needed = 1 + int(np.floor(np.log10(nframes)))
    max_range = lennard_jones_simulation.get_box_width()
    for i in range(nframes):
        if verbose:
            print(i)
        for j in range(steps_per_frame):
            lennard_jones_simulation.step()
        data = lennard_jones_simulation.get_center_data(plottable = True)
        write_plot(*data, i, digits_needed, fileformat, image_dir,
                   [0, max_range], [0, max_range])

if __name__ == "__main__":
    T = 0
    generate_initial_data = False
    initial_data = [[10.447833222655293, 14.568590866348524],
                    [10.955711564419916, 13.705105881991935],
                    [12.455736695446042, 14.458385725680376],
                    [18.327575746576574, 18.260198743422801],
                    [15.033609708034794, 10.380903859316414],
                    [10.595908032629549, 18.863411627249342],
                    [12.874226554506739, 16.125643959679394],
                    [18.971859968748678, 11.72376059966542],
                    [10.236525971008868, 10.911828527773878],
                    [9.5804995335412872, 12.298415733547381],
                    [14.574158608003229, 12.147325253724805],
                    [11.474358880880841, 14.529641153447949],
                    [14.463317201050796, 19.560116155635601],
                    [13.558901569928709, 17.659957790553531],
                    [18.154306290814006, 11.20440747441822],
                    [17.541229247487255, 15.511344424805865],
                    [15.506727977076945, 17.657606085479763],
                    [16.952808427899818, 9.9885769167518994],
                    [18.923776403675529, 10.71038774760903],
                    [14.878740202954047, 9.4225991999147638],
                    [9.8016434348856549, 18.186695093786795],
                    [14.492414575509949, 15.85393986987688],
                    [16.780956447899154, 19.220823338363669],
                    [17.255766700533055, 18.324602966718938],
                    [19.572406829212088, 17.123527148615626],
                    [12.935805298114087, 13.497876902339586],
                    [14.114859129741077, 14.935680680025262],
                    [11.408848242074576, 19.630368845568729],
                    [17.931279373456043, 10.19088293565432],
                    [18.170464639820029, 9.1861379010140212],
                    [13.263848507954778, 10.742425114720026],
                    [13.554777020597312, 12.684834885410845],
                    [18.917277005572082, 14.45885953475344],
                    [15.946216528854888, 14.891197344422851],
                    [16.400378952516188, 11.611263171183108],
                    [16.020361826828726, 16.762758392536348],
                    [11.271836271759557, 12.656556841199885],
                    [9.3304052823993917, 11.36581528758712],
                    [19.78706613534542, 11.126796571700742],
                    [19.010969069007015, 16.257578759380742],
                    [18.010091782992973, 16.330282218974315],
                    [19.660173931487513, 10.108407182200439],
                    [16.224905680735933, 20.530048986072071],
                    [12.03865675432308, 15.428526687303638],
                    [19.509925951387984, 15.238952871169253],
                    [14.087902901947288, 10.084713128043576],
                    [20.656847979674783, 10.52967944943838],
                    [16.092243222467204, 18.472709283350337],
                    [10.643544642450271, 11.806400425915198],
                    [14.000653636790295, 13.534807930481083],
                    [14.037252356851196, 18.533475518110844],
                    [12.553098790803014, 11.526825780903247],
                    [18.66345647834234, 17.348508418114005],
                    [16.539882468793497, 15.86781880999389],
                    [17.203333611095161, 10.943406141369286],
                    [8.479464548718715, 12.149006634980065],
                    [17.003381612099751, 16.702190398021674],
                    [20.422509866913479, 17.826236833811098],
                    [15.401548686196762, 11.504837898264487],
                    [15.016734122448433, 15.074720103388319],
                    [15.588431847880511, 15.850727025136326],
                    [10.936652515450824, 15.455252345952864],
                    [10.116558217281137, 17.245404652401152],
                    [13.415692392493106, 14.333719803336329],
                    [13.909538023306167, 16.774316877938503],
                    [14.514466013777659, 17.671457250016026],
                    [17.756710479124184, 19.210753703526912],
                    [20.609875415579904, 16.899652532954409],
                    [14.924928728746226, 16.772759133204303],
                    [13.322871989216152, 19.175613871172171],
                    [11.671680097661882, 11.811808494733011],
                    [16.050078726439175, 10.734168715484193],
                    [17.687147268736194, 17.372729231789492],
                    [11.263006219859257, 20.557711780317451],
                    [12.270683663371901, 20.247634958344754],
                    [12.339247227444973, 19.206097371375666],
                    [9.7365854707487429, 19.24492249766514],
                    [12.693758599587287, 18.233782070016986],
                    [16.63386157182315, 17.606372413184257],
                    [12.587870978955866, 17.231004550842837],
                    [19.95771548064462, 16.189751834803495],
                    [16.852906786357497, 14.651834190270334],
                    [10.805455260975116, 17.840779909825645],
                    [17.891402867136236, 14.614768614539672],
                    [11.514822045853574, 18.510113948935317],
                    [11.973205141595093, 13.436604723807248],
                    [15.019546920072127, 18.458729955752784],
                    [13.050560506254211, 15.202265341553963],
                    [13.51684718597922, 11.741319855554623],
                    [20.965429269996541, 15.99023118676328],
                    [14.326049288990864, 11.171185696506955],
                    [15.842758473960268, 9.7314440757329912],
                    [19.994809178816233, 14.354098298613339],
                    [19.411285411805974, 18.090239855707456],
                    [15.601230605712946, 19.380522956264382],
                    [15.283694168701846, 20.314484962567825],
                    [18.500090169527944, 15.351112578402647],
                    [10.183207618947272, 13.056614218275612],
                    [10.487767225962056, 20.002465459942382],
                    [12.574348939394364, 12.552897828800294]]
    if len(sys.argv) > 1:
        T = float(sys.argv[1])
    movie = False
    fileformat = "png"
    if len(sys.argv) > 2:
        if sys.argv[2] == "movie":
            print("Making movie instead of displaying to screen")
            movie = True
            image_dir = sys.argv[3]
            nframes = int(sys.argv[4])
        elif sys.argv[2] == "pdf":
            print("Making pdf frames instead of displaying to screen")
            movie = True
            fileformat = "pdf"
            image_dir = sys.argv[3]
            nframes = int(sys.argv[4])
    if len(sys.argv) > 5:
        print("Generating new initial data. Your mileage may vary")
        print("Getting initial data (solid phase)")
        initial_data = Molec.get_initial_data(100, 30, 8, 14, 1000)
        print(len(initial_data))
        print("Done")
    sim = Molec.LennardJones(initial_data, m = 1., rmin = 1., epsilon = 0.8,
                             T = T, alpha = 0.2, h = 1e-4, box_width = 30)
    #import pdb; pdb.set_trace()
    try:
        if movie:
            molecule_gif(sim, image_dir, nframes, fileformat = fileformat)
        else:
            viewer = MolecularViewer(sim)
            viewer.show()
    finally:
        sim.deconstruct()
