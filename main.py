# play.py
# first attempt at a HTM implementation

import math
import random
import time
import os



#########################
# Scalar encoder function
#########################

class Encoder():
    """
    Provides encode(), which turns an int into a SDR. Holds a bunch of hard-coded
    constants.
    """

    # Record high and low temperature for NYC
    minval = -30
    maxval = 120

    # Length of output
    W = 11     # Number of 1's (must be odd)
    N = 400    # Total number of bits (must be > w). Must be a perfect square for now.

    # Derived quantities
    inputrange = float(maxval - minval)
    halfwidth = (W - 1)/2
    padding = halfwidth
    resolution = inputrange / (N - 2 * halfwidth)
    

    def encode(self, input):
        """
        Encodes an int in a SDR and puts the encoded value in a list.
        @param input: Data to encode, an int.
        @returns: list of bits of length W.
        """
    
        # Sanity checks
        if input == None: return [None]
    
        if input is not None and not isinstance(input, int):
            raise TypeError("Expected a scalar input but got input of type %s" % type(input))
    
        if input < self.minval:
            raise Exception('input (%s) less than minvalue (%s)' % (str(input), str(self.minval)))
    
        if input > self.maxval:
            raise Exception('input (%s) greater than maxval (%s - %s)' % (str(input), str(self.maxval)))
    
        if type(input) is float and math.isnan(input):
            input = None
    
        output = [0 for i in range(self.N)]
    
        # Compute the center bin. We use the first bit to be set in the output as the index
        centerbin = int(((input - self.minval) + self.resolution/2) / self.resolution) + self.padding
        minbin = centerbin - self.halfwidth
        bucketid = minbin
    
        if bucketid is None:
            # None is returned for missing value
            for i in range(self.N): output[i] = 0
    
        else:
            # The bucket index is the index of the first bit to set in the output
            for i in range(self.N): output[i] = 0
            minbin = int(bucketid)
            maxbin = int(minbin + 2 * self.halfwidth)
            assert minbin >= 0
            assert maxbin < self.N
            for i in range(minbin, maxbin + 1): output[i] = 1
    
        return output



##################
# Helper functions
##################

def getBitAt(input, pos):
    """
    Treats input as a 2D array and returns the bit at position (x,y) in the array.
    @param input: a SDR (list).
    @param pos: (x, y) pos in a 2D array.
    @returns: a bit.
    """
    return input[int(math.sqrt(len(input))) * pos[0] + pos[1]]



##############
# HTM classes
##############

class Synapse():
    """
    Segment of a synapse that connects to another cell. Keeps track of its presynpatic
    cell's position, permanence value and whether or not the Synapse is active.
    """
    
    threshold = 0.6
    """Threshold value for permanence. If self.permanence > threshold, then this synapse is valid."""

    permdelta = 0.01
    """Amount by which permanence value Increments or decrements when learning."""


    def __init__(self, inputpos, perm = 0.0):
        self.permanence = perm
        self.inputpos = inputpos


    @property
    def valid(self):
        """Whether or not this synapse is valid."""
        return self.permanence > self.threshold



class Dendrite():
    """
    Stems from cells (sometimes shared by a whole column). Keeps track of current Synapses
    and of cells that might form synapses while learning.
    """
    
    npotential = 20
    """The number of potential synapses for this Dendrite."""


    def __init__(self):
        def randpos(maxval):
            return (random.randint(0, maxval), random.randint(0, maxval))

        potential = [randpos(int(math.sqrt(Region.ncolumns)) - 1) for i in range(self.npotential)]
        self.synapses = [Synapse(p, random.random()) for p in potential]



class Cell():
    """
    A single computational unit. Keeps track of its dendrite segments, which determine
    where a Cell gets its input from and its active/inactive/predictive state.
    """

    
    def __init__(self, proximal):
        """
        A Column must send its shared Dendrite to each of its cells.
        @param proximal: a Dendrite, shared with other Cells in the same Column.
        """
        
        # Dendrite segments
        self.distal = []
        self.proximal = proximal
        
        # This is not a property b/c Columns tell their Cells when to activate
        self.state = "Inactive"



class Column():
    """
    An array of Cells. All the cells in one column share a single proximal Dendrite
    segment.
    """

    ncells = 4
    """Number of cells in this Column."""


    def __init__(self):
        self.proximal = Dendrite()
        self.cells = [Cell(self.proximal) for c in range(self.ncells)]
        self.boost = 1.0


    @property
    def state(self):
        """A Column is active whenever at least one of its Cells is active."""
        if "Active" in [c.state for c in self.cells]:
            return "Active"
        else:
            return "Inactive"


    @state.setter
    def state(self, state):
        # When a Column changes state, it must signal all of its Cells
        if state == "Active":
            predictive = [c for c in self.cells if c.state == "Predictive"]
            if predictive:
                # If there's at least one predictive cell, activate only that one
                for pred in predictive :
                    pred.state = "Active"
                for npred in list(set(self.cells) - set(predictive)):
                    npred.state = "Inactive"
            else:
                # else, activate every cell (surprise input!)
                for c in self.cells:
                    c.state = "Active"
                
        if state == "Inactive":
            for c in self.cells:
                c.state = "Inactive"



class Region():
    """A Region is an array of Columns."""

    ncolumns = 400 # needs to be a perfect square for now!
    """Number of columns in a Region."""
    
    density = 2./100
    """The sparsity of the active columns."""


    def __init__(self):
        self.columns = [Column() for c in range(self.ncolumns)]
        """A list of Columns. Call Region.columns2D for a 2D array."""


    @property
    def columns2D(self):
        """Return the Columns as a 2D array."""
        cols2D = []
        side = int(math.sqrt(self.ncolumns))
        for x in range(side):
            cols2D.append([])
            for y in range(side):
                cols2D[x].append(self.columns[x*side + y])

        return cols2D


    def process(self, input):
        """Process some sinput data (spatial pooler)."""
        
        # For any given input, determine how many valid synapses on each column are
        # connected to active input bits
        numvalid = {}

        for col in self.columns:
            numvalid[col] = sum([1 for s in col.proximal.synapses
                                     if s.valid
                                     and getBitAt(input, s.inputpos) == 1])
        
            # The number of active synapses is multiplied by a "boosting" factor which is
            # dynamically determined by how often a column is active relative to its neighbors.
            numvalid[col] *= col.boost

        print("Valid synapses: " + str(numvalid[self.columns[-1]]))
         
        # The columns with the highest activations after boosting disable all but a fixed
        # percentage of the columns within an inhibition radius. The inhibition radius is
        # itself dynamically determined by the spread (or "fan-out") of input bits. There
        # is now a sparse set of active columns.

        # inhibition code here

        ### CHANGE ME
        srtd = sorted(self.columns, key=(lambda c: numvalid[c]))
        index = Encoder.N - int(Encoder.N * self.density)
        inactivecols = srtd[:index]
        activecols = srtd[index:]
        for col in activecols:
            col.state = "Active"
        for col in inactivecols:
            col.state = "Inactive"
        ### CHANGE ME
        

        # For each active column, adjust the permanence of all the potential synapses.
        # These changes may validate inactive synapses, and vice-versa.
        for col in activecols:
            # Permanence of synapses with active input is increased.
            activesyn = [s for s in col.proximal.synapses
                             if getBitAt(input, s.inputpos) == 1]
            for s in activesyn:
                s.permanence = min(s.permanence + Synapse.permdelta, 1.0)
                s.valid   # Force recomputation of state

            # Permanence of synapses with inactive input bits are decreased.
            inactivesyn = list(set(col.proximal.synapses) - set(activesyn))
            for s in inactivesyn:
                s.permanence = max(0.0, s.permanence - Synapse.permdelta)
                s.valid   # Force recomputation of state

        # return numvalid


    def _prettyprint(self):
        """Prints a visualization of active/predictive/inactive columns."""
        # print_binary_matrix(list(map(lambda r: (list(map(lambda c: c.state == "Active", r))), self.columns2D)))
        pretty_print_matrix(list(map(lambda r: (list(map(lambda c: sum([cell.state == "Active" and 1 or 0 for cell in c.cells]), r))), self.columns2D)),
                                    {4: "■", 3: "X", 2: "x", 1: "o", 0: "·"})



##################
# Helper functions
##################

def pretty_print_matrix(mat, rules):
    for x in range(20):
        s = ""
        for y in range(20):
            s += rules[mat[x][y]]
        print(s)


def roll_array(data):
    mat = []
    side = int(math.sqrt(len(data)))
    for x in range(side):
        mat.append([])
        for y in range(side):
            mat[x].append(data[x*side + y])

    return mat


def log(msg):
    filename = "logs/main.txt"
    os.system("touch " + filename)
    os.system("echo '" + msg + "' > " + filename)



########
# MAIN
########
    

if __name__ == "__main__":
    log("\n\nStarting at " + str(time.time()))
    random.seed()
    
    enc = Encoder()
    r = Region()
    while(True):
        x = random.randint(enc.minval, enc.maxval - 1)

        data = enc.encode(x)
        print("Input: " + str(x))
        print("\n")
        
        r.process(data)
        r._prettyprint()

        print("\n")
        for tu in [(s.inputpos, round(s.permanence, 2)) for s in r.columns[-1].proximal.synapses]:
            print(tu)
        
        time.sleep(1)
        os.system("clear")
