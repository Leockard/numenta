# play.py
# first attempt at a HTM implementation

import math
import random
import time
import os



#########################
# Scalar encoder function
#########################

def encode(input):
    """
    Encodes an int in a SDR and puts the encoded value in a list.  @param input: Data to
    encode. This should be validated by the encoder.  @returns: list of same length W
    """
    # Record high and low temperature for NYC
    minval = -30
    maxval = 120

    # Length of output
    W = 11     # number of 1's (must be odd)
    N = 400    # total number of bits (must be > w)

    # Derived quantities
    inputrange = float(maxval - minval)
    resolution = inputrange / N
    halfwidth = (W - 1)/2
    padding = halfwidth

    # Sanity checks
    if input == None: return [None]

    if input is not None and not isinstance(input, int):
        raise TypeError("Expected a scalar input but got input of type %s" % type(input))

    if input < minval:
        raise Exception('input (%s) less than minvalue (%s)' % (str(input), str(minval)))

    if input > maxval:
        raise Exception('input (%s) greater than maxval (%s - %s)' % (str(input), str(maxval)))

    if type(input) is float and math.isnan(input):
        input = None


    output = [0 for i in range(N)]

    # Compute the center bin. We use the first bit to be set in the output as the index
    centerbin = int(((input - minval) + resolution/2) / resolution) + padding
    minbin = centerbin - halfwidth
    bucketid = minbin

    if bucketid is None:
        # None is returned for missing value
        for i in range(N): output[i] = 0

    else:
        # The bucket index is the index of the first bit to set in the output
        for i in range(N): output[i] = 0
        minbin = int(bucketid)
        maxbin = int(minbin + 2 * halfwidth)
        assert minbin >= 0
        assert maxbin < N
        for i in range(minbin, maxbin + 1): output[i] = 1

    return output


##################
# Helper functions
##################

def getBitAt(input, pos):
    return input[20 * pos[0] + pos[1]]



##############
# HTM classes
##############

class Synapse():

    threshold = 0.6

    def __init__(self, inputpos, perm = 0.0):
        self.permanence = perm
        self.inputpos = inputpos
        self.valid = perm > self.threshold


class Dendrite():

    npotential = 20

    def __init__(self):
        self.potential = [(random.randint(0, int(math.sqrt(Region.ncolumns)) - 1), random.randint(0, int(math.sqrt(Region.ncolumns)) - 1)) for i in range(self.npotential)]
        self.synapses = [Synapse(p, random.random()) for p in self.potential]


class Cell():
    
    def __init__(self, proximal):
        self.distal = []
        self.proximal = proximal
        
        # This is not a property because Columns must tell their Cells whether or not
        # they're active
        self.state = "Inactive"   


class Column():

    def __init__(self):
        self.proximal = Dendrite()
        self.cells = [Cell(self.proximal) for c in range(4)]
        self.boost = 1.0


    @property
    def state(self):
        if "Active" in [c.state for c in self.cells]:
            return "Active"
        else:
            return "Inactive"

    @state.setter
    def state(self, state):
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

    ncolumns = 400 # needs to be a perfect square for now!

    def __init__(self):
        self.columns = [Column() for c in range(self.ncolumns)]


    @property
    def columns2D(self):
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
         
        # The columns with the highest activations after boosting disable all but a fixed
        # percentage of the columns within an inhibition radius. The inhibition radius is
        # itself dynamically determined by the spread (or "fan-out") of input bits. There
        # is now a sparse set of active columns.

        # inhibition code here

        ### CHANGE ME
        states = {}
        index = int(400 - math.floor(len(input) * (2./100)))
        srtd = sorted(self.columns, key=(lambda c: numvalid[c]))
        for col in srtd[:index]:
            col.state = "Inactive"
        for col in srtd[index+1:]:
            col.state = "Active"

        
        # For each of the active columns, we adjust the permanence values of all the potential
        # synapses. The permanence values of synapses aligned with active input bits are
        # increased. The permanence values of synapses aligned with inactive input bits are
        # decreased. The changes made to permanence values may change some synapses from being
        # valid to not valid, and vice-versa.

        # learning code here


        # When a column becomes active, it looks at all its cells. If one or more cells
        # are already in the predictive state, only those cells become active.

        # for col in self.columns: col.state = states[col]


        # return numvalid


    def _prettyprint(self):
        """Prints a visualization of active/predictive/inactive columns."""
        for x in range(20):
            s = ""
            
            for y in range(20):
                if self.columns2D[x][y].state == "Active":
                    s += "■"
                else:
                    s += "."
            print(s)



if __name__ == "__main__":
    r = Region()
    while(True):
        r.process(encode(random.randint(-30, 110)))
        r._prettyprint()
        print("\n")
        time.sleep(1)
        os.system("clear")
