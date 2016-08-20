"""
Attempt at a HTM implementation directly from the paper.
"""

import math
import random
import time
import os


###########################################################################
# Scalar encoder function
###########################################################################

class Encoder:
    """
    Provides encode(), which turns an intteger into a SDR.  Holds a bunch of
    hard-coded constants.
    """

    # Record high and low temperature for NYC
    minval = -30
    maxval = 120

    # Length of output
    W = 11     # Number of 1's (must be odd)
    N = 400    # Total number of bits (must be > w). Perfect square for now.

    # Derived quantities
    inputrange = float(maxval - minval)
    halfwidth = (W - 1)/2
    padding = halfwidth
    resolution = inputrange / (N - 2 * halfwidth)


    def encode(self, integer):
        """Encodes a single integer in a SDR.

        @param integer: Data to encode, an int.
        @returns: list of bits of length Encoder.W.
        """
        assert isinstance(integer, int), \
            'Expected an int but got {}'.format(type(integer))
        assert integer >= self.minval, \
            'Input {} less than minvalue {}'.format(integer, self.minval)
        assert integer < self.maxval, \
            'Input {} greater than maxval {}'.format(integer, self.maxval)

        output = [0 for _ in range(self.N)]

        # Compute the indices where the 1s start and end
        centerbin = integer - self.minval + self.resolution//2
        centerbin = centerbin // self.resolution + self.padding
        minbin = round(centerbin - self.halfwidth)
        maxbin = round(minbin + 2 * self.halfwidth)

        assert minbin >= 0
        assert maxbin < self.N
        for i in range(minbin, maxbin + 1):
            output[i] = 1

        return output



##################
# Helper functions
##################

def get_bit_at(sparse, x, y):
    """
    Returns the bit at position (x,y) in the array.

    @param sparse: an unrolled SDR
    @param x, y: (x, y) pos in the rolled up SDR
    @returns: a bit
    """
    side = int(math.sqrt(len(sparse)))
    return sparse[side * x + y]



##############
# HTM classes
##############

class Synapse():
    """
    Segment of a synapse that connects to another cell.  Keeps track of its
    presynpatic cell's position, permanence value and whether or not the
    Synapse is active.
    """

    THRESHOLD = 0.6
    """If self.permanence > threshold, then this synapse is valid."""

    PERMDELTA = 0.01
    """Amount by which permanence value changes when learning."""


    def __init__(self, inputpos, perm=0.0):
        self.permanence = perm
        self.inputpos = inputpos

    def is_valid(self):
        """Whether or not this synapse is valid."""
        return self.permanence > Synapse.THRESHOLD

    def increase_perm(self):
        """Increase permanence value."""
        self.permanence = min(self.permanence + Synapse.PERMDELTA, 1.0)

    def decrease_perm(self):
        """Decrease permanence value."""
        self.permanence = max(0.0, self.permanence - Synapse.PERMDELTA)


class Dendrite():
    """
    Stems from cells (sometimes shared by a whole column).  Keeps track of
    current Synapses and of cells that might form synapses while learning.
    """

    NPOTENTIAL = 20
    """The number of potential synapses for this Dendrite."""


    def __init__(self):
        def randpos():
            """Return a random (x, y) where x, y are valid column indices."""
            maxval = round(math.sqrt(Region.NCOLUMNS)) - 1
            return (random.randint(0, maxval), random.randint(0, maxval))

        potential = (randpos() for _ in range(Dendrite.NPOTENTIAL))
        self.synapses = [Synapse(pos, random.random()) for pos in potential]


class Cell():
    """
    A single computational unit.  Keeps track of its dendrite segments,
    which determine where a Cell gets its input from, and its state.
    """

    def __init__(self, proximal):
        """
        A Column must send its shared Dendrite to each of its cells.

        @param proximal: a Dendrite, shared among Cells in the same Column.
        """

        self.distal = []
        self.proximal = proximal

        # This is not a property because Columns tell their Cells when to
        # activate
        self.state = Column.INACTIVE


class Column():
    """
    An array of Cells.  All the cells in one column share a single proximal
    Dendrite segment.
    """

    NCELLS = 4
    """Number of cells in this Column."""

    INACTIVE = 0
    ACTIVE = 1
    PREDICTIVE = 2

    def __init__(self):
        self.proximal = Dendrite()
        self.cells = [Cell(self.proximal) for _ in range(self.NCELLS)]
        self.boost = 1.0

    def num_active_cells(self):
        """Return the current number of active cells."""
        return sum([cell.state == Column.ACTIVE for cell in self.cells])

    @property
    def state(self):
        """A Column is active whenever at least one of its Cells is active."""
        if Column.ACTIVE in [c.state for c in self.cells]:
            return Column.ACTIVE
        else:
            return Column.INACTIVE

    @state.setter
    def state(self, state):
        """Set the state of each Cell in this Column."""
        if state == Column.ACTIVE:
            predictive = {c for c in self.cells if c.state == Column.PREDICTIVE}

            if predictive:
                for cell in predictive:
                    cell.state = Column.ACTIVE
                for cell in set(self.cells) - predictive:
                    cell.state = Column.INACTIVE

            else:
                for cell in self.cells:
                    cell.state = Column.ACTIVE

        elif state == Column.INACTIVE:
            for cell in self.cells:
                cell.state = Column.INACTIVE



class Region():
    """A Region is an array of Columns."""

    NCOLUMNS = 400 # needs to be a perfect square for now!
    """Number of columns in a Region."""

    density = 2./100
    """The sparsity of the active columns."""


    def __init__(self):
        self.columns = [Column() for _ in range(Region.NCOLUMNS)]
        """A list of Columns. Call Region.column_matrix for a 2D array."""

    @property
    def column_matrix(self):
        """Return the Columns as a 2D array."""
        side = int(math.sqrt(Region.NCOLUMNS))

        return [[self.columns[row * side + col] for col in range(side)]
                for row in range(side)]


    def process(self, sparse):
        """Process some input data (runs spatial and temporal pooler).

        sparse: SDR of the input data.
        """
        ### Spatial pooler
        ### input: the Sparse Ddistributed Representation of inpute data
        ### output: the list of active columns, activecols

        # Determine how many valid synapses on each column are connected to
        # active input bits
        numvalid = {}

        for col in self.columns:
            numvalid[col] = len([syn for syn in col.proximal.synapses
                                 if syn.is_valid() and
                                 get_bit_at(sparse, *syn.inputpos) == 1])

            # The number of active synapses is multiplied by a "boosting"
            # factor which is dynamically determined by how often a column
            # is active relative to its neighbors.
            numvalid[col] *= col.boost

        # The columns with the highest activations after boosting disable
        # all but a fixed percentage of the columns within an inhibition
        # radius.  The inhibition radius is determined by the spread of
        # input bits.  There is now a sparse set of active columns.

        # inhibition code here


        ### CHANGE ME
        srtd = sorted(self.columns, key=numvalid.get)
        index = Encoder.N - round(Encoder.N * self.density)

        # Here, we compute which columns are the "winners" (soon to be
        # activated).  We do not set col.state yet, because that will force
        # every cell in each column to update its state too.  We need to do
        # that after the permanence values have been recomputed.
        inactivecols = srtd[:index]
        activecols = srtd[index:]
        ### CHANGE ME


        # For each active column, adjust the permanence of all the
        # potential synapses.  These changes may validate inactive
        # synapses, and vice-versa.
        for col in activecols:
            activesyn = {syn for syn in col.proximal.synapses
                         if get_bit_at(sparse, *syn.inputpos) == 1}

            for synapse in activesyn:
                synapse.increase_perm()

            inactivesyn = set(col.proximal.synapses) - activesyn
            for synapse in inactivesyn:
                synapse.decrease_perm()


        ### Temporal pooler
        ### input: the list of active columns, activecols
        ### output: None

        # Once the permanence values have been updated, we change the state
        # of each column, which in turn forces every cell to update its
        # state.
        for col in activecols:
            col.state = Column.ACTIVE
        for col in inactivecols:
            col.state = Column.INACTIVE


        # Count how many active synapses are connected to active cells.  If
        # there are enoguh of them, mark the segment as active.

        # For every dendrite segment (proximal AND distal?) on every cell
        # in the region, count how many established (active?) synapses are
        # connected to active cells.  If the number exceeds a threshold,
        # that dendrite segment is marked as active.  Cells with active
        # dendrite segments are put in the predictive state unless they are
        # already active.  Cells with no active dendrites and not active
        # due to bottom-up input become or remain inactive.  The collection
        # of cells now in the predictive state is the prediction of the
        # region.

        return


    def __str__(self):
        """Return a string representation of all columns."""
        rules = {4: "■", 3: "X", 2: "x", 1: "o", 0: "·"}
        strings = [''.join([rules[col.num_active_cells()] for col in row])
                   for row in self.column_matrix]

        return '\n'.join(strings)


################################################
# MAIN
################################################

def main():
    """Run one HTM region learning over randomly-generated data."""
    random.seed()

    encoder = Encoder()
    region = Region()
    while True:
        value = random.randint(encoder.minval, encoder.maxval - 1)
        data = encoder.encode(value)

        region.process(data)
        print(region)

        time.sleep(1)
        os.system("clear")


if __name__ == "__main__":
    main()
