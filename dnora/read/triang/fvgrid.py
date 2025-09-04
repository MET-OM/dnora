import numpy as np
#from matplotlib.tri.triangulation import Triangulation
import matplotlib.tri as triangulation # CubicTriInterpolator, triangulation

# This function is taken dircetly from the PyFVCOM package
# https://github.com/pwcazenave/pyfvcom
def read_sms_mesh(mesh, nodestrings=False):
    """
    Reads in the SMS unstructured grid format. Also creates IDs for output to
    MIKE unstructured grid format.

    Parameters
    ----------
    mesh : str
        Full path to an SMS unstructured grid (.2dm) file.
    nodestrings : bool, optional
        Set to True to return the IDs of the node strings as a dictionary.

    Returns
    -------
    triangle : ndarray
        Integer array of shape (nele, 3). Each triangle is composed of
        three points and this contains the three node numbers (stored in
        nodes) which refer to the coordinates in X and Y (see below). Values
        are python-indexed.
    nodes : ndarray
        Integer number assigned to each node.
    X, Y, Z : ndarray
        Coordinates of each grid node and any associated Z value.
    types : ndarray
        Classification for each node string based on the number of node
        strings + 2. This is mainly for use if converting from SMS .2dm
        grid format to DHI MIKE21 .mesh format since the latter requires
        unique IDs for each boundary (with 0 and 1 reserved for land and
        sea nodes).
    nodestrings : list, optional (see nodestrings above)
        Optional list of lists containing the node IDs (python-indexed) of the
        node strings in the SMS grid.

    """

    fileRead = open(mesh, 'r')
    lines = fileRead.readlines()
    fileRead.close()

    triangles = []
    nodes = []
    types = []
    nodeStrings = []
    nstring = []
    x = []
    y = []
    z = []

    # MIKE unstructured grids allocate their boundaries with a type ID flag.
    # Although this function is not necessarily always the precursor to writing
    # a MIKE unstructured grid, we can create IDs based on the number of node
    # strings in the SMS grid. MIKE starts counting open boundaries from 2 (1
    # and 0 are land and sea nodes, respectively).
    typeCount = 2
    positive_last = True

    for line in lines:
        line = line.strip()
        if line.startswith('E3T'):
            ttt = line.split()
            t1 = int(ttt[2]) - 1
            t2 = int(ttt[3]) - 1
            t3 = int(ttt[4]) - 1
            triangles.append([t1, t2, t3])
        elif line.startswith('ND '):
            xy = line.split()
            x.append(float(xy[2]))
            y.append(float(xy[3]))
            z.append(float(xy[4]))
            nodes.append(int(xy[1]))
            # Although MIKE keeps zero and one reserved for normal nodes and
            # land nodes, SMS doesn't. This means it's not straightforward
            # to determine this information from the SMS file alone. It would
            # require finding nodes which are edge nodes and assigning their
            # ID to one. All other nodes would be zero until they were
            # overwritten when examining the node strings below.
            types.append(0)
        elif line.startswith('NS '):
            allTypes = line.split(' ')

            for nodeID in allTypes[2:]:
                types[np.abs(int(nodeID)) - 1] = typeCount
                if int(nodeID) > 0:
                    if positive_last:
                        # The first number after a negative number counts the nodestring id, and
                        # should not be used as a node id
                        nstring.append(int(nodeID) - 1)

                    positive_last = True
                else:
                    nstring.append(np.abs(int(nodeID)) - 1)
                    nodeStrings.append(nstring)
                    nstring = []
                    positive_last = False

                # Count the number of node strings, and output that to types.
                # Nodes in the node strings are stored in nodeStrings.
                if int(nodeID) < 0:
                    typeCount += 1

    # Convert to numpy arrays.
    triangle = np.asarray(triangles)
    nodes = np.asarray(nodes)
    types = np.asarray(types)
    X = np.asarray(x)
    Y = np.asarray(y)
    Z = np.asarray(z)
    if nodestrings:
        return triangle, nodes, X, Y, Z, types, nodeStrings
    else:
        return triangle, nodes, X, Y, Z, types


def sigma_tanh(nlev, dl, du):
    '''Generate a tanh sigma coordinate distribution'''

    dist = np.zeros(nlev)

    for k in np.arange(nlev-1):
        x1 = dl + du
        x1 = x1 * (nlev - 2 - k) / (nlev - 1)
        x1 = x1 - dl
        x1 = np.tanh(x1)
        x2 = np.tanh(dl)
        x3 = x2 + np.tanh(du)
        dist[k+1] = (x1 + x2) / x3 - 1.0

    return dist

def find_connected_nodes(n, triangles):
    """
    Return the IDs of the nodes surrounding node number `n'.

    Parameters
    ----------
    n : int
        Node ID around which to find the connected nodes.
    triangles : ndarray
        Triangulation matrix to find the connected nodes. Shape is [nele, 3].

    Returns
    -------
    surroundingidx : ndarray
        Indices of the surrounding nodes.

    See Also
    --------
    PyFVCOM.grid.find_connected_elements().

    """

    eidx = np.max((np.abs(triangles - n) == 0), axis=1)
    surroundingidx = np.unique(triangles[eidx][triangles[eidx] != n])

    return surroundingidx

def trigradient(x, y, z, t=None):
    """
    Returns the gradient of `z' defined on the irregular mesh with Delaunay
    triangulation `t'. `dx' corresponds to the partial derivative dZ/dX,
    and `dy' corresponds to the partial derivative dZ/dY.

    Parameters
    ----------
    x, y, z : array_like
        Horizontal (`x' and `y') positions and vertical position (`z').
    t : array_like, optional
        Connectivity table for the grid. If omitted, one will be calculated
        automatically.

    Returns
    -------
    dx, dy : ndarray
        `dx' corresponds to the partial derivative dZ/dX, and `dy'
        corresponds to the partial derivative dZ/dY.

    Example
    -------

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from PyFVCOM.grid import trigradient
    >>> from matplotlib.tri.triangulation import Triangulation
    >>> x, y = np.meshgrid(np.arange(-2, 2, 0.1), np.arange(-2, 2, 0.1))
    >>> x[1:2:-1, :] = x[1:2:-1, :] + 0.1 / 2
    >>> tt = Triangulation(x.ravel(), y.ravel())
    >>> z = x * np.exp(-x**2 - y**2)
    >>> dx, dy = trigradient(x.ravel(), y.ravel(), z.ravel())
    >>> dzdx = (1 / x - 2 * x) * z
    >>> dzdy = -2 * y * z
    >>> plt.figure(1)
    >>> plt.quiver(x.ravel(), y.ravel(), dzdx.ravel(), dzdy.ravel(),
    >>>            color='r', label='Exact')
    >>> plt.quiver(x.ravel(), y.ravel(), dx, dy,
    >>>            color='k', label='trigradient', alpha=0.5)
    >>> tp = plt.tripcolor(x.ravel(), y.ravel(), tt.triangles, z.ravel(),
    >>>                    zorder=0)
    >>> plt.colorbar(tp)
    >>> plt.legend()

    Notes
    -----
    Lifted from:
        http://matplotlib.org/examples/pylab_examples/trigradient_demo.html

    """

    if np.any(t):
        tt = triangulation.Triangulation(x.ravel(), y.ravel())
    else:
        tt = triangulation.Triangulation(x.ravel(), y.ravel(), t)

    tci = triangulation.CubicTriInterpolator(tt, z.ravel())
    # Gradient requested here at the mesh nodes but could be anywhere else:
    dx, dy = tci.gradient(tt.x, tt.y)

    return dx, dy

def get_attached_unique_nodes(this_node, trinodes):
    """
    Find the nodes on the boundary connected to `this_node'.

    Parameters
    ----------
    this_node : int
        Node ID.
    trinodes : ndarray
        Triangulation table for an unstructured grid.

    Returns
    -------
    connected_nodes : ndarray
        IDs of the nodes connected to `this_node' on the boundary. If `this_node' is not on the boundary,
        `connected_nodes' is empty.

    """

    all_trinodes = trinodes[(trinodes[:, 0] == this_node) | (trinodes[:, 1] == this_node) | (trinodes[:, 2] == this_node), :]
    u, c = np.unique(all_trinodes, return_counts=True)

    return u[c == 1]


def grid_metrics(tri, noisy=False):
    """
    Calculate unstructured grid metrics (most of FVCOM's tge.F).

    Parameters
    ----------
    tri : ndarray
        Triangulation table for the grid.
    noisy : bool
        Set to True to enable verbose output (default = False)

    Returns
    -------
    ntve : ndarray
        The number of neighboring elements of each grid node
    nbve : ndarray
        nbve(i,1->ntve(i)) = ntve elements containing node i
    nbe : ndarray
        Indices of tri for the elements connected to each element in the domain. To visualise:
            plt.plot(x[tri[1000, :], y[tri[1000, :], 'ro')
            plt.plot(x[tri[nbe[1000], :]] and y[tri[nbe[1000], :]], 'k.')
        plots the 999th element nodes with the nodes of the surrounding elements too.
    isbce : ndarray
        Flag if element is on the boundary (True = yes, False = no)
    isonb : ndarray
        Flag if node is on the boundary (True = yes, False = no)

    Notes
    -----
    This is more or less a direct translation from FORTRAN (FVCOM's tge.F).

    """

    m = len(np.unique(tri.ravel()))

    # Allocate all our arrays. Use masked by default arrays so we only use valid indices.
    isonb = np.zeros(m).astype(bool)
    ntve = np.zeros(m, dtype=int)
    nbe = np.ma.array(np.zeros(tri.shape, dtype=int), mask=True)
    nbve = np.ma.array(np.zeros((m, 10), dtype=int), mask=True)
    # Number of elements connected to each node (ntve) and the IDs of the elements connected to each node (nbve).
    if noisy:
        print('Counting neighbouring nodes and elements')
    for i, (n1, n2, n3) in enumerate(tri):
        nbve[tri[i, 0], ntve[n1]] = i
        nbve[tri[i, 1], ntve[n2]] = i
        nbve[tri[i, 2], ntve[n3]] = i
        # Only increment the counters afterwards as Python indexes from 0.
        ntve[n1] += 1
        ntve[n2] += 1
        ntve[n3] += 1

    if noisy:
        print('Getting neighbouring elements for each element')

    # Get the element IDs connected to each element.
    for i, (n1, n2, n3) in enumerate(tri):
        for j1 in range(ntve[n1]):
            for j2 in range(ntve[n2]):
                if nbve[n1, j1] == nbve[n2, j2] and nbve[n1, j1] != i:
                    nbe[i, 2] = nbve[n1, j1]
        for j2 in range(ntve[n2]):
            for j3 in range(ntve[n3]):
                if nbve[n2, j2] == nbve[n3, j3] and nbve[n2, j2] != i:
                    nbe[i, 0] = nbve[n2, j2]
        for j1 in range(ntve[n1]):
            for j3 in range(ntve[n3]):
                if nbve[n1, j1] == nbve[n3, j3] and nbve[n1, j1] != i:
                    nbe[i, 1] = nbve[n3, j3]

    if noisy:
        print('Getting boundary element IDs')
    isbce = np.max(nbe.mask, axis=1)

    if noisy:
        print('Getting boundary node IDs')

    # Get the boundary node IDs. Holy nested list comprehensions, Batman!
    boundary_element_node_ids = np.unique(tri[isbce, :]).ravel()
    boundary_nodes = []
    for i in boundary_element_node_ids:
        current_nodes = get_attached_unique_nodes(i, tri)
        if np.any(current_nodes):
            boundary_nodes += current_nodes.tolist()
    boundary_nodes = np.unique(boundary_nodes)
    # Make a boolean of that.
    isonb[boundary_nodes] = True

    return ntve, nbve, nbe, isbce, isonb

def smoothfield(fieldin, tri, SmoothPts, Niter = 1, SmoothFactor = 1.0):
    ''' '''
    #if SmoothPts.all(0):
    #    SmoothPts = range(0,len(fieldin))
    for k in range(1,Niter+1):
        field = np.copy(fieldin)
        for r in SmoothPts:
            nnodes = find_connected_nodes(r, tri)
            fave = np.average(np.append(fieldin[nnodes],fieldin[r]))
            field[r] = SmoothFactor*fave + (1 - SmoothFactor) * fieldin[r]

    return field
