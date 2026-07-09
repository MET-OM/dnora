import numpy as np
def reorganize_triangulation_with_edges_first(tri, lon, lat, edge_nodes):
    other_nodes = np.setdiff1d(np.unique(tri), edge_nodes)
    # Re-organize nodes so that boundary nodes are first
    new_tri = np.copy(tri)
    new_lon = np.copy(lon)
    new_lat = np.copy(lat)
    for ind, node in enumerate(edge_nodes):
        new_tri[tri == node] = ind
        new_lon[ind] = lon[node]
        new_lat[ind] = lat[node]
    for ind, node in enumerate(other_nodes):
        new_tri[tri == node] = ind + len(edge_nodes)
        new_lon[ind + len(edge_nodes)] = lon[node]
        new_lat[ind + len(edge_nodes)] = lat[node]
    return new_tri, new_lon, new_lat



def reindex_triangulation_when_nodes_removed(tri, kept_nodes):
    """If nodes are removed the triangulation needs to be reindexed. 
    
    kept_nodes is a list of nodes that remain, with index using the old triangulation"""
    N = tri.shape[0]
    
    keep = np.zeros(N, dtype=bool)
    keep[kept_nodes] = True

    # Index mapping
    old_to_new = -np.ones(N, dtype=int)
    old_to_new[keep] = np.arange(keep.sum())

    tri_reindexed = old_to_new[tri]

    return tri_reindexed 

def force_nodes_as_first_in_triangulation(lon, lat, tri, bnd_nodes, nodes_to_be_first):
    """Reorders the index values in old_inds to new_inds and keeps consistens lon, lat and trianulation values"""
    

    order = np.ones(len(lon)).astype(int)
    order[0:len(nodes_to_be_first)] = nodes_to_be_first
    order[len(nodes_to_be_first):] = np.array(list(set(range(len(lon))).difference(set(nodes_to_be_first))))

    lon, lat = lon[order], lat[order]
    

    N = order.size
    inv = np.empty(N, dtype=int)
    inv[order] = np.arange(N)        # inverse permutation: old -> new

    tri = inv[tri]
    bnd_nodes = inv[bnd_nodes]   

    return lon, lat, tri, bnd_nodes


def get_boundary_edges(tri):
    """Finds boundary edges in triangulation"""
    # Extract all edges from triangles
    edges = np.vstack([
        tri[:, [0, 1]],  # First edge of each triangle
        tri[:, [1, 2]],  # Second edge of each triangle
        tri[:, [2, 0]]   # Third edge of each triangle
    ])

    # Sort the edges to ensure consistent ordering (e.g., [1, 2] is the same as [2, 1])
    edges = np.sort(edges, axis=1)

    # Count the occurrences of each edge
    from collections import Counter

    edge_counts = Counter(map(tuple, edges))

    # Boundary edges are those that occur only once
    boundary_edges = [edge for edge, count in edge_counts.items() if count == 1]
    return boundary_edges

def get_boundary_nodes(edges):
    """Gets all the nodes that form the edges in a consequtive order"""
    ind = 0
    shorelines = [] # Shorelines of boundary and all islands
    bnd_nodes = []
    edge_array = np.array(edges)
    n0, n1 = edges[ind]
    bnd_nodes.append(n0)
    bnd_nodes.append(n1)
    edge_array[ind,:] = -1
    next_node = n1
    # print(f"Node {n0} connects to {n1}")
    # print(f"Next node: {next_node}")
    while np.max(edge_array) > -1:
        try:
            ind = np.argwhere(edge_array == next_node)[:,0][0]
        except IndexError:
            shorelines.append(bnd_nodes)
            bnd_nodes = []
            ind = np.where(np.max(edge_array,axis=1)>-1)[0][0]
            
        n0, n1 = edges[ind]
        
        if not bnd_nodes or n0 == bnd_nodes[-1]:
            next_node = n1
        else:
            next_node = n0
        #print(f"Node {n0} connects to {n1}")
        #print(f"Next node: {next_node}")
        bnd_nodes.append(next_node)
        edge_array[ind,:] = -1

    # Assume that the longest consequitve shoreline is boundary/land and rest are islands
    shore_ind = np.argmax([len(bnd) for bnd in shorelines])
    return shorelines[shore_ind][:-1]


def find_nodes_with_three_neighbours(tri) -> list[int]:
    """
    Returns nodes that have exactly three distinct neighboring nodes.
    tri: array of shape (T, 3) with zero-based node indices.
    """
    tri = np.asarray(tri)
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]

    # Build undirected edges from each triangle
    edges = np.vstack([
        np.stack([a, b], axis=1),
        np.stack([b, c], axis=1),
        np.stack([c, a], axis=1),
    ])

    # Make undirected (sorted endpoints) and drop duplicates
    edges = np.sort(edges, axis=1)
    edges = np.unique(edges, axis=0)

    # Degree: each edge contributes one to both endpoints
    n_nodes = tri.max() + 1
    deg = np.bincount(
        np.concatenate([edges[:, 0], edges[:, 1]]),
        minlength=n_nodes
    )
    return np.where(deg == 3)[0].tolist()


def triangle_angles_xy(x, y):
    """
    x, y: arrays of shape (N, 3), where row i is [Ax, Bx, Cx] and [Ay, By, Cy].
    Returns angles of shape (N, 3) in order (angle at A, at B, at C).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.shape != y.shape or x.ndim != 2 or x.shape[1] != 3:
        raise ValueError("x and y must be arrays of shape (N, 3)")

    # Vectors at each vertex
    ABx = x[:, 1] - x[:, 0]; ABy = y[:, 1] - y[:, 0]
    ACx = x[:, 2] - x[:, 0]; ACy = y[:, 2] - y[:, 0]

    BAx = x[:, 0] - x[:, 1]; BAy = y[:, 0] - y[:, 1]
    BCx = x[:, 2] - x[:, 1]; BCy = y[:, 2] - y[:, 1]

    CAx = x[:, 0] - x[:, 2]; CAy = y[:, 0] - y[:, 2]
    CBx = x[:, 1] - x[:, 2]; CBy = y[:, 1] - y[:, 2]

    def ang(ux, uy, vx, vy):
        dot = ux * vx + uy * vy
        cross = ux * vy - uy * vx  # 2D cross-product z-component
        a = np.arctan2(np.abs(cross), dot)
        return np.degrees(a)

    A = ang(ABx, ABy, ACx, ACy)  # angle at A
    B = ang(BAx, BAy, BCx, BCy)  # angle at B
    C = ang(CAx, CAy, CBx, CBy)  # angle at C

    return np.stack([A, B, C], axis=1)


def outlier_angle_indices_median(tri_angles):
    """
    tri_angles: array (N, 3)
    Returns rows, cols as above.
    """
    A = np.asarray(tri_angles, dtype=float)

    med = np.median(A, axis=1, keepdims=True)
    dev = np.abs(A - med)
    cols = np.argmax(dev, axis=1)
    return cols

def length_of_path(x, y):
    """Calculates the cartesian path length of a list of points"""
    dx = np.diff(x)
    dy = np.diff(y)
    return np.hypot(dx, dy).sum()
