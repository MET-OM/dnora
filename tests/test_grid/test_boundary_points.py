import unittest
import dnora
from dnora.grid import Grid
from dnora.grid.mask import Edges
from dnora.utils.grid import (
    identify_boundary_edges,
    create_ordered_boundary_list,
    get_coords_for_boundary_edges,
)
from dnora.executer.inputfile.swan_functions import create_swan_segment_coords
import numpy as np
from copy import copy


class IdentifyBnd(unittest.TestCase):
    def test_all_edges(self):
        grid = Grid(lon=(0, 1), lat=(0, 1))
        grid.set_spacing(nx=4, ny=4)

        grid.set_boundary_points(Edges(edges=["N", "W", "S", "E"]))

        self.assertEqual(
            set(identify_boundary_edges(grid.boundary_mask())),
            set(["N", "W", "S", "E"]),
        )

    def test_one_missing(self):
        grid = Grid(lon=(0, 1), lat=(0, 1))
        grid.set_spacing(nx=4, ny=5)

        full_list = ["N", "W", "S", "E"]

        for not_this in full_list:
            incomplete_list = copy(full_list)
            incomplete_list.remove(not_this)

            grid.set_boundary_points(Edges(edges=incomplete_list))
            self.assertEqual(
                set(identify_boundary_edges(grid.boundary_mask())), set(incomplete_list)
            )


class OrderBnd(unittest.TestCase):
    def test_order_full_list(self):
        edges = ["S", "E"]
        self.assertEqual(create_ordered_boundary_list(edges), ["E", "S"])

        edges = ["N", "E"]
        self.assertEqual(create_ordered_boundary_list(edges), ["N", "E"])

        edges = ["W", "E", "S"]
        self.assertEqual(create_ordered_boundary_list(edges), ["E", "S", "W"])

        edges = ["W", "E", "N"]
        self.assertEqual(create_ordered_boundary_list(edges), ["W", "N", "E"])

        edges = ["W", "E"]
        self.assertEqual(create_ordered_boundary_list(edges), [])


class BoundaryCoordinates(unittest.TestCase):
    def test_single_edge(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        edges = ["N"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([3, 3])), True)

        edges = ["S"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([1, 0])), True)
        self.assertEqual(np.array_equal(lats, np.array([2, 2])), True)

        edges = ["W"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 0])), True)
        self.assertEqual(np.array_equal(lats, np.array([2, 3])), True)

        edges = ["E"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([1, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([3, 2])), True)

    def test_two_edge(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        edges = ["N", "E"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 1, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([3, 3, 2])), True)

        edges = ["S", "W"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([1, 0, 0])), True)
        self.assertEqual(np.array_equal(lats, np.array([2, 2, 3])), True)

        edges = ["W", "N"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 0, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([2, 3, 3])), True)

        edges = ["E", "S"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([1, 1, 0])), True)
        self.assertEqual(np.array_equal(lats, np.array([3, 2, 2])), True)

    def test_three_edge(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        edges = ["N", "E", "S"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 1, 1, 0])), True)
        self.assertEqual(np.array_equal(lats, np.array([3, 3, 2, 2])), True)

        edges = ["S", "W", "N"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([1, 0, 0, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([2, 2, 3, 3])), True)

        edges = ["W", "N", "E"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 0, 1, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([2, 3, 3, 2])), True)

        edges = ["E", "S", "W"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([1, 1, 0, 0])), True)
        self.assertEqual(np.array_equal(lats, np.array([3, 2, 2, 3])), True)

    def test_four_edge(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        edges = ["N", "E", "S", "W"]
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 1, 1, 0, 0])), True)
        self.assertEqual(np.array_equal(lats, np.array([3, 3, 2, 2, 3])), True)

    def test_invalid_edge(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        edges = []
        lons, lats = get_coords_for_boundary_edges(
            edges, grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([])), True)
        self.assertEqual(np.array_equal(lats, np.array([])), True)


class FullPipeline(unittest.TestCase):
    def test_single_edge(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        grid.set_boundary_points(Edges(edges=["N"]))
        lons, lats = create_swan_segment_coords(
            grid.boundary_mask(), grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([3, 3])), True)

    def test_two_edges(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        grid.set_boundary_points(Edges(edges=["N", "W"]))
        lons, lats = create_swan_segment_coords(
            grid.boundary_mask(), grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 0, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([2, 3, 3])), True)

    def test_three_edges(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        grid.set_boundary_points(Edges(edges=["N", "W", "E"]))
        lons, lats = create_swan_segment_coords(
            grid.boundary_mask(), grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([0, 0, 1, 1])), True)
        self.assertEqual(np.array_equal(lats, np.array([2, 3, 3, 2])), True)

    def test_no_edges(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)

        lons, lats = create_swan_segment_coords(
            grid.boundary_mask(empty=True), grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([])), True)
        self.assertEqual(np.array_equal(lats, np.array([])), True)

    def test_invalid_edges(self):
        grid = Grid(lon=(0, 1), lat=(2, 3))
        grid.set_spacing(nx=4, ny=5)
        grid.set_boundary_points(Edges(edges=["W", "E"]))

        lons, lats = create_swan_segment_coords(
            grid.boundary_mask(Empty=True), grid.edges("lon"), grid.edges("lat")
        )
        self.assertEqual(np.array_equal(lons, np.array([])), True)
        self.assertEqual(np.array_equal(lats, np.array([])), True)


if __name__ == "__main__":
    unittest.main()
