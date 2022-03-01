# coding: utf-8
"""
"""


import genice2.lattices
from genice2.cell import cellvectors
from logging import getLogger
from math import sin, pi, cos, floor, gcd
import numpy as np


def usage():
    logger = getLogger()
    logger.info(__doc__)


desc = {
    "ref": {
    },
    "usage": usage(),
    "brief": ""
}

# ベクトルは座標の対で表現する。inflation中は、細分しかしないので、
# 頂点の重なりは気にしない。さいごにえられたタイルを整数化する。
# 必要ならそこからさらに細分を行う。
# 手でタイルを貼るのと同じ感覚で、辺の横にこのタイル、という指定ができると良い。
# 新たなタイルをくっつける時は辺の向きに対し左につけるものとする。

# edgeをclassにする。(あとで拡張しやすいように)

# edgeを拡張する。houseとtriangleで埋める。

# houseとtriangleを、パリティ付き菱形と六角に変換する。まだ周期境界は考えない
# parityはaugmentの辺の性質である。菱形の場合はパリティが変更される。
# 辺の辞書をつくり、辺のオーナーをさがし、辞書をてがかりに再帰的にパリティを決定する。

# 最後に、decorationを行う。これは簡単。

from math import *
theta = pi / 5.0  # 36 degree
thetah = pi / 10.0  # 18 degree
omega = (sqrt(5.0) + 1.0) / 2.0
omega1 = omega + 1.0
# longer edge of the triangle (shorter is the unit)
sin36 = sin(theta)
cos36 = cos(theta)
sin18 = sin(thetah)
cos18 = cos(thetah)
L = 2.0 * sin36


class Tile():
    """
    A tile is a polygon defined by a list of edges.
    """

    def __init__(self, edges):
        # edges is a tuple of edges
        self.edges = edges
        self.parity = None

    def __getitem__(self, n):
        return self.edges[n]

    def __len__(self):
        return len(self.edges)


class Edge():
    """
    Turtle graphics
    """

    def __init__(self, s, d):
        self.start = s
        self.delta = d
        self.augment = None

    def proceed(self):
        self.start = (self.start[0] + self.delta[0],
                      self.start[1] + self.delta[1])

    def offset(self, delta):
        self.start = (self.start[0] + delta[0],
                      self.start[1] + delta[1])

    def rotate(self, angle):
        # rotate around start
        angle *= pi / 180
        sn = sin(angle)
        cs = cos(angle)
        self.delta = (self.delta[0] * cs - self.delta[1] * sn,
                      self.delta[0] * sn + self.delta[1] * cs)

    def copy(self):
        return Edge(self.start, self.delta)

    def triangle(self):
        e0 = self.copy()
        e1 = self.copy()
        e1.proceed()
        e1.delta = (e1.delta[0] / L, e1.delta[1] / L)
        e1.rotate(126.)
        e2 = e1.copy()
        e2.proceed()
        e2.rotate(108.)
        self.augment = Tile((e0, e1, e2))

    def house(self):
        e0 = self.copy()
        e1 = self.copy()
        e1.proceed()
        e1.delta = (e1.delta[0] / L, e1.delta[1] / L)
        e1.rotate(90.)
        e2 = e1.copy()
        e2.proceed()
        e2.rotate(36.)
        e3 = e2.copy()
        e3.proceed()
        e3.rotate(108.)
        e4 = e3.copy()
        e4.proceed()
        e4.rotate(36.)
        self.augment = Tile((e0, e1, e2, e3, e4))


def invert(self):
    start = (self.start[0] + self.delta[0],
             self.start[1] + self.delta[1])
    delta = (-self.delta[0], -self.delta[1])
    return Edge(start, delta)


def diff(fro, to):
    start = (fro.start[0] + fro.delta[0],
             fro.start[1] + fro.delta[1])
    delta = (to.start[0] + to.delta[0] - start[0],
             to.start[1] + to.delta[1] - start[1])
    return Edge(start, delta)
# 位置を指定されると、周の辺を返す。内部はどうする?


def decagon_at(pos, edgelength, dir):
    fro = Edge((0, 0), (L * edgelength / 2.0 / sin(pi / 10), 0))
    to = fro.copy()
    fro.rotate(+18. + dir * 36.)
    to.rotate(-18. + dir * 36.)
    edge = diff(fro, to)
    edge.offset(pos)
    return polygon_along(edge, 10)


def pentagon_at(pos, edgelength, dir):
    fro = Edge((0, 0), (edgelength, 0))
    to = fro.copy()
    fro.rotate(+72. + dir * 36.)
    to.rotate(dir * 36.)
    edge = diff(fro, to)
    edge.offset(pos)
    return polygon_along(edge, 5)


#edge is (startpos, delta)
def polygon_along(edge, n):
    e = invert(edge)
    # start downward, clockwise
    edges = []
    for i in range(n):
        f = e.copy()
        if n == 5 or (n == 10 and i % 2 == 0):
            f.triangle()
        elif n == 10 and i % 2 == 1:
            f.house()
        edges.append(f)
        e.proceed()
        e.rotate(360.0 / n)
    return Tile(edges)

# "which" edgeof the boomelang will be shared.


def boomelang_along(edge, which):
    turn = (2, 2, 1, 1, 1, 1, 2, 2, -2)
    which -= 9
    e = invert(edge)
    # start downward, clockwise
    edges = []
    for i in range(which, 9):
        if i >= 0:
            edges.append(e.copy())
        e.proceed()
        e.rotate(turn[i] * 36.0)
    for i in (0, 1, 2, 4, 6, 7, 8):
        edges[i].triangle()
    for i in (3, 5):
        edges[i].house()
    return Tile(edges)


def drawpoly(polygon, label=""):
    # to check overlaps
    cx = cy = 0
    for edge in polygon.edges:
        cx += edge.start[0]
        cy += edge.start[1]
    cx /= len(polygon)
    cy /= len(polygon)
    beginpath(polygon[0].start[0], polygon[0].start[1])
    for i in range(1, len(polygon)):
        lineto(polygon[i].start[0], polygon[i].start[1])
    endpath()
    for edge in polygon:
        if edge.augment is not None:
            drawpoly(edge.augment)
    fill(0)
    text(label, cx, cy)


def D(pos, edgelength, dir):
    polygons = []
    decagon = decagon_at(pos, edgelength, dir)
    polygons.append(decagon)
    for i in range(10):
        edge = decagon[i - 1]
        pentagon = polygon_along(edge, 5)
        polygons.append(pentagon)

        edge = pentagon[3]
        boomelang = boomelang_along(edge, 6)
        polygons.append(boomelang)

        edge = boomelang[2]
        pentagon = polygon_along(edge, 5)
        polygons.append(pentagon)

        edge = boomelang[3]
        pentagon = polygon_along(edge, 5)
        polygons.append(pentagon)
    return polygons


def P(pos, edgelength, dir):
    polygons = []
    pentagon = pentagon_at(pos, edgelength, dir)
    polygons.append(pentagon)
    for i in range(5):
        edge = pentagon[i - 1]
        boomelang = boomelang_along(edge, 6)
        polygons.append(boomelang)

        edge = boomelang[0]
        penta2 = polygon_along(edge, 5)
        polygons.append(penta2)

        edge = boomelang[1]
        boom2 = boomelang_along(edge, 4)
        polygons.append(boom2)

        for j in (0, 2, 3, 6):
            edge = boom2[j]
            penta2 = polygon_along(edge, 5)
            polygons.append(penta2)
    return polygons


def X(pos, edgelength, dir):
    polygons = []
    r = L * edgelength
    dx, dy = -r * sin(dir * pi / 5), r * cos(dir * pi / 5)
    x, y = pos[0] - dx / 2, pos[1] - dy / 2
    pentagon = polygon_along(Edge((x, y), (dx, dy)), 5)
    polygons.append(pentagon)

    edge = pentagon[0]
    boom1 = boomelang_along(edge, 2)
    polygons.append(boom1)

    # right wing
    edge = boom1[3]
    penta = polygon_along(edge, 5)
    polygons.append(penta)

    edge = boom1[4]
    boom2 = boomelang_along(edge, 7)
    polygons.append(boom2)

    edge = boom2[1]
    boom3 = boomelang_along(edge, 3)
    polygons.append(boom3)

    edge = boom3[7]
    boom4 = boomelang_along(edge, 4)
    polygons.append(boom4)

    edge = boom4[1]
    deca2 = polygon_along(edge, 10)
    polygons.append(deca2)

    edge = boom4[7]
    deca3 = polygon_along(edge, 10)
    polygons.append(deca3)
    for i in (6, 8, -1):
        polygons.append(polygon_along(deca3[i], 5))
    edge = deca3[7]
    boom6 = boomelang_along(edge, 7)
    polygons.append(boom6)

    edge = boom3[1]
    deca4 = polygon_along(edge, 10)
    polygons.append(deca4)
    for i in (1, 2, 3):
        polygons.append(polygon_along(deca4[i], 5))

    # center decagon
    edge = boom1[1]
    deca1 = polygon_along(edge, 10)
    polygons.append(deca1)
    for i in (1, 2, 7, 8):
        polygons.append(polygon_along(deca1[i], 5))

    # left wing
    edge = boom1[7]
    boom5 = boomelang_along(edge, 4)
    polygons.append(boom5)

    edge = boom5[1]
    deca5 = polygon_along(edge, 10)
    polygons.append(deca5)

    edge = boom5[7]
    deca6 = polygon_along(edge, 10)
    polygons.append(deca6)
    for i in (8, -1):
        polygons.append(polygon_along(deca6[i], 5))

    edge = boom1[5]
    penta = polygon_along(edge, 5)
    polygons.append(penta)

    edge = penta[3]
    deca7 = polygon_along(edge, 10)
    polygons.append(deca7)
    for i in (5, 7, 8, -1):
        polygons.append(polygon_along(deca7[i], 5))
    edge = deca7[6]
    boom7 = boomelang_along(edge, 7)
    polygons.append(boom7)
    return polygons


def Y(pos, edgelength, dir):
    polygons = []
    r = edgelength * cos36
    x0, y0 = -r, 0
    dx, dy = L * edgelength * cos18, -L * edgelength * sin18
    cs = cos(dir * pi / 5)
    sn = sin(dir * pi / 5)
    dx, dy = dx * cs - dy * sn, dx * sn + dy * cs
    x0, y0 = x0 * cs - y0 * sn, x0 * sn + y0 * cs
    x, y = pos[0] + dx + x0, pos[1] + dy + y0
    boom0 = boomelang_along(Edge((x, y), (-dx, -dy)), 0)
    polygons.append(boom0)
    for i in (5, 6):
        polygons.append(polygon_along(boom0[i], 5))

    # right wing
    edge = boom0[7]
    deca1 = polygon_along(edge, 10)
    polygons.append(deca1)
    for i in (3, 4, 8, -1):
        polygons.append(polygon_along(deca1[i], 5))

    edge = boom0[1]
    boom1 = boomelang_along(edge, 4)
    polygons.append(boom1)

    edge = boom1[7]
    deca2 = polygon_along(edge, 10)
    polygons.append(deca2)
    for i in (7, 8, -1):
        polygons.append(polygon_along(deca2[i], 5))

    edge = boom1[1]
    boom2 = boomelang_along(edge, 4)
    polygons.append(boom2)

    edge = boom2[1]
    deca3 = polygon_along(edge, 10)
    polygons.append(deca3)
    penta = polygon_along(deca3[8], 5)
    polygons.append(penta)

    edge = penta[3]
    deca6 = polygon_along(edge, 10)
    polygons.append(deca6)
    for i in (6, 7, 8, -1):
        polygons.append(polygon_along(deca6[i], 5))

    edge = boom2[7]
    deca4 = polygon_along(edge, 10)
    polygons.append(deca4)
    polygons.append(polygon_along(deca4[-1], 5))
    polygons.append(polygon_along(deca4[4], 5))
    polygons.append(boomelang_along(deca4[3], 1))

    # left wing
    edge = boom0[4]
    boom3 = boomelang_along(edge, 1)
    polygons.append(boom3)
    penta = polygon_along(boom3[6], 5)
    polygons.append(penta)

    edge = penta[2]
    deca8 = polygon_along(edge, 10)
    polygons.append(deca8)
    for i in (7, 8, -1):
        polygons.append(polygon_along(deca8[i], 5))
    penta = polygon_along(deca8[6], 5)
    polygons.append(penta)

    edge = penta[3]
    deca9 = polygon_along(edge, 10)
    polygons.append(deca9)
    polygons.append(polygon_along(deca9[7], 5))

    edge = deca9[8]
    boom8 = boomelang_along(edge, 7)
    polygons.append(boom8)

    edge = boom3[7]
    boom4 = boomelang_along(edge, 4)
    polygons.append(boom4)

    edge = boom4[1]
    boom5 = boomelang_along(edge, 3)
    polygons.append(boom5)
    for i in (4, 5, 6):
        polygons.append(polygon_along(boom5[i], 5))

    edge = boom4[7]
    boom6 = boomelang_along(edge, 4)
    polygons.append(boom6)
    polygons.append(polygon_along(boom6[0], 5))

    edge = boom5[7]
    boom7 = boomelang_along(edge, 4)
    polygons.append(boom7)
    for i in (6,):
        polygons.append(polygon_along(boom7[i], 5))

    edge = boom7[1]
    deca5 = polygon_along(edge, 10)
    polygons.append(deca5)
    for i in (1, 2):
        polygons.append(polygon_along(deca5[i], 5))
    penta = polygon_along(deca5[3], 5)
    polygons.append(penta)
    edge = penta[2]
    deca7 = polygon_along(edge, 10)
    polygons.append(deca7)
    polygons.append(polygon_along(deca7[8], 5))
    polygons.append(boomelang_along(deca7[-1], 7))

    return polygons


def test():
    pass
    # P((200,200),15,0)
    # X((200,200),15,0)
    # Y((200,200),15,0)
    # D((200,200),15,0)


def test2(e):
    t = omega**6
    x = t * e
    cx = WIDTH / 2
    cy = HEIGHT / 2
    polygons = []
    polygons += D((cx - x / 2, cy - x / 2 * L), x / t, 0)
    polygons += P((cx + x / 2, cy - x / 2 * L), x / t, 0)
    polygons += D((cx - x / 2, cy + x / 2 * L), x / t, 0)
    polygons += P((cx + x / 2, cy + x / 2 * L), x / t, 0)
    polygons += X((cx - x / 2, cy), x / t, 0)
    polygons += Y((cx + x / 2, cy), x / t, 0)
    polygons += P((cx - x / 2 - x * cos36, cy), x / t, 0)
    polygons += D((cx + x / 2 + x * cos36, cy), x / t, 0)
    # X((250-x/2,250),x/18.0,0)
    return polygons


def test3():
    x = 150
    t = omega**6
    D((250, 250 - x / 2 * L), x / t, 0)
    D((250, 250 + x / 2 * L), x / t, 0)
    X((250, 250), x / t, 0)
    P((250 - x * cos36, 250), x / t, 0)
    P((250 + x * cos36, 250), x / t, 5)
    # X((250-x/2,250),x/18.0,0)
    Y((250 - x * cos36, 250 - x / 2 * L), x / t, 0)


def test4(x):
    polygons = []
    pentagon = pentagon_at((250, 250), x, 0)
    polygons.append(pentagon)

    edge = pentagon[3]
    boom1 = boomelang_along(edge, 6)
    polygons.append(boom1)
    polygons.append(polygon_along(boom1[4], 5))

    edge = pentagon[0]
    boom2 = boomelang_along(edge, 8)
    polygons.append(boom2)
    return polygons


def test5():
    x = 150
    decagon = decagon_at((250, 250), x, 0)
    polygons.append(decagon)


def test6(x):
    polygons = []
    pentagon = pentagon_at((250, 250), x, 0)
    polygons.append(pentagon)
    boom1 = boomelang_along(pentagon[0], 4)
    polygons.append(boom1)
    boom2 = boomelang_along(boom1[0], 1)
    polygons.append(boom2)
    penta1 = polygon_along(pentagon[1], 5)
    polygons.append(penta1)
    return polygons


def test7(x):
    polygons = []
    pentagon = pentagon_at((250, 250), x, 0)
    polygons.append(pentagon)
    boom2 = boomelang_along(pentagon[0], 4)
    polygons.append(boom2)
    boom1 = boomelang_along(pentagon[1], 0)
    polygons.append(boom1)
    penta1 = polygon_along(boom2[0], 5)
    polygons.append(penta1)
    return polygons


def test8(x):
    polygons = []
    pentagon = pentagon_at((250, 250), x, 0)
    polygons.append(pentagon)
    boom2 = boomelang_along(pentagon[0], 3)
    polygons.append(boom2)
    penta1 = polygon_along(pentagon[1], 5)
    polygons.append(penta1)
    deca1 = polygon_along(pentagon[2], 10)
    polygons.append(deca1)
    boom1 = boomelang_along(pentagon[3], 0)
    polygons.append(boom1)
    penta2 = polygon_along(boom1[6], 5)
    polygons.append(penta2)
    penta3 = polygon_along(boom2[6], 5)
    polygons.append(penta3)
    penta4 = polygon_along(deca1[7], 5)
    polygons.append(penta4)
    penta5 = polygon_along(deca1[8], 5)
    polygons.append(penta5)
    boom3 = boomelang_along(deca1[6], 1)
    polygons.append(boom3)
    penta6 = polygon_along(boom1[5], 5)
    polygons.append(penta6)
    return polygons


def test9(x):
    polygons = []
    decagon = decagon_at((250, 250), x, 0)
    polygons.append(decagon)
    boom1 = boomelang_along(decagon[1], 1)
    polygons.append(boom1)
    penta1 = polygon_along(decagon[0], 5)
    polygons.append(penta1)
    penta2 = polygon_along(boom1[3], 5)
    polygons.append(penta2)
    boom2 = boomelang_along(boom1[4], 7)
    polygons.append(boom2)
    penta3 = polygon_along(decagon[2], 5)
    polygons.append(penta3)
    penta4 = polygon_along(boom1[5], 5)
    polygons.append(penta4)
    boom3 = boomelang_along(boom1[7], 4)
    polygons.append(boom3)
    penta5 = polygon_along(boom3[5], 5)
    polygons.append(penta5)
    penta6 = polygon_along(boom3[6], 5)
    polygons.append(penta6)
    boom4 = boomelang_along(boom2[4], 1)
    polygons.append(boom4)
    deca2 = polygon_along(boom4[7], 10)
    polygons.append(deca2)
    penta7 = polygon_along(boom2[5], 5)
    polygons.append(penta7)
    penta8 = polygon_along(boom2[2], 5)
    polygons.append(penta8)
    penta9 = polygon_along(boom2[3], 5)
    polygons.append(penta9)
    penta10 = polygon_along(boom4[3], 5)
    polygons.append(penta10)
    penta11 = polygon_along(boom4[5], 5)
    polygons.append(penta11)
    penta12 = polygon_along(boom4[6], 5)
    polygons.append(penta12)
    return polygons


def ears_boomelang(boom, ears=""):
    tiles = []
    if "A" in ears:
        tiles.append(polygon_along(boom[2], 5))
    if "B" in ears or "L" in ears:
        tiles.append(polygon_along(boom[3], 5))
    if "C" in ears:
        tiles.append(polygon_along(boom[4], 5))
    if "D" in ears or "R" in ears:
        tiles.append(polygon_along(boom[5], 5))
    if "E" in ears:
        tiles.append(polygon_along(boom[6], 5))
    return tiles


def test10(x):
    polygons = []
    decagon = decagon_at((250, 250), x, 0)
    polygons.append(decagon)
    boom1 = boomelang_along(decagon[1], 7)
    polygons.append(boom1)
    polygons += ears_boomelang(boom1, "RL")
    penta1 = polygon_along(boom1[6], 5)
    polygons.append(penta1)
    boom2 = boomelang_along(penta1[2], 3)
    polygons.append(boom2)
    polygons += ears_boomelang(boom2, "A")
    boom3 = boomelang_along(boom2[5], 7)
    polygons.append(boom3)
    polygons += ears_boomelang(boom3, "LC")
    boom4 = boomelang_along(boom3[5], 7)
    polygons.append(boom4)
    polygons += ears_boomelang(boom4, "LRE")
    deca1 = polygon_along(boom4[1], 10)
    polygons.append(deca1)
    boom5 = boomelang_along(deca1[3], 1)
    polygons.append(boom5)
    polygons += ears_boomelang(boom5, "R")
    boom6 = boomelang_along(boom5[7], 4)
    polygons.append(boom6)
    polygons += ears_boomelang(boom6, "ABD")
    deca2 = polygon_along(boom6[1], 10)
    polygons.append(deca2)
    penta2 = polygon_along(deca2[5], 5)
    polygons.append(penta2)
    penta3 = polygon_along(deca2[6], 5)
    polygons.append(penta3)
    boom7 = boomelang_along(boom6[7], 4)
    polygons.append(boom7)
    polygons += ears_boomelang(boom7, "ABDE")
    return polygons


def gridify(grid, pos):
    # pos: tuple of two floats
    # grid: dict of known grid points
    ix = int(floor(pos[0] * 1000.0 + 0.5))
    iy = int(floor(pos[1] * 1000.0 + 0.5))
    if (ix, iy) in grid:
        return grid[(ix, iy)]
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            grid[(ix + dx, iy + dy)] = (ix, iy)
    return ix, iy


def extract_augments(polygons):
    # pickup the augment polygons
    tiles = []
    for poly in polygons:
        for edge in poly.edges:
            tiles.append(edge.augment)
    return tiles


def EdgeOwners(tiles):
    edgeowners = dict()
    grid1 = dict()
    grid2 = dict()
    for tile in tiles:
        if tile is None:
            continue
        iedges = []
        iedgesr = []
        for edge in tile.edges:
            start = gridify(grid1, edge.start)
            delta = gridify(grid2, edge.delta)
            iedge = Edge(start, delta)
            iedges.append(iedge)
            rev = invert(edge)
            start = gridify(grid1, rev.start)
            delta = gridify(grid2, rev.delta)
            iedge = Edge(start, delta)
            iedgesr.append(iedge)
        itile = Tile(tuple(iedges))
        for j in range(len(iedges)):
            edge = iedges[j]
            rev = iedgesr[j]
            if (edge.start, edge.delta) in edgeowners:
                edgeowners[(edge.start, edge.delta)][0] = itile
            else:
                edgeowners[(edge.start, edge.delta)] = [None] * 2
                edgeowners[(edge.start, edge.delta)][0] = itile
            if (rev.start, rev.delta) in edgeowners:
                edgeowners[(rev.start, rev.delta)][1] = itile
            else:
                edgeowners[(rev.start, rev.delta)] = [None] * 2
                edgeowners[(rev.start, rev.delta)][1] = itile
    return edgeowners


def set_parity2(queue, parity, edgeowners):
    while True:
        if len(queue) == 0:
            return
        itile, p = queue.pop(0)
        if itile is None:
            continue
        if itile in parity:
            if parity[itile] != p:
                drawitile(itile, 0.333)
                # print itile,len(itile.edges),parity[itile]
            continue
        else:
            break
    parity[itile] = p
    if len(itile.edges) == 3:
        e = itile.edges[0]
        e0 = (e.start, e.delta)
        nei0 = edgeowners[e0][1]
        if nei0 is not None:
            if len(nei0.edges) == 3:
                queue.append((nei0, -p))
            else:
                queue.append((nei0, p))
        for i in (1, 2):
            e = itile.edges[i]
            e0 = (e.start, e.delta)
            nei0 = edgeowners[e0][1]
            queue.append((nei0, p))
    else:  # house tile
        for e in itile.edges:
            e0 = (e.start, e.delta)
            nei0 = edgeowners[e0][1]
            queue.append((nei0, p))


def decorate_triangle(tile, parity):
    atoms = []
    c = tile.edges[0].start
    b = tile.edges[1].start
    a = tile.edges[2].start
    atoms.append((a, 0))
    atoms.append((b, 0))
    atoms.append((c, 0))
    cx = (b[0] + c[0]) / 2
    cy = (b[1] + c[1]) / 2
    cx = (a[0] + cx * 3) / 4
    cy = (a[1] + cy * 3) / 4
    atoms.append(((cx, cy), -parity))
    cx = (a[0] + b[0]) / 2
    cy = (a[1] + b[1]) / 2
    atoms.append(((cx, cy), parity))
    cx = (a[0] + c[0]) / 2
    cy = (a[1] + c[1]) / 2
    atoms.append(((cx, cy), parity))
    return atoms


def decorate_rectangle(tile, parity):
    atoms = []
    d = tile.edges[0].start
    a = tile.edges[1].start
    b = tile.edges[2].start
    c = tile.edges[3].start
    atoms.append((a, 0))
    atoms.append((b, 0))
    atoms.append((c, 0))
    atoms.append((d, 0))
    cx = (a[0] * 3 + b[0] * 3 + c[0] * 1 + d[0] * 1) / 8
    cy = (a[1] * 3 + b[1] * 3 + c[1] * 1 + d[1] * 1) / 8
    atoms.append(((cx, cy), -parity))
    cx = (a[0] * 1 + b[0] * 1 + c[0] * 3 + d[0] * 3) / 8
    cy = (a[1] * 1 + b[1] * 1 + c[1] * 3 + d[1] * 3) / 8
    atoms.append(((cx, cy), -parity))
    cx = (a[0] * 1 + b[0] * 4 + c[0] * 4 + d[0] * 1) / 10
    cy = (a[1] * 1 + b[1] * 4 + c[1] * 4 + d[1] * 1) / 10
    atoms.append(((cx, cy), parity))
    cx = (a[0] * 4 + b[0] * 1 + c[0] * 1 + d[0] * 4) / 10
    cy = (a[1] * 4 + b[1] * 1 + c[1] * 1 + d[1] * 4) / 10
    atoms.append(((cx, cy), parity))
    return atoms


def decorate_house(tile, parity):
    triangle = Tile(tile.edges[4:5] + tile.edges[2:4])
    rectangle = Tile(tile.edges[0:3] + tile.edges[4:5])
    a1 = decorate_triangle(triangle, parity)
    a2 = decorate_rectangle(rectangle, parity)
    return a1 + a2


def drawitile(polygon, p, label=""):
    # fill(random(),1,1,0.5)
    # to check overlaps
    if len(polygon.edges) == 5:
        return decorate_house(polygon, p)
    else:
        return decorate_triangle(polygon, p)


def decorate(polygons):
    # pick up triangles and houses
    tiles = extract_augments(polygons)
    edgeowners = EdgeOwners(tiles)

    # first edge
    itile = None
    for edge in edgeowners:
        if edgeowners[edge][0] is not None and edgeowners[edge][1] is not None:
            # first tile
            itile = edgeowners[edge][0]
            break
    parity = dict()
    #set_parity(parity, itile, 1, edgeowners, 6)
    queue = [(itile, 1)]
    while len(queue) > 0:
        set_parity2(queue, parity, edgeowners)

    atoms = []
    for itile in parity:
        p = parity[itile]
        atoms += drawitile(itile, p)
    return atoms


def igridify(grid, pos):
    # pos: tuple of two floats
    # grid: dict of known grid points
    ix, iy = pos
    if (ix, iy) in grid:
        return grid[(ix, iy)]
    for dx in range(-1, 2):
        for dy in range(-1, 2):
            grid[(ix + dx, iy + dy)] = (ix, iy)
    return ix, iy


class Lattice(genice2.lattices.Lattice):
    def __init__(self, **kwargs):
        logger = getLogger()
        # global sides, rows, bondlen, density, cell, waters, coord

        a = 6
        b = 0
        for k, v in kwargs.items():
            if k == "size":
                arg = int(v)
            elif v is True:
                # unlabeled option
                arg = int(k)
            else:
                sys.exit(1)


        edgelength = 40

        if arg == 4:
            # test4
            w = int(edgelength * sqrt(5) * omega**1 * 1000 + 0.5)
            h = int(omega**3 * edgelength * L * 1000 + 0.5)
            polygons = test4(edgelength)  # small approximant #1
        elif arg == 6:
            # test6
            w = int(edgelength * sqrt(5) * omega**2 * 1000 + 0.5)
            h = int(omega**2 * edgelength * L * 1000 + 0.5)
            polygons = test6(edgelength)  # small approximant #2
        elif arg == 8:
            # test8
            w = int(edgelength * sqrt(5) * omega**3 * 1000 + 0.5)
            h = int(omega**3 * edgelength * L * 1000 + 0.5)
            polygons = test8(edgelength)  # small approximant #3
        elif arg == 9:
            # test9
            w = int(edgelength * sqrt(5) * omega**3 * 1000 + 0.5)
            h = int(omega**4 * edgelength * L * 1000 + 0.5)
            polygons = test9(edgelength)  # small approximant #4 designed by matto
        elif arg == 10:
            # test10
            w = int(edgelength * sqrt(5) * omega**4 * 1000 + 0.5)
            h = int(omega**4 * edgelength * L * 1000 + 0.5)
            polygons = test10(edgelength)  # small approximant #4 designed by matto

        # polygons = test2(edgelength) #huge
        # polygons = test7(edgelength) #small approximant #3
        #polygons = P((250,250),20,0)
        # for poly in polygons:
        #    drawpoly(poly)

        atoms = decorate(polygons)

        # set Z coord
        depth = edgelength * 0.34    # A-X layer distance

        coord = dict()
        grid = dict()
        for atom in atoms:
            x, y = atom[0]
            x = (x + w) % w
            y = (y + h) % h
            x, y = igridify(grid, (x, y))
            layer = atom[1]
            z = layer * depth
            if layer == 0:
                coord[(x, y, z + 2 * depth)] = 1
                coord[(x, y, z + 6 * depth)] = 1
            coord[(x, y, z)] = 1
            coord[(x, y, z + 4 * depth)] = 1

        # import yaplotlib as yp
        import numpy as np
        box3 = np.array((w * 1e-3, h * 1e-3, depth * 4 * 2))
        # cell matrix
        self.cell = np.diag(box3)
        # relative coord
        coord = np.array([[x / 1000, y / 1000, z] for x, y, z in coord]) / box3

        from genice.FrankKasper import toWater
        self.waters = np.array([w for w in toWater(coord, self.cell)])
        self.coord = "relative"
        self.density = 0.8
        self.bondlen = 14
