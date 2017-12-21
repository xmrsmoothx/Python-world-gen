import math
import numpy as np
import random
from scipy.spatial import Voronoi


def lerp(t, a, b):
    return a + t * (b - a)


def smoothcurve(t):
    return t * t * (3. - 2. * t)


def lengthDirX(length, angle):
    radian_angle = math.radians(angle)
    return length * math.cos(radian_angle)


def lengthDirY(length, angle):
    radian_angle = math.radians(angle)
    return length * math.sin(radian_angle)


def A(dx, dy):
    return math.degrees(math.atan2(dy, dx))


def drawCircle(drawer, x, y, radius, color):
    x1 = x - radius
    x2 = x + radius
    y1 = y - radius
    y2 = y + radius
    drawer.ellipse([(x1, y1), (x2, y2)], color)


def drawTrapezoid(drawer, x1, y1, x2, y2, r1, r2, color):
    directAngle = A(x2 - x1, y2 - y1)
    pAngle = directAngle - 90
    pAngle2 = pAngle - 180
    p1 = (x1 + lengthDirX(r1, pAngle), y1 + lengthDirY(r1, pAngle))
    p2 = (x1 + lengthDirX(r1, pAngle2), y1 + lengthDirY(r1, pAngle2))
    p3 = (x2 + lengthDirX(r2, pAngle2), y2 + lengthDirY(r2, pAngle2))
    p4 = (x2 + lengthDirX(r2, pAngle), y2 + lengthDirY(r2, pAngle))
    drawer.polygon([p1, p2, p3, p4], color, color)


def centroidnp(arr):
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    return sum_x / length, sum_y / length


def unitCos(x):
    return (math.cos(x * math.pi) + 1) / 2


def distMod(x, maxDist):
    return unitCos(x / maxDist)


def stick(val, minimum, maximum):
    if abs(val - minimum) < abs(val - maximum):
        return minimum
    else:
        return maximum


def clamp(x, minimum, maximum):
    if x < minimum:
        return minimum
    elif x > maximum:
        return maximum
    else:
        return x


def strDivider(length):
    n = ""
    for l in range(length):
        n += "_"
    return n


def relaxLloyd(pts, strength):
    for i in range(strength):
        vor = Voronoi(pts)
        newpts = []
        for idx in range(len(vor.points)):
            pt = vor.points[idx, :]
            region = vor.regions[vor.point_region[idx]]
            if -1 in region:
                newpts.append(pt)
            else:
                vxs = np.asarray([vor.vertices[i, :] for i in region])
                newpt = centroidnp(vxs)
                newpts.append(newpt)
        pts = np.array(newpts)
