import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import matplotlib.animation as animation
from matplotlib import style
import math
from itertools import combinations
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
line, = ax.plot(0, 0)
x_data=[]
y_data=[]


# returns a list which contain the elements of the Farey sequence
def farey_sequence(n, descending, x):
    seq = [[0, 1]]
    (a, b, c, d) = (0, 1, 1, n)
    if descending:
        (a, c) = (1, n-1)
        print(a, b)

    while(c <= n and not descending) or (a > 0 and descending):
        k = (n + b) // d
        (a, b, c, d) = (c, d, k * c - a, k * d - b)
        seq.append([a, b])

    # if x is greater than 1 we must translate the sequence
    if x > 1:
        original_seq = seq.copy()
        original_seq.remove(original_seq[0])
        for i in range(math.ceil(x)):
            for j in original_seq:
                seq.append([j[1] * i + j[0], j[1]])

    return seq

# draws the axis and calls the necessary geodesic and Farey drawing functions
def draw(x, p):
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    theta = intersection(x, p)
    theta_start = theta[2]
    theta_end = theta[3]
    #geodesic((x - p, 0), p, 'pink', theta_start, ax)
    special_case_geodesic((x - p, 0), p, theta_start, theta_end)
    special_case_geodesic((abs(x - p) + 1, 0), p, theta_start + np.radians(90), theta_end + np.radians(90))

    # control the depth to which the Farey diagram is drawn
    depth = 8
    draw_farey(farey_sequence(depth, False, x))

    plt.xlim(0, math.ceil(x)+.01)
    plt.ylim(0, math.ceil(x)+.01)

    #plt.savefig('fareyGauss4.1.png', dpi=200)
    plt.show()


# draws a geodesic arc
def geodesic(center, radius, color, angle, ax):
    wedge = Wedge(center, radius, angle, angle + 180, fill=False, ec=color)
    ax.add_artist(wedge)
    return wedge


def special_case_geodesic(center, radius, theta_start, theta_end):
    n = 64
    t = np.linspace(theta_start, theta_end, n+1)
    x = center[0] + radius*np.cos(t)
    y = radius*np.sin(t)

    plt.plot(x, y)


# draws the Farey sequence
def draw_farey(seq):
    for i in range(seq[len(seq) - 1][0]+1):
        plt.axvline(x=i, color='blue')

    for i in combinations(seq, 2):
        if not farey_neighbors(i[0], i[1]):
            continue
        if i[0][1] == 0:
            p1 = 0
            p2 = i[1][0] / i[1][1]

        elif i[1][1] == 0:
            p1 = i[0][0] / i[0][1]
            p2 = 0
        else:
            p1 = i[0][0] / i[0][1]
            p2 = i[1][0] / i[1][1]
        midpoint = (p1 + p2) / 2
        radius = abs(p1 - p2)/2
        geodesic((midpoint, 0), radius, 'blue', 0, ax)


# returns true if two numbers are farey neighbors
def farey_neighbors(p1, p2):
    if abs(p1[0]*p2[1] - p1[1]*p2[0]) == 1:
        return True
    else:
        return False


def gauss_map(x):
    return (1/x) % 1


# if intersection, check whether it is with x = 1 or a curved geodesic
def intersection(x, p):
    # let g be a tuple which contains x and p
    # f is the semi circle with radius .5 centered at (.5, 0)
    # if x < 1, we know that the intersection occurs with the semicircle geodesic
    radius = p
    if x < 1:
        d = abs(.5 - (x - radius))
        print('d: ' + str(d))
        # we need to find the point of intersection
        xpos = (radius**2 - .25 + d**2) / (2 * d) - abs(x-radius)
        print('x: ' + str(xpos))
        print('sqrt( ' + str(radius**2 - (xpos + abs((x - p)))**2))
        ypos = math.sqrt(radius**2 - (xpos + abs((x - p)))**2)
        print('y: ' + str(ypos))
        theta_end = np.arcsin(ypos/radius)
        theta_start = np.arccos(abs(x-radius)/radius)
        print('theta start: ' + str(np.degrees(theta_start)))
        print('theta end: ' + str(np.degrees(theta_end)))

        return [xpos, ypos, theta_start, theta_end]


if __name__ == '__main__':
    print(1/math.sqrt(2))
    x = float(input('Please input a value for x: '))
    p = float(input('Input a value for p such that ceil(x) >= p >= x/2: '))

    while p < x/2 or p > math.ceil(x):
        p = float(input('Please input a valid p such that ceil(x) >= p >= x/2: '))

    draw(x, p)



