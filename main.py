import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
import math
from itertools import combinations

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)


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


def draw(x, p):
    ax.spines['left'].set_position('zero')
    ax.spines['bottom'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    geodesic((x-p, 0), p, 'orange', 0, ax)
    draw_farey(farey_sequence(10, False, x))

    plt.xlim(0, math.ceil(x))
    plt.ylim(0, math.ceil(x))
    plt.show()


def geodesic(center, radius, color, angle, ax):
    wedge = Wedge(center, radius, angle, angle + 180, fill=False, ec=color)
    ax.add_artist(wedge)
    return wedge


def draw_farey(seq):
    for i in range(seq[len(seq) - 1][0]*2):
        plt.axvline(x=i/2, color='blue')

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


def farey_neighbors(p1, p2):
    if abs(p1[0]*p2[1] - p1[1]*p2[0]) == 1:
        return True
    else:
        return False


if __name__ == '__main__':
    x = float(input('Please input a value for x: '))
    p = float(input('Input a value for p such that ceil(x) >= p >= x/2: '))

    while p < x/2 or p > math.ceil(x):
        p = float(input('Please input a valid p such that ceil(x) >= p >= x/2: '))

    draw(x, p)
