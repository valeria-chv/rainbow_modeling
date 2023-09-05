from math import cos, pi, sin
from math import tan, atan, asin


def arctg(x):
    return atan(x)*180/3.1415


def tg(x):
    return tan(x*3.1415/180)


def sin_(x):
    return sin(x*3.1415/180)

def cos_(x):
    return cos(x*3.1415/180)


def arcsin(x):
    return asin(x)*180/3.1415


def sign(x):
    if x >= 0: return 1
    else: return -1


def turn_x(x, y, z, alpha):
    alpha = alpha * pi / 180
    buf = y
    y = cos(alpha) * y - sin(alpha) * z
    z = cos(alpha) * z + sin(alpha) * buf
    return x, y, z


def turn_y(x, y, z, alpha):
    alpha = alpha * pi / 180
    buf = x
    x = cos(alpha) * x - sin(alpha) * z
    z = cos(alpha) * z + sin(alpha) * buf
    return x, y, z


def turn_z(x, y, z, alpha):
    alpha = alpha * pi / 180
    buf = x
    x = cos(alpha) * x - sin(alpha) * y
    y = cos(alpha) * y + sin(alpha) * buf
    return x, y, z


def transform(x, y, z, alpha_x, alpha_y, alpha_z, scale_k, w, h):
    x, y, z = turn_x(x, y, z, alpha_x)
    x, y, z = turn_y(x, y, z, alpha_y)
    x, y, z = turn_z(x, y, z, alpha_z)
    x = x * scale_k + w / 2
    y = y * scale_k + h / 2
    return round(x), round(y), round(z)
