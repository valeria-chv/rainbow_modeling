'''
50, 120 -> -50, 350
-130, -80 -> -230, 150
'''

import sys
import win2
import time

from func import funcs
# from draw_hor import float_horizon
from my_math import transform, arcsin, tg, sign, cos_, arctg, sin_


from PyQt5 import QtWidgets, QtCore, QtGui
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtGui import QImage, QPixmap
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor


wind = None


def pixelDraw(image, x, y, colorPen):
    image.setPixelColor(round(x), round(y), colorPen)


# алгоритм Брезенхема для окружности
def drawCircleBrez(image, xc, yc, rad, colorPen):
    xc, yc = round(xc), round(yc)

    x, y = 0, round(rad)
    deltaI = 2 * (1 - rad)  # deltaI = (x+1)^2 + (y-1)^2 - R^2
    limit = round(rad / (2)**(1/2))

    while y >= limit:
        pixelDraw(image, xc + x, yc + y, colorPen)
        pixelDraw(image, xc + x, yc - y, colorPen)
        pixelDraw(image, xc - x, yc + y, colorPen)
        pixelDraw(image, xc - x, yc - y, colorPen)

        pixelDraw(image, xc + y, yc - x, colorPen)
        pixelDraw(image, xc + y, yc + x, colorPen)
        pixelDraw(image, xc - y, yc - x, colorPen)
        pixelDraw(image, xc - y, yc + x, colorPen)

        if deltaI < 0:  # диаг. пиксель лежит внутри окружности
            d1 = deltaI + deltaI + y + y - 1
            x += 1
            if d1 <= 0:  # горизонтальный шаг
                deltaI += x + x + 1  # уже увеличенный x
            else:  # диагональный шаг
                y -= 1
                deltaI += 2 * (x - y + 1)
        elif deltaI >= 0:  # диаг. пиксель на/ вне окружности
            # диагональный шаг
            x += 1
            y -= 1
            deltaI += 2 * (x - y + 1)


# алгоритм Брезенхема с целочисленными данными
def brezIntDraw(image, xs, ys, xe, ye, colorPen):
    xs = round(xs)
    ys = round(ys)
    swap = False
    dx = (xe - xs)
    dy = (ye - ys)
    xSign = sign(dx)
    ySign = sign(dy)
    dx = abs(dx)
    dy = abs(dy)
    if dy > dx:
        dx, dy = dy, dx
        swap = True
    e = 2 * dy - dx
    i = 1
    while i <= dx + 1:

        pixelDraw(image, xs, ys, colorPen)
        if swap:
            if e >= 0:
                xs += xSign
                e -= 2 * dx
            ys += ySign
        else:
            if e >= 0:
                ys += ySign
                e -= 2 * dx
            xs += xSign
        e += 2 * dy
        i += 1



def comparing_now_point_and_array(xs, ys, zs, xe, ye, ze, array):
    len_z = len(array)

    if xs > xe: sign_x = -1
    else: sign_x = 1
    if ys > ye: sign_y = -1
    else: sign_y = 1
    if zs > ze: sign_z = -1
    else: sign_z = 1

    if sign_z == 1:
        i1 = 0
    else:
        i1 = len_z - 1
    j1 = 0; k1 = 0

    while (i1 > 0 and sign_z == -1 and array[i1][j1][k1][2] > round(zs) or i1 < len_z - 1 and sign_z == 1 and array[i1][j1][k1][2] < round(zs)):
        i1 += sign_z

    if array[i1][j1][k1][2] == round(zs):
        len_y = len(array[i1])
        if sign_y == 1:
            j1 = 0
        else:
            j1 = len_y - 1

        while (j1 > 0 and sign_y == -1 and array[i1][j1][k1][1] > round(ys) or j1 < len_y - 1 and sign_y == 1 and array[i1][j1][k1][1] < round(ys)):
            j1 += sign_y

        if array[i1][j1][k1][1] == round(ys):
            len_x = len(array[i1][j1])
            if sign_x == 1:
                k1 = 0
            else:
                k1 = len_x - 1

            while (k1 > 0 and sign_x == -1 and array[i1][j1][k1][0] > round(xs) or k1 < len_x - 1 and sign_x == 1 and array[i1][j1][k1][0] < round(xs)):
                k1 += sign_x

            if array[i1][j1][k1][0] == round(xs):
                return 1
    return 0


# пересечение луча и точки графически (и массива)
def intersection_of_beam_and_point(image, color, xs, ys, zs, xe, ye, ze, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h, array):
    xs = round(xs)
    ys = round(ys)
    zs = round(zs)
    if (abs(xe - xs) >= abs(ye - ys) and abs(xe - xs) >= abs(ze - zs)):
        l_points = abs((xe - xs))
    elif (abs(ye - ys) >= abs(xe - xs) and abs(ye - ys) >= abs(ze - zs)):
        l_points = abs((ye - ys))
    else:
        l_points = abs((ze - zs))
    dx = (xe - xs) / l_points
    dy = (ye - ys) / l_points
    dz = (ze - zs) / l_points
    i = 1

    while i <= l_points+1:
        key = comparing_now_point_and_array(xs, ys, zs, xe, ye, ze, array)
        if key:
            return 1, round(xs), round(ys), round(zs)

        x0_beam_2d, y0_beam_2d, z_buf = transform(round(xs), round(ys), round(zs), alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h)
        pixelDraw(image, round(x0_beam_2d), round(y0_beam_2d), color)
        xs += dx
        ys += dy
        zs += dz
        i += 1

    return 0, round(xs), round(ys), round(zs)

'''
# пересечение луча и точки (x, y, z)
def intersection_of_beam_and_point(x0_beam, y0_beam, z0_beam, alpha_x, alpha_y, alpha_z, x, y, z):
    if alpha_z < 0: alpha_z += 360
    if x0_beam == x: alpha_z_ = 0
    else: alpha_z_ = arctg((y - y0_beam)/(x - x0_beam))  # угол входящего в каплю луча
    if (x0_beam - x < 0 and y0_beam - y > 0): alpha_z_ = 360 + alpha_z_
    if (x0_beam - x > 0 and y0_beam - y < 0): alpha_z_ = 180 + alpha_z_
    if (x0_beam - x > 0 and y0_beam - y > 0): alpha_z_ = 180 + alpha_z_

    if alpha_y < 0: alpha_y += 360
    if z0_beam == z: alpha_y_ = 0
    else: alpha_y_ = arctg((x0_beam - x)/(z0_beam - z))
    if (z0_beam - z < 0 and x0_beam - x > 0): alpha_y_ = 360 + alpha_y_
    if (z0_beam - z > 0 and x0_beam - x < 0): alpha_y_ = 180 + alpha_y_
    if (z0_beam - z > 0 and x0_beam - x > 0): alpha_y_ = 180 + alpha_y_

    if alpha_x < 0: alpha_x += 360
    if y0_beam == y: alpha_x_ = 0
    else: alpha_x_ = arctg((z0_beam - z)/(y0_beam - y))
    if (y0_beam - y < 0 and z0_beam - z > 0): alpha_x_ = 360 + alpha_x_
    if (y0_beam - y > 0 and z0_beam - z < 0): alpha_x_ = 180 + alpha_x_
    if (y0_beam - y > 0 and z0_beam - z > 0): alpha_x_ = 180 + alpha_x_

    if abs(alpha_z - alpha_z_) > 1 or abs(alpha_y - alpha_y_) > 1 or abs(alpha_x - alpha_x_) > 1:
        #print(x0_beam, y0_beam, z0_beam, alpha_x, alpha_x_, alpha_y, alpha_y_, alpha_z, alpha_z_, x, y, z)
        return 0
    else:
        return 1
    ''''''
    if abs((x - x0_beam) * cos_(alpha_y) - (z - z0_beam) * cos_(alpha_z) < 0.05):
        print(1)
        return 1
    else:
        #print(2)
        return 0'''


def draw_after_intersection_of_beam_and_drop(image, x0_beam, y0_beam, z0_beam, alpha_x, alpha_y, alpha_z, Q, x_drop, y_drop, z_drop, x_eye, z_eye, color, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h, beta_x, beta_y, beta_z):
    x_beam_beg = x_drop
    y_beam_beg = y_drop
    z_beam_beg = z_drop

    b = y_drop - tg(beta_z)*x_drop
    x_beam_end = x_eye
    y_beam_end = tg(beta_z)*x_beam_end + b
    z_beam_end = z_eye

    x0_beam_2d, y0_beam_2d, z_buf = transform(x0_beam, y0_beam, z0_beam, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h)
    x_beam_beg_2d, y_beam_beg_2d, z_buf = transform(x_beam_beg, y_beam_beg, z_beam_beg, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h)
    x_beam_end_2d, y_beam_end_2d, z_buf = transform(x_beam_end, y_beam_end, z_beam_end, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h)

    brezIntDraw(image, x0_beam_2d, y0_beam_2d, x_beam_beg_2d, y_beam_beg_2d, color)
    brezIntDraw(image, x_beam_beg_2d, y_beam_beg_2d, x_beam_end_2d, y_beam_end_2d, color)

# отрисовка при пересечении луча и капли
def find_beta_after_intersection_of_beam_and_drop(alpha_x, alpha_y, alpha_z, Q):
    R = 1
    n = 4/3
    E = Q/R

    gamma = 4*arcsin(E/n) - 2*arcsin(E)
    #print("gamma: ", gamma)
    # угол выходящего из капли луча
    beta_z = alpha_z + 180 - gamma
    beta_y = alpha_y
    beta_x = alpha_x

    return beta_x, beta_y, beta_z


def find_o_(x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, x_eye, y_eye, z_eye, beta_x, beta_y, beta_z):
    x_o_ = (y_eye + y_light)/(tg(beta_z) + tg(alpha_z))*(-1)
    y_o_ = tg(beta_z)*x_o_ + y_eye
    z_o_ = z_eye
    R_rainbow = tg(alpha_z+beta_z)*((x_o_ - x_eye)**2 + (y_o_ - y_eye)**2)**(1/2)

    x_rainbow = x_o_
    y_rainbow = tg(beta_z)*x_rainbow + y_eye
    z_rainbow = (R_rainbow**2 - y_rainbow**2)**(1/2)

    return x_o_, y_o_, z_o_, R_rainbow, x_rainbow, y_rainbow, z_rainbow


# def func_for_all_drops(image, scene, x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, Q, x_eye, y_eye, z_eye, color_res, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h):
#    func_for_one_drop(image, scene, x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, Q, x_eye, y_eye, z_eye, color_res, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h, x_drop, y_drop, z_drop, R_drop)


def norm_angles(alpha_x, alpha_y, alpha_z):


    return alpha_x, alpha_y, alpha_z


def find_coord_wiht_one_of_them(x0, y0, z0, alpha_x, alpha_y, alpha_z, coord, key):
    #alpha_x_, alpha_y_, alpha_z_ = norm_angles(alpha_x, alpha_y, alpha_z)
    #print(key, x0, y0, z0, alpha_x, alpha_y, alpha_z, coord)
    if key == 1:
        # coord, 0, 0 -> углы: 90,0,90
        if alpha_z >= 90: alpha_z_ = alpha_z - 90
        else: alpha_z_ = 360 + alpha_z - 90
        alpha_y_ = alpha_y
        if alpha_x >= 90: alpha_x_ = alpha_x - 90
        else: alpha_x_ = 360 + alpha_x - 90
        # coord, 0, 0
        coord -= x0

        x = coord * cos_(alpha_z_)
        y = coord * sin_(alpha_z_)
        z = 0

        x = x * cos_(alpha_y_) + z * sin_(alpha_y_)
        z = -1 * x * sin_(alpha_y_) + z * cos_(alpha_y_)

        y = y * cos_(alpha_x_) - z * sin_(alpha_x_)
        z = y * sin_(alpha_x_) + z * cos_(alpha_x_)

        z *= coord/x
        y *= coord/x
        x *= coord/x
        x += x0; y += y0; z += z0
        #print("x, y, z: ", x, y, z)
    if key == 2:
        # 0, coord, 0 -> улгы: 90, 90, 0
        if alpha_z >= 90: alpha_z_ = alpha_z - 90
        else: alpha_z_ = 360 + alpha_z - 90
        if alpha_y >= 90: alpha_y_ = alpha_y - 90
        else: alpha_y_ = 360 + alpha_y - 90
        alpha_x_ = alpha_x
        # 0, coord, 0
        coord -= y0

        x = -1 * coord * sin_(alpha_z_)
        y = coord * cos_(alpha_z_)
        z = 0

        x = x * cos_(alpha_y_) + z * sin_(alpha_y_)
        z = -1 * x * sin_(alpha_y_) + z * cos_(alpha_y_)

        y = y * cos_(alpha_x_) - z * sin_(alpha_x_)
        z = y * sin_(alpha_x_) + z * cos_(alpha_x_)

        z *= coord/y
        x *= coord/y
        y *= coord/y

        x += x0; y += y0; z += z0

    if key == 3:
        # 0, 0, coord -> улгы: 0, 90, 90
        if alpha_x >= 90: alpha_x_ = alpha_x - 90
        else: alpha_x_ = 360 + alpha_x - 90
        if alpha_y >= 90: alpha_y_ = alpha_y - 90
        else: alpha_y_ = 360 + alpha_y - 90
        alpha_z_ = alpha_z
        # 0, 0, coord
        coord -= z0

        x = 0
        y = 0
        z = coord

        x = x * cos_(alpha_y_) + z * sin_(alpha_y_)
        z = -1 * x * sin_(alpha_y_) + z * cos_(alpha_y_)

        y = y * cos_(alpha_x_) - z * sin_(alpha_x_)
        z = y * sin_(alpha_x_) + z * cos_(alpha_x_)

        z *= coord/z
        x *= coord/z
        y *= coord/z

        x += x0; y += y0; z += z0

    return x, y, z


def func_for_one_drop(image, scene, x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, Q, x_eye, y_eye, z_eye, color_res, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h, x_drop_max, y_drop_max, z_drop_max, R_drop, array_drops, x_eye_min, y_eye_min, z_eye_min, array_eyes):
    p = QPixmap()

    dx_drop = x_drop_max - x_light
    dy_drop = y_drop_max - y_light
    dz_drop = z_drop_max - z_light

    if dx_drop > dy_drop and dx_drop > dz_drop:
        xe, ye, ze = find_coord_wiht_one_of_them(x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, x_drop_max, 1)
    elif dy_drop > dx_drop and dy_drop > dz_drop:
        xe, ye, ze = find_coord_wiht_one_of_them(x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, y_drop_max, 2)
    else:
        xe, ye, ze = find_coord_wiht_one_of_them(x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, z_drop_max, 3)
    #print(x_light, y_light, z_light, xe, ye, ze)
    key, x_drop, y_drop, z_drop = intersection_of_beam_and_point(image, QtCore.Qt.red, x_light, y_light, z_light, xe, ye, ze, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h, array_drops)
    if key:
        print(1)
        beta_x, beta_y, beta_z = find_beta_after_intersection_of_beam_and_drop(alpha_x, alpha_y, alpha_z, Q)

        dx_eye = x_eye_min - x_drop
        dy_eye = y_eye_min - y_drop
        dz_eye = z_eye_min - z_drop

        if dx_eye < dy_eye and dx_eye < dz_eye:
            xe, ye, ze = find_coord_wiht_one_of_them(x_drop, y_drop, z_drop, beta_x, beta_y, beta_z, x_eye_min, 1)
        elif dy_eye < dx_eye and dy_eye < dz_eye:
            xe, ye, ze = find_coord_wiht_one_of_them(x_drop, y_drop, z_drop, beta_x, beta_y, beta_z, y_eye_min, 2)
        else:
            xe, ye, ze = find_coord_wiht_one_of_them(x_drop, y_drop, z_drop, beta_x, beta_y, beta_z, z_eye_min, 3)
        #print(alpha_x, alpha_y, alpha_z, beta_x, beta_y, beta_z)
        print(22, x_drop, y_drop, z_drop, xe, ye, ze)
        #print(array_eyes)
        key, x_drop, y_drop, z_drop = intersection_of_beam_and_point(image, QtCore.Qt.blue, x_drop, y_drop, z_drop, xe, ye, ze, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h, array_eyes)
        #draw_after_intersection_of_beam_and_drop(image, x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, Q, x_drop, y_drop, z_drop, x_eye, z_eye, color_res, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h, beta_x, beta_y, beta_z)
        if key:
            print(3)
            #draw_after_intersection_of_beam_and_drop(image, x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, Q, x_drop, y_drop, z_drop, x_eye, z_eye, color_res, alpha_x_all, alpha_y_all, alpha_z_all, scale, w, h, beta_x, beta_y, beta_z)
            #x_o_, y_o_, z_o_, R_rainbow, x_rainbow, y_rainbow, z_rainbow = find_o_(x_light, y_light, z_light, alpha_x, alpha_y, alpha_z, x_eye, y_eye, z_eye, beta_x, beta_y, beta_z)

            #print(1)
            #print("x_o, y_o, z_o:", x_o_, y_o_, z_o_)
            #print("x_rainbow, y, z, R:", x_rainbow, y_rainbow, z_rainbow, R_rainbow)

            #xc, yc, zc = transform(x_o_, y_o_, z_o_, alpha_x + alpha_x_all, alpha_y + alpha_y_all, alpha_z + alpha_z_all, scale, w, h)
            #print("xc, yc: ", xc, yc)
            #drawCircleBrez(image, xc, yc, 10, color_res)

            #xc, yc, zc = transform(x_rainbow, y_rainbow, z_rainbow, alpha_x + alpha_x_all, alpha_y + alpha_y_all, alpha_z + alpha_z_all, scale, w, h)
            #drawCircleBrez(image, xc, yc, 7, color_res)

            #p = QPixmap()
            #p.convertFromImage(image)
            #scene.addPixmap(p)
            #self.show()


    p.convertFromImage(image)
    scene.addPixmap(p)
        #self.show()


def drawing_axes(self):
    xs, ys, z_buf = transform(0, 0, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    xe, ye, z_buf = transform(1000, 0, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xs, ys, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(990, 10, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(990, -10, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(990, 0, 10, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(990, 0, -10, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)

    #xs, ys, z_buf = transform(0, 0, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    xe, ye, z_buf = transform(0, -1000, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xs, ys, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(10, -990, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(-10, -990, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(0, -990, 10, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(0, -990, -10, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)

    #xs, ys, z_buf = transform(0, 0, 0, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    xe, ye, z_buf = transform(0, 0, 1000, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xs, ys, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(0, 10, 990, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(0, -10, 990, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(10, 0, 990, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)
    xc, yc, z_buf = transform(-10, 0, 990, self.alpha_x_all, self.alpha_y_all, self.alpha_z_all, self.scale_k, self.w, self.h)
    brezIntDraw(self.image, xc, yc, xe, ye, QtCore.Qt.black)


class Visual(QtWidgets.QMainWindow, win2.Ui_MainWindow):
    def __init__(self):
        super().__init__()
        self.setupUi(self)

        self.scene = QtWidgets.QGraphicsScene()
        self.graphicsView.setScene(self.scene)
        self.h = self.graphicsView.height()
        self.w = self.graphicsView.width()
        self.scene.setSceneRect(0, 0, self.w-2, self.h-2)
        self.image = QImage(self.w, self.h, QImage.Format_RGB30)
        self.image.fill(Qt.white)
        self.pen_res = QtGui.QPen(QtCore.Qt.red)
        self.pen_res.setWidth(0)

        self.scale_k = 0.2
        self.x_begin = -100
        self.x_end = 100
        self.z_begin = -100
        self.z_end = 100
        self.x_step = 0.5
        self.z_step = 0.5
        self.number_func = 0
        self.alpha_x_all = 0
        self.alpha_y_all = 0
        self.alpha_z_all = 0
        self.color_back = QtCore.Qt.white
        self.flag_scale_or_not = False
        self.color_res = QtCore.Qt.red

        self.radioButton_func_1.clicked.connect(self.set_number_func_0)
        self.radioButton_func_2.clicked.connect(self.set_number_func_1)
        self.radioButton_func_3.clicked.connect(self.set_number_func_2)
        self.radioButton_func_4.clicked.connect(self.set_number_func_3)
        self.radioButton_func_5.clicked.connect(self.set_number_func_4)
        self.radioButton_func_6.clicked.connect(self.set_number_func_5)

        self.radioButtonBlack_res.clicked.connect(self.set_black_res)
        self.radioButtonBlue_res.clicked.connect(self.set_blue_res)
        self.radioButtonGreen_res.clicked.connect(self.set_green_res)
        self.radioButtonRed_res.clicked.connect(self.set_red_res)
        self.radioButtonWhite_res.clicked.connect(self.set_white_res)
        self.radioButtonYellow_res.clicked.connect(self.set_yellow_res)

        self.pushButton_k.clicked.connect(self.scale)
        self.pushButton_x_turn.clicked.connect(self.turn_x1)
        self.pushButton_y_turn.clicked.connect(self.turn_y1)
        self.pushButton_z_turn.clicked.connect(self.turn_z1)
        self.pushButton_clean.clicked.connect(self.clean_screen)

        self.pushButton_RES.clicked.connect(self.draw_res)

        # rainbow

        # eyes
        self.eye_z0 = -100
        self.eye_zend = 100
        self.eye_y0 = 500
        self.eye_yend = 700
        self.eye_x0 = 10
        self.eye_xend = 11

        self.eye_step_z = 1
        self.eye_step_y = 1
        self.eye_step_x = 1

        self.array_eyes = []
        self.create_array_eyes()
        # eyes

        # drops
        self.R_drop = 1

        self.drop_z0 = -100
        self.drop_zend = 100
        self.drop_y0 = -150
        self.drop_yend = 150
        self.drop_x0 = 600
        self.drop_xend = 700

        self.drop_step_z = 5
        self.drop_step_y = 5
        self.drop_step_x = 5

        self.array_drops = []
        self.create_array_drops()
        # drops

        self.x_light = 0
        self.y_light = 0
        self.z_light = 0
        self.alpha_x = 170
        self.alpha_y = 0
        self.alpha_z = 0

        self.Q_red = 0.86
        self.Q_orange = 0.86
        self.Q_yellow = 0.86
        self.Q_green = 0.86
        self.Q_cyan = 0.86
        self.Q_blue = 0.86
        self.Q_purple = 0.86

        self.n_water = 4/3

    def create_array_eyes(self):
        #x_eye = self.eye_x0
        for z_eye in range(self.eye_z0, self.eye_zend, self.eye_step_z):
            array_eyes_z = []
            for y_eye in range(self.eye_y0, self.eye_yend, self.eye_step_y):
                array_eyes_y = []
                for x_eye in range(self.eye_x0, self.eye_xend, self.eye_step_x):
                    array_eyes_y.append([x_eye, y_eye, z_eye])
                array_eyes_z.append(array_eyes_y)
            self.array_eyes.append(array_eyes_z)

    def create_array_drops(self):
        for z in range(self.drop_z0, self.drop_zend, self.drop_step_z):
            arr_drops_z = []
            for y in range(self.drop_y0, self.drop_yend, self.drop_step_y):
                arr_drops_y = []
                for x in range(self.drop_x0, self.drop_xend, self.drop_step_x):
                    arr_drops_y.append([x, y, z])
                arr_drops_z.append(arr_drops_y)
            self.array_drops.append(arr_drops_z)

        # rainbow

    def set_number_func_0(self):
        self.number_func = 0

    def set_number_func_1(self):
        self.number_func = 1

    def set_number_func_2(self):
        self.number_func = 2

    def set_number_func_3(self):
        self.number_func = 3

    def set_number_func_4(self):
        self.number_func = 4

    def set_number_func_5(self):
        self.number_func = 5

    def clean_screen(self):
        self.scene.clear()
        self.trans_matrix = [[int(i == j) for i in range(4)] for j in range(4)]
        self.flag_scale_or_not = False
        self.alpha_x = 0
        self.alpha_y = 0
        self.alpha_z = 0

    def set_black_res(self):
        self.pen_res.setColor(QtCore.Qt.black)
        self.color_res = QtCore.Qt.black

    def set_white_res(self):
        self.pen_res.setColor(QtCore.Qt.white)
        self.color_res = QtCore.Qt.white

    def set_blue_res(self):
        self.pen_res.setColor(QtCore.Qt.blue)
        self.color_res = QtCore.Qt.blue

    def set_red_res(self):
        self.pen_res.setColor(QtCore.Qt.red)
        self.color_res = QtCore.Qt.red

    def set_green_res(self):
        self.pen_res.setColor(QtCore.Qt.green)
        self.color_res = QtCore.Qt.green

    def set_yellow_res(self):
        self.pen_res.setColor(QtCore.Qt.yellow)
        self.color_res = QtCore.Qt.yellow

    # Поворот вокруг Х
    def turn_x1(self):
        try:
            alpha = float(self.lineEdit_x_turn.text())
        except Exception:
            QMessageBox.warning(wind, "Внимание!",
                                "Неверное значение угла поворота X")
            return
        self.alpha_x_all += alpha
        self.alpha_y_all += 0
        self.alpha_z_all += 0
        self.draw_res()

    # Поворот вокруг У
    def turn_y1(self):
        try:
            alpha = float(self.lineEdit_y_turn.text())
        except Exception:
            QMessageBox.warning(wind, "Внимание!",
                                "Неверное значение угла поворота Y")
            return
        self.alpha_x_all += 0
        self.alpha_y_all += alpha
        self.alpha_z_all += 0
        self.draw_res()

    # Поворот вокруг Z
    def turn_z1(self):
        try:
            alpha = float(self.lineEdit_z_turn.text())
        except Exception:
            QMessageBox.warning(wind, "Внимание!",
                                "Неверное значение угла поворота Z")
            return
        self.alpha_x_all += 0
        self.alpha_y_all += 0
        self.alpha_z_all += alpha
        self.draw_res()

    # Масштабирование
    def scale(self):
        try:
            sc_k = float(self.lineEdit_k.text())
        except Exception:
            QMessageBox.warning(wind, "Внимание!",
                                "Неверное значение масштаба")
            return

        self.scale_k = sc_k
        self.flag_scale_or_not = True
        self.draw_res()

    # Чтение значений начала конца и шага по оси X и Z
    def read_x_z_value(self):
        # print(self.flag_scale_or_not)
        try:
            x_b = float(self.lineEdit_x_begin.text())
            x_e = float(self.lineEdit_x_end.text())
            x_s = float(self.lineEdit_x_step.text())

            z_b = float(self.lineEdit_z_begin.text())
            z_e = float(self.lineEdit_z_end.text())
            z_s = float(self.lineEdit_z_step.text())

        except Exception:
            QMessageBox.warning(wind, "Внимание!",
                                "Неверное значение параметров X и Z")
            return

        if self.flag_scale_or_not is False:
            f = funcs[self.number_func]
            y2 = max(f(x_b, z_b), f(x_b, z_e), f(x_e, z_b), f(x_e, z_e))
            y1 = min(f(x_b, z_b), f(x_b, z_e), f(x_e, z_b), f(x_e, z_e))
            #print(y1, y2, "y1 and y2")
            dy = y2 - y1
            if abs(y2 - y1) < 1e-6:
                dy = 1

            k1 = int(self.h / dy) - 7
            k2 = int(self.w / (x_e - x_b)) - 7

            #print("SCALE START", k1, k2)
            self.scale_k = min(k1, k2)
        self.x_begin = int(x_b)
        self.x_end = int(x_e)
        self.x_step = x_s

        self.z_begin = int(z_b)
        self.z_end = int(z_e)
        self.z_step = z_s

    def draw_res(self):
        self.scene.clear()
        self.image.fill(QtCore.Qt.white)

        p = QPixmap()

        # отрисовка осей
        drawing_axes(self)
        # отрисовка осей

        # отрисовка глаз (1)
        for z_array_eye in self.array_eyes:
            for y_array_eye in z_array_eye:
                for eye in y_array_eye:
                    x_eye, y_eye, z_buf = transform(eye[0], eye[1], eye[2],\
                                                    self.alpha_x_all, self.alpha_y_all, self.alpha_z_all,\
                                                    self.scale_k, self.w, self.h)
                    drawCircleBrez(self.image, x_eye, y_eye, 1, QtCore.Qt.blue)
        # отрисовка глаз

        # отрисовка капель (2)
        for z_array_drops in self.array_drops:
            for y_array_drops in z_array_drops:
                for drop in y_array_drops:
                    x_drop_2d, y_drop_2d, z_buf = transform(drop[0], drop[1], drop[2],\
                                                    self.alpha_x_all, self.alpha_y_all, self.alpha_z_all,\
                                                    self.scale_k, self.w, self.h)
                    drawCircleBrez(self.image, x_drop_2d, y_drop_2d, 1, QtCore.Qt.green)
        # отрисовка капель

        p.convertFromImage(self.image)
        self.scene.addPixmap(p)

        q = 0.86
        alpha_x = 90
        while q < 0.865:
            for alpha_z in range(80, 101, 5):
                for alpha_y in range(-10, 11, 5):
                    # функция для всех капель
                    func_for_one_drop(self.image, self.scene, self.x_light, self.y_light, self.z_light, alpha_x, alpha_y, alpha_z, q, eye[0], eye[1], eye[2], self.color_res, self.alpha_x_all, self.alpha_y_all, \
                            self.alpha_z_all, self.scale_k, self.w, self.h, max(self.drop_x0, self.drop_xend), max(self.drop_y0, self.drop_yend), max(self.drop_z0, self.drop_zend), self.R_drop, self.array_drops, \
                            min(self.eye_x0, self.eye_xend), max(self.eye_y0, self.eye_yend), max(self.eye_z0, self.eye_zend), self.array_eyes)


            q += 0.005

        #p = QPixmap()
        #p.convertFromImage(self.image)
        #self.scene.addPixmap(p)


def main():
    global wind
    app = QtWidgets.QApplication(sys.argv)
    wind = Visual()
    wind.show()
    app.exec_()


if __name__ == "__main__":
    main()
