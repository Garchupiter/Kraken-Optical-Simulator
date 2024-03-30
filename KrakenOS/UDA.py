#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  3 22:44:18 2023

@author: joelherreravazquez
"""

import numpy as np
import matplotlib.path as mpl_path
# import matplotlib.pyplot as plt

import pyvista as pv
# from matplotlib.path import Path

#############################################################################


def interpolar_perimetro(px, py, num_puntos_interpolados):
    # Combinar las listas de coordenadas en un solo arreglo numpy
    coordenadas = np.column_stack((px, py))

    # Inicializar listas para las coordenadas interpoladas
    x_interpolado = []
    y_interpolado = []

    # Iterar a través de los puntos originales y agregarlos a las listas interpoladas
    for i in range(len(px) - 1):
        x1, y1 = coordenadas[i]
        x2, y2 = coordenadas[i + 1]

        # Interpolar puntos entre (x1, y1) y (x2, y2)
        for j in range(num_puntos_interpolados + 1):
            alpha = j / num_puntos_interpolados
            x_interp = x1 + alpha * (x2 - x1)
            y_interp = y1 + alpha * (y2 - y1)

            x_interpolado.append(x_interp)
            y_interpolado.append(y_interp)

    # Agregar el último punto para cerrar el polígono
    x_interpolado.append(px[-1])
    y_interpolado.append(py[-1])

    return x_interpolado, y_interpolado


#############################################################################


def calcular_centroide(vertices):
    # las coordenadas de los vértices son un arreglo numpy
    vertices = np.array(vertices)

    # Calcular la suma de las coordenadas de los vértices
    suma_coord = np.sum(vertices, axis=0)

    # Calcular el promedio dividiendo por la cantidad de vértices
    centroide = suma_coord / len(vertices)

    return centroide


#############################################################################


class UDA():
    """Hit_Solver.
    """

    def __init__(self, x_poligono, y_poligono):

        # Crear un objeto de ruta a partir de las coordenadas del polígono
        self.ruta_poligono = mpl_path.Path(np.column_stack((x_poligono, y_poligono)))

        ppx = []
        ppy = []

        ppxx = []
        ppyy = []

        filas = 100
        columnas = 100


        # Crear coordenadas x e y utilizando np.linspace
        x = np.linspace(np.min(x_poligono), np.max(x_poligono), columnas)
        y = np.linspace(np.min(y_poligono), np.max(y_poligono), filas)

        # Crear el grid combinando las coordenadas x e y utilizando np.meshgrid
        xx, yy = np.meshgrid(x, y)

        px = xx.ravel()
        py = yy.ravel()


        # Verificar si los puntos están dentro del polígono
        puntos_dentro = self.ruta_poligono.contains_points(np.column_stack((px, py)))

        # Imprimir los resultados
        for i, punto_dentro in enumerate(puntos_dentro):
            if punto_dentro:
                ppx.append(px[i])
                ppy.append(py[i])
            #     print(f'El punto ({px[i]}, {py[i]}) está dentro del polígono.')
            else:
                ppxx.append(px[i])
                ppyy.append(py[i])
            #     print(f'El punto ({px[i]}, {py[i]}) está fuera del polígono.')

        # Visualización del polígono y los puntos

        num_puntos_interpolados = 100
        x_poligono, y_poligono = interpolar_perimetro(x_poligono, y_poligono, num_puntos_interpolados)

        ppx = ppx + x_poligono
        ppy = ppy + y_poligono

        ppx = np.asarray(ppx)
        ppy = np.asarray(ppy)
        ppz = np.zeros_like(ppy)

        ppxx = np.asarray(ppxx)
        ppyy = np.asarray(ppyy)

        puntos = np.column_stack((ppx, ppy, ppz))

        nube_puntos = pv.PolyData(puntos)

        # Genera la malla a partir de la nube de puntos utilizando Delaunay 3D
        malla = nube_puntos.delaunay_2d()

        # Índice de la celda que deseas eliminar

        N = malla.n_cells
        X = []
        Y = []
                
        # print("Atributos y métodos:", dir(malla))
        
        # # Obtener solo los atributos
        # print("Atributos:", vars(malla))
        
        # Obtener solo los métodos
        # metodos = [attr for attr in dir(malla) if callable(getattr(malla, attr))]
        # print("Métodos:", metodos)
        
        for i in range(N):
            try:
                celda = malla.get_cell(i)
            except:
                celda = malla.extract_cells(i)
            puntos_celda = celda.points
            [x,y,z] = calcular_centroide(puntos_celda)
            X.append(x)
            Y.append(y)

        X = np.asarray(X)
        Y = np.asarray(Y)

        puntos_dentro = self.ruta_poligono.contains_points(np.column_stack((X, Y)))

        indices_celdas_a_mantener = []
        for i in range(N):
            if puntos_dentro[i] == True:
                indices_celdas_a_mantener.append(i)


        # Crea una nueva malla que contiene solo las celdas que deseas mantener
        malla_actualizada = malla.extract_cells(indices_celdas_a_mantener)
        try:
            self.UDA_Surf = malla_actualizada.clean()
        except:
            self.UDA_Surf = malla_actualizada



    def Hit(self, x, y):
        stat = self.ruta_poligono.contains_points([[x, y]])
        return stat[0]


###############################################################################





# aa = 2
# bb = 2
# # Coordenadas del polígono cerrado (x, y)
# px = [1, 3, 2.5, 4, 2, 1,  .1+aa, .3+aa, .4+aa, .2+aa, .1+aa]
# py = [1, 1, 2.5, 3, 4, 1,  .1+bb, .1+bb, .3+bb, .4+bb, .1+bb]



# Poligono = UDA(px, py)

# objeto = Poligono.UDA_Surf
# # Visualiza el polígono y su perímetro
# p = pv.Plotter()

# p.add_mesh(objeto, color="red", line_width=5, label="Perímetro Exterior")

# p.show()
