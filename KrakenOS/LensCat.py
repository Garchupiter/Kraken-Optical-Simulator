'''
I modified the code copied from below
https://github.com/quartiq/rayopt/blob/master/rayopt/zemax.py
https://github.com/cihologramas/pyoptools/blob/master/pyoptools/raytrace/library/library.py
'''

from pickle import DICT
from struct import Struct
from typing import Dict, List
import numpy as np
from .SurfClass import surf
from .SurfTools import *
import re

# below parameters are unknown or not essential 
zmf_category = ["GCAT",  # glass catalog names
                "OPDX",  # opd
                "RAIM",  # ray aiming
                "CONF",  # configurations
                "ENPD", "PUPD",  # pupil
                "EFFL",  # focal lengths
                "VERS",  # version
                "MODE",  # mode
                "NOTE",  # note
                "TYPE",  # surface type
                "HIDE",  # surface hide
                "MIRR",  # surface is mirror
                "PARM",  # aspheric parameters
                "SQAP",  # square aperture?
                "XDAT", "YDAT",  # xy toroidal data
                "OBNA",  # object na
                "PKUP",  # pickup
                "MAZH", "CLAP", "PPAR", "VPAR", "EDGE", "VCON",
                "UDAD", "USAP", "TOLE", "PFIL", "TCED", "FNUM",
                "TOL", "MNUM", "MOFF", "FTYP", "SDMA", "GFAC",
                "PUSH", "PICB", "ROPD", "PWAV", "POLS", "GLRS",
                "BLNK", "COFN", "NSCD", "GSTD", "DMFS", "ISNA",
                "VDSZ", "ENVD", "ZVDX", "ZVDY", "ZVCX", "ZVCY",
                "ZVAN", "XFLN", "YFLN", "VDXN", "VDYN", "VCXN",
                "VCYN", "VANN", "FWGT", "FWGN", "WWGT", "WWGN",
                "WAVN", "WAVM", "XFLD", "YFLD", "MNCA", "MNEA",
                "MNCG", "MNEG", "MXCA", "MXCG", "RGLA", "TRAC",
                "FLAP", "TCMM", "FLOA", "PMAG", "TOTR", "SLAB",
                "POPS", "COMM", "PZUP", "LANG", "FIMP"]

def zmf2dict(file_list: List[str]) -> Dict:
    '''
    parsing .zmf data 
    Función que lee una librería de Zemax (archivo con terminación zmf), y genera un diccionario con las descripciones
    de cada componente. La llave es la referencia de cada componente

    Parameters
    ----------
    file_list : List[str]
        list of lens catalog file path

    Returns
    -------
    cat_data : Dict
        catalog data
    '''
    cat_data = {}
    for fn in file_list:    
        f = open(fn, "rb")
        head = Struct("<I")
        lens = Struct("<100sIIIIIIIdd")
        shapes = "?EBPM"
        (version,) = head.unpack(f.read(head.size))
        assert version in (1001,)
        while True:
            li = f.read(lens.size)
            if len(li) != lens.size:
                if len(li) > 0:
                    print(f, "additional data", repr(li))
                break
            li = list(lens.unpack(li))
            li[0] = li[0].decode("latin1").strip("\0")
            li[3] = shapes[li[3]]
            description = f.read(li[7])
            assert len(description) == li[7]
            description = zmf_obfuscate(description, li[8], li[9])
            description = description.decode("latin1")
            assert description.startswith("VERS {:06d}\n".format(li[1]))
            cat_data[li[0]] = zmf_parsing(description)
    return cat_data

def zmf_obfuscate(data: str, a, b):
    iv = np.cos(6 * a + 3 * b)
    iv = np.cos(655 * (np.pi / 180) * iv) + iv
    p = np.arange(len(data))
    k = 13.2 * (iv + np.sin(17 * (p + 3))) * (p + 1)
    k = (int(("{:.8e}".format(_))[4:7]) for _ in k)
    data = np.fromstring(data, np.uint8)
    data ^= np.fromiter(k, np.uint8, len(data))
    return data.tostring()

def zmf_parsing(data: str):
    lens_info = {}
    for line in data.splitlines():
        if not line.strip():
            continue
        line = line.strip().split(" ", 1)
        cmd = line[0]
        args = len(line) == 2 and line[1] or ""
        if cmd == "UNIT":
            unit = {"MM": 1, # basic scale mm
                    "INCH": 25.4, # in mm
                    "IN": 25.4, # in mm
                    }[args.split()[0]]
        elif cmd == "NAME":
            lens_info['description'] = args.strip("\"")
        elif cmd == "SURF":
            current_surf = f'SUFR {args}'
            lens_info[current_surf] = {}
        elif cmd == "CURV":
            if float(args.split()[0])*unit == 0:
                lens_info[current_surf]['Rc'] = 0
            else:
                lens_info[current_surf]['Rc'] = 1 / (float(args.split()[0])*unit)

        elif cmd == "DISZ":
            if float(args)==float('INFINITY'):
                lens_info[current_surf]['Thickness'] = 0
            else:
                lens_info[current_surf]['Thickness'] = float(args)*unit
        elif cmd == "GLAS":
            if '___BLANK' in args:
                # for unknown material
                # reference for the zemax glass data format
                # https://github.com/nzhagen/zemaxglass/blob/master/ZemaxGlass_user_manual.pdf
                # NM <glass_name> <dispersion_formula_number> <MIL> <N(d)> <V(d)> <Exclude_sub> <status> <melf_freq> 
                # format
                # ___BLANK 1 0 <refractive index> <Abbe number> 0 0 0 0 0 0 
                # ex) ___BLANK 1 0 1.52216 5.88E+1 0 0 0 0 0 0 
                args = args.replace(' ', ',')
                lens_info[current_surf]['Glass'] = args
                args = args.split(',')

                lens_info[current_surf]['Refractive_index'] = float(args[3])
                lens_info[current_surf]['Abbe_num'] = float(args[4])

            else:
                args = args.split()
                lens_info[current_surf]['Glass'] = args[0]

        elif cmd == "DIAM":
            lens_info[current_surf]['Diameter'] = float(args.split()[0])*2*unit
        elif cmd == "STOP":
            lens_info[current_surf]['STOP'] = True
        elif cmd == "WAVL":
            lens_info['wavelengths'] = [float(i)*1e-6 for i in args.split() if i]
        elif cmd == "COAT":
            lens_info[current_surf]['coating'] = args.split()[0]
        elif cmd == "CONI":
            lens_info[current_surf]['conic'] = float(args.split()[0])
        elif cmd == "PARM":
            i, j = args.split()
            i = int(i) - 1
            j = float(j)
            if i < 0:
                if j:
                    # print("aspheric 0 degree not supported", cmd, args)
                    pass
                continue

            if lens_info[current_surf].get('aspherics') is None:
                lens_info[current_surf]['aspherics'] = []
            while len(lens_info[current_surf]['aspherics']) <= i:
                lens_info[current_surf]['aspherics'].append(0.)
            
            lens_info[current_surf]['aspherics'][i] = j
        else:
            pass
    return lens_info    
    
def surflist2dict(surflist: List) -> Dict:
    '''
    Convert surface list to dictionary form

    Parameters
    ----------
    surflist : List
        surface list

    Returns
    -------
    Dict
        surface list in dictionary form
    '''
    surfdict = {}

    for idx, surf in enumerate(surflist):
        surfdict[f'SUFR {idx}'] = {'Rc' : surf.__dict__.get('Rc', 0),
                                    'Thickness' : surf.__dict__.get('Thickness', 0),
                                    'Glass' : surf.__dict__.get('Glass', 'AIR'),
                                    'Diameter' : surf.__dict__.get('Diameter', 0),
                                    'conic' : surf.__dict__.get('k', 0),
                                    'aspherics' : surf.__dict__.get('AspherData', [0]*200)}
    return surfdict

def cat2surf(cat_dict: DICT, Thickness: float=0, Glass: str='AIR', inverse: bool=False,
            DespX=0.0, DespY=0.0, DespZ=0.0, 
            TiltX=0.0, TiltY=0.0, TiltZ=0.0, AxisMove=0.0):
    '''
    Convert lens catalog dict to surface list

    Parameters
    ----------
    cat_dict : DICT
        Parsed lens catalog dict
    Thickness : float, optional
        Last thickness which is space before next surface, by default 0
    Glass : str, optional
        refractive index or glass type of last spacinc, by default 'AIR'
    inverse : bool, optional
        determine whether reverse or not a lens block, by default False
    Desp[X, Y, Z] : float, optional
        Displacement in [X, Y, Z] axis, respectively, by default 0.0
    Tilt[X, Y, Z] : float, optional
        Tiliting in [X, Y, Z] axis, respectively, by default 0.0
    AxisMove : float, optional
        Defines what will happen to the optical axis after a coordinate transformation. 
        If the value is 0, the transformation is only carried out to the surface in question. 
        If the value is 1 then the transformation also affects the optical axis. 
        Therefore, the other surfaces will follow the transformation. 
        If the value is different, for example 2, then the optical axis will be affected twice. 
        , by default 0.0

    Returns
    -------
    _type_
        surface list
    '''
    # get surface
    surf_name = [surface for surface in cat_dict.keys() if ('SUFR' in surface)]

    # find maximum diameter
    max_diameter = max([cat_dict[surf].get('Diameter', 0) for surf in surf_name])

    # exclude surface whose diameter is different with maximum diameter
    surf_name = [surface for surface in surf_name if cat_dict[surface].get('Diameter') not in [None, 0]]

    surf_list = []
    if inverse:
        surf_name = surf_name[::-1]
        
        for idx, surface in enumerate(surf_name):
            if int(AxisMove)==0:
                sf = surf(DespX=DespX, DespY=DespY, DespZ=DespZ, AxisMove=AxisMove)
            elif int(AxisMove)==1:
                if idx==0:
                    sf = surf(DespX=DespX, DespY=DespY, DespZ=DespZ,
                        TiltX=TiltX, TiltY=TiltY, TiltZ=TiltZ,
                            Order=0, AxisMove=AxisMove)
                else:
                    sf = surf()

            current_surf = cat_dict[surface]
            
            # current value => surface properties
            sf.Rc = -current_surf.get('Rc', 0)
            sf.Diameter = current_surf.get('Diameter', max_diameter) # default 25 mm
            sf.k = current_surf.get('conic', 0)
            sf.AspherData = -1*np.array(current_surf.get('aspherics', [0]*200))
            
            # next value => properties of material
            if idx!=len(surf_name)-1:
                next_surf = cat_dict[surf_name[idx+1]]
                sf.Thickness = next_surf.get('Thickness', 0)
                sf.Glass = next_surf.get('Glass', 'AIR')

            else:
                sf.Thickness = 0
                sf.Glass = Glass
            surf_list.append(sf)
            
    else:
        for idx, surface in enumerate(surf_name):
            if int(AxisMove)==0:
                sf = surf(DespX=DespX, DespY=DespY, DespZ=DespZ, AxisMove=AxisMove)
            elif int(AxisMove)==1:
                if idx==0:
                    sf = surf(DespX=DespX, DespY=DespY, DespZ=DespZ,
                        TiltX=TiltX, TiltY=TiltY, TiltZ=TiltZ,
                            Order=0, AxisMove=AxisMove)
                else:
                    sf = surf()

            current_surf = cat_dict[surface]
            sf.Rc = current_surf.get('Rc', 0)
            sf.Thickness = current_surf.get('Thickness', 0)

            sf.Diameter = current_surf.get('Diameter', max_diameter)
            sf.k = current_surf.get('conic', 0)
            sf.AspherData = np.array(current_surf.get('aspherics', [0]*200))

            sf.Glass = current_surf.get('Glass', 'AIR')
            
            if idx==len(surf_name)-1:
                sf.Thickness = 0
                sf.Glass = Glass
                
            surf_list.append(sf)

    if int(AxisMove)!=0:
        surf_list[-1].AxisMove = AxisMove

    return surf_list


# def zmx_parse(data):
#     """Función que lee e interpreta un archivo zmx de zemax.
#     - Solo funciona para lentes esfericas y dobletes.
#     - No tiene ne cuenta recubrimientos antireflectivos
#     - Asume que las medidas están en milimetros. Hay que arreglar esto
#     """

#     lines = data.splitlines()

#     # Interpretar el encabezado
#     while True:
#         line = lines.pop(0)

#         if lines[0].startswith("SURF"):
#             break

#     # Separar las superficies en una lista de diccionarios
#     surflist = []
#     for line in lines:
#         if line.startswith("SURF"):
#             surflist.append(dict())
#             continue
#         line = line.lstrip()
#         code = line[:4]
#         data = line[5:].split()
#         data = [convert(d) for d in data]
#         surflist[-1][code] = data

#     # Eliminar el plano objeto y el plano imagen

#     surflist.pop(0)
#     surflist.pop()

#     # Identificar el tipo de lentes a partir de el numero de superficies
#     # validas

#     ns = len(surflist)
#     return surflist

# def convert(d):
#     try:
#         return int(d)
#     except ValueError:
#         try:
#             return float(d)
#         except ValueError:
#             return d

# def zmx_read(fn):

#     with open(fn, 'r', encoding='utf16') as f:
#         data = f.read()
#     return zmx_parse(data)