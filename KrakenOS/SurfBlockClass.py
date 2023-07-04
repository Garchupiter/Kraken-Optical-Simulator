from typing import Dict, List, Union
from .LensCat import *

class SurfBlock():
    '''
    grouping surface list to use it just like a single optical component
    '''
    count = 0

    def __init__(self, surflist: Union[List, Dict], 
                    name: str, 
                    LastThickness: float=0.0, 
                    LastGlass: str='AIR',
                    inverse: bool=False,
                    **kwargs):
        '''
        init SurfBlock class

        Parameters
        ----------
        surflist : Union[List, Dict]
            surface list to be grouped
        name : str
            optical components name which is also assigned to the last surface
        LastThickness : float, optional
            define spacing between last surface of SurfBlock and next surface to be aligned, by default 0.0
        LastGlass : str, optional
            define material between last surface of SurfBlock and next surface to be aligned, by default 'AIR'
        inverse : bool, optional
            if True, reverse the optical components, by default False
        '''

        self.kwargs = kwargs
        self.DespX = self.kwargs.get('DespX', 0)
        self.DespY = self.kwargs.get('DespY', 0)
        self.DespZ = self.kwargs.get('DespZ', 0)
        self.TiltX = self.kwargs.get('TiltX', 0)
        self.TiltY = self.kwargs.get('TiltY', 0)
        self.TiltZ = self.kwargs.get('TiltZ', 0)
        self.AxisMove = self.kwargs.get('AxisMove', 0)
        self.LastGlass = LastGlass
        self.LastThickness = LastThickness
        self.surflist = surflist
        self.inverse = inverse
        self.name = name

        self.gen_surflist()

    def __call__(self):
        return self.surfblock
    
    def __repr__(self):
        return self.name

    def gen_surflist(self):
        if type(self.surflist)==list:
            # convert surface list to dict
            self.surflist = surflist2dict(self.surflist)

        self.surfblock = cat2surf(cat_dict=self.surflist, Thickness=self.LastThickness, 
                                    Glass=self.LastGlass, inverse=self.inverse,
                                    DespX=self.DespX, DespY=self.DespY, DespZ=self.DespZ,
                                    TiltX=self.TiltX, TiltY=self.TiltY, TiltZ=self.TiltZ,
                                    AxisMove=self.AxisMove)
        
        # assign name to last surface
        self.surfblock[-1].Name = self.name

def alignment(lens_list: List, Distances: Dict={}):
    surf_align = []

    for idx, surfblock in enumerate(lens_list):
        if surfblock.__class__.__name__=='SurfBlock':

            if surfblock()[-1].Name in Distances.keys():
                surfblock()[-1].Thickness = Distances.get(surfblock()[-1].Name, 0)
            elif idx in Distances.keys():
                surfblock()[-1].Thickness = Distances.get(idx, 0)

            
            surf_align += surfblock()
        
        elif surfblock.__class__.__name__=='surf':
            if surfblock.Name in Distances.keys():
                surfblock.Thickness = Distances.get(surfblock.Name, 0)
            elif idx in Distances.keys():
                surfblock.Thickness = Distances.get(idx, 0)
        
            surf_align.append(surfblock)
    return surf_align


