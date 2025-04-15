import pkg_resources


import KrakenOS as Kos


import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import pickle
import math
import itertools as IT
import sys

class PSFMap:
    def __init__(self, Sys, fieldxs=[0, 20], fieldys=[0, 20], surf=7, aperVal=1.80, stepx=1, stepy=1, Wx = 0.589, aperType="STOP"):
        self.Surf = surf # Defining M1 as enter pupil diameter
        self.AperVal = aperVal
        self.AperType = aperType
        self.W = 0.486
        self.Sys = Sys  # System
        self.RayosList = []
        self.PupList = []
        self.PSFList = []
        self.PSFPupList = []
        self.coordList = []
        self.Wx = Wx
        self.ray_n = 0
        self.stepx, self.stepy = stepx, stepy
        self.fieldxs = fieldxs
        self.fieldys = fieldys
        self.POZList = []
        self.CR_pSourceList = []
        self.CR_ray_index = 0
        self.SpotMapRays = []
        self.Sys.IgnoreVignetting()
        self.FieldXList = []
        self.FieldYList = []
    
    def wide_angle_ray_trace(self, Rayos, x, y, z, L, M, N):
        '''
        wide_angle_ray_trace: Traces the rays for the given x, y, z, L, M, N
        Parameters:
            x: x values
            y: y values
            z: z values
            L: L values
            M: M values
            '''
        pSource_0 = [x[0], y[0], z[0]]
        dCos=[L[0], M[0], N[0]]
        self.Sys.Trace(pSource_0, dCos, self.Wx)
        Rayos.push()
        return Rayos

    def calculate_wide_angle_pupil(self, fieldx, fieldy, sampling):
        '''
        calculate_wide_angle_pupil: Calculates the pupil for the given fieldx and fieldy
        Parameters:
            fieldx: fieldx value
            fieldy: fieldy value
        '''
        print("Create pupil with values: ", "Surf", self.Surf, "Wx", self.Wx,
              "AperType", self.AperType, "AperVal", self.AperVal, "FieldX", fieldx, "FieldY", fieldy)
        Pup = Kos.PupilCalc(self.Sys, self.Surf, self.Wx, self.AperType, self.AperVal)
        Pup.Samp=sampling
        Pup.FieldType = "angle"
        Pup.Ptype = "chief"
        self.POZList.append(Kos.precise_det_exit_pupil(self.Sys, Pup, self.Wx))
        self.FieldXList.append(fieldx)
        self.FieldYList.append(fieldy)
        PSFFieldX = 0
        Pup.FieldX = fieldx
        Pup.FieldY = fieldy

        Pup.Pattern()

        X, Y, Z, L, M, N = Pup.Pattern2FieldPlus()
        self.CR_pSourceList.append([X[0], Y[0], Z[0]])
        Rayos = Kos.raykeeper(self.Sys)

        Pup.Ptype = "rtheta"

        Pup.rad = 1
        Pup.theta = 0
        Pup.Pattern()
        X0, Y0, Z0, L0, M0, N0 = Pup.Pattern2FieldPlus()

        Pup.theta = 90
        Pup.Pattern()
        X1, Y1, Z1, L1, M1, N1 = Pup.Pattern2FieldPlus()

        Pup.theta = 180
        Pup.Pattern()
        X2, Y2, Z2, L2, M2, N2 = Pup.Pattern2FieldPlus()


        Pup.theta = 270
        Pup.Pattern()
        X3, Y3, Z3, L3, M3, N3 = Pup.Pattern2FieldPlus()

        Rayos = self.wide_angle_ray_trace(Rayos, X, Y, Z, L, M, N)
        Rayos = self.wide_angle_ray_trace(Rayos, X0, Y0, Z0, L0, M0, N0)
        Rayos = self.wide_angle_ray_trace(Rayos, X1, Y1, Z1, L1, M1, N1)
        Rayos = self.wide_angle_ray_trace(Rayos, X2, Y2, Z2, L2, M2, N2)
        CRRay = self.wide_angle_ray_trace(Rayos, X3, Y3, Z3, L3, M3, N3)

        self.RayosList.append(CRRay)
        self.ray_n = self.ray_n + 1
        self.PupList.append(Pup)
    
    def calc_wide_angle_pupil_list(self, sampling=3, spotdiagram=False):
        '''
        calc_wide_angle_pupil_list: Calculates the pupil list for the given fieldxs and fieldys

        Parameters:
            sampling: sampling for the pupil
            spotdiagram: if spotdiagram is True, it will generate the spot diagram

        Expected format for fieldxs and fieldys:
        fieldxs = [-FOV, ..., 0, ..., +FOV]
        fieldys = [-FOV, ..., 0, ..., +FOV]

        or

        step = #fixed step size
        fieldxs = np.arange(-FOV, +FOV, step)
        fieldys = np.arange(-FOV, +FOV, step)
        '''
        for fieldx in self.fieldxs:
            for fieldy in self.fieldys:
                self.calculate_wide_angle_pupil(fieldx, fieldy, sampling)
                if spotdiagram == True:
                    self.generate_spot_map_pupil(fieldx, fieldy) #choose to calculate spotdiagram for each fieldx, fieldy
    
    def calc_wide_angle_pupil_list_scanline(self, sampling=3, spotdiagram=False):
        '''
        calc_wide_angle_pupil_list: Calculates the pupil list for the given fieldxs and fieldys using scanline method

        Parameters:
            sampling: sampling for the pupil
            spotdiagram: if spotdiagram is True, it will generate the spot diagram

        Expected format for fieldxs and fieldys:
        fieldxs = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]
        fieldys = [y1, y2, y3, y4, y5, y6, y7, y8, y9, y10]

        Data Type: A List or a numpy array of floats or integers

        fieldxs and fiedys should be read from a CSV file.
        Below is an example format with the following Assumptions:
        1. The CSV file has two columns: FieldX and FieldY

        FieldX | FieldY
        ----------------
        x1     | y1
        x2     | y2
        x3     | y3
        .     | .
        .     | .
        .     | .
        xn     | yn
        '''

        #For scanline method, we need to calculate the fieldx and fieldy values
        #for each fieldx and fieldy simultaneously
        for fieldx, fieldy in zip(self.fieldxs, self.fieldys):
            self.calculate_wide_angle_pupil(fieldx, fieldy, sampling)
            if spotdiagram == True:
                self.generate_spot_map_pupil(fieldx, fieldy)

    def plot_spot_map(self, show=1, save=1):
        '''
        plot_spot_map: Plots the spot map for the given fieldxs and fieldys
        
        Parameters:
            show: if show is 0, it will not show the plot
            save: if save is 0, it will not save the plot
        '''

        f, ax = plt.subplots(figsize=(15,15))
        for i in range(0, self.ray_n):
            X,Y,Z,L,M,N = self.RayosList[i].pick(-1)
            if self.Wx == 0.486:
                ax.plot(X,Y, 'x', c="b")
            if self.Wx == 0.531:
                ax.plot(X,Y, 's', c="g")
            if self.Wx == 0.656:
                ax.plot(X,Y, '^', c="r")
        plt.xlabel('X Axis', axes=ax)
        plt.ylabel('Y Axis', axes=ax)
        plt.title('Spot Diagram')
        plt.axis('square')
        if show == 1:
            plt.show()
        if save == 1:
            f.savefig('SpotDiagram' + str(self.Wx) + '.png')
            pickle.dump(f, open('SpotDiagram' + str(self.Wx) + '.pkl', 'wb'))
        return f
    
    def generate_spot_map_pupil(self, fieldx, fieldy, sampling=10, show=1, save=1):
        SptRay = Kos.raykeeper(self.Sys)
        print("Generate Spot Diagram at (FieldX: ", fieldx, "FieldY: ", fieldy, ')')
        Pup = Kos.PupilCalc(self.Sys, self.Surf, self.Wx, self.AperType, self.AperVal)
        Pup.Samp=sampling
        Pup.FieldType = "angle"
        Pup.Ptype = "hexapolar"
        Pup.FieldX = fieldx
        Pup.FieldY = fieldy
        Pup.Pattern()
        X,Y,Z,L,M,N=Pup.Pattern2FieldPlus()

        for i in range(0,len(X)):
            pSource_0 = [X[i], Y[i], Z[i]]
            dCos=[L[i], M[i], N[i]]
            self.Sys.Trace(pSource_0, dCos, self.Wx)
            SptRay.push()
        
        self.SpotMapRays.append(SptRay)
        self.plot_wide_angle_spot_map(show, save)

    def plot_wide_angle_spot_map(self, show=1, save=1):
        f, ax = plt.subplots(figsize=(15,15))
        for i in range(0, self.ray_n):
            try:
                X,Y,Z,L,M,N = self.SpotMapRays[i].pick(-1)
                if self.Wx == 0.486:
                    ax.plot(X,Y, 'x', c="b")
                if self.Wx == 0.531:
                    ax.plot(X,Y, 's', c="g")
                if self.Wx == 0.656:
                    ax.plot(X,Y, '^', c="r")
            except Exception as e:
                print("***Exception*** ", e)
                continue
        plt.xlabel('X Axis', axes=ax)
        plt.ylabel('Y Axis', axes=ax)
        plt.title('Spot Diagram')
        plt.axis('square')
        if show == 1:
            plt.show()
        if save == 1:
            f.savefig('SpotDiagram' + str(self.Wx) + '.png')
            pickle.dump(f, open('SpotDiagram' + str(self.Wx) + '.pkl', 'wb'))
        return f

    def generate_psf_map(self, pixels=144, PupilSample=4, plot=0, sqr=0, wideAngle=False, NC=28, cmap=plt.cm.jet, downsample=False):
        '''
        generate_psf_map: Generates a PSF map for the given fieldxs and fieldys
        Parameters:
            pixels: number of pixels for the PSF
            PupilSample: Pupil sample for the PSF map
            plot: if plot is 0, it will not plot the PSF map
        '''
        pixel_res = pixels

        if plot != 0:
            fig = plt.figure(figsize=(15,15))
            columns = int(np.sqrt(len(self.PupList)))
            rows = int(len(self.PupList)/columns)
        for i in range(len(self.PupList)):
            pixels = pixel_res # Reset the pixels to the original value at the beginning of each iteration
            Pup_i = self.PupList[i]
            fieldx_i = self.FieldXList[i]
            fieldy_i = self.FieldYList[i]
            x_shift = abs(fieldx_i)/self.stepx
            y_shift = abs(fieldy_i)/self.stepy

            print("i: ", i, "FieldX: ", fieldx_i, "FieldY: ", fieldy_i, "x_shift: ", x_shift, "y_shift: ", y_shift)
            if wideAngle == False:
                X, Y, Z, P2V = Kos.Phase(Pup_i)
                A = np.ones(NC)
                Zcoef, Mat, RMS2Chief, RMS2Centroid, FITTINGERROR = Kos.Zernike_Fitting(X, Y, Z, A)
            else:
                #print(self.RayosList[i].pick(0))
                NPx, NPy, NPz, L, M, N, Pup_i = Kos.RequestRays(self.RayosList[i], Pup_i)
                P2V, Wi, ERX, ERY, SYSTEM, RAYS = Kos.CreateRefSphere(self.RayosList[i],
                                                                        Pup_i,
                                                                        self.Sys,
                                                                        self.POZList[i],
                                                                        NPx,
                                                                        NPy,
                                                                        NPz,
                                                                        L,
                                                                        M,
                                                                        N,
                                                                        self.Wx,
                                                                        self.CR_pSourceList[i])

                A = np.ones(NC)
                Zcoef, Mat, RMS2Chief, RMS2Centroid, FITTINGERROR = Kos.Zernike_Fitting(ERY, ERX, Wi, A)

            COEF = Zcoef
            Focal = Pup_i.EFFL
            #Diameter was experimentally determined compared OpticStudio 
            if self.Wx == 0.486:
                Diameter = 2.5 * Pup_i.RadPupInp
            elif self.Wx == 0.531:
                Diameter = 2.80 * Pup_i.RadPupInp
            else:
                Diameter = 3.45 * Pup_i.RadPupInp # DJ - red wavelength creates a smaller PSF diameter than the others.
            #if group plots are disabled, plot the um plots of the PSFs
            if plot == 0:
                I, param = Kos.psf(COEF, Focal, Diameter, self.Wx, PupilSample=PupilSample, pixels=pixels, plot=0, sqr=sqr)
            elif plot == 1:
                I, param = Kos.psf(COEF, Focal, Diameter, self.Wx, PupilSample=PupilSample, pixels=pixels, plot=0, sqr=sqr)
            else:
                I, param = Kos.psf(COEF, Focal, Diameter, self.Wx, PupilSample=PupilSample, pixels=pixels, plot=1, sqr=sqr, cmap=cmap)
            
            if downsample==True:
                I = Image.fromarray(I).resize((9, 9), Image.BILINEAR)
                I = np.array(I)
                pixels = 9 # We downsample the pixel resolution to 9x9

            if fieldx_i == 0 and fieldy_i == 0:
                x_min = -pixels/2
                x_max = pixels/2
                y_min = -pixels/2
                y_max = pixels/2
            elif fieldx_i == 0 and fieldy_i > 0:
                x_min = -pixels/2
                x_max = pixels/2
                y_min = -pixels/2 + (pixels*y_shift)
                y_max = pixels/2 + (pixels*y_shift)
            elif fieldx_i == 0 and fieldy_i < 0:
                x_min = -pixels/2
                x_max = pixels/2
                y_min = -pixels/2 - (pixels*y_shift)
                y_max = pixels/2 - (pixels*y_shift)
            elif fieldx_i > 0 and fieldy_i == 0:
                x_min = -pixels/2 + (pixels*x_shift)
                x_max = pixels/2 + (pixels*x_shift)
                y_min = -pixels/2
                y_max = pixels/2
            elif fieldx_i > 0 and fieldy_i > 0:
                x_min = -pixels/2 + (pixels*x_shift)
                x_max = pixels/2 + (pixels*x_shift)
                y_min = -pixels/2 + (pixels*y_shift)
                y_max = pixels/2 + (pixels*y_shift)
            elif fieldx_i > 0 and fieldy_i < 0:
                x_min = -pixels/2 + (pixels*x_shift)
                x_max = pixels/2 + (pixels*x_shift)
                y_min = -pixels/2 - (pixels*y_shift)
                y_max = pixels/2 - (pixels*y_shift)
            elif fieldx_i < 0 and fieldy_i == 0:
                x_min = -pixels/2 - (pixels*x_shift)
                x_max = pixels/2 - (pixels*x_shift)
                y_min = -pixels/2
                y_max = pixels/2
            elif fieldx_i < 0 and fieldy_i > 0:
                x_min = -pixels/2 - (pixels*x_shift)
                x_max = pixels/2 - (pixels*x_shift)
                y_min = -pixels/2 + (pixels*y_shift)
                y_max = pixels/2 + (pixels*y_shift)
            elif fieldx_i < 0 and fieldy_i < 0:
                x_min = -pixels/2 - (pixels*x_shift)
                x_max = pixels/2 - (pixels*x_shift)
                y_min = -pixels/2 - (pixels*y_shift)
                y_max = pixels/2 - (pixels*y_shift)
            
            self.PSFList.append(I)
            self.PSFPupList.append(Pup_i)
            self.coordList.append([x_min, x_max, y_min, y_max])
            if plot == 1:
                fig.add_subplot(rows, columns, (i+1))
                fig.tight_layout()
                im = plt.imshow(I, extent=param, cmap=cmap)
                plt.title('(' + str(fieldx_i) + ','+ str(fieldy_i) + ')')
                plt.ylabel('V[μm]')
                plt.xlabel('U[μm]')
                plt.plot()
        numpyPSFs = np.array(self.PSFList)
        np.save('numpyPSFs'+str(self.Wx)+'.npy', numpyPSFs)
        if plot == 1:
            fig.subplots_adjust(right=0.8)
            cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
            fig.colorbar(im, cbar_ax)
            return fig, columns, rows, numpyPSFs
        else:
            return numpyPSFs
        
    def merge_psf_map(self):
        '''
        merge_psf_map: Merges the PSF map for the given fieldxs and fieldys
        '''
        xmin = 0
        xmax = 0
        ymin = 0
        ymax = 0
        for x_min_i, x_max_i, y_min_i, y_max_i in self.coordList:
            if xmin > x_min_i:
                xmin = x_min_i
            if xmax < x_max_i:
                xmax = x_max_i
            if ymin > y_min_i:
                ymin = y_min_i
            if ymax < y_max_i:
                ymax = y_max_i
        im = Image.new('RGB', (int(xmax-xmin), int(ymax-ymin)))

        print("xmin: ", xmin, "xmax: ", xmax, "ymin: ", ymin, "ymax: ", ymax)

        for i in range(len(self.PSFList)):
            x_min_i, x_max_i, y_min_i, y_max_i = self.coordList[i]
            fieldx_i = self.PSFPupList[i].FieldX
            fieldy_i = self.PSFPupList[i].FieldY
            x_min_i = int(x_min_i + abs(xmin))
            x_max_i = int(x_max_i + abs(xmax))
            y_min_i = int(y_min_i + abs(ymin))
            y_max_i = int(y_max_i + abs(ymax))
            I = self.PSFList[i]
            print("i: ", i, "FieldX: ", fieldx_i, "FieldY: ", fieldy_i, "x_min: ", x_min_i, "x_max: ", x_max_i, "y_min: ", y_min_i, "y_max: ", y_max_i)
            I = I / np.max(I)
            I = np.uint8(255 * I)
            I = np.stack((I, I, I), axis=-1)
            I = Image.fromarray(I)

            if fieldx_i > 0 and fieldy_i > 0 and abs(fieldx_i) == abs(fieldy_i):
                I = I.rotate(-90)
            elif fieldx_i > 0 and fieldy_i < 0 and abs(fieldx_i) == abs(fieldy_i):
                I = I.rotate(90)
            elif fieldx_i < 0 and fieldy_i > 0 and abs(fieldx_i) == abs(fieldy_i):
                I = I.rotate(90)
            elif fieldx_i < 0 and fieldy_i < 0 and abs(fieldx_i) == abs(fieldy_i):
                I = I.rotate(-90)
            elif fieldx_i == 0 and fieldy_i > 0:
                I = I.rotate(180)
            elif fieldx_i == 0 and fieldy_i < 0:
                I = I.rotate(180)
            else:
                I = I.transpose(Image.FLIP_TOP_BOTTOM)
            print(I.size)
            x_sz, y_sz = I.size
            x_max_i = x_min_i + x_sz
            y_max_i = y_min_i + y_sz
            im.paste(I, (x_min_i, y_min_i, x_max_i, y_max_i))
        out = im.transpose(Image.FLIP_TOP_BOTTOM)
        return out
    
    def crop_to_res(self, im_in, new_width, new_height):
        '''
        crop_to_res: Crops the PSF map to the resolution
        '''
        width, height = im_in.size   # Get dimensions


        left = (width - new_width)/2
        top = (height - new_height)/2
        right = (width + new_width)/2
        bottom = (height + new_height)/2

        # Crop the center of the image
        return im_in.crop((left, top, right, bottom))
        
