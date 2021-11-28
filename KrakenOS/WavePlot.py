
import numpy as np
import matplotlib.pyplot as plt
import KrakenOS as Kos

def ZernikeDataImage2Plot(datos, Type='interferogram'):
    """ZernikeDataImage2Plot.

    Parameters
    ----------
    datos :
        datos
    Type :
        Type
    """
    if (Type == 'interferogram'):


        datos = np.sin(((2 * np.pi) * datos))

    if (Type == 'phase'):
        datos = datos
    img = (datos / np.max(datos))
    img = (img * 255)
    img = img.astype(int)
    img = np.asarray(img)
    plt.figure(987)
    plt.imshow(img, cmap='gray')
    plt.show()
    return 0

def WavefrontData2Image(z_coeff, res=323):
    """WavefrontData2Image.

    Parameters
    ----------
    z_coeff :
        z_coeff
    res :
        res
    """
    TamImag = int(res)
    r = (TamImag / 2.0)
    ARRAY_ZERNIKE = np.zeros((TamImag, TamImag))
    H = []
    K = []
    X = []
    Y = []
    for h in range(0, TamImag):
        for k in range(0, TamImag):
            x = ((h - r) / r)
            y = ((k - r) / r)
            RP = np.sqrt(((x * x) + (y * y)))
            if (RP <= 1):
                H.append(h)
                K.append(k)
                X.append(x)
                Y.append(y)
    H = np.asarray(H)
    K = np.asarray(K)
    X = np.asarray(X)
    Y = np.asarray(Y)
    Z = Kos.Wavefront_Zernike_Phase(X, Y, z_coeff)
    ARRAY_ZERNIKE[(H, K)] = Z
    return ARRAY_ZERNIKE

