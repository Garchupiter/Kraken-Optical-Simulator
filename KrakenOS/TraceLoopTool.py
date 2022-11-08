def TraceLoop(x, y, z, L, M, N, W, Container, clean = 1):
    """TraceLoop.

    Parameters
    ----------
    x :
        x
    y :
        y
    z :
        z
    L :
        L
    M :
        M
    N :
        N
    W :
        W
    Container :
        Container
    """
    if clean ==1:
        Container.clean()
    System = Container.SYSTEM
    for i in range(0, len(x)):
        pSource_0 = [x[i], y[i], z[i]]
        dCos = [L[i], M[i], N[i]]
        System.Trace(pSource_0, dCos, W)
        Container.push()
    return 0

def NsTraceLoop(x, y, z, L, M, N, W, Container, clean = 1):
    """NsTraceLoop.

    Parameters
    ----------
    x :
        x
    y :
        y
    z :
        z
    L :
        L
    M :
        M
    N :
        N
    W :
        W
    Container :
        Container
    """
    if clean ==1:
        Container.clean()
    System = Container.SYSTEM
    for i in range(0, len(x)):
        pSource_0 = [x[i], y[i], z[i]]
        dCos = [L[i], M[i], N[i]]
        System.NsTrace(pSource_0, dCos, W)
        Container.push()
    return 0

