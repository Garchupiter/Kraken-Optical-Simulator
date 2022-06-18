rotate(a=[0,-90,0])
    {
    linear_extrude(height = 10, center=true)
        {    
        polygon(points=[[0,10],[10,0],[0,-10],[0,10]]);
        }
    }