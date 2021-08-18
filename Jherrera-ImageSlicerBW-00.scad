


/////////////////////////////////////////////////
module Corte(A,B){


    r=0.3;
    g=0.1;
    b=0.2;
    

    ax=-atan2((B/2.0),A);
    ay=0;
    az=0;
        
    
    
    rotate(a=[ax,ay,az])
    {
        color( c = [r, g, b], alpha = 0.9 )
        {
            translate([-A/2.0,-B,0])
            {
                cube([2*A,B,A*1.1]);
            }
        }
    }

}
/////////////////////////////////////////////////
module Base(A,B){
    r=0.1;
    g=0.1;
    b=0.3;
    
    ax=0;
    ay=0;
    az=0;
            
    rotate(a=[ax,ay,az])
    {
        color( c = [r, g, b], alpha = 0.9 )
        {
            translate([0,0,0])
            {
                cube([A,B,A]);
            }
        }
    }

}
/////////////////////////////////////////////////
module Nucleo(A,B){
    difference(){
        Base(A,B);
        Corte(A,B);
    }
}
/////////////////////////////////////////////////
module Corte_Superior(A,B){
    rotate(a=[90,0,0])
    {
        linear_extrude(height = B*3, center=true)
        {
            h=1.0;
            polygon(points=[[0-h,A+h],[A+h,A+h],[A+h,0-h],[0-h,A+h]]);
        }
    }
}
///////////////////////////////////////////////////////////////////
module inferior(A,B){
    difference(){
        Nucleo(A,B);
        Corte_Superior(A,B);
    }
}

////////////////////////////////////////////////////////////////////////

module superior(A,B,P){
    difference(){
        translate([0,0,A-0.01])
            {
            rotate(a=[0,45,0])
                {
                cube([sqrt(2*A*A),B,P]);
                }
            }
            translate([0,0,A+P/2.0])
            {
                    cube([sqrt(2*A*A),3*B,P], center=true);
            }
        }
}

///////////////////////////////////////////////////////////////////////
module completo(A,B,P){
    translate([0,-B/2,0])
    {
    union()
    {
        inferior(A,B);
        superior(A,B,P);
    }
    }
}

////////////////////////////////////////////////////////////////////////
module Prisma(pm,B,A){
translate([0,0,A-pm])
    {
        rotate(a=[90,0,0])
        {
            linear_extrude(height = B, center=true)
            {
                polygon(points=[[0,pm],[pm,pm],[pm,0],[0,pm]]);
            }
        }
    }
}
///////////////////////////////////////////////////////////////////

A=10;
B=2;
P=0.5;
pm=1.0;

translate([0,0,-A])
    {
    union()
        {
        completo(A,B,P);
        Prisma(pm,B, A);
        }
    }