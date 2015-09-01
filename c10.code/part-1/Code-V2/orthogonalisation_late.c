#include "global.h"

double plijk(int r,int i,int j,int k)
{
    int p1,p2,p3;
    find_perm_late(r,&p1,&p2,&p3);
    
    double b1,b2,b3,b4,b5,b6;
    b1 = get_basis_late(i,p1)*get_basis_late(j,p2)*get_basis_late(k,p3);
    b2 = get_basis_late(i,p2)*get_basis_late(j,p3)*get_basis_late(k,p1);
    b3 = get_basis_late(i,p3)*get_basis_late(j,p1)*get_basis_late(k,p2);
    b4 = get_basis_late(i,p3)*get_basis_late(j,p2)*get_basis_late(k,p1);
    b5 = get_basis_late(i,p2)*get_basis_late(j,p1)*get_basis_late(k,p3);
    b6 = get_basis_late(i,p1)*get_basis_late(j,p3)*get_basis_late(k,p2);
    
    double result = (b1+b2+b3+b4+b5+b6)/(6.0);
    return result;
}
 
