#include "global.h"

double plijk(int r,int i,int j,int k){
            
    int p1,p2,p3;
    find_perm_late(r,&p1,&p2,&p3);
    int pmax_late = get_pmax_late();
            
    double* restrict basis_late_flat = get_basis_late_array();
    double (*restrict basis_late)[pmax_late+1] = (double (*restrict)[pmax_late+   1]) basis_late_flat;
        
    double b1,b2,b3,b4,b5,b6;
    b1 = basis_late[i][p1]*basis_late[j][p2]*basis_late[k][p3];
    b2 = basis_late[i][p2]*basis_late[j][p3]*basis_late[k][p1];
    b3 = basis_late[i][p3]*basis_late[j][p1]*basis_late[k][p2];
    b4 = basis_late[i][p3]*basis_late[j][p2]*basis_late[k][p1];
    b5 = basis_late[i][p2]*basis_late[j][p1]*basis_late[k][p3];
    b6 = basis_late[i][p1]*basis_late[j][p3]*basis_late[k][p2];

    double result = (b1+b2+b3+b4+b5+b6)/(6.0);

    return result;
} 
