#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include "global.h"

void calculate_gamma_3D(int n, int i, double *mvec) {

    int j,k,m,t1,t2,t3;
    double x,y,z;
    int terms = get_terms_prim();
    int lsize = get_lmax()+1;
    double s1,s2,s3;

    int xsize = get_b_xsize();
    double *xvec = create_vector(xsize);
    get_b_xvec(xvec);
    double xmax = xvec[xsize-1];

    for(m=0;m<terms;m++) mvec[m] = 0.0;

    s1 = pow(2.0*i+1.0,1.0/3.0)*(get_cl(i)+get_noise(i)/(get_beam(i)*get_beam(i)));
    for(j=i;j<lsize;j++){
        s2 = pow(2.0*j+1.0,1.0/3.0)*(get_cl(j)+get_noise(j)/(get_beam(j)*get_beam(j)));
        t1 = i+j;
        t2 = j%2;
        t3 = t1%2;
        if(t1>get_lmax())t1=get_lmax();
        if(t2==0&&t3==0){
            for(k=j;k<t1+1;k+=2){
                x = calculate_xint(i,j,k,n,xsize,xvec);
                s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl(k)+get_noise(k)/(get_beam(k)*get_beam(k)));
                z = permsix(i,j,k)*calculate_geometric(i,j,k)/sqrt(s1 * s2 * s3);
                for(m=0;m<terms;m++){
                    y = plijk(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }else if(t2==0&&t3==1){
            for(k=j+1;k<t1+1;k+=2){
                x = calculate_xint(i,j,k,n,xsize,xvec);
                s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl(k)+get_noise(k)/(get_beam(k)*get_beam(k)));
                z = permsix(i,j,k)*calculate_geometric(i,j,k)/sqrt(s1 * s2 * s3);
                for(m=0;m<terms;m++){
                    y = plijk(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }else if(t2==1&&t3==0){
            for(k=j+1;k<t1+1;k+=2){
                x = calculate_xint(i,j,k,n,xsize,xvec);
                s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl(k)+get_noise(k)/(get_beam(k)*get_beam(k)));
                z = permsix(i,j,k)*calculate_geometric(i,j,k)/sqrt(s1 * s2 * s3);
                for(m=0;m<terms;m++){
                    y = plijk(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }else if(t2==1&&t3==1){
            for(k=j;k<t1+1;k+=2){
                x = calculate_xint(i,j,k,n,xsize,xvec);
                s3 = pow(2.0*k+1.0,1.0/3.0)*(get_cl(k)+get_noise(k)/(get_beam(k)*get_beam(k)));
                z = permsix(i,j,k)*calculate_geometric(i,j,k)/sqrt(s1 * s2 * s3);
                for(m=0;m<terms;m++){
                    y = plijk(m,i,j,k);
                    mvec[m] += x*y*z;
                }
            }
        }
    }

    for(m=0;m<terms;m++) mvec[m] *= 6.0*deltaphi*deltaphi;
    return;
}
