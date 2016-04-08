#include <math.h>
#include "global.h"

_OFFLOADABLE double calculate_geometric(int l1, int l2, int l3)
{
	double result;
	double s1,s3,s6,s13,s23,s33,s16,s26,s36;
	int sign;
	double fact,p1,p2,p3;
	int test;

	test = (l1+l2+l3)%2;

	if(!test && l1+l2>=l3 && l2+l3>=l1 && l3+l1>=l2)
	{
		s1 = l1+l2+l3+1.0;
		s3 = l1+l2+l3+1.0/3.0;
		s6 = l1+l2+l3+1.0/6.0;
		s13 = l2+l3-l1+1.0/3.0;
		s23 = l3+l1-l2+1.0/3.0;
		s33 = l1+l2-l3+1.0/3.0;
		s16 = l2+l3-l1+1.0/6.0;
		s26 = l3+l1-l2+1.0/6.0;
		s36 = l1+l2-l3+1.0/6.0;

		fact = M_1_PI * M_1_PI / 2.0;
		p1 = fact*((2.0*l1+1.0)*(2.0*l2+1.0)*(2.0*l3+1.0)) / s1;
		p2 = s3 / (s13*s23*s33);
		p3 = sqrt(s16*s26*s36/s6);

		result = p1*p2*p3;
	} 
	else 
	{
		result = 0.0;
	}

	return result;
}
