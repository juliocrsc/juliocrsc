#include <julio.hpp>
// x^3 - 9x + 3
double Julio::funcaoSecantNewton(double x) {
	return pow(x, 3) - 9 * x + 3;
}

// 3*x^2 - 9
 double Julio::derivada(double x) {
		return 3 * pow(x, 2) - 9;
}

double Julio::metodoNewton(double x0, double e1, double e2){
	double x1=0;
	if(abs(funcaoSecantNewton(x0)) < e1){
		return x0;
	}else{
		int k=0;
		x1 = x0 - (funcaoSecantNewton(x0)/derivada(x0));
		while(abs(funcaoSecantNewton(x1))  >e1 || abs(x1-x0) > e2  ){
			x1 = x0 - (funcaoSecantNewton(x0)/derivada(x0));
			x0= x1;
			k++;
		}		
	}
	return x1;
}

double Julio::secant(double x0, double x1, double e1, double e2){
	double x2 = 0.0;
	if(abs(funcaoSecantNewton(x0)) < e1)
		return x0;
	else{
		int k=0;		
		do{
			x2 = x1 - (funcaoSecantNewton(x1)/ (funcaoSecantNewton(x1) - funcaoSecantNewton(x0)));
			x0 = x1;
			x1 = x2;
		}while(abs(funcaoSecantNewton(x2)) > e1 || abs(x2 - x1) > e2);		
	}
	return x2;
}