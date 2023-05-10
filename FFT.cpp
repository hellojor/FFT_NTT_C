/*
 * FFT_org.cpp
 *
 * Description
 * This progrom wanna to show the FFT and IFFT algorithm
 * Using DIF-FFT (Centleman-Sande) and DIT-FFT (Cooley-Tukey) to implement FFT and IFFT respectively
 *
 * 
 * Using "g++ FFT_org.cpp -o FFT_org.out" to compile the cpp file
 * and using "./FFT_org.out" to run the program
 *
 * History
 * 2023/05/10	jorjor	First release
 * */

#include <iostream>
#include <complex>
#include <cmath>

using namespace std;

typedef complex<double> Complex;

enum {
	normal = 0,
	inverse
};

void BFU_CT(Complex* arr, int i, int j, Complex w) {
	// DIT-FFT
	// Cooley-Tukey butterfly unit
    Complex temp1 = arr[i];
    Complex temp2 = arr[j];
    arr[i] = temp1 + w * temp2;
    arr[j] = temp1 - w * temp2;
}

void BFU_GS(Complex* arr, int i, int j, Complex w) {
	// DIF-FFT
	// Gentleman-Sande butterfly unit
    Complex temp1 = arr[i];
    Complex temp2 = arr[j];
    arr[i] = temp1 + temp2;
    arr[j] = (temp1 - temp2) * w;
}

int reverse(int num,int len) { // reverse bit
    int out = 0;
    for (int i = 0; i < len; i++){
        if (num & 1) {
            out += pow(2, len - i - 1);
        }
        num >>= 1;
    }
    return out;
}

void print(double* arr, int len) {
    cout << "[ ";
    for (int i = 0; i < len; i++) {
        cout << arr[i];
		if (i != len - 1) {
            cout << ", ";
        }
    }
    cout << " ]" << endl;
}

void print_complex(Complex* arr, int len) {
    cout << "[ ";
    for (int i = 0; i < len; i++) {
        cout << "(" << arr[i].real() << " + " << arr[i].imag() << "i)";
        if (i != len - 1) {
            cout << ", ";
        }
    }
    cout << " ]" << endl;
}

Complex W(int m, int n, bool stat) {
	// acos(-1) = pi
    Complex w;
    w.real(cos(2*acos(-1)*m/n));
	if (stat == inverse) {
		w.imag(sin(2*acos(-1)*m/n));
	}
	else {
		w.imag(sin(-2*acos(-1)*m/n));
	}
    return w;
}

void right_rotate(double* arr, int len){
	double *temp = new double[len];
	
	// copy
	for (int i = 0; i < len; i++) {
		temp[i] = arr[i];
	}

	arr[0] = temp[len - 1];
	for (int i = 0; i < len - 1; i++) {
		arr[i+1] = temp[i];
	}

	delete temp;
}

double* convolution(double x1[], double x2[], int len){
	double *a = new double[len];
	double *b = new double[len];
	double *out = new double[len];
	// reverse b only
	for (int i = 0; i < len; i++) {
		a[i] = x1[i];
		b[i] = x2[len - i - 1];
	}
	
	// fixed a
	// rotate b
	int sum;
	for (int step = 0; step < len; step++) {
		// right rotate
		right_rotate(b, len);

		sum = 0;
		for (int i = 0; i < len; i++) {
			sum += (a[i] * b[i]);
		}
		out[step] = sum;
	}
	
	delete a;
	delete b;
	
	return out;
}

int main() {
    int n = 8; // 4
    double x1[] = {1, 2, 2, 0, 1, 2, 2, 0};
    double x2[] = {1, 2, 3, 4, 5, 6, 7, 8};
    //double x1[] = {1, 2, 2, 0};
    //double x2[] = {1, 2, 3, 4};
    
    cout << "***** Original array *****" << endl;
    cout << "x1: "; print(x1, n);
    cout << "x2: "; print(x2, n); cout << endl;

    
	Complex* x1_complex;
    Complex* x2_complex;
    x1_complex = new Complex[n];
    x2_complex = new Complex[n];

    for (int i = 0; i < n; i++) { // preprocessing
        x1_complex[i].real(x1[i]);
        x1_complex[i].imag(0);

        x2_complex[i].real(x2[i]);
        x2_complex[i].imag(0);
    }

    cout << "***** After Preprocessing *****" << endl;
    cout << "x1: "; print_complex(x1_complex, n);
    cout << "x2: "; print_complex(x2_complex, n); cout << endl;
	
	// FFT
    for (int step = log2(n); step >= 1; step--) {								// 
        cout << step << endl;
		for (int idx = 0; idx < (n/pow(2, step)); idx++) {					// 
            for (int distance = 0; distance < pow(2, step - 1); distance++) {	//
				int i = idx * pow(2, step) + distance;
				int j = idx * pow(2, step) + distance + pow(2, step - 1);
				cout << i << ' ' << j << "(" << distance << ", " << pow(2, step) << ") " << endl;
				Complex w = W(distance, pow(2, step), normal);
                BFU_GS(x1_complex, i, j, w);
                BFU_GS(x2_complex, i, j, w);
            }
        }
    }

    cout << "***** After FFT *****" << endl;
    cout << "x1: "; print_complex(x1_complex, n);
    cout << "x2: "; print_complex(x2_complex, n); cout << endl;
	
    Complex* X_multi = new Complex[n];

	// convolution
    for (int i = 0; i < n; i++) {
        X_multi[i] = x1_complex[i] * x2_complex[i];
    }

    cout << "***** After Convolution *****" << endl;
    cout << "X = x1 * x2 : "; print_complex(X_multi, n); cout << endl;
	
	// reverse
	Complex* X_complex = new Complex[n];

	for (int i = 0; i < n; i++) {
		X_complex[i] = X_multi[i];
	} 

	// IFFT
    for (int step = 1; step <= log2(n); step++) {								// 
        for (int idx = 0; idx < (n/pow(2, step)); idx++) {					// 
            for (int distance = 0; distance < pow(2, step - 1); distance++) {	//
				int i = idx * pow(2, step) + distance;
				int j = idx * pow(2, step) + distance + pow(2, step - 1);
				Complex w = W(distance, pow(2, step), inverse);
                BFU_CT(X_complex, i, j, w);
            }
        }
    }
	for (int i = 0; i < n; i++) {
		X_complex[i] /= n;
	}


    cout << "***** After IFFT *****" << endl;
    cout << "x: "; print_complex(X_complex, n); cout << endl;

    cout << "***** Original Convolution *****" << endl;
	double * org_conv = convolution(x1, x2, n);
    cout << "org_conv: "; print(org_conv, n); cout << endl;

	delete x1_complex;
	delete x2_complex;
	delete X_multi;
	delete X_complex;
    
	return 0;
}

