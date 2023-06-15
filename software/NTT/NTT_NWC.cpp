/*
 * NNT_NWC.cpp
 *
 * Description
 * This progrom wanna to show the in-place FNTT and IFNTT algorithm, which F means Finite Field
 * using Negative Wrapped Convolution (NWC)
 * using DIF (Gentleman-Sande) and DIT (Cooley-Tukey) to implement NTT and INTT respectively
 * without preprocessing (bit reverse)
 *
 * paper reference : A Hardware Accelerator for Polynomial Multiplication Operation of CRYSTALS-KYBER PQC Scheme
 * link : https://ieeexplore.ieee.org/document/9474139
 * 
 * Using "g++ NTT_NWC.cpp -o NTT_NWC.out" to compile the cpp file
 * and using "./NTT_NWC.out" to run the program
 *
 * History
 * 2023/06/15	jorjor	First release
 * */

#include <iostream>
#include <cmath>

#define q 3329
#define n 256

using namespace std;

int wn[n] = {0};
int wn_inv[n] = {0};
int wq[n] = {0};

int gcd(int a, int b) {
	if (b == 0) return a;
	else		return gcd(a, a%b);
}

int quickmod(int a, int b) {
	// a ** b % q
	int ans = 1;
	while (b != 0) {
		if (b & 1) ans = (ans * a) % q;
		a = (a * a) % q;
		b >>= 1;
	}
	return ans;
}

int InverseMod(int a) {
	for (int b = 2; b < q; b++) {
		if ((a * b) % q == 1){
			return b;
		}
	}
	return -1;
}

int DIV2(int a) {
	return (a >> 1) + (a & 1) * ((q + 1) / 2);
}

void print(int* arr) {
    cout << "[ ";
    for (int i = 0; i < n; i++) {
        cout << arr[i];
		if (i != n - 1) {
            cout << ", ";
        }
    }
    cout << " ]" << endl;
}

int bitreverse(int num, int len) {
	int result = 0;

	for (int i = len - 1; i >= 0; i--) {
		result += ((num & 1) * pow(2, i));
		num >>= 1;
	}

	return result;
}

int modq(int num){
	int modnum = num % q;
	if (num < 0){
		return modnum + q;
	}
	else {
		return modnum;
	}
}

void BFU_CT(int *arr, int i, int j, int wn) {
	// DIT-FFT
	// Cooley Tukey algorithm
	int temp1 = arr[i];
	int temp2 = arr[j];
	arr[i] = (temp1 + wn * temp2) % q;
	arr[j] = (temp1 - wn * temp2) % q;

	arr[i] = (arr[i] < 0) ? arr[i] + q : arr[i];
	arr[j] = (arr[j] < 0) ? arr[j] + q : arr[j];
}

void BFU_GS(int *arr, int i, int j, int wn) {
	// DIF-FFT
	// Gentleman Sande algorithm
	int temp1 = arr[i];
	int temp2 = arr[j];
	arr[i] = (temp1 + temp2) % q;
	arr[j] = ((temp1 - temp2) * wn) % q;

	arr[i] = (arr[i] < 0) ? arr[i] + q : arr[i];
	arr[j] = (arr[j] < 0) ? arr[j] + q : arr[j];
}

void naive_polynomial_multiplication(int *x1, int *x2, int *arr) {
	int temp[2*n] = {0};
	for (int i = 0; i < 2*n; i++) {
		temp[i] = 0;
	}
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			temp[i + j] += modq(x1[i] * x2[j]);
		}
	}

	for (int i = 0; i < n; i++) {
		arr[i] = modq(temp[i] - temp[i+n]);
	}
}	

void NTT(int x_ntt[n]){
 	int k = 1;
	for (int i = 1; i <= log2(n) - 1; i++) {
		int m = pow(2, log2(n) - i);
		for (int s = 0; s < n; s += 2*m) {
			for (int j = s; j < s + m; j++) {
				int A = x_ntt[j];
				int B = x_ntt[j + m];
				int W = wn[k];
				int T = modq(W * B);
				int E = modq(A + T);
				int O = modq(A - T);
				x_ntt[j] = E;
				x_ntt[j + m] = O;
			}
			k++;
		}
	}
}

void INTT(int x_intt[n]){
 	int k = 0;
	for (int i = log2(n) - 1; i >= 1; i--) {
		// cout << "stage " << log2(n) - i << endl;
		int m = pow(2, log2(n) - i);
		for (int s = 0; s < n; s += 2*m) {
			for (int j = s; j < s + m; j++) {
				// cout << j << "," << j + m << " | ";
				int A = x_intt[j];
				int B = x_intt[j + m];
				int W = wn_inv[k];
				int E = modq(A + B);
				int O = modq((A - B) * W);
				x_intt[j] = DIV2(E);
				x_intt[j + m] = DIV2(O);
			}
			// cout << endl;
			k++;
		}
	}

}
void PWM(int *out, int a[256], int b[256]) {
	
	int a0, a1;
	int b0, b1;

	for (int i = 0; i < n / 2; i++) {
		a0 = a[2*i];
		a1 = a[2*i+1];
		b0 = b[2*i];
		b1 = b[2*i+1];

		out[2*i] = modq(modq(a0 * b0) + modq(a1 * b1) * wq[i]);
		out[2*i+1] = modq(a0 * b1 + a1 * b0);
	}
}

int main(){

	/* set seed to 0 */
	srand(0);

	int x1[n] = {0};
	int x2[n] = {0};
	for (int i = 0; i < 256; i++) x1[i] = i;
	for (int i = 0; i < 256; i++) x2[i] = i;
	

    cout << "***** Original array *****" << endl;
    cout << "x1: "; print(x1);
    cout << "x2: "; print(x2); cout << endl;

	/* build the array of w */
	int w = 17;
	
	wn[0] = 1;
	for (int i = 1; i < n; i++){
		wn[i] = (wn[i-1] * w) % q;
	}

	/* build the array of winv */
	int winv = InverseMod(w);
	
	wn_inv[0] = 1;
	for (int i = 1; i < n; i++){
		wn_inv[i] = (wn_inv[i-1] * winv) % q;
	}
    
	for (int i = 0; i < n; i++){
		wq[i] = wn[2*bitreverse(i, 7)+1];
	}

	/* bitreverse */
	int temp[n] = {0};
	for (int i = 0; i < n; i++) temp[i] = wn[i];
	for (int i = 0; i < n; i++){
		wn[i] = temp[bitreverse(i, 7)];
	}

	for (int i = 0; i < n; i++) temp[i] = wn_inv[i];
	for (int i = 0; i < n; i++){
		wn_inv[i] = temp[bitreverse(i, 7)+1];
	}
	cout << "***** Wn array *****" << endl;
    cout << "wn: "; print(wn);
	
	cout << "***** Wn_inv array *****" << endl;
    cout << "wn_inv: "; print(wn_inv); cout << endl;
	
	cout << "***** Wq array *****" << endl;
    cout << "wq: "; print(wq); cout << endl;
	
	int x1_ntt[n] = {0};
	int x2_ntt[n] = {0};

	for (int i = 0; i < n; i++) {
		x1_ntt[i] = x1[i];
		x2_ntt[i] = x2[i];
	}

	/* NTT */
	
	NTT(x1_ntt);
	NTT(x2_ntt);
	
	cout << "***** After NTT *****" << endl;
    cout << "x1_ntt: "; print(x1_ntt);
    cout << "x2_ntt: "; print(x2_ntt); cout << endl;

	/* Coefficient-Wise Multiplication */
	int X_CWM[n] = {0};
	
	PWM(X_CWM, x1_ntt, x2_ntt) ;

	cout << "***** After CWM *****" << endl;
    cout << "X_CWM: "; print(X_CWM); cout << endl;
	
	// INTT
	
	INTT(X_CWM);
	
    cout << "***** After INTT *****" << endl;
    cout << "x_intt: "; print(X_CWM); cout << endl;
	
	
	int naive_result[n] = {0};
	naive_polynomial_multiplication(x1, x2, naive_result);

    cout << "***** Naive polynomial multiplication *****" << endl;
    cout << "naive_result: "; print(naive_result); cout << endl;

	return 0;
}
