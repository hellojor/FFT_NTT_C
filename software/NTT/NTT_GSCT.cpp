/*
 * NNT_GSCT.cpp
 *
 * Description
 * This progrom wanna to show the NTT and INTT algorithm
 * using DIF (Gentleman-Sande) and DIT (Cooley-Tukey) to implement NTT and INTT respectively
 * without preprocessing (bir reverse)
 * 
 * Using "g++ NTT_GSCT.cpp -o NTT_GSCT.out" to compile the cpp file
 * and using "./NTT_GSCT.out" to run the program
 *
 * History
 * 2023/05/11	jorjor	First release
 * */

#include <iostream>
#include <cmath>

#define q 3329

using namespace std;

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

int findw(int n) {
	// find primitive root
	bool CONTINUE = true;
	int pr;
	for (pr = 2; pr < q; pr++){
		for (int i = 1; i < q; i++) {
			if (i == q - 1) {
				CONTINUE = false;
			}
			else if (quickmod(pr, i) % q == 1) {
				break;
			}
		}
		if (CONTINUE == false) {
			break;
		}
	}

	int w = quickmod(pr, (int)(q - 1)/n);
	return w;
}

int InverseMod(int a) {
	for (int b = 2; b < q; b++) {
		if ((a * b) % q == 1){
			return b;
		}
	}
	return -1;
}

void DIV2(int *a) {
	int temp = (*a) ;
	(*a) = (temp >> 1) + (temp & 1) * ((q + 1) / 2);
}

void print(int* arr, int len) {
    cout << "[ ";
    for (int i = 0; i < len; i++) {
        cout << arr[i];
		if (i != len - 1) {
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

void naive_polynomial_multiplication(int *x1, int *x2, int *arr, int n) {
	for (int i = 0; i < n; i++) {
		arr[i] = 0;
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			arr[(i+j)%n] = (arr[(i+j)%n] + (x1[i] * x2[j])) % q;
		}
	}
}	

int main(){

	// set seed to 0
	srand(0);

	int n = 8;
	//int x1[] = {100, 100, 100, 100};
	//int x2[] = {100, 200, 300, 400};
	int x1[] = {100, 200, 200, 0, 100, 200, 500, 0};
	int x2[] = {100, 200, 300, 400, 500, 600, 700, 800};
	
	/*int *x1 = new int[n];
	int *x2 = new int[n];

	for (int i = 0 ; i < n; i++){
		x1[i] = rand() % q;
		x2[i] = rand() % q;
	}*/

    cout << "***** Original array *****" << endl;
    cout << "x1: "; print(x1, n);
    cout << "x2: "; print(x2, n); cout << endl;

	// build the array of w
	int w = findw(n);
	int winv = InverseMod(w);
	
	int *wn = new int[n];
	wn[0] = 1;
	for (int i = 1; i < n; i++){
		wn[i] = (wn[i-1] * w) % q;
	}

	// build the array of w
	int *wn_inv = new int[n];
	wn_inv[0] = 1;
	for (int i = 1; i < n; i++){
		wn_inv[i] = (wn_inv[i-1] * winv) % q;
	}
    
    cout << "***** Wn array *****" << endl;
    cout << "wn: "; print(wn, n);
	cout << "***** Wn_inv array *****" << endl;
    cout << "wn_inv: "; print(wn_inv, n); cout << endl;
	
	int *x1_ntt = new int[n];
	int *x2_ntt = new int[n];

	for (int i = 0; i < n; i++) {
		x1_ntt[i] = x1[i];
		x2_ntt[i] = x2[i];
	}
	// NTT
	for (int step = log2(n); step >= 1; step--) {
		for (int idx = 0; idx < (n/pow(2, step)); idx++) {
			for (int distance = 0; distance < pow(2, step - 1); distance++) {
				int i = idx * pow(2, step) + distance;
				int j = idx * pow(2, step) + distance + pow(2, step - 1);
				BFU_GS(x1_ntt, i, j, wn[distance * int(n/pow(2, step))]);
				BFU_GS(x2_ntt, i, j, wn[distance * int(n/pow(2, step))]);
			}
		}
	}

    cout << "***** After NTT *****" << endl;
    cout << "x1_ntt: "; print(x1_ntt, n);
    cout << "x2_ntt: "; print(x2_ntt, n); cout << endl;

	int *X_multi = new int[n];

	for (int i = 0; i < n; i++) {
		X_multi[i] = (x1_ntt[i] * x2_ntt[i]) % q;
	}
    cout << "***** Convolution *****" << endl;
    cout << "X_multi: "; print(X_multi, n); cout << endl;
	
	int *X_intt = new int[n];
	
	for (int i = 0; i < n; i++) {
		X_intt[i] = X_multi[i];
	}


	// INTT
	for (int step = 1; step <= log2(n); step++) {
		for (int idx = 0; idx < (n/pow(2, step)); idx++) {
			for (int distance = 0; distance < pow(2, step - 1); distance++) {
				int i = idx * pow(2, step) + distance;
				int j = idx * pow(2, step) + distance + pow(2, step - 1);
				BFU_CT(X_intt, i, j, wn_inv[distance * int(n/pow(2, step))]);
			}
		}
	}

	int rv = quickmod(n, q - 2);
	for (int i = 0; i < n; i++) {
		X_intt[i] = (X_intt[i] * rv) % q;
	}

    cout << "***** After INTT *****" << endl;
    cout << "X_intt: "; print(X_intt, n); cout << endl;
	
	int *naive_result = new int[n];
	naive_polynomial_multiplication(x1, x2, naive_result, n);

    cout << "***** Naive polynomial multiplication *****" << endl;
    cout << "naive_result: "; print(naive_result, n); cout << endl;

# if 0
	
	int *x1_intt = new int[n];
	int *x2_intt = new int[n];

	// reverse
	for(int i = 0; i < n; i++) {
		x1_intt[i] = x1_ntt[i];
		x2_intt[i] = x2_ntt[i];
	}


	// INTT
	for (int step = 1; step <= log2(n); step++) {
		for (int idx = 0; idx < (n/pow(2, step)); idx++) {
			for (int distance = 0; distance < pow(2, step - 1); distance++) {
				int i = idx * pow(2, step) + distance;
				int j = idx * pow(2, step) + distance + pow(2, step - 1);
				BFU_CT(x1_intt, i, j, wn_inv[distance * int(n/pow(2, step))]);
				BFU_CT(x2_intt, i, j, wn_inv[distance * int(n/pow(2, step))]);
			}
		}
	}

	rv = quickmod(n, q - 2);
	for (int i = 0; i < n; i++) {
		x1_intt[i] = (x1_intt[i] * rv) % q;
		x2_intt[i] = (x2_intt[i] * rv) % q;
	}

    cout << "***** After INTT *****" << endl;
    cout << "x1_intt: "; print(x1_intt, n);
    cout << "x2_intt: "; print(x2_intt, n); cout << endl;

	delete x1_intt;
	delete x2_intt;
#endif
	
	delete x1_ntt;
	delete x2_ntt;
	delete X_multi;
	delete X_intt;
	delete naive_result;

	return 0;
}
