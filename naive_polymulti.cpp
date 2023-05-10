/*
 * naive polynomial multiplication.cpp
 *
 * Description
 * This progrom wanna to implement the naive polynomial algorithm
 *
 * 
 * Using "g++ naive_polymulti.cpp -o naive_polymulti.out" to compile the cpp file
 * and using "./naive_polymulti.out" to run the program
 *
 * History
 * 2023/05/10	jorjor	First release
 * */

#include <iostream>

using namespace std;

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

int main(){
	
	int n = 4; // 8
	//double x1[] = {1, 2, 2, 0, 1, 2, 2, 0};
	//double x2[] = {1, 2, 3, 4, 5, 6, 7, 8};
	double x1[] = {1, 2, 2, 0};
	double x2[] = {1, 2, 3, 4};

    cout << "***** Original array *****" << endl;
    cout << "x1: "; print(x1, n);
    cout << "x2: "; print(x2, n); cout << endl;

	double *ans = new double[n];
	// set zero
	for (int i = 0; i < n; i++) {
		ans[i] = 0;
	}
	
	/*
	double *buff = new double[2*n];
	// set zero
	for (int i = 0; i < n * 2; i++) {
		buff[i] = 0;
	}
	
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			buff[i+j] += x1[i] * x2[j];
		}
	}

	for (int i = 2 * n - 1; i >= 0; i--) {
		cout << buff[i] << ' ';
	}
	cout << endl;
	*/

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			ans[(i+j)%n] += x1[i] * x2[j];
		}
	}
	
	// show answer
    cout << "***** Result *****" << endl;
    cout << "result: "; print(ans, n);
	
	//delete buff;
	delete ans;

	return 0;
}
