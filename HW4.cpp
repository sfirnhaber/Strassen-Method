// HW4.cpp
// Sam Firnhaber

#include <iostream>
#include <fstream>
using namespace std;

const string IN_FILE = "input.txt";
const string OUT_FILE = "output.txt";
ofstream outFile;

//Reads the from the input file one 2D array worth of values
void fillArray(ifstream& file, int**& matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> matrix[i][j];
        }
    }
}

void outputArray(int**& matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            outFile << matrix[i][j] << " ";
        }
        outFile << "\n";
    }
}

//Prints out array
void printArray(int**& matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

//Prints the S values used in the Strassen Method
void printSValues(int s[], int n) {
    cout << "S values: " << endl;
    for (int i = 0; i < n; i++) {   
        cout << "S" << (i + 1) << ": " << s[i] << endl;
    }
    cout << endl;
}

//Allocates memory for a 2D array to be used
void allocateArray(int**& matrix, int n) {
    matrix = new int* [n];
    for (int i = 0; i < n; i++)
        matrix[i] = new int[n];
}

//Deallocated the memory needed for a 2D array
void deallocateArray(int**& matrix, int n) {
    for (int i = 0; i < n; i++)
        delete[] matrix[i];
    delete[] matrix;
}

//Adds two 20 matrices together
void addMatrices(int**& A, int**& B, int**& C, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] + B[i][j];
}

//Subtracts two 20 matrices together
void subtractMatrices(int**& A, int**& B, int**& C, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
}

//Multiply two 20 matrices using Strassen's Method
void strassenMultiply(int**& A, int**& B, int**& C) {

    //S values
    int s[10];
    s[0] = B[0][1] - B[1][1];
    s[1] = A[0][0] + A[0][1];
    s[2] = A[1][0] + A[1][1];
    s[3] = B[1][0] - B[0][0];
    s[4] = A[0][0] + A[1][1];
    s[5] = B[0][0] + B[1][1];
    s[6] = A[0][1] - A[1][1];
    s[7] = B[1][0] + B[1][1];
    s[8] = A[0][0] - A[1][0];
    s[9] = B[0][0] + B[0][1];

    //print S values
    printSValues(s, 10);

    //P values
    int p[7];
    p[0] = A[0][0] * s[0];
    p[1] = s[1] * B[1][1];
    p[2] = s[2] * B[0][0];
    p[3] = A[1][1] * s[3];
    p[4] = s[4] * s[5];
    p[5] = s[6] * s[7];
    p[6] = s[8] * s[9];

    //C values
    C[0][0] = p[4] + p[3] - p[1] + p[5];
    C[0][1] = p[0] + p[1];
    C[1][0] = p[2] + p[3]; 
    C[1][1] = p[4] + p[0] - p[2] - p[6];

}

void StrassenMethodR(int**& A, int**& B, int**& C, int n) {

    if (n == 2)
        strassenMultiply(A, B, C);
    else {

        int m = n / 2;
        int** a11; int** a12; int** a21; int** a22;
        int** b11; int** b12; int** b21; int** b22;
        int** c11; int** c12; int** c21; int** c22;
        int** p1; int** p2; int** p3; int** p4; int** p5; int** p6; int** p7;
        int** s1; int** s2; int** s3; int** s4; int** s5; int** s6; int** s7; int** s8; int** s9; int** s10;
        int** aMatrix; int** bMatrix;

        //Massive memory allocation
        allocateArray(a11, m); allocateArray(a12, m); allocateArray(a21, m); allocateArray(a22, m);
        allocateArray(b11, m); allocateArray(b12, m); allocateArray(b21, m); allocateArray(b22, m);
        allocateArray(c11, m); allocateArray(c12, m); allocateArray(c21, m); allocateArray(c22, m);
        allocateArray(p1, m); allocateArray(p2, m); allocateArray(p3, m); allocateArray(p4, m); allocateArray(p5, m); allocateArray(p6, m); allocateArray(p7, m);
        allocateArray(s1, m); allocateArray(s2, m); allocateArray(s3, m); allocateArray(s4, m); allocateArray(s5, m); allocateArray(s6, m); allocateArray(s7, m); allocateArray(s8, m); allocateArray(s9, m); allocateArray(s10, m);
        allocateArray(aMatrix, m); allocateArray(bMatrix, m);

        //Submatrices
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                a11[i][j] = A[i][j];
                a12[i][j] = A[i][j + m];
                a21[i][j] = A[i + m][j];
                a22[i][j] = A[i + m][j + m];
                b11[i][j] = B[i][j];
                b12[i][j] = B[i][j + m];
                b21[i][j] = B[i + m][j];
                b22[i][j] = B[i + m][j + m];
            }
        }

        //S values
        subtractMatrices(b12, b22, s1, m);
        addMatrices(a11, a12, s2, m);
        addMatrices(a21, a22, s3, m);
        subtractMatrices(b21, b11, s4, m);
        addMatrices(a11, a22, s5, m);
        addMatrices(b11, b22, s6, m);
        subtractMatrices(a12, a22, s7, m);
        addMatrices(b21, b22, s8, m);
        subtractMatrices(a11, a21, s9, m);
        addMatrices(b11, b12, s10, m);

        //P values
        StrassenMethodR(a11, s1, p1, m);
        StrassenMethodR(s2, b22, p2, m);
        StrassenMethodR(s3, b11, p3, m);
        StrassenMethodR(a22, s4, p4, m);
        StrassenMethodR(s5, s6, p5, m);
        StrassenMethodR(s7, s8, p6, m);
        StrassenMethodR(s9, s10, p7, m);

        //C values
        //C12
        addMatrices(p1, p2, c12, m);
        //C21
        addMatrices(p3, p4, c21, m);
        //C11
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                c11[i][j] = 0;
        addMatrices(c11, p5, c11, m);
        addMatrices(c11, p4, c11, m);
        subtractMatrices(c11, p2, c11, m);
        addMatrices(c11, p6, c11, m);
        //C22
        for (int i = 0; i < m; i++)
            for (int j = 0; j < m; j++)
                c22[i][j] = 0;
        addMatrices(c22, p5, c22, m);
        addMatrices(c22, p1, c22, m);
        subtractMatrices(c22, p3, c22, m);
        subtractMatrices(c22, p7, c22, m);

        //Creating C
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < m; j++) {
                C[i][j] = c11[i][j];
                C[i][j + m] = c12[i][j];
                C[i + m][j] = c21[i][j];
                C[i + m][j + m] = c22[i][j];
            }
        }

        //Massive memory deallocation
        deallocateArray(a11, m); deallocateArray(a12, m); deallocateArray(a21, m); deallocateArray(a22, m);
        deallocateArray(b11, m); deallocateArray(b12, m); deallocateArray(b21, m); deallocateArray(b22, m);
        deallocateArray(c11, m); deallocateArray(c12, m); deallocateArray(c21, m); deallocateArray(c22, m);
        deallocateArray(p1, m); deallocateArray(p2, m); deallocateArray(p3, m); deallocateArray(p4, m); deallocateArray(p5, m); deallocateArray(p6, m); deallocateArray(p7, m);
        deallocateArray(s1, m); deallocateArray(s2, m); deallocateArray(s3, m); deallocateArray(s4, m); deallocateArray(s5, m); deallocateArray(s6, m); deallocateArray(s7, m); deallocateArray(s8, m); deallocateArray(s9, m); deallocateArray(s10, m);
        deallocateArray(aMatrix, m); deallocateArray(bMatrix, m);
    }

}

void StrassenMethod(int**& A, int**& B, int**& C, int n) {
    int** Crecursive;
    allocateArray(Crecursive, n);

    StrassenMethodR(A, B, Crecursive, n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            C[i][j] = Crecursive[i][j];
        }
    }

    deallocateArray(Crecursive, n);
}

int main() {

    int n;

    ifstream file;
    file.open(IN_FILE);
    if (!file.is_open()) //Error opening file
        return -1;

    //Reads input file data into arrays
    file >> n;
    int** matrixA;
    int** matrixB;
    int** matrixC;
    allocateArray(matrixA, n);
    allocateArray(matrixB, n);
    allocateArray(matrixC, n);
    fillArray(file, matrixA, n);
    fillArray(file, matrixB, n);
    file.close();

    //Computes the Strassen Method
    outFile.open(OUT_FILE);
    if (n > 1)
        StrassenMethod(matrixA, matrixB, matrixC, n);
    else
        matrixC[0][0] = matrixA[0][0] * matrixB[0][0];
    printArray(matrixC, n);

    //Prints new matrix to output file
    outputArray(matrixC, n);
    outFile.close();

    //Clears memory used in arrays
    deallocateArray(matrixA, n);
    deallocateArray(matrixB, n);
    deallocateArray(matrixC, n);

    return 0;

}
