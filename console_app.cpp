#include <bits/stdc++.h>

#define EPSILON 0.000001

using namespace std;

//some demo non linear equations

//  demo equation : 4x²-32x+63=0
// demo formate   : n(degree) a b c....(co-efficient)
//  demo input1   : 2 4 -32 63

//  demo equation : 8x³-108x²+478x-693=0
// demo formate   : n(degree) a b c....(co-efficient)
//  demo input1   : 3 8 -108 478 -693

//  demo equation : 16x⁴-320x³+2360x²-7600x+9009=0
// demo formate   : n(degree) a b c....(co-efficient)
//  demo input1   : 4 16 -320 2360 -7600 9009

// int no_of_iteration = 0;

struct ab_pair
{
    float a;
    float b;
};

vector<float> taking_input_co_efficient_arr(int power_of_equation)
{
    vector<float> co_efficient_arr(power_of_equation + 1, 0);
    for (int i = 0; i < power_of_equation + 1; i++)
    {
        cin >> co_efficient_arr[i];
    }
    return co_efficient_arr;
}

float equation_value(float x, int power_of_equation, vector<float> co_efficient_arr)
{
    float value = 0;
    float power = power_of_equation;
    // cout<<"check1:"<<co_efficient_arr.size()<<endl;
    for (int i = 0; i < co_efficient_arr.size(); i++)
    {
        value = value + pow(x, power) * co_efficient_arr[i];
        // cout << "check2:" << value << endl;

        power--;
    }

    return value;
}

float derivative_equation_value(float x, int power_of_equation, vector<float> co_efficient_arr)
{
    int co_power = power_of_equation;
    float value = 0;
    int power = power_of_equation - 1;

    for (int i = 0; i < co_efficient_arr.size() - 1; i++)
    {
        value += pow(x, power) * co_efficient_arr[i] * co_power;
        power--;
        co_power--;
    }
    return value;
}

vector<ab_pair> make_pairab(int power_of_equation, vector<float> co_efficient_arr)
{

    float range_a = 0;
    float range_b = 0;

    float c = sqrt((pow((co_efficient_arr[1] / co_efficient_arr[0]), 2) - 2 * (co_efficient_arr[2] / co_efficient_arr[0])));

    range_a = floor(-1 * c);
    range_b = ceil(1 * c);

    vector<ab_pair> pairab;

    pairab.reserve(power_of_equation);

    float temp_a = range_a;
    float temp_b = range_a + 1;

    while (temp_b < range_b + 1)
    {
        if (equation_value(temp_a, power_of_equation, co_efficient_arr) * equation_value(temp_b, power_of_equation, co_efficient_arr) < 0)
        {
            pairab.push_back({temp_a, temp_b});
        }
        temp_a = temp_b;
        temp_b = temp_b + 1;
    }

    return pairab;
}

float bisection_method(float temp_a, float temp_b, int power_of_equation, vector<float> co_efficient_arr)
{
    // int temp_no_of_iteration = 1;
    float mid = temp_a + (temp_b - temp_a) / 2;
    float prev_mid = mid;

    while (fabs(equation_value(mid, power_of_equation, co_efficient_arr)) > 0.000001) // To check floating point tolerance fabs function is necessary
    {
        // temp_no_of_iteration++;
        if (equation_value(temp_a, power_of_equation, co_efficient_arr) * equation_value(mid, power_of_equation, co_efficient_arr) < 0)
        {
            temp_b = mid;
        }
        else
        {
            temp_a = mid;
        }

        mid = temp_a + (temp_b - temp_a) / 2;
        cout << fabs((mid - prev_mid)) << endl;
        prev_mid = mid;
    }

    // no_of_iteration += temp_no_of_iteration;

    return mid;
}

float false_position_method(float temp_a, float temp_b, float power_of_equation, vector<float> co_efficient_arr)
{
    // int temp_no_of_iteration = 1;
    float ftemp_b = equation_value(temp_b, power_of_equation, co_efficient_arr);
    float ftemp_a = equation_value(temp_a, power_of_equation, co_efficient_arr);
    float mid = (temp_a * ftemp_b - temp_b * ftemp_a) / (ftemp_b - ftemp_a);
    float prev_mid = mid;

    while (fabs(equation_value(mid, power_of_equation, co_efficient_arr)) > 0.000001) // flaoting point tolarance check e fabs deya lagbe
    {
        // temp_no_of_iteration++;
        if (equation_value(temp_a, power_of_equation, co_efficient_arr) * equation_value(mid, power_of_equation, co_efficient_arr) < 0)
        {
            temp_b = mid;
        }
        else
        {
            temp_a = mid;
        }

        ftemp_b = equation_value(temp_b, power_of_equation, co_efficient_arr);
        ftemp_a = equation_value(temp_a, power_of_equation, co_efficient_arr);

        mid = (temp_a * ftemp_b - temp_b * ftemp_a) / (ftemp_b - ftemp_a);

        // cout<<fabs((mid-prev_mid))<<endl;

        prev_mid = mid;
    }

    // no_of_iteration += temp_no_of_iteration;

    return mid;
}

float scant_method(float x, float y, int power_of_equation, vector<float> co_efficient_arr)
{
    float f_x = equation_value(x, power_of_equation, co_efficient_arr);
    float f_y = equation_value(y, power_of_equation, co_efficient_arr);
    float z;

    // Ensuring that f(x) and f(y) are not equal initially
    if (fabs(f_x - f_y) < EPSILON)
    {
        cout << "Warning: Initial points too close or identical function values." << endl;
        return x;
    }

    while (fabs(y - x) > EPSILON)
    {
        z = y - f_y * (y - x) / (f_y - f_x);
        x = y;
        f_x = f_y;
        y = z;
        f_y = equation_value(z, power_of_equation, co_efficient_arr);
    }
    return z;
}

float newton_raphson_method(float x, int power_of_equation, vector<float> co_efficient_arr)
{
    float h = equation_value(x, power_of_equation, co_efficient_arr) / derivative_equation_value(x, power_of_equation, co_efficient_arr);
    int iteration = 0;
    float y = 0;

    while (fabs(h) > EPSILON)
    {

        y = x;
        x = y - h;
        h = equation_value(x, power_of_equation, co_efficient_arr) / derivative_equation_value(x, power_of_equation, co_efficient_arr);
        // iteration++;
        cout << fabs(x - y) << endl;
    }

    // cout << "No of iteration :" << iteration << endl;
    return x;
}

void printing_roots(int power_of_equation, vector<float> result)
{
    for (int i = 0; i < power_of_equation; i++)
    {
        cout << "Root" << i + 1 << " : " << result[i] << endl;
    }
}

vector<float> function_call(int power_of_equation, vector<ab_pair> pairab, vector<float> co_efficient_arr, int method_no)
{
    vector<float> result;
    result.resize(power_of_equation);

    if (method_no == 1)
    {
        for (int i = 0; i < power_of_equation; i++)
        {
            result[i] = (bisection_method(pairab[i].a, pairab[i].b, power_of_equation, co_efficient_arr));
        }
        cout << "Solved by Bisection method" << endl;
    }
    else if (method_no == 2)
    {
        for (int i = 0; i < power_of_equation; i++)
        {
            result[i] = (false_position_method(pairab[i].a, pairab[i].b, power_of_equation, co_efficient_arr));
        }
        cout << "Solved by False Position method" << endl;
    }

    else if (method_no == 3)
    {

        for (int i = 0; i < power_of_equation; i++)
        {
            result[i] = newton_raphson_method(((pairab[i].a + pairab[i].b) / 2.0), power_of_equation, co_efficient_arr);
        }
        cout << "Solved by Newton-Raphson method" << endl;
    }

    else if (method_no == 4)
    {
        for (int i = 0; i < power_of_equation; i++)
        {
            result[i] = (scant_method(pairab[i].a + 0.5, pairab[i].b + 0.5, power_of_equation, co_efficient_arr));
            // 0.5 is added for random initial value
        }
        cout << "Solved by Scant method" << endl;
    }
    else
    {
        cout << "Select a valid method" << endl;
    }

    return result;
}

// Function to perform LU Decomposition
void luDecomposition(const vector<vector<double>>& mat, int n, vector<vector<double>>& lower, vector<vector<double>>& upper) {
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += lower[i][j] * upper[j][k];
            }
            upper[i][k] = mat[i][k] - sum;
        }

        for (int k = i; k < n; k++) {
            if (i == k) {
                lower[i][i] = 1;
            } else {
                double sum = 0;
                for (int j = 0; j < i; j++) {
                    sum += lower[k][j] * upper[j][i];
                }
                lower[k][i] = (mat[k][i] - sum) / upper[i][i];
            }
        }
    }
}

// Forward substitution to solve L*y = b
void forwardSubstitution(const vector<vector<double>>& L, vector<double>& y, const vector<double>& b, int n) {
    for (int i = 0; i < n; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum);
    }
}

// Backward substitution to solve U*x = y
void backwardSubstitution(const vector<vector<double>>& U, vector<double>& x, const vector<double>& y, int n) {
    for (int i = n - 1; i >= 0; i--) {
        double sum = 0;
        for (int j = i + 1; j < n; j++) {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
}

// Function to display a vector
void displayVector(const vector<double>& vec) {
    for (double val : vec) {
        cout << setw(10) << setprecision(5) << val << " ";
    }
    cout << endl;
}

// Function to display a matrix
void displayMatrix(const vector<vector<double>>& matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << setw(10) << setprecision(5) << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// Function to calculate the inverse of a matrix using LU decomposition
vector<vector<double>> matrixInverse(const vector<vector<double>>& lower, const vector<vector<double>>& upper, int n) {
    vector<vector<double>> inverse(n, vector<double>(n));

    for (int i = 0; i < n; i++) {
        vector<double> b(n, 0);
        b[i] = 1;  // Setting the i-th element to 1 for the identity matrix

        vector<double> y(n, 0);
        vector<double> x(n, 0);

        // Solving for each column of the inverse
        forwardSubstitution(lower, y, b, n);
        backwardSubstitution(upper, x, y, n);

        // Assign the solution to the corresponding column of the inverse matrix
        for (int j = 0; j < n; j++) {
            inverse[j][i] = x[j];
        }
    }

    return inverse;
}
// Runge-Kutta 4th order method for a linear ODE
void rungeKuttaLinear(double y0, double t0, double t_end, double h) {
    double t = t0;
    double y = y0;

    cout << "t\t y (Linear ODE)" << endl;

    // Linear differential equation: dy/dt = ay + b
    double a = 1.0; // Coefficient
    double b = 0.0; // Constant term

    while (t <= t_end) {
        cout << t << "\t " << y << endl;

        double k1 = h * (a * y + b);
        double k2 = h * (a * (y + k1 / 2) + b);
        double k3 = h * (a * (y + k2 / 2) + b);
        double k4 = h * (a * (y + k3) + b);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6; // Update y
        t += h; // Increment time
    }
}

// Runge-Kutta 4th order method for a trigonometric ODE
void rungeKuttaTrigonometric(double y0, double t0, double t_end, double h) {
    double t = t0;
    double y = y0;

    cout << "\nSolving Trigonometric ODE:" << endl;
    cout << "t\t y (Trigonometric ODE)" << endl;

    // Trigonometric differential equation: dy/dt = sin(t)
    while (t <= t_end) {
        cout << t << "\t " << y << endl;

        double k1 = h * sin(t);
        double k2 = h * sin(t + h / 2);
        double k3 = h * sin(t + h / 2);
        double k4 = h * sin(t + h);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6; // Update y
        t += h; // Increment time
    }
}

int main()

{
    cout<<"1.Nonlinear equation:"<<endl;
    cout<<"2.LU factorization,Runge kutta,Matrix inversion:"<<endl;
    cout<<"Enter your choice:";
    int choice;
    if(choice==1)
    {
        // taking no of degree for equation
    int power_of_equation = 0;
    cin >> power_of_equation;

    // declaring co efficient vector and taking input from function
    vector<float> co_efficient_arr = taking_input_co_efficient_arr(power_of_equation);

    // declaring pair/range(a,b) vector and initializing it by function
    vector<ab_pair> pairab = make_pairab(power_of_equation, co_efficient_arr);

    // which method(out of 4) to call
    int method_no = 0;
    cout << "Method : ";
    cin >> method_no;

    // declaring result vector and initializing it with calling specific function
    vector<float> result = function_call(power_of_equation, pairab, co_efficient_arr, method_no);

    // printing result vector by function
    printing_roots(power_of_equation, result);
    }
    else if(choice==2)
    {
         cout<<"1. LU factorization"<<endl
                <<"2.Matrix inversion"<<endl
        <<"3.Runge Kutta Method"<<endl;


          int select;
            cout<<"Select what type of equation you want to solve:";
            cin>>select;

            if(select==1)
            {
                int n;
    cout << "Enter the number of variables: ";
    cin >> n;
    vector<vector<double>> mat(n, vector<double>(n));

    // Taking user input for the matrix A
    cout << "Enter the elements of the matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "Element [" << i << "][" << j << "]: ";
            cin >> mat[i][j];
        }
    }

    // Declaring lower and upper matrices
    vector<vector<double>> lower(n, vector<double>(n, 0));
    vector<vector<double>> upper(n, vector<double>(n, 0));

    // Performing LU Decomposition
    luDecomposition(mat, n, lower, upper);

    vector<double> b(n);
            cout << "Enter the elements of the vector b:" << endl;
            for (int i = 0; i < n; i++) {
                cout << "Element [" << i << "]: ";
                cin >> b[i];
            }

            vector<double> y(n, 0); // Solution to L*y = b
            vector<double> x(n, 0); // Solution to U*x = y

            forwardSubstitution(lower, y, b, n);
            backwardSubstitution(upper, x, y, n);

            cout << "\nSolution vector x:" << endl;
            displayVector(x);


            }
            else if(select==2)
            {
                          int n;
    cout << "Enter the number of variables: ";
    cin >> n;
    vector<vector<double>> mat(n, vector<double>(n));

    // Taking user input for the matrix A
    cout << "Enter the elements of the matrix A:" << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << "Element [" << i << "][" << j << "]: ";
            cin >> mat[i][j];
        }
    }

    // Declaring lower and upper matrices
    vector<vector<double>> lower(n, vector<double>(n, 0));
    vector<vector<double>> upper(n, vector<double>(n, 0));

    // Performing LU Decomposition
    luDecomposition(mat, n, lower, upper);
                vector<vector<double>> inverse = matrixInverse(lower, upper, n);
            cout << "\nInverse of Matrix A:" << endl;
            displayMatrix(inverse, n);

            }
            else if (select==3)
            {
                double y0;      // Initial condition
    double t0 ;     // Initial time
    double t_end ;  // End time
    double h ;      // Step size
    cout<<"Enter the value of Initial Condition,Initial time,End time,step size"<<endl;
    cin>>y0>>t0>>t_end>>h;

    int choice;
        cout<<"1.solve Linear equation:"<<endl
        <<"2.Tigonometric(sine) equation"<<endl;
        cout<<"Enter your choice:";
        cin>>choice;

if(choice==1)rungeKuttaLinear(y0, t0, t_end, h);        // Solve linear ODE
    else if (choice==2)rungeKuttaTrigonometric(y0, t0, t_end, h); // Solve trigonometric ODE
else cout<<"No solution choose"<<endl;
    }


}
else
    cout<<"No choice selected"<<endl;

}
