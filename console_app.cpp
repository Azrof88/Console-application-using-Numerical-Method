#include <bits/stdc++.h>

#define EPSILON 0.000001

using namespace std;

void jacobi_iteration(const vector<vector<double>> &a, const vector<double> &b, vector<double> &x_init, int max_iterations, double tolerance)
{
    int n = b.size();
    vector<double> x = x_init;
    vector<double> x_new(n);
    vector<double> mae_list;

    cout << setw(10) << "Iteration" << setw(10) << "xi" << setw(10) << "yi" << setw(10) << "zi"
         << setw(10) << "x(i+1)" << setw(10) << "y(i+1)" << setw(10) << "z(i+1)"
         << setw(10) << "ex" << setw(10) << "ey" << setw(10) << "ez"
         << setw(10) << "mae" << setw(10) << "mse" << setw(10) << "rmse" << endl;

    for (int iteration = 0; iteration < max_iterations; iteration++)
    {
        for (int i = 0; i < n; i++)
        {
            double sum_ax = 0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    sum_ax += a[i][j] * x[j];
                }
            }
            x_new[i] = (b[i] - sum_ax) / a[i][i];
        }

        double ex = (n > 0) ? abs(x_new[0] - x[0]) : 0;
        double ey = (n > 1) ? abs(x_new[1] - x[1]) : 0;
        double ez = (n > 2) ? abs(x_new[2] - x[2]) : 0;

        double mae = (ex + ey + ez) / n;
        double mse = (ex * ex + ey * ey + ez * ez) / n;
        double rmse = sqrt(mse);

        cout << setw(10) << iteration + 1 << setw(10) << fixed << setprecision(4) << x[0]
             << setw(10) << x[1] << setw(10) << x[2]
             << setw(10) << x_new[0] << setw(10) << x_new[1] << setw(10) << x_new[2]
             << setw(10) << ex << setw(10) << ey << setw(10) << ez
             << setw(10) << mae << setw(10) << mse << setw(10) << rmse << endl;

        if (sqrt(pow(x_new[0] - x[0], 2) + pow(x_new[1] - x[1], 2) + pow(x_new[2] - x[2], 2)) < tolerance)
        {
            cout << "Converged after " << iteration + 1 << " iterations." << endl;
            break;
        }

        x = x_new;
        mae_list.push_back(mae);
    }
}

void gauss_seidel(const vector<vector<double>> &a, const vector<double> &b, int max_iterations, double tolerance)
{
    int n = b.size();
    vector<double> x(n, 0);
    vector<double> mae_list;

    cout << setw(10) << "Iteration" << setw(10) << "xi" << setw(10) << "yi" << setw(10) << "zi"
         << setw(10) << "x(i+1)" << setw(10) << "y(i+1)" << setw(10) << "z(i+1)"
         << setw(10) << "ex" << setw(10) << "ey" << setw(10) << "ez"
         << setw(10) << "mae" << setw(10) << "mse" << setw(10) << "rmse" << endl;

    for (int iteration = 0; iteration < max_iterations; iteration++)
    {
        vector<double> x_new = x;

        for (int i = 0; i < n; i++)
        {
            double sum_ax = 0;
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    sum_ax += a[i][j] * x_new[j];
                }
            }
            x_new[i] = (b[i] - sum_ax) / a[i][i];
        }

        double ex = (n > 0) ? abs(x_new[0] - x[0]) : 0;
        double ey = (n > 1) ? abs(x_new[1] - x[1]) : 0;
        double ez = (n > 2) ? abs(x_new[2] - x[2]) : 0;

        double mae = (ex + ey + ez) / n;
        double mse = (ex * ex + ey * ey + ez * ez) / n;
        double rmse = sqrt(mse);

        cout << setw(10) << iteration + 1 << setw(10) << fixed << setprecision(4) << x[0]
             << setw(10) << x[1] << setw(10) << x[2]
             << setw(10) << x_new[0] << setw(10) << x_new[1] << setw(10) << x_new[2]
             << setw(10) << ex << setw(10) << ey << setw(10) << ez
             << setw(10) << mae << setw(10) << mse << setw(10) << rmse << endl;

        if (sqrt(pow(x_new[0] - x[0], 2) + pow(x_new[1] - x[1], 2) + pow(x_new[2] - x[2], 2)) < tolerance)
        {
            cout << "Converged after " << iteration + 1 << " iterations." << endl;
            break;
        }

        x = x_new;
        mae_list.push_back(mae);
    }
}

void calculate_errors(const vector<double> &true_values, const vector<double> &estimated_values, double &mse, double &mae, double &rmse)
{
    double sum_squared_error = 0;
    double sum_absolute_error = 0;
    int n = true_values.size();

    for (int i = 0; i < n; ++i)
    {
        double error = true_values[i] - estimated_values[i];
        sum_squared_error += error * error;
        sum_absolute_error += abs(error);
    }

    mse = sum_squared_error / n;
    mae = sum_absolute_error / n;
    rmse = sqrt(mse);
}

void printMatrix(const vector<vector<double>> &matrix)
{
    for (const auto &row : matrix)
    {
        for (double val : row)
        {
            cout << setw(10) << setprecision(4) << val << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void rearrangeRows(vector<vector<double>> &matrix, int col)
{
    int n = matrix.size();
    int flag = 0;

    for (int row = 0; row < n; ++row)
    {
        if (matrix[row][col] == 0)
        {
            for (int l = 0; l < col; l++)
            {

                if (matrix[row][l] < 0 || matrix[row][l] > 0)
                {
                    flag = 1;

                    break;
                }
            }

            if (flag == 0)
            {
                for (int k = row + 1; k < n; ++k)
                {
                    if (matrix[k][col] != 0)
                    {

                        swap(matrix[row], matrix[k]);
                        break;
                    }
                }
            }
        }
    }
}

vector<vector<double>> gaussElimination(vector<vector<double>> &matrix)
{
    int n = matrix.size();
    int m = matrix[0].size();

    for (int col = 0; col < n - 1; ++col)
    {

        rearrangeRows(matrix, col);
        for (int row = col + 1; row < n; ++row)
        {
            if (matrix[row][col] != 0)
            {
                double multiplier = matrix[row][col] / matrix[col][col];
                for (int j = col; j < m; ++j)
                {
                    matrix[row][j] -= multiplier * matrix[col][j];
                }
            }
        }
    }
    return matrix;
}

vector<vector<double>> gaussJordan(vector<vector<double>> &matrix)
{
    int n = matrix.size();
    int m = matrix[0].size();

    for (int col = 0; col < n; ++col)
    {

        rearrangeRows(matrix, col);

        for (int row = 0; row < n; ++row)
        {
            if (row != col)
            {
                double multiplier = matrix[row][col] / matrix[col][col];
                for (int j = 0; j < m; ++j)
                {
                    matrix[row][j] -= multiplier * matrix[col][j];
                }
            }
        }
    }

    return matrix;
}

struct ab_pair
{
    float a;
    float b;
};

vector<float> taking_input_co_efficient_arr(int power_of_equation)
{
    vector<float> co_efficient_arr(power_of_equation + 1, 0);
    cout<<"Enter coefficients of equation : ";
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
        //cout << "Warning: Initial points too close or identical function values." << endl;
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
void luDecomposition(const vector<vector<double>> &mat, int n, vector<vector<double>> &lower, vector<vector<double>> &upper)
{
    for (int i = 0; i < n; i++)
    {
        for (int k = i; k < n; k++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
            {
                sum += lower[i][j] * upper[j][k];
            }
            upper[i][k] = mat[i][k] - sum;
        }

        for (int k = i; k < n; k++)
        {
            if (i == k)
            {
                lower[i][i] = 1;
            }
            else
            {
                double sum = 0;
                for (int j = 0; j < i; j++)
                {
                    sum += lower[k][j] * upper[j][i];
                }
                lower[k][i] = (mat[k][i] - sum) / upper[i][i];
            }
        }
    }
}

// Forward substitution to solve L*y = b
void forwardSubstitution(const vector<vector<double>> &L, vector<double> &y, const vector<double> &b, int n)
{
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < i; j++)
        {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum);
    }
}

// Backward substitution to solve U*x = y
void backwardSubstitution(const vector<vector<double>> &U, vector<double> &x, const vector<double> &y, int n)
{
    for (int i = n - 1; i >= 0; i--)
    {
        double sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum += U[i][j] * x[j];
        }
        x[i] = (y[i] - sum) / U[i][i];
    }
}

// Function to display a vector
void displayVector(const vector<double> &vec)
{
    for (double val : vec)
    {
        cout << setw(10) << setprecision(5) << val << " ";
    }
    cout << endl;
}

// Function to display a matrix
void displayMatrix(const vector<vector<double>> &matrix, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << setw(10) << setprecision(5) << matrix[i][j] << " ";
        }
        cout << endl;
    }
}

// Function to calculate the inverse of a matrix using LU decomposition
vector<vector<double>> matrixInverse(const vector<vector<double>> &lower, const vector<vector<double>> &upper, int n)
{
    vector<vector<double>> inverse(n, vector<double>(n));

    for (int i = 0; i < n; i++)
    {
        vector<double> b(n, 0);
        b[i] = 1; // Setting the i-th element to 1 for the identity matrix

        vector<double> y(n, 0);
        vector<double> x(n, 0);

        // Solving for each column of the inverse
        forwardSubstitution(lower, y, b, n);
        backwardSubstitution(upper, x, y, n);

        // Assign the solution to the corresponding column of the inverse matrix
        for (int j = 0; j < n; j++)
        {
            inverse[j][i] = x[j];
        }
    }

    return inverse;
}
// Runge-Kutta 4th order method for a linear ODE
void rungeKuttaLinear(double y0, double t0, double t_end, double h)
{
    double t = t0;
    double y = y0;

    cout << "t\t y (Linear ODE)" << endl;

    // Linear differential equation: dy/dt = ay + b
    double a = 1.0; // Coefficient
    double b = 0.0; // Constant term

    while (t <= t_end)
    {
        cout << t << "\t " << y << endl;

        double k1 = h * (a * y + b);
        double k2 = h * (a * (y + k1 / 2) + b);
        double k3 = h * (a * (y + k2 / 2) + b);
        double k4 = h * (a * (y + k3) + b);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6; // Update y
        t += h;                               // Increment time
    }
}

// Runge-Kutta 4th order method for a trigonometric ODE
void rungeKuttaTrigonometric(double y0, double t0, double t_end, double h)
{
    double t = t0;
    double y = y0;

    cout << "\nSolving Trigonometric ODE:" << endl;
    cout << "t\t y (Trigonometric ODE)" << endl;

    // Trigonometric differential equation: dy/dt = sin(t)
    while (t <= t_end)
    {
        cout << t << "\t " << y << endl;

        double k1 = h * sin(t);
        double k2 = h * sin(t + h / 2);
        double k3 = h * sin(t + h / 2);
        double k4 = h * sin(t + h);

        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6; // Update y
        t += h;                               // Increment time
    }
}

void slnoflinear()
{
    cout << "1. Jacobi iterative method" << endl
         << "2. Gauss-Seidel iterative method" << endl
         << "3. Gauss elimination" << endl
         << "4. Gausss-Jordan elimination" << endl
         << "5. LU factorization" << endl;

    int select;
    cout << "Select what type of equation you want to solve:";
    cin >> select;

    if (select == 1)
    {
        int n;
        cout << "Enter the number of equations (and variables): ";
        cin >> n;

        vector<vector<double>> a(n, vector<double>(n));
        vector<double> b(n);

        cout << "Enter the coefficients of the equations (" << n << "x" << n + 1 << " matrix):" << endl;
        for (int i = 0; i < n; i++)
        {
            cout << "Row " << i + 1 << " (coefficients followed by constant): ";
            for (int j = 0; j < n; j++)
            {
                cin >> a[i][j];
            }
            cin >> b[i];
        }

        int max_iterations;
        double tolerance;
        cout << "Enter the maximum number of iterations: ";
        cin >> max_iterations;
        cout << "Enter the tolerance for convergence: ";
        cin >> tolerance;

        vector<double> x_init(n, 0.0);

        jacobi_iteration(a, b, x_init, max_iterations, tolerance);

        vector<double> true_values(n, 1.0);
        vector<double> final_solution_jacobi(n, 0.0);

        double mse_jacobi, mae_jacobi, rmse_jacobi;
        calculate_errors(true_values, final_solution_jacobi, mse_jacobi, mae_jacobi, rmse_jacobi);

        cout << "Jacobi Method: MSE: " << mse_jacobi << ", MAE: " << mae_jacobi << ", RMSE: " << rmse_jacobi << endl;
    }

    else if (select == 2)
    {
        int n;
        cout << "Enter the number of equations (and variables): ";
        cin >> n;

        vector<vector<double>> a(n, vector<double>(n));
        vector<double> b(n);

        cout << "Enter the coefficients of the equations (" << n << "x" << n + 1 << " matrix):" << endl;
        for (int i = 0; i < n; i++)
        {
            cout << "Row " << i + 1 << " (coefficients followed by constant): ";
            for (int j = 0; j < n; j++)
            {
                cin >> a[i][j];
            }
            cin >> b[i];
        }

        int max_iterations;
        double tolerance;
        cout << "Enter the maximum number of iterations: ";
        cin >> max_iterations;
        cout << "Enter the tolerance for convergence: ";
        cin >> tolerance;

        vector<double> x_init(n, 0.0);
        gauss_seidel(a, b, max_iterations, tolerance);

        vector<double> true_values(n, 1.0);
        vector<double> final_solution_jacobi(n, 0.0);
        vector<double> final_solution_gs(n, 0.0);

        double mse_gs, mae_gs, rmse_gs;
        calculate_errors(true_values, final_solution_gs, mse_gs, mae_gs, rmse_gs);

        cout << "Gauss-Seidel Method: MSE: " << mse_gs << ", MAE: " << mae_gs << ", RMSE: " << rmse_gs << endl;
    }

    else if (select == 3)
    {
        int n;
        cout << "Enter the number of equations (and variables): ";
        cin >> n;
        vector<vector<double>> matrix(n, vector<double>(n + 1));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= n; j++)
            {
                cin >> matrix[i][j];
            }
        }
        cout << "Original Matrix:\n";
        printMatrix(matrix);

        // Perform Gauss Elimination
        matrix = gaussElimination(matrix);
        cout << "\nGauss Elimination Matrix\n";
        printMatrix(matrix);
    }
    else if (select == 4)
    {
        int n;
        cout << "Enter the number of equations (and variables): ";
        cin >> n;
        vector<vector<double>> matrix(n, vector<double>(n + 1));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= n; j++)
            {
                cin >> matrix[i][j];
            }
        }
        cout << "Original Matrix:\n";
        printMatrix(matrix);

        // Perform Gauss Elimination
        matrix = gaussElimination(matrix);

        matrix = gaussJordan(matrix);
        cout << "\nGauss Jordan Elimination matrix\n";
        printMatrix(matrix);
        cout << endl;
    }
    else if (select == 5)
    {
        int n;
        cout << "Enter the number of variables: ";
        cin >> n;
        vector<vector<double>> mat(n, vector<double>(n));

        // Taking user input for the matrix A
        cout << "Enter the elements of the matrix A:" << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
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
        for (int i = 0; i < n; i++)
        {
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
    else
        cout << "NO method selected";
}

void slnofnonlinear()
{
    cout << "1.Bisection Method" << endl
         << "2.False Position Method" << endl
         << "3.Newton-Raphson Method" << endl
         << "4.Secant Method" << endl;

    int select;

    cout<<"Select a method : ";

    cin>>select;

    cout << "Enter no of degree for equation : ";
    int power_of_equation = 0;
    cin >> power_of_equation;

    // declaring co efficient vector and taking input from function
    vector<float> co_efficient_arr = taking_input_co_efficient_arr(power_of_equation);

    // declaring pair/range(a,b) vector and initializing it by function
    vector<ab_pair> pairab = make_pairab(power_of_equation, co_efficient_arr);



    // declaring result vector and initializing it with calling specific function
    vector<float> result = function_call(power_of_equation, pairab, co_efficient_arr, select);

    // printing result vector by function
    printing_roots(power_of_equation, result);
}

void slnofdiff()
{
    double y0;
    double t0;
    double t_end;
    double h;
    cout << "Enter the value of Initial Condition, Initial time, End time, step size" << endl;
    cin >> y0 >> t0 >> t_end >> h;

    int choice;
    cout << "1. Solve Linear equation:" << endl
         << "2. Trigonometric (sine) equation" << endl;
    cout << "Enter your choice: ";
    cin >> choice;

    if (choice == 1)
    {
        rungeKuttaLinear(y0, t0, t_end, h);
    }
    else if (choice == 2)
    {
        rungeKuttaTrigonometric(y0, t0, t_end, h);
    }
    else
    {
        cout << "No solution choose" << endl;
    }
}

void inversematrix()
{
    int n;
    cout << "Enter the number of variables: ";
    cin >> n;
    vector<vector<double>> mat(n, vector<double>(n));

    // Taking user input for the matrix A
    cout << "Enter the elements of the matrix A:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
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

int main()
{
    char c;
    while (1)
    {
        cout << "want to solve a equation(y/n):" << endl;

        cout << "Enter your choice :";
        cin >> c;
        if (c == 'y')
        {
            cout << "1. Solution of Linear Equations" << endl
                 << "2. Solution of Non-linear Equations" << endl
                 << "3. Solution of Differential Equations" << endl
                 << "4.Matrix inversion" << endl;
            int select;
            cout << "Select what type of equation you want to solve:";
            cin >> select;

            if (select == 1)
            {
                slnoflinear();
                cout << endl;
            }
            else if (select == 2)
            {
                slnofnonlinear();
                cout << endl;
            }
            else if (select == 3)
            {
                slnofdiff();
                cout << endl;
            }
            else if (select == 4)
            {
                inversematrix();
            }
            else
                cout << "Nothing is seleceted" << endl;
        }
        else
        {
            cout << "No equation want to solve" << endl;
            break;
        }
    }

    return 0;
}
