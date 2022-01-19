//Do not edit the code below (unless you know what you are doing)

#include"solution.h"

int solution::f_calls = 0;
int solution::g_calls = 0;
int solution::H_calls = 0;

solution::solution(double L)
{
	x = L;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(const matrix& A)
{
	x = A;
	g = NAN;
	H = NAN;
	y = NAN;
}

solution::solution(int n, double* A)
{
	x = matrix(n, A);
	g = NAN;
	H = NAN;
	y = NAN;
}

int get_dim(const solution& A)
{
	return get_len(A.x);
}

void solution::clear_calls()
{
	f_calls = 0;
	g_calls = 0;
	H_calls = 0;
}

ostream& operator<<(ostream& S, const solution& A)
{
	S << "x = " << A.x << endl;
	S << "y = " << A.y << endl;
	S << "f_calls = " << solution::f_calls << endl;
	if (solution::g_calls > 0)
		S << "g_calls = " << solution::g_calls << endl;
	if (solution::H_calls)
		S << "H_calls = " << solution::H_calls << endl;
	return S;
}

//You can edit the following code

void solution::fit_fun(matrix* ud, matrix* ad)
{
	++f_calls;
#if LAB_NO==2 && (LAB_PART==1 || LAB_PART==2)
	y = -cos(0.1 * x()) * exp(-pow(0.1 * x() - 2 * 3.14, 2)) + 0.002 * pow(0.1 * x(), 2);
#elif LAB_NO==2 && LAB_PART==3
	matrix Y0 = matrix(3, new double[3]{ 5,1,10 });
	matrix* Y = solve_ode(0, 1, 1000, Y0, ud, &x);
	int n = get_len(Y[0]);
	double max = Y[1](0, 2);
	for (int i = 1; i < n; ++i)
		if (max < Y[1](i, 2))
			max = Y[1](i, 2);
	*ud = Y[1];
	y = abs(max - 50);
#elif LAB_NO==3 && (LAB_PART==1 || LAB_PART==2)
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
#elif LAB_NO==3 && LAB_PART==3
	matrix Y0(2, 1);
	matrix* Y = solve_ode(0, 0.1, 100, Y0, ud, &x);
	double a_ref = 3.14, omg_ref = 0;
	int n = get_len(Y[0]);
	y = 0;
	for (int i = 0; i < n; ++i) {
		y = y + 10 * pow(a_ref - Y[1](i, 0), 2) + pow(omg_ref - Y[1](i, 1), 2) + pow(x(0) * (a_ref - Y[1](i, 0)) + x(1) * (omg_ref - Y[1](i, 1)), 2);
	}
	y = y * 0.1;
	//(*ud) = Y[1];
#elif LAB_NO==4 && LAB_PART==1
	double arg = 3.14 * sqrt(pow(x(0) / 3.14, 2) + pow(x(1) / 3.14, 2));
	y = sin(arg) / arg;
	if ((*ad)(1) > 1) //KARA ZEWNETRZNA
	{
		if (-x(0) + 1 > 0)
			y = y + (*ad)(0) * pow(-x(0) + 1, 2);
		if (-x(1) + 1 > 0)
			y = y + (*ad)(0) * pow(-x(1) + 1, 2);
		if (norm(x) - (*ud)(0) > 0)
			y = y + (*ad)(0) * pow(norm(x) - (*ud)(0), 2);
	}
	else {
		if (-x(0) + 1 > 0)
			y = 1e10;
		else
			y = y - (*ad)(0) / (-x(0) + 1);
		if (-x(1) + 1 > 0)
			y = 1e10;
		else
		{
			y = y - (*ad)(0) / (-x(1) + 1);
		}
		if (norm(x) - (*ud)(0) > 0)
			y = 1e10;
		else
			y = y - (*ad)(0) / (norm(x) - (*ud)(0));
	}
#elif LAB_NO==4 && LAB_PART==2
	matrix Y0(4, new double[4]{ 0,x(0), 100, 0 });
	matrix omega(x(1));
	matrix* Y = solve_ode(0, 0.01, 7, Y0, &omega);
	int n = get_len(Y[0]);
	int i0 = 0, i50 = 0;
	for (int i = 1; i < n; ++i)
	{
		if (abs(Y[1](i, 2) - 50) < abs(Y[1](i50, 2) - 50))
			i50 = i;
		if (abs(Y[1](i, 2)) < abs(Y[1](i0, 2)))
			i0 = i;
	}
	y = -Y[1](i0, 0);
	if (abs(x(0)) - 10 > 0)
		y = y + (*ad)(0) * pow(abs(x(0)) - 10, 2);
	if (abs(x(1) - 20 > 0))
		y = y + (*ad)(0) * pow(abs(x(1)) - 20, 2);
	if (abs(Y[1](i50, 0) - 5) - 1 > 0)
		y = y + (*ad)(0) * pow(abs(Y[1](i50, 0) - 5) - 1, 2);

#elif LAB_NO==5 && (LAB_PART==1 || LAB_PART==2)
	if (ad == nullptr) { //liczymy f(x)
		y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
	}
	else //g(alfa)
	{
		solution temp;
		temp.x = ad[0] + x * ad[1]; //ad[0] x, x = alfa, ad[1] = d
		temp.fit_fun(ud);
		y = temp.y;
		--f_calls;
	}
#elif LAB_NO==5 && LAB_PART==3
	int m = 100, n = get_len(x); //x - theta
	static matrix X(n, m), Y(1, m);
	if (f_calls == 1) {
		ifstream S("XData.txt");
		S >> X;
		S.close();
		S.open("YData.txt");
		S >> Y;
		S.close();
		cout << X << endl << endl;
		cout << Y << endl << endl;
	}
	double h;
	y = 0;
	for (int i = 0; i < m; ++i) {
		h = (trans(x) * X[i])(); //pusty nawias okragly pozwoli na wyjecie double z macierzy
		h = 1 / (1 + exp(-h));
		y = y - Y(0, i) * log(h) - (1 - Y(0, i)) * log(1 - h);

	}
	y = y / m;
#elif LAB_NO==6 && LAB_PART==1
	if (ad == nullptr) {
		y = matrix(2, 1);
		y(0) = ud[0]() * (pow(x(0) - 2, 2) + pow(x(1) - 2, 2));
		y(1) = 1 / ud[0]() * (pow(x(0) + 2, 2) + pow(x(1) + 2, 2));
	}
	else { //y = g(alfa)
		solution temp;
		temp.x = ad[0] + x * ad[1]; //x = alfa - d³ugoœæ kroku
		temp.fit_fun(ud);
		y = ud[1]() * temp.y(0) + (1 - ud[1]()) * temp.y(1);
		--f_calls;
	}
#elif LAB_NO==6 && LAB_PART==2
	if (ad == nullptr) {
		y = matrix(3, 1);
		double ro = 7800, P = 1e3, E = 207e9;
		y(0) = ro * x(0) * 3.14 * pow(x(1), 2) / 4; //masa
		y(1) = 64 * P * pow(x(0), 3) / (3 * E * 3.14 * pow(x(1), 4)); //ugiêcie
		y(2) = 32 * P * x(0) / (3.14 * pow(x(1), 3)); //naprê¿enie
	}
	else { //g(alfa)
		solution T;
		T.x = ad[0] + x * ad[1];
		T.fit_fun(ud);
		matrix yn(2, 1);
		yn(0) = (T.y(0) - 0.12) / (3.06 - 0.12);
		yn(1) = (T.y(1) - 4.2e-5) / (0.026 - 4.2e-5);
		y = (*ud)() * yn(0) + (1 - (*ud)()) * yn(1);
		//y = (*ud)() * T.y(0) + (1 - (*ud)()) * T.y(1);

		double c = 1e10;
		if (T.x(0) < 0.2) {
			y = y + c * pow(0.2 - T.x(0), 2);
		}
		if (T.x(0) > 1) {
			y = y + c * pow(T.x(0) - 1, 2);
		}
		if (T.x(1) < 0.01) {
			y = y + c * pow(0.01 - T.x(1), 2);
		}
		if (T.x(1) > 0.05) {
			y = y + c * pow(0.05 - T.x(1), 2);
		}
		if (T.y(1) > 0.05) {
			y = y + c * pow(T.y(1) - 0.05, 2);
		}
		if (T.y(2) > 3e8) {
			y = y + c * pow(T.y(2) - 3e8, 2);
		}
		--f_calls;
	}

#elif LAB_NO==7 && LAB_PART==1
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * 3.14 * x(0)) - cos(2.5 * 3.14 * x(1)) + 2;
#elif LAB_NO==7 && LAB_PART==2
	int N = 1001; //N liczba punktów
	static matrix X(N, 2);
	if (solution::f_calls == 1) {
		ifstream S("polozenia.txt");
		S >> X;
		//cout << X;
		S.close();
	}
	matrix Y0(4, new double[4]{ 0,0,0,0 });
	matrix* Y = solve_ode(0, 0.1, 100, Y0, &x[0]);	//w x(0) s¹ wartoœci b1, b2, w x[1] s¹ odchylenia standardowe
	y = 0;
	for (int i = 0; i < N; ++i) {
		y = y + abs(X(i, 0) - Y[1](i, 0)) + abs(X(i, 1) - Y[1](i, 2)); //X - polozenie ciezarkow z pliku
		//cout << y << endl;
	}
	//cout << y;
	y = y / (2 * N); //2N bo policzylismy 2N bledow, N punktow dla kazdego 2 roznice

#endif
}

void solution::grad(matrix* ud, matrix* ad)
{
	++g_calls;
#if LAB_NO==5 && (LAB_PART==1 || LAB_PART==2)
	g = matrix(2, 1);
	g(0) = 10 * x(0) + 8 * x(1) - 34;
	g(1) = 8 * x(0) + 10 * x(1) - 38;
#elif LAB_NO==5 && LAB_PART==3
	int m = 100, n = get_len(x);
	static matrix X(n, m), Y(1, m);
	double h;
	g = matrix(n, 1);
	if (g_calls == 1) {
		ifstream S("XData.txt");
		S >> X;
		S.close();
		S.open("YData.txt");
		S >> Y;
		S.close();
		cout << X << endl << endl;
		cout << Y << endl << endl;
	}
	for (int j = 0; j < n; ++j) {
		for (int i = 0; i < m; ++i)
		{
			h = (trans(x) * X[i])(); //pusty nawias okragly pozwoli na wyjecie double z macierzy
			h = 1 / (1 + exp(-h));
			g(j) = g(j) + X(j, i) * (h - Y(0, i));
		}
		g(j) = g(j) / m;
}
#endif
}

void solution::hess(matrix* ud, matrix* ad)
{
	++H_calls;
#if LAB_NO==5 && (LAB_PART==1 || LAB_PART==2)
	H = matrix(2, 2);
	H(0, 0) = H(1, 1) = 10;
	H(0, 1) = H(1, 0) = 8;
#endif
}