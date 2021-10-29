/***************************************************
Code written for the optimization exercises purposes
by Lukasz Sztangret, PhD
Department of Applied Computer Science and Modelling
AGH University of Science and Technology
***************************************************/

#include"opt_alg.h"

int main() {
    try {
        cout << "LAB NUMBER " << LAB_NO << endl;
        cout << "LAB PART " << LAB_PART << endl << endl;
#if LAB_NO == 0

#elif LAB_NO == 1 && LAB_PART == 1
        double t0 = 0, dt = 0.1, tend = 50;
        matrix Y0 = matrix(2, new double[2]{0, 0}); //pierwsza wartosc x1 od 0, druga wartosc predkosc od 0(??)
        matrix *Y = solve_ode(t0, dt, tend,
                              Y0); //funkcja solve_ode zwraca dwie macierze. Pierwsza to czas, druga rozwi�zania w kroku czasowym
        matrix out = hcat(Y[0], Y[1]);
        ofstream sout("wyniki.csv");
        sout << out;
        sout.close();


#elif LAB_NO == 1 && LAB_PART == 2

#elif LAB_NO==2 && LAB_PART==1
        double x0;
        double alpha = 0.5 ,d = 1, epsilon = 1e-5, gamma = 1e-200;
        int Nmax = 1000;

        matrix tab(300,12);
        for (int i = 1; i <= 3; i++) {
            for (int j = 0; j < 100; j++) {
                x0 = 200.0*rand_mat(1, 1)()-100.0;
                double* p = expansion(x0, d, alpha, Nmax);
                tab(j*i, 0) = x0;
                tab(j*i,1) = p[0];
                tab(j*i,2) = p[1];
                tab(j*i, 3) = solution::f_calls;
                solution::clear_calls();

                solution opt_F = fib(p[0], p[1], epsilon);
                tab(j + i, 4) = opt_F.x();
                tab(j + i, 5) = opt_F.y();
                tab(j + i, 6) = solution::f_calls;
                solution::clear_calls();

                solution opt_L = lag(p[0], p[1], epsilon, gamma, Nmax);
                tab(j + i, 7) = opt_L.x();
                tab(j + i, 8) = opt_L.y();
                tab(j + i, 9) = solution::f_calls;
                solution::clear_calls();
            }

            d += 0, 5;
        }
        ofstream sout("wyniki.csv");
        double *p = expansion(x0, d, alpha, Nmax);

        sout << tab;
        sout.close();


#elif LAB_NO==2 && LAB_PART==2
        
        double x0 = -20, d = 1, alpha = 2, epsilon = 1e-5, gamma = 1e-200;
        int Nmax = 1000;

        matrix tab(300, 12);
            for (int j = 0; j < 100; j++) {
               
                solution opt_F = fib(-100, 100, epsilon);
                tab(j + i, 4) = opt_F.x();
                tab(j + i, 5) = opt_F.y();
                tab(j + i, 6) = solution::f_calls;
                solution::clear_calls();

                solution opt_L = lag(-100, 100, epsilon, gamma, Nmax);
                tab(j + i, 7) = opt_L.x();
                tab(j + i, 8) = opt_L.y();
                tab(j + i, 9) = solution::f_calls;
                solution::clear_calls();
            }

            d += 0, 5;
        }
        ofstream sout("wyniki1.csv");
        double* p = expansion(x0, d, alpha, Nmax);


        matrix ab_F(1, 1, 200);

        solution opt_F = fib(-100, 100, epsilon, &ab_F);
        std::cout << ab_F;

        solution::clear_calls();


        matrix ab_L(1, 1, 200);

        solution::clear_calls();

        solution opt_L = lag(-100, 100, epsilon, gamma, Nmax, &ab_L);
        std::cout << ab_L;*/


#elif LAB_NO==2 && LAB_PART==3

        //ustalany zmienne
        solution test(0.001);
        test.fit_fun();


        solution opt_f = fib();
        solution::clear_calls();
        solution opt_l = lag();
        //fib
        //lab
        matrix Y0 = matrix(3, new double[3], { 5,1,10 });

        matrix *Y_F = matrix(0,1,1000 Y0, nullptr, %opt_f.x);
        matrix* Y_L = matrix(0, 1, 1000 Y0, nullptr, % opt_l.x);
        //matrixt (5 1 10)
        // solve_ode (dla lab i fab)
        matrix SYM(get_len(Y_F[0], 7));
        sym
        cout << SYM;

#elif LAB_NO==3 && LAB_PART==1

#elif LAB_NO==3 && LAB_PART==2

#elif LAB_NO==3 && LAB_PART==3

#elif LAB_NO==4 && LAB_PART==1

#elif LAB_NO==4 && LAB_PART==2

#elif LAB_NO==5 && LAB_PART==1

#elif LAB_NO==5 && LAB_PART==2

#elif LAB_NO==5 && LAB_PART==3

#elif LAB_NO==6 && LAB_PART==1

#elif LAB_NO==6 && LAB_PART==2

#elif LAB_NO==7 && LAB_PART==1

#elif LAB_NO==7 && LAB_PART==2

#endif
    }
    catch (char *EX_INFO) {
        cout << EX_INFO << endl;
    }
    system("pause");
    return 0;
}
