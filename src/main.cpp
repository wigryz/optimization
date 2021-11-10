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

#elif LAB_NO == 2 && LAB_PART == 1
        double x0;
        double alpha = 1.0, d = 1, epsilon = 1e-5, gamma = 1e-200;
        int Nmax = 1000;

        double diffBA[3] = {0.0, 0.0, 0.0};
        int l_wywolan[3] = {0,0,0};

        double min_glob_fibx[3] = {0.0, 0.0, 0.0};
        double min_lok_fibx[3] = {0.0, 0.0, 0.0};
        double min_glob_fiby[3] = {0.0, 0.0, 0.0};
        double min_lok_fiby[3] = {0.0, 0.0, 0.0};
        int l_fib[3] = {0,0,0};


        double min_glob_lagx[3] = {0.0, 0.0, 0.0};
        double min_lok_lagx[3] = {0.0, 0.0, 0.0};
        double min_glob_lagy[3] = {0.0, 0.0, 0.0};
        double min_lok_lagy[3] = {0.0, 0.0, 0.0};
        int l_lag[3] = {0,0,0};

        matrix tab(300, 12);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 100; j++) {
                x0 = 200.0 * rand_mat(1, 1)() - 100.0; //losowanko z przedzialu od -100 do 100
                double *p = expansion(x0, d, alpha, Nmax); //metoda exp
                tab(i * 100 + j, 0) = x0; //poczatkowy (losujemy)
                tab(i * 100 + j, 1) = p[0]; //przedzial obliczony metoda exp
                tab(i * 100 + j, 2) = p[1];
                tab(i * 100 + j, 3) = solution::f_calls; //ilosc wywołań
                diffBA[i] += p[1] - p[0];
                l_wywolan[i] += solution::f_calls;
                solution::clear_calls(); //sprzątanko

                solution opt_F = fib(p[0], p[1], epsilon);
                tab(i * 100 + j, 4) = opt_F.x();
                tab(i * 100 + j, 5) = opt_F.y();
                tab(i * 100 + j, 6) = solution::f_calls;
                tab(i * 100 + j, 7) = opt_F.y() < 0.5 ? 1 : 0; //
//                std::cout << i << " => " << opt_F << std::endl;
                min_glob_fibx[i] += opt_F.x();
                min_lok_fibx[i] = opt_F.y();
                min_glob_fiby[i] = opt_F.x();
                min_lok_fiby[i] = opt_F.y();
                l_fib[i] += solution::f_calls;
                solution::clear_calls();

                solution opt_L = lag(p[0], p[1], epsilon, gamma, Nmax);
                tab(i * 100 + j, 8) = opt_L.x();
                tab(i * 100 + j, 9) = opt_L.y();
                tab(i * 100 + j, 10) = solution::f_calls;

                min_glob_lagx[i] += opt_L.x();
                min_lok_lagx[i] = opt_L.y();
                min_glob_lagy[i] = opt_L.x();
                min_lok_lagy[i] = opt_L.y();
                l_fib[i] += solution::f_calls;
                tab(i * 100 + j, 11) = opt_F.y() < 0.5 ? 1 : 0; //liczba wywolan funkcji celu
                solution::clear_calls();
            }

            alpha += 0.5;
        }
        ofstream sout("wyniki.csv");

        sout << tab;
        for(int i=0 ; i < 3 ; i++) {
            std::cout << i << ". -> b-a eksp: " << diffBA[i] / 300.0 << " calls: " << (double) l_wywolan[i] / 300.0 << std::endl;
            std::cout << i << ". fib gl -> x*: " << min_glob_fibx[i] / 300.0 << " y*: " << min_glob_fiby[i] / 300.0 << " f_calls " << (double)l_fib[i]/300.0 << std::endl;
            std::cout << i << ". lab gl -> x*: " << min_glob_lagx[i] / 300.0 << " y*: " << min_glob_lagy[i] / 300.0 << " f_calls " << (double)l_lag[i]/300.0 << std::endl;
        }
        sout.close();


#elif LAB_NO == 2 && LAB_PART == 2

        double x0 = -20, d = 1, alpha = 2, epsilon = 1e-5, gamma = 1e-200;
        int Nmax = 1000;

        matrix tab(300, 12);
            for (int j = 0; j < 100; j++) {
               
                solution opt_F = fib(-100, 100, epsilon);
                tab(j, 0) = opt_F.x();
                tab(j , 1) = opt_F.y();
                tab(j , 2) = solution::f_calls;
                solution :: clear_calls();
                      
                solution opt_L = lag(-100, 100, epsilon, gamma, Nmax);
                tab(j , 4) = opt_L.x();
                tab(j , 5) = opt_L.y();
                tab(j , 6) = solution::f_calls;
                solution::clear_calls();
            }
        ofstream sout("wyniki2.csv");
        sout << tab;
        sout.close();


#elif LAB_NO==2 && LAB_PART==3

        //ustalany zmienne
        double alpha = 0.5, d = 1, epsilon = 1e-5, gamma = 1e-200;
        int Nmax = 1000;

        solution test(0.001);
        test.fit_fun();


        solution opt_f = fib(-100, 100, epsilon); //w sumie nie wiem jaki ma tu byc przedzial
        solution::clear_calls();
        solution opt_l = lag(-100, 100, epsilon, gamma, Nmax);
        //fib
   
        matrix Y0 = matrix(3, new double[3]{ 5,1,10 });

        //matrix *Y_F = matrix(0,1,1000, Y0, nullptr, &opt_f.x);   // tak przepisalam, ale chyba źle xd
        //matrix* Y_L = matrix(0, 1, 1000,Y0, nullptr, &opt_l.x);


        matrix* Y_F = solve_ode(0, 1, 1000, Y0, nullptr, &opt_f.x);
        matrix* Y_L = solve_ode(0, 1, 1000, Y0, nullptr, &opt_l.x);
        matrix SYM(get_len(Y_F[0]),7);
        SYM(0, 0) = matrix Y_F; //I tu jakos te wartosci trzeba wpisac do tego Mecierza/Tablicy SYM (7 wartosci)


        ofstream sout("wyniki1.csv");
        sout << SYM;

#elif LAB_NO==3 && LAB_PART==1
        srand(time(NULL));
        ofstream sout("wynikiPart1.csv");
        double s = 0.1, alphaHJ = 0.5, alphaR = 3, beta = 0.5, epsilon = 1e-3;
        int Nmax = 5000;

        for (int h = 0; h < 3; h++) {
            double s = (double)rand() / RAND_MAX;
            matrix s0(2, 1, s);

            cout << "Dlugosc kroku dla iteracji " << h << ".: " << s << "\n";
            for (int i = 0; i < 100; i++) {
                matrix x0 = 2 * rand_mat(2, 1) - 1;

                solution optHJ = HJ(x0, s, alphaHJ, epsilon, Nmax);
                sout << x0(0, 0) << ";" << x0(1, 0) << ";" << optHJ.x(0) << ";" << optHJ.x(1) << ";" << optHJ.y(0) << ";"
                    << solution::f_calls << ";" << ((optHJ.y(0) < 0.000001) ? 1 : 0) << ";";
                solution::clear_calls();

                solution optR = Rosen(x0, s0, alphaR, beta, epsilon, Nmax);
                sout << optR.x(0) << ";" << optR.x(1) << ";" << optR.y(0) << ";"
                    << solution::f_calls << ";" << ((optHJ.y(0) < 0.000001) ? 1 : 0) << "\n";
                solution::clear_calls();
            }
        }
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
    return 0;
}
