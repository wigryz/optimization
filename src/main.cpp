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
        double t0 = 0, s = 0.1, tend = 50;
        matrix Y0 = matrix(2, new double[2]{0, 0}); //pierwsza wartosc x1 od 0, druga wartosc predkosc od 0(??)
        matrix *Y1 = solve_ode(t0, s, tend,
                              Y0); //funkcja solve_ode zwraca dwie macierze. Pierwsza to czas, druga rozwi�zania w kroku czasowym
        matrix out = hcat(Y1[0], Y1[1]);
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


#elif LAB_NO == 2 && LAB_PART == 3

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
        matrix symulacja(get_len(Y_F[0]),7);
        symulacja(0, 0) = matrix Y_F; //I tu jakos te wartosci trzeba wpisac do tego Mecierza/Tablicy symulacja (7 wartosci)


        ofstream sout("wyniki1.csv");
        sout << symulacja;

#elif LAB_NO == 3 && LAB_PART == 1
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
#elif LAB_NO == 3 && LAB_PART == 2
        ofstream sout("wynikiPart2.csv");
        double s = 0.1, alphaHJ = 0.5, alphaR = 2, beta = 0.5, epsilon = 1e-3;
        int Nmax = 5000;
        matrix x0 = 2 * rand_mat(2, 1) - 1;
        matrix s0(2, 1, s);

        matrix XSHJ = trans(x0);
        matrix XSRosen = trans(x0);

        solution optHJ = HJ(x0, s, alphaHJ, epsilon, Nmax, &XSHJ);
        cout << XSHJ << endl << endl;

        solution::clear_calls();
        solution optR = Rosen(x0, s0, alphaR, beta, epsilon, Nmax, &XSRosen);
        cout << XSRosen << endl << endl;
        solution::clear_calls();
#elif LAB_NO == 3 && LAB_PART == 3
        double s = 0.1, alphaHJ = 0.5, alphaR = 2, beta = 0.5, epsilon = 1e-3;
        int Nmax = 5000;
        matrix x0 = 10 * rand_mat(2, 1);
        matrix s0(2, 1, s);

        solution optHJ = HJ(x0, s, alphaHJ, epsilon, Nmax);
        cout << "HJ:\n" << optHJ << endl;

        solution::clear_calls();
        solution optR = Rosen(x0, s0, alphaR, beta, epsilon, Nmax);
        cout << "Rosen:\n" << optR << endl;
        solution::clear_calls();

        matrix Y0(2, 1);
        matrix ud;
        matrix* Y1 = solve_ode(0, 0.1, 100, Y0, &ud, &optHJ.x);
        matrix ud2;
        matrix* Y2 = solve_ode(0, 0.1, 100, Y0, &ud2, &optR.x);
        ofstream out("simHJ.csv");
        out << hcat(Y1[0], Y1[1]);
        out.close();
        out.open("simR.csv");
        out << hcat(Y2[0], Y2[1]);
        out.close();

#elif LAB_NO == 4 && LAB_PART == 1
        matrix x0, wyniki(300, 12);
        double a[3] = { 4, 4.4934, 5 }, c_ex = 1, c_in = 10, dc_ex = 2, dc_in = 0.5, epsilon = 1e-10;
        int Nmax = 10000;
        solution opt;
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 100; j++) {
                do
                    x0 = 5 * rand_mat(2, 1) + 1;
                while (norm(x0) > a[i]);
                wyniki(100 * i + j, 0) = x0(0);
                wyniki(100 * i + j, 1) = x0(1);
                matrix a1 = matrix(a[i]);
                opt = pen(x0, c_ex, dc_ex, epsilon, Nmax, &a1);
                wyniki(100 * i + j, 2) = opt.x(0);
                wyniki(100 * i + j, 3) = opt.x(1);
                wyniki(100 * i + j, 4) = norm(opt.x);
                wyniki(100 * i + j, 5) = opt.y();
                wyniki(100 * i + j, 6) = opt.f_calls;
                solution::clear_calls();
                opt = pen(x0, c_in, dc_in, epsilon, Nmax, &a1);
                wyniki(100 * i + j, 7) = opt.x(0);
                wyniki(100 * i + j, 8) = opt.x(1);
                wyniki(100 * i + j, 9) = norm(opt.x);
                wyniki(100 * i + j, 10) = opt.y();
                wyniki(100 * i + j, 11) = opt.f_calls;
                solution::clear_calls();
            }
            ofstream out("wyniki_tab_1.csv");
            out << wyniki;
            out.close();
        }
#elif LAB_NO == 4 && LAB_PART == 2
        matrix x0(2, 1);
        double a, c = 10, dc = 2, epsilon = 1e-4;
        int Nmax = 5000;
        x0(0) = 20 * rand_mat(1, 1)() - 10;
        x0(1) = 40 * rand_mat(1, 1)() - 20;
        cout << x0 << endl << endl;
        solution opt = pen(x0, c, dc, epsilon, Nmax);
        opt.y = -opt.y;
        cout << opt << endl;
        solution::clear_calls();
        matrix Y0(4, new double[4]{0, opt.x(0), 100, 0});
        matrix xd = matrix(opt.x(1));
        matrix *Y = solve_ode(0, 0.01, 7, Y0, &xd);
        matrix symulacja = Y[0];
        symulacja = hcat(symulacja, Y[1][0]);
        symulacja = hcat(symulacja, Y[1][2]);
        ofstream out("wyniki_sym.csv");
        out << symulacja;
        out.close();
#elif LAB_NO == 5 && LAB_PART == 1

        ofstream sout("wynikiPart1.csv");
        matrix x0, wyniki(100, 18);

        for (int j = 0; j < 100; j++) {
         
            matrix x0 = 20 * rand_mat(2, 1) - 10;
            
            wyniki(j, 0) = x0(0);   //punkty
            wyniki(j, 1) = x0(1); 

            double h0 = 0.12, epsilon = 1e-5;
            int Nmax = 10000;

            solution optSD, optCG, optN;
            optSD = SD(x0, h0, epsilon, Nmax);
            wyniki(j,2) = optSD.x(0);
            wyniki(j, 3) = optSD.x(1);
            wyniki(j, 4) = optSD.y();
            wyniki(j, 5) = optSD.f_calls;
            wyniki(j, 6) = optSD.g_calls;


            //cout << optSD << endl << endl;
            solution::clear_calls();
            optCG = CG(x0, h0, epsilon, Nmax);
            wyniki(j, 7) = optCG.x(0);
            wyniki(j, 8) = optCG.x(1);
            wyniki(j, 9) = optCG.y();
            wyniki(j, 10) = optCG.f_calls;
            wyniki(j, 11) = optCG.g_calls;

            //cout << optCG << endl << endl;
            solution::clear_calls();
            optN = Newton(x0, h0, epsilon, Nmax);
            wyniki(j, 12) = optN.x(0);
            wyniki(j, 13) = optN.x(1);
            wyniki(j, 14) = optN.y();
            wyniki(j, 15) = optN.f_calls;
            wyniki(j, 16) = optN.g_calls;
            wyniki(j, 17) = optN.H_calls;
            //cout << optN << endl << endl;
            solution::clear_calls();


        }
        sout << wyniki;
        sout.close();




#elif LAB_NO == 5 && LAB_PART == 2
        ofstream sout("wynikiPart2.csv");
        matrix wyniki(100, 18);
        matrix x0 = 20 * rand_mat(2, 1) - 10;
        for (int j = 0; j < 100; j++) {

            double h0 = 0.05, epsilon = 1e-5;
            int Nmax = 10000;
            solution optSD, optCG, optN;
            optSD = SD(x0, h0, epsilon, Nmax); // Stw�rz matrix ud pusty
            wyniki(j, 2) = optSD.x(0);
            wyniki(j, 3) = optSD.x(1);
            wyniki(j, 4) = optSD.y();
            wyniki(j, 5) = optSD.f_calls;
            wyniki(j, 6) = optSD.g_calls;
            //cout << optSD << endl << endl;
            solution::clear_calls();

            optCG = CG(x0, h0, epsilon, Nmax);
            wyniki(j, 7) = optCG.x(0);
            wyniki(j, 8) = optCG.x(1);
            wyniki(j, 9) = optCG.y();
            wyniki(j, 10) = optCG.f_calls;
            wyniki(j, 11) = optCG.g_calls;
            //cout << optCG << endl << endl;
            solution::clear_calls();

            optN = Newton(x0, h0, epsilon, Nmax);
            wyniki(j, 12) = optN.x(0);
            wyniki(j, 13) = optN.x(1);
            wyniki(j, 14) = optN.y();
            wyniki(j, 15) = optN.f_calls;
            wyniki(j, 16) = optN.g_calls;
            wyniki(j, 17) = optN.H_calls;
            //cout << optN << endl << endl;
            solution::clear_calls();
        }

        sout << wyniki;
        sout.close();

#elif LAB_NO == 5 && LAB_PART == 3
        matrix x0(3, new double[3]{ -1, 0.1, 0.1 });
        solution test(x0);
        test.fit_fun();
        test.grad();
        cout << test<<endl;
        cout << test.g << endl;

#elif LAB_NO == 6 && LAB_PART == 1

#elif LAB_NO == 6 && LAB_PART == 2

        // N - liczbaw wymiar�w,
        // limits - zakres rozwiazan poczatkowych,
        // mi - wielkosc populacji poczatkowej,
        // lambda - wielkosc pop potomnej,
        // sigma0 - odchylenie standardowe zmiany X
#elif LAB_NO == 7 && LAB_PART == 1
        int N = 2, Nmax = 5000, mi = 20, lambda = 40;
        double sigmaValue = 0.01;
        ofstream sout("tabela1_k6.csv");
        for (int i = 0; i < 5; i++) {
            sout << "sigma = " << sigmaValue << "\n";
            for (int j = 0; j < 100; j++) {
                double epsilon = 1e-3;
                matrix limits(2, 2), sigma0(2, 1);
                limits(0, 0) = limits(1, 0) = -5;
                limits(0, 1) = limits(1, 1) = 5;
                sigma0(0) = sigma0(1) = sigmaValue;
                solution optEA = EA(N, limits, mi, lambda, sigma0, epsilon, Nmax);
                sout << optEA.x(0) << "\t" << optEA.x(1) << "\t" << optEA.y(0) << "\t" << solution::f_calls << "\n";
                solution::clear_calls();
            }
            sigmaValue *= 10.0;
        }
        sout.close();
#elif LAB_NO == 7 && LAB_PART == 2
        //parametry poczatkowe
        int N = 2, Nmax = 5000, mi = 20, lambda = 40;
        double epsilon = 1e-3;
        matrix limits(2, 2), sigma0(2, 1);
        limits(0, 0) = limits(1, 0) = 0.1;
        limits(0, 1) = limits(1, 1) = 3;
        sigma0(0) = sigma0(1) = 10;
        //optymalizacja
        solution optEA = EA(N, limits, mi, lambda, sigma0, epsilon, Nmax);
        ofstream sout("symulacja.csv");
        sout << optEA.x(0, 0) << ";" << optEA.x(1, 0) << ";" << optEA.y(0) << ";" << solution::f_calls << endl;
        solution::clear_calls();
        matrix Y0(4, new double[4]{0, 0, 0, 0 });
        matrix ud(2, 1);
        ud(0) = optEA.x(0, 0);
        ud(1) = optEA.x(1, 0);
        //symulacja
        matrix* R = solve_ode(0, 0.1, 100, Y0, &ud);
        matrix result = hcat(R[1][0], R[1][2]);
        sout << result << endl;
        sout.close();
#endif
    }
    catch (char *EX_INFO) {
        cout << EX_INFO << endl;
    }
    return 0;
}
