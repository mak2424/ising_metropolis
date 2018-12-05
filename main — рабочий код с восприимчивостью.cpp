#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <random>

#define N 32
#define NumMC 50000*N*N
#define HEATING 1
#define J 1 //('1'-ferro, '-1'-antiferro)

int My = 0;

using namespace std;

//1D массив
void Init_1D_Array(int *array_1D); //заполнение одномерного массива
double Full_Energy_Square_Ising_1D_array(int *array_1D); //функция подсчета энергии системы (4 соседа)
double Full_Energy_Square_Ising_recalculation_1D_array(int *array_1D, double old_E, int x); //функция пересчета энергии системы (4 соседа)

//2D массив
void Init_2D_Array(int **array_2D); //заполнение двумерного массива
double Full_Energy_Square_Ising_2D_array(int **array_2D); //функция подсчета энергии системы (4 соседа)
double Full_Energy_Square_Ising_recalculation_2D_array(int **array_2D, double old_E, int x, int y); //функция пересчета энергии системы (4 соседа)


int main()
{
    time_t s_time = time(NULL);
    //int mas[N][N];
    //выделение места в памяти под динамический одномерный массив:
    int *array_1D = new int [N*N];
    //выделение места в памяти под динамический двумерный массив:
    //int **array_2D=new int* [N];
    //for(int i=0; i<N; i++)
        //array_2D[i]=new int [N];
    
    Init_1D_Array(array_1D);
    //Init_2D_Array(array_2D);
    
    int a;//,b;
    double r, ver;
    double E_aver = 0, E_aver2 = 0;
    float T = 0.0001;
    double C, hi;
    
    ofstream C_data("C.txt");
    ofstream hi_data("hi.txt");
    ofstream E_data("E.txt");
    
    default_random_engine generator; //тип генератора случайных чисел
    generator.seed(time(NULL)); //стартовое число
    uniform_int_distribution<int> distribution_int(0,N*N-1); //целочисленное равномерное распределение
    uniform_real_distribution<double> distribution_real(0,1); //вещественное равномерное распределение
    
    //cout << "E_1D = " << Full_Energy_Square_Ising_1D_array(array_1D) << endl;
    //cout << "E_2D = " << Full_Energy_Square_Ising_2D_array(array_2D) << endl;
    //system("pause");
    
    double E_prev = Full_Energy_Square_Ising_1D_array(array_1D);
    //double E_prev = Full_Energy_Square_Ising_2D_array(array_2D);
    cout << "E_start = " << E_prev << endl << endl;
    double E_new;
    
    double Mpr_aver = 0, Mpr2_aver = 0;
    
    for(T = 0.0001; T<4.01; T+=0.01)
    {
        cout << "T = " << T << endl;
        
        E_aver = 0;
        E_aver2 = 0;
        Mpr_aver = 0;
        Mpr2_aver = 0;
        
        //прогрев
        if(HEATING == 1){
            cout << "Heating...\n";
            
            for(unsigned int MCS = 0; MCS<NumMC; MCS++){
                //Metropolis
                //--start
                a = distribution_int(generator);
                //b = distribution_int(generator);
                
                //array_2D[a][b] *= (-1);
                //E_new = Full_Energy_Square_Ising_recalculation_2D_array(array_2D, E_prev, a, b);
                array_1D[a] *= (-1);
                E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_prev, a);
                
                if(E_new >= E_prev){
                    r = distribution_real(generator);
                    ver = exp(-((E_new-E_prev)/T));
                    if (r>ver){
                        //array_2D[a][b] *= (-1);
                        array_1D[a] *= (-1);
                        //E_new = E_prev;
                        E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_new, a);
                    }
                }
                E_prev = E_new;
                //--finish
            }
        }
        
        cout << "Monte Carlo...\n";
        for(unsigned int MCS = 0; MCS<NumMC; MCS++)
        {
            //Metropolis
            //--start
            a = distribution_int(generator);
            //b = distribution_int(generator);
            
            //array_2D[a][b] *= (-1);
            //E_new = Full_Energy_Square_Ising_recalculation_2D_array(array_2D, E_prev, a, b);
            array_1D[a] *= (-1);
            E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_prev, a);
            
            if(E_new >= E_prev)
            {
                r = distribution_real(generator);
                ver = exp(-(E_new-E_prev)/T);
                if (r>ver)
                {
                    //array_2D[a][b] *= (-1);
                    array_1D[a] *= (-1);
                    E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_new, a);
                    //E_new = E_prev;
                }
            }
            E_prev = E_new;
            //--finish
            
            //средние
            E_aver += (E_new - E_aver) / (MCS + 1);
            E_aver2 += (E_new*E_new - E_aver2) / (MCS + 1);
            
//            My = 0;
//            ///ПЕРЕНЕС В ПОДСЧЕТ ЭНЕРГИИ ДЛЯ 1D 
//            for(int i = 0; i<N*N; i++){
//                My += array_1D[i];
//            }
            Mpr_aver += (abs(My) - Mpr_aver) / (double)(MCS + 1);
            Mpr2_aver += ((My*My) - Mpr2_aver) / (double)(MCS + 1);
        }
        
        //cout << "Mpr_aver = " << Mpr_aver << endl << endl;
        
        C = ((E_aver2 - E_aver*E_aver)/(T*T))/(N*N);
        cout << "C = " << C << endl;
        C_data << T << "\t" << C << endl;
        
        hi = ((Mpr2_aver - Mpr_aver*Mpr_aver)/T) / (double)(N*N);
        cout << "hi = " << hi << endl << endl;
        hi_data << T << "\t" << hi << endl;
        
        
        E_data << T << "\t" << E_aver << endl;
    }
    
    cout << "\nE_prev = " << E_prev << endl;
    //cout << "Full Energy = " << Full_Energy_Square_Ising_1D_array(array_1D) << endl;
    //cout << "Full Energy = " << Full_Energy_Square_Ising_2D_array(array_2D) << endl;
    
    C_data.close();
    hi_data.close();
    E_data.close();
    
    // удаление двумерного динамического массива
    //for (int i = 0; i < N; i++)
        //delete [] array_2D[i];
    // удаление одномерного динамического массива
    delete [] array_1D;
    
    cout << "\n" << time(NULL) - s_time << " sec.\n";
    cout << "\nfinish\n";
    return 0;
}



///--------------------//
/// 
//заполнение одномерного массива
void Init_1D_Array(int *array_1D){
    
    default_random_engine generator; //тип генератора случайных чисел
    generator.seed(time(NULL)); //стартовое число
    uniform_int_distribution<int> distribution_int(-1,1); //целочисленное равномерное распределение
    
    int rand_n;
    
    for(int i=0; i<N*N; i++){
        rand_n = distribution_int(generator);
        
        while(rand_n == 0)
            rand_n = distribution_int(generator);
        
        array_1D[i]=1;
        //array_1D[i]=rand_n;
        //printf("%3d", mas[i]);
    }
    //cout << endl;
}

//заполнение двумерного массива
void Init_2D_Array(int **array_2D){

    default_random_engine generator; //тип генератора случайных чисел
    generator.seed(time(NULL)); //стартовое число
    uniform_int_distribution<int> distribution_int(-1,1); //целочисленное равномерное распределение
    
    int rand_n;
    
    for(int i=0; i<N; i++){
        for(int j=0; j<N; j++){
            rand_n = distribution_int(generator);
            
            while(rand_n == 0)
                rand_n = distribution_int(generator);
            
            array_2D[i][j]=rand_n;
            //printf("%3d", mas[i][j]);
        }
        //cout << endl;
    }
    //cout << endl;
}

//функция подсчета энергии системы для 2D массива (4 соседа)
double Full_Energy_Square_Ising_2D_array(int **array_2D)
{
    double E_sys = 0;
    int k[4],l[4];
    
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            if(i==0)
            {
                k[0] = i;
                k[1] = i+1;
                k[2] = i;
                k[3] = N-1;
            }
            else if(i==N-1)
            {
                k[0] = i;
                k[1] = 0;
                k[2] = i;
                k[3] = i-1;
            }
            else if (i>0 && i<N-1)
            {
                k[0] = i;
                k[1] = i+1;
                k[2] = i;
                k[3] = i-1;
            }
            
            if(j==0)
            {
                l[0] = j+1;
                l[1] = j;
                l[2] = N-1;
                l[3] = j;
            }
            else if(j==N-1)
            {
                l[0] = 0;
                l[1] = j;
                l[2] = j-1;
                l[3] = j;
            }
            else if (j>0 && j<N-1)
            {
                l[0] = j+1;
                l[1] = j;
                l[2] = j-1;
                l[3] = j;
            }
            
            for(int d=0; d<4; d++)
                E_sys += -J*(array_2D[i][j] * array_2D[k[d]][l[d]]);
        }
    }
    return E_sys/2;
}


//функция пересчета энергии системы для 2D массива (4 соседа)
double Full_Energy_Square_Ising_recalculation_2D_array(int **array_2D, double old_E, int x, int y)
{
    double E_sys = 0;
    int k[4];
    int l[4];
    int i=x;
    int j=y;
    
    if(i==0)
    {
        k[0] = i;
        k[1] = i+1;
        k[2] = i;
        k[3] = N-1;
    }
    else if(i==N-1)
    {
        k[0] = i;
        k[1] = 0;
        k[2] = i;
        k[3] = i-1;
    }
    else if (i>0 && i<N-1)
    {
        k[0] = i;
        k[1] = i+1;
        k[2] = i;
        k[3] = i-1;
    }
    
    if(j==0)
    {
        l[0] = j+1;
        l[1] = j;
        l[2] = N-1;
        l[3] = j;
    }
    else if(j==N-1)
    {
        l[0] = 0;
        l[1] = j;
        l[2] = j-1;
        l[3] = j;
    }
    else if (j>0 && j<N-1)
    {
        l[0] = j+1;
        l[1] = j;
        l[2] = j-1;
        l[3] = j;
    }
    
    for(int d=0; d<4; d++)
        E_sys += -J*(array_2D[i][j] * array_2D[k[d]][l[d]]);
    
    return old_E+2*E_sys;
}

//функция подсчета энергии системы для 1D массива (4 соседа)
double Full_Energy_Square_Ising_1D_array(int *array_1D)
{
    double E_sys = 0;
    int k[4],l[4];
    
    for(int i=0; i<N; i++)
    {
        for(int j=0; j<N; j++)
        {
            if(i==0)
            {
                k[0] = i;
                k[1] = i+1;
                k[2] = i;
                k[3] = N-1;
            }
            else if(i==N-1)
            {
                k[0] = i;
                k[1] = 0;
                k[2] = i;
                k[3] = i-1;
            }
            else if (i>0 && i<N-1)
            {
                k[0] = i;
                k[1] = i+1;
                k[2] = i;
                k[3] = i-1;
            }
            
            if(j==0)
            {
                l[0] = j+1;
                l[1] = j;
                l[2] = N-1;
                l[3] = j;
            }
            else if(j==N-1)
            {
                l[0] = 0;
                l[1] = j;
                l[2] = j-1;
                l[3] = j;
            }
            else if (j>0 && j<N-1)
            {
                l[0] = j+1;
                l[1] = j;
                l[2] = j-1;
                l[3] = j;
            }
            
            for(int d=0; d<4; d++)
                E_sys += -J*(array_1D[i*N+j] * array_1D[k[d]*N+l[d]]);
        }
    }
    
    My = 0;
    for(int i = 0; i<N*N; i++){
        My += array_1D[i];
    }
    
    return E_sys/2;
}


//функция пересчета энергии системы для 1D массива (4 соседа)
double Full_Energy_Square_Ising_recalculation_1D_array(int *array_1D, double old_E, int x)
{
    double E_sys = 0;
    int k[4];
    int l[4];
    int j=x%N;
    int i=(x-j)/N;
    
    if(i==0)
    {
        k[0] = i;
        k[1] = i+1;
        k[2] = i;
        k[3] = N-1;
    }
    else if(i==N-1)
    {
        k[0] = i;
        k[1] = 0;
        k[2] = i;
        k[3] = i-1;
    }
    else if (i>0 && i<N-1)
    {
        k[0] = i;
        k[1] = i+1;
        k[2] = i;
        k[3] = i-1;
    }
    
    if(j==0)
    {
        l[0] = j+1;
        l[1] = j;
        l[2] = N-1;
        l[3] = j;
    }
    else if(j==N-1)
    {
        l[0] = 0;
        l[1] = j;
        l[2] = j-1;
        l[3] = j;
    }
    else if (j>0 && j<N-1)
    {
        l[0] = j+1;
        l[1] = j;
        l[2] = j-1;
        l[3] = j;
    }
    
    for(int d=0; d<4; d++)
        E_sys += -J*(array_1D[i*N+j] * array_1D[k[d]*N+l[d]]);
    
    My += 2*array_1D[x];
    
    return old_E+2*E_sys;
}
