#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <random>

#define N 5
#define NumMC 100000*N
#define HEATING 1
#define J 1 //('1'-ferro, '-1'-antiferro)

using namespace std;

//1D array
void Init_1D_Array(int *array_1D); //заполнение одномерного массива
double Full_Energy_Square_Ising_1D_array(int *array_1D); //функция подсчета энергии системы (2 соседа)
double Full_Energy_Square_Ising_recalculation_1D_array(int *array_1D, double old_E, int x); //функция пересчета энергии системы (2 соседа)

int main()
{
    time_t s_time = time(NULL);
    //выделение места в памяти под динамический одномерный массив:
    int *array_1D = new int [N];
    
    Init_1D_Array(array_1D);
    
    int a;
    double r, ver;
    double E_aver = 0, E_aver2 = 0;
    float T = 0.0001;
    double C;
    
    ofstream C_data("C.txt");
    ofstream E_data("E.txt");
    
    default_random_engine generator; //тип генератора случайных чисел
    generator.seed(time(NULL)); //стартовое число
    uniform_int_distribution<int> distribution_int(0,N-1); //целочисленное равномерное распределение
    uniform_real_distribution<double> distribution_real(0,1); //вещественное равномерное распределение
    
    //cout << "E_1D = " << Full_Energy_Square_Ising_1D_array(array_1D) << endl;
    //system("pause");
    
    double E_prev = Full_Energy_Square_Ising_1D_array(array_1D);
    cout << "E_start = " << E_prev << endl << endl;
    double E_new;
    
    for(T = 0.0001; T<4.01; T+=0.1)
    {
        cout << "T = " << T << endl;
        
        E_aver = 0;
        E_aver2 = 0;
        
        //прогрев
        if(HEATING == 1){
            cout << "Heating...\n";
            
            for(unsigned int MCS = 0; MCS<NumMC; MCS++){
                //Metropolis
                //--start
                a = distribution_int(generator);
                
                array_1D[a] *= (-1);
                E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_prev, a);
                
                if(E_new >= E_prev){
                    r = distribution_real(generator);
                    ver = exp(-((E_new-E_prev)/T));
                    if (r>ver){
                        array_1D[a] *= (-1);
                        E_new = E_prev;
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
            
            array_1D[a] *= (-1);
            E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_prev, a);
            
            if(E_new >= E_prev)
            {
                r = distribution_real(generator);
                ver = exp(-(E_new-E_prev)/T);
                if (r>ver)
                {
                    array_1D[a] *= (-1);
                    E_new = E_prev;
                }
            }
            E_prev = E_new;
            //--finish
            
            //средние
            E_aver += (E_new - E_aver) / (MCS + 1);
            E_aver2 += (E_new*E_new - E_aver2) / (MCS + 1);
        }
        
        C = ((E_aver2 - E_aver*E_aver)/(T*T))/(N);
        cout << "C = " << C << endl << endl;
        C_data << T << "\t" << C << endl;
        
        E_data << T << "\t" << E_aver << endl;
    }
    
    cout << "\nE_prev = " << E_prev << endl;
    //cout << "Full Energy = " << Full_Energy_Square_Ising_1D_array(array_1D) << endl;
    
    C_data.close();
    E_data.close();
    
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
    
    for(int i=0; i<N; i++){
        rand_n = distribution_int(generator);
        
        while(rand_n == 0)
            rand_n = distribution_int(generator);
        
        array_1D[i]=rand_n;
        printf("%3d", array_1D[i]);
    }
    cout << endl;
}

//функция подсчета энергии системы для 1D массива (2 соседа)
double Full_Energy_Square_Ising_1D_array(int *array_1D)
{
    double E_sys = 0;
    int k[2];
    
    for(int i=0; i<N; i++)
    {
        if(i==0)
        {
            k[0] = i+1;
            k[1] = N-1;
        }
        else if(i==N-1)
        {
            k[0] = 0;
            k[1] = i-1;
        }
        else if (i>0 && i<N-1)
        {
            k[0] = i+1;
            k[1] = i-1;
        }
        
        for(int d=0; d<2; d++)
            E_sys += -J*(array_1D[i] * array_1D[k[d]]);
    }
    return E_sys/2;
}


//функция пересчета энергии системы для 1D массива (2 соседа)
double Full_Energy_Square_Ising_recalculation_1D_array(int *array_1D, double old_E, int x)
{
    double E_sys = 0;
    int k[2];
    int i=x;
    
    if(i==0)
    {
        k[0] = i+1;
        k[1] = N-1;
    }
    else if(i==N-1)
    {
        k[0] = 0;
        k[1] = i-1;
    }
    else if (i>0 && i<N-1)
    {
        k[0] = i+1;
        k[1] = i-1;
    }
    
    for(int d=0; d<2; d++)
        E_sys += -J*(array_1D[i] * array_1D[k[d]]);
    
    return old_E+2*E_sys;
}
