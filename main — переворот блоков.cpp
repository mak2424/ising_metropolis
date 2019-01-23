/** 
 * Список функций:
 * - Init_1D_Array - заполнение одномерного массива
 * - Full_Energy_Square_Ising_1D_array - функция подсчета энергии системы (4 соседа) для одномерного массива
 * - Full_Energy_Square_Ising_recalculation_1D_array - функция пересчета энергии системы (4 соседа) для одномерного массива
 * - Double_Spin_flip_ver1 - переворот двух спинов (версия 1) для одномерного массива
 * - Double_Spin_flip - переворот двух спинов (версия 2) для одномерного массива
 * - Block_Spin_flip - переворот блока спинов для одномерного массива
 * 
 * - Init_2D_Array - заполнение двумерного массива
 * - Full_Energy_Square_Ising_2D_array - функция подсчета энергии системы (4 соседа) для двумерного массива
 * - Full_Energy_Square_Ising_recalculation_2D_array - функция пересчета энергии системы (4 соседа) для двумерного массива
 * 
 * - Neighbors_searching - определение соседей одного спина
 * - Neighbors_saving - функция формирование массива соседей каждого спина
*/

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <fstream>
#include <random>
#include <map>
#include <vector>

typedef unsigned long long ulong;

#define N 10 // MAX 360
//#define NumMC 1*N*N
#define NumMC 1
#define HEATING 0
#define MONTE_CARLO 0
#define J 1 //('1'-ferro, '-1'-antiferro)

std::default_random_engine generator; //тип генератора случайных чисел

double E_aver = 0;
double E_aver2 = 0;
int My = 0;

using namespace std;

//1D массив
void Init_1D_Array(int *array_1D); //заполнение одномерного массива
double Full_Energy_Square_Ising_1D_array(int *array_1D); //функция подсчета энергии системы (4 соседа)
double Full_Energy_Square_Ising_recalculation_1D_array(int *array_1D, double old_E, int x); //функция пересчета энергии системы (4 соседа)

double Double_Spin_flip_ver1(int *array_1D, double old_E, int x, float T); //переворот двух спинов версия 1
double Double_Spin_flip(int *array_1D, double old_E, int x, float T); //переворот двух спинов
double Block_Spin_flip(int *array_1D, double old_E, float T, int **neighbors_array); //переворот блока спинов


//2D массив
void Init_2D_Array(int **array_2D); //заполнение двумерного массива
double Full_Energy_Square_Ising_2D_array(int **array_2D); //функция подсчета энергии системы (4 соседа)
double Full_Energy_Square_Ising_recalculation_2D_array(int **array_2D, double old_E, int x, int y); //функция пересчета энергии системы (4 соседа)

//общие
void Neighbors_searching(int *k, int *l, int i, int j); //определение соседей одного спина
void Neighbors_saving(int **neighbors_array); //функция формирование массива соседей каждого спина


int main()
{
    time_t s_time = time(NULL);
    generator.seed(1); //стартовое число
    srand(1);
    //int mas[N][N];
    
    //выделение места в памяти под динамический одномерный массив спинов:
    int *array_1D = new int [N*N];
    
    //выделение места в памяти под динамический двумерный массив соседей:
    int **neighbors_array = new int* [N*N];
    for (int i = 0; i < N*N; ++i)
        neighbors_array[i] = new int [4];
    
    //выделение места в памяти под динамический двумерный массив:
    //int **array_2D=new int* [N];
    //for(int i=0; i<N; i++)
        //array_2D[i]=new int [N];
    
    //заполняем массив спинов
    Init_1D_Array(array_1D);
    //Init_2D_Array(array_2D);
    
    //заполняем массив соседей
    Neighbors_saving(neighbors_array);
    
    int x;//,b;
    ///double r, ver;
    float T = 0.0001;
    double C, hi;
    
    ofstream C_data("C.txt");
    ofstream hi_data("hi.txt");
    ofstream E_data("E.txt");
    
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
        ///Block_Spin_flip(array_1D, E_prev, T);
        
        cout << "T = " << T << endl;
        
        E_aver = 0;
        E_aver2 = 0;
        Mpr_aver = 0;
        Mpr2_aver = 0;
        
        //прогрев
        if(HEATING)
        {
            cout << "Heating...\n";
            
            for(ulong MCS = 0; MCS<NumMC; MCS++)
            {
                //Metropolis
                //--start
                x = distribution_int(generator);
                //b = distribution_int(generator);
                
                //array_2D[a][b] *= (-1);
                //E_new = Full_Energy_Square_Ising_recalculation_2D_array(array_2D, E_prev, a, b);
                E_new = Double_Spin_flip(array_1D, E_prev, x, T);
                E_prev = E_new;
                
                if(E_new < (-2*N*N) || E_new > (2*N*N))
                {
                    cout<<"ERROR2 E = "<<E_new<<endl;
                    return 0;
                }
                
                /**
                array_1D[x] *= (-1);
                E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_prev, x);
                
                if(E_new >= E_prev){
                    r = distribution_real(generator);
                    ver = exp(-((E_new-E_prev)/T));
                    if (r>ver){
                        //array_2D[a][b] *= (-1);
                        array_1D[x] *= (-1);
                        //E_new = E_prev;
                        E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_new, x);
                    }
                }
                E_prev = E_new;*/
                //--finish
            }
        }
        
        if(MONTE_CARLO)
        {
            cout << "Monte Carlo...\n";
            for(ulong MCS = 0; MCS<NumMC; MCS++)
            {
                //Metropolis
                //--start
                x = distribution_int(generator);
                //b = distribution_int(generator);
                
                //array_2D[a][b] *= (-1);
                //E_new = Full_Energy_Square_Ising_recalculation_2D_array(array_2D, E_prev, a, b);
                
                ///E_new = Double_Spin_flip(array_1D, E_prev, x, T);
                E_new = Block_Spin_flip(array_1D, E_prev, T, neighbors_array);
                E_prev = E_new;
                
                ///cout << "E_new = "<<E_new<<endl;
                ///cout << "My = "<<My<<endl;
                ///cout<<"E_new(real) = "<<Full_Energy_Square_Ising_1D_array(array_1D)<<endl;
                ///cout << "My(real) = "<<My<<endl<<endl;
                
                /**
                for(int i=0;i<N*N;++i){
                    cout<<array_1D[i];
                    if((i+1)%N==0) cout<<endl;
                }
                //*/
                
                if(E_new < (-2*N*N) || E_new > (2*N*N))
                {
                    cout<<"ERROR2 E = "<<E_new<<endl;
                    return 0;
                }
                
                //system("pause");
                
                
                /**
                array_1D[x] *= (-1);
                E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_prev, x);
                
                if(E_new >= E_prev)
                {
                    r = distribution_real(generator);
                    ver = exp(-(E_new-E_prev)/T);
                    if (r>ver)
                    {
                        //array_2D[a][b] *= (-1);
                        array_1D[x] *= (-1);
                        E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_new, x);
                        //E_new = E_prev;
                    }
                }
                E_prev = E_new;
                //*/
                //--finish
                
                //средние
                //E_aver += (E_new - E_aver) / (MCS + 1);
                //E_aver2 += (E_new*E_new - E_aver2) / (MCS + 1);
                
                Mpr_aver += (abs(My) - Mpr_aver) / (double)(MCS + 1);
                Mpr2_aver += ((My*My) - Mpr2_aver) / (double)(MCS + 1);
                
            }
        }
        
        //cout << "Mpr_aver = " << Mpr_aver << endl << endl;
        
        cout << "E_aver = " << E_aver << endl << endl;
        ///system("pause");
        
        C = ((E_aver2 - E_aver*E_aver)/(T*T))/(N*N);
        cout << "C = " << C << endl;
        C_data << T << "\t" << C << endl;
        
        hi = ((Mpr2_aver - Mpr_aver*Mpr_aver)/T) / (double)(N*N);
        cout << "hi = " << hi << endl << endl;
        hi_data << T << "\t" << hi << endl;
        
        
        E_data << T << "\t" << E_aver << endl;
    }
    
    //cout << "\nE_prev = " << E_prev << endl;
    //cout << "Full Energy = " << Full_Energy_Square_Ising_1D_array(array_1D) << endl;
    //cout << "Full Energy = " << Full_Energy_Square_Ising_2D_array(array_2D) << endl;
    
    C_data.close();
    hi_data.close();
    E_data.close();
    
    
    /** Очистка памяти */
    
    // удаление двумерного динамического массива
    //for (int i = 0; i < N; i++)
        //delete [] array_2D[i];
    
    // удаление одномерного динамического массива
    delete [] array_1D;
    
    // удаление массива соседей
    for (int i = 0; i < N*N; ++i)
        delete [] neighbors_array[i];
    
    
    cout << "\n" << time(NULL) - s_time << " sec.\n";
    cout << "\nfinish\n";
    return 0;
}



///--------------------//
/// 
//заполнение одномерного массива
void Init_1D_Array(int *array_1D){
    //default_random_engine generator; //тип генератора случайных чисел
    //generator.seed(time(NULL)); //стартовое число
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
    //default_random_engine generator; //тип генератора случайных чисел
    //generator.seed(time(NULL)); //стартовое число
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


//определение соседей
void Neighbors_searching(int *k,int *l, int i, int j)
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
}

//функция формирование массива соседей каждого спина
void Neighbors_saving(int **neighbors_array)
{
    /** -------------- сохраняем соседей для каждого спина -------------->>>> */
    int k[4],l[4]; //координаты 4x соседей спина
    int i, j;
    
    for(int ii=0; ii<N*N; ++ii)
    {
        j=ii%N;
        i=(ii-j)/N;
        Neighbors_searching(k,l,i,j); //функция поиска соседей
        ///cout << ii << ": ";
        for(int neighbor_num=0; neighbor_num<4; ++neighbor_num)
        {
            neighbors_array[ii][neighbor_num] = k[neighbor_num]*N+l[neighbor_num];
            ///cout << neighbors[ii][jj] << ", ";
        }
        ///cout << endl;
    }
    /** 
 * на выходе получаем массив neighbors[][] с соседями каждого спина
 * Схема соседей:
 *      3
 *      |
 *  2 --S-- 0
 *      |
 *      1
*/
    /** <<<<----------------------------------------------------------------- */
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
            Neighbors_searching(k,l,i,j);
            
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
    
    Neighbors_searching(k,l,i,j);
    
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
            Neighbors_searching(k,l,i,j);
            
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
    
    Neighbors_searching(k,l,i,j);
    
    for(int d=0; d<4; d++)
        E_sys += -J*(array_1D[i*N+j] * array_1D[k[d]*N+l[d]]);
    
    My += 2*array_1D[x];
    
    return old_E+2*E_sys;
}

//переворот двух спинов версия 1
double Double_Spin_flip_ver1(int *array_1D, double old_E, int x, float T)
{
    ///default_random_engine generator; //тип генератора случайных чисел
    ///generator.seed(time(NULL)); //стартовое число
    
    //srand(time(NULL));
    int local_random_neighbor = rand()%4;
    ///cout<<endl<<"local_random_neighbor = "<<local_random_neighbor<<endl<<endl;
    int global_random_neighbor_num = 0;
    int state[4][5];
    
    double E_sys = 0;
    int k[4],l[4]; //координаты соседей
    int i, j;
    
    int configuration[2]={1};
    int min[2];
    
    double P12, P13, P14;
    double r;
    
    for(int state_num=0; state_num<4; state_num++)
    {
        E_sys = 0;
        j=x%N;
        i=(x-j)/N;
        
        for(int c=0; c<2; c++)
        {
            Neighbors_searching(k,l,i,j);
            
            for(int d=0; d<4; d++)
            {
                /**
                cout<<"spin= "<<i*N+j<<"("<<array_1D[i*N+j]<<
                      "), neigh= "<<k[d]*N+l[d]<<"("<<array_1D[k[d]*N+l[d]]<<")"<<
                      endl;*/
                
                if(c!=0 || d!=local_random_neighbor)
                {
                    E_sys += -J*(array_1D[i*N+j] * array_1D[k[d]*N+l[d]]);
                    //cout<<"E = "<<E_sys<<endl;
                    //system("pause");
                }
            }
            
            if(c==0){
                i = k[local_random_neighbor];
                j = l[local_random_neighbor];
                global_random_neighbor_num = i*N+j;
                ///cout<<"global_random_neighbor_num = "<<global_random_neighbor_num<<endl;
            }
        }
        
        state[state_num][0] = E_sys;
        state[state_num][1] = x;
        state[state_num][2] = array_1D[x];
        state[state_num][3] = global_random_neighbor_num;
        state[state_num][4] = array_1D[global_random_neighbor_num];
        
        configuration[0] = array_1D[x];
        configuration[1] = array_1D[global_random_neighbor_num];
        /**
        cout<<"configuration: <" << configuration[0] << configuration[1]<<">" << endl;
        
        cout << "E = "<<state[state_num][0]<<
                ", spin = "<<state[state_num][1]<<
                ", array_1D[x] = "<<state[state_num][2]<<
                ", global_random_neighbor_num = "<<state[state_num][3]<<
                ", array_1D[global_random_neighbor_num] = "<<state[state_num][4]<<
                endl<<endl;
        */
        for(int t = 0; t<2; t++)
        {
            if(configuration[t] == 1)
            {
                configuration[t] = -1;
                break;
            }
            else configuration[t] = 1;
        }
        array_1D[x] = configuration[0];
        array_1D[global_random_neighbor_num] = configuration[1];
        
    }
    
    min[0]=state[1][0];
    min[1]=1;
    
    for(int state_num=2; state_num<4; state_num++)
    {
        if(min[0]>state[state_num][0])
        {
            min[0] = state[state_num][0];//E
            min[1] = state_num;//номер состояния
        }
    }
    
    ///cout<<"Emin = "<<min[0]<<", conf = "<<min[1]<<endl;
    
    ///cout<<"x = "<<x<<", neig = "<<global_random_neighbor_num<<endl;
    
    if(min[0]<state[0][0])
    {
        /**
        cout<<
               "Emin < E1"<<endl<<
               "P=1"<<endl;*/
        
        array_1D[x] = state[min[1]][2];
        array_1D[global_random_neighbor_num] = state[min[1]][4];
        
        My += -(state[0][2]+state[0][4])+(array_1D[x] + array_1D[global_random_neighbor_num]);
        old_E += -state[0][0]+state[min[1]][0];
    }
    else
    {
        //T=10;//УБРАТЬ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        P12=exp((double)(state[0][0]-state[1][0])/T);
        P13=exp((double)(state[0][0]-state[2][0])/T);
        P14=exp((double)(state[0][0]-state[3][0])/T);
        /**
        cout<<
               "Emin > E1"<<endl<<
               "P12=exp((E1-E2)/T) = "<<P12<<endl<<
               "P13=exp((E1-E3)/T) = "<<P13<<endl<<
               "P14=exp((E1-E4)/T) = "<<P14<<endl<<
               endl;*/
        
        uniform_real_distribution<double> distribution_real(0,3); //вещественное равномерное распределение
        
        r = distribution_real(generator);
        ///cout<<"P12+P13+P14= "<<P12+P13+P14<<endl;
        ///cout<<"r= "<<r<<endl;
        
        if(r<3-(P12+P13+P14))
        //if(r<1-P12)
        {
            ///cout<<"r<3-(P12+P13+P14)"<<endl;
            
            array_1D[x] = state[0][2];
            array_1D[global_random_neighbor_num] = state[0][4];
        }
        else if(r<3-(P13+P14))
        //if(r<1)
        {
            ///cout<<"r<3-(P13+P14)"<<endl;
            
            array_1D[x] = state[1][2];
            array_1D[global_random_neighbor_num] = state[1][4];
            
            My += -(state[0][2]+state[0][4])+(array_1D[x] + array_1D[global_random_neighbor_num]);
            old_E += -state[0][0]+state[1][0];
        }
        else //if(r<3-P14)
        if(r<1+P13)
        {
            ///cout<<"r<3-P14"<<endl;
            
            array_1D[x] = state[2][2];
            array_1D[global_random_neighbor_num] = state[2][4];
            
            My += -(state[0][2]+state[0][4])+(array_1D[x] + array_1D[global_random_neighbor_num]);
            old_E += -state[0][0]+state[2][0];
        }
        else if(r<=3)
        //if(r<1+P13+P14)
        {
            ///cout<<"r<=3"<<endl;
            
            array_1D[x] = state[3][2];
            array_1D[global_random_neighbor_num] = state[3][4];
            
            My += -(state[0][2]+state[0][4])+(array_1D[x] + array_1D[global_random_neighbor_num]);
            old_E += -state[0][0]+state[3][0];
        }
        
    }
    
    ///cout<<"state: "<<array_1D[x]<<array_1D[global_random_neighbor_num]<<endl<<endl;
    ///system("pause");
    
    return old_E;
}


//переворот двух спинов, версия 2
double Double_Spin_flip(int *array_1D, double old_E, int x, float T)
{
    int num_of_random_spins = 2;//количество случайных спинов
    int states_amount = 1<<num_of_random_spins;//количество конфигураций
    int local_random_neighbor = rand()%4;//локальный номер случайного соседа
    ///cout<<endl<<"local_random_neighbor = "<<local_random_neighbor<<endl<<endl;
    int global_random_neighbor_num;//глобальный номер случайного соседа
    
    int state[states_amount][5]; //массив состояний выбранных спинов, 
                                 //каждое состояние хранит: 
                                 //[0] энергия выбранных спинов, с учетом соседей
                                 //[1] номер случайного спина, 
                                 //[2] значение случайного спина, 
                                 //[3] глобальный номер случайного соседа, 
                                 //[4] значение случайного соседа
    
    double E_spins = 0; //энергия выбранных спинов, с учетом соседей
    
    int configuration[num_of_random_spins]={1};//варианты конфигураций выбранных спинов
    int min[2];//в массиве хранится значение минимальной энергии и номер конфигурации
    
    double Z = 0;//статсумма
    double Z_state[states_amount];//массив статсум по конфигурациям
    double P[states_amount];//массив вероятностей энергий
    double r;//число от 0 до 1
    
    int k[4],l[4]; //координаты 4x соседей спина
    int i, j;
    
    //проходим по всем состояниям выбранных спинов
    for(int state_num=0; state_num<states_amount; state_num++)
    {
        E_spins = 0;
        j=x%N;
        i=(x-j)/N;
        
        //проходим по выбранным спинам
        for(int random_spin=0; random_spin<num_of_random_spins; random_spin++)
        {
            Neighbors_searching(k,l,i,j);
            
            for(int d=0; d<4; d++)//проходим по соседям, у каждого спина 4 соседа
            {
                //**
                cout<<"spin= "<<i*N+j<<"("<<array_1D[i*N+j]<<
                      "), neigh= "<<k[d]*N+l[d]<<"("<<array_1D[k[d]*N+l[d]]<<")"<<
                      endl;
                //*/
                
                if(random_spin!=0 || d!=local_random_neighbor)
                {
                    E_spins += -J*(array_1D[i*N+j] * array_1D[k[d]*N+l[d]]);//получаем энергию спинов, с учетом соседей
                    //cout<<"E = "<<E_sys<<endl;
                    //system("pause");
                }
            }
            
            if(random_spin<num_of_random_spins-1){
                i = k[local_random_neighbor];
                j = l[local_random_neighbor];
                global_random_neighbor_num = i*N+j;
                cout<<"global_random_neighbor_num = "<<global_random_neighbor_num<<endl;
            }
        }
        
        state[state_num][0] = E_spins; //энергия выбранных спинов, с учетом соседей
        state[state_num][1] = x; //номер случайного спина
        state[state_num][2] = array_1D[x]; //значение случайного спина
        state[state_num][3] = global_random_neighbor_num; //глобальный номер случайного соседа
        state[state_num][4] = array_1D[global_random_neighbor_num]; //значение случайного соседа
        
        configuration[0] = array_1D[x];//исправить на несколько спинов!!!!!!!!!!!!!!!!!
        configuration[1] = array_1D[global_random_neighbor_num];//исправить на несколько спинов!!!!!!!!!!!!
        
        /**
        cout<<"configuration: <" << configuration[0] << configuration[1]<<">" << endl;
        
        cout << "E = "<<state[state_num][0]<<
                ", spin = "<<state[state_num][1]<<
                ", array_1D[x] = "<<state[state_num][2]<<
                ", global_random_neighbor_num = "<<state[state_num][3]<<
                ", array_1D[global_random_neighbor_num] = "<<state[state_num][4]<<
                endl<<endl;
        //*/
        
        //переход на следующую конфигурацию
        for(int t = 0; t<num_of_random_spins; t++)
        {
            if(configuration[t] == 1)
            {
                configuration[t] = -1;
                break;
            }
            else configuration[t] = 1;
        }
        
        array_1D[x] = configuration[0];//исправить на несколько спинов!!!!!!!!!!!!!!!
        array_1D[global_random_neighbor_num] = configuration[1];//исправить на несколько спинов!!!!!!!!!!!!!!
        
    }
    
    //поиск конфигурации с минимальной энергией
    min[0]=state[0][0];
    min[1]=1;
    
    for(int state_num=1; state_num<states_amount; state_num++)
    {
        if(min[0]>state[state_num][0])
        {
            min[0] = state[state_num][0];//E
            min[1] = state_num;//номер состояния
        }
    }
    
    ///cout<<"Emin = "<<min[0]<<", conf = "<<min[1]<<endl;
    ///cout<<"x = "<<x<<", neig = "<<global_random_neighbor_num<<endl;
    
    //считаем статсумму
    Z = 0;
    for(int state_num=0; state_num<states_amount; state_num++){
        Z_state[state_num] = exp((double)(-(state[state_num][0]-min[0])/T));
        Z += Z_state[state_num];
    }
    
    ///cout << "Z = " << Z << endl;
    
    //считаем вероятности
    for(int state_num=0; state_num<states_amount; state_num++)
    {
        P[state_num]= Z_state[state_num]/Z;
        ///cout << "P"<<state_num<<"= " << P[state_num] << endl;
        ///E_aver+= P[state_num]*(old_E-state[0][0]+state[state_num][0]);
    }
    
    
    /**
    cout<<"E1_tot= "<<old_E<<endl;
    cout<<"E2_tot= "<<(old_E-state[0][0]+state[1][0])<<endl;
    cout<<"E3_tot= "<<(old_E-state[0][0]+state[2][0])<<endl;
    cout<<"E4_tot= "<<(old_E-state[0][0]+state[3][0])<<endl;
    
    cout<<"P1*E1_tot= "<<P1*old_E<<endl;
    cout<<"P2*E2_tot= "<<P2*(old_E-state[0][0]+state[1][0])<<endl;
    cout<<"P3*E3_tot= "<<P3*(old_E-state[0][0]+state[2][0])<<endl;
    cout<<"P4*E4_tot= "<<P4*(old_E-state[0][0]+state[3][0])<<endl;
    
    cout<<"E_aver= "<<E_aver<<endl<<endl;
    
    cout<<"E1= "<<state[0][0]<<endl;
    cout<<"E2= "<<state[1][0]<<endl;
    cout<<"E3= "<<state[2][0]<<endl;
    cout<<"E4= "<<state[3][0]<<endl;
    cout<<"Emin= "<<min[0]<<endl;
    
    cout<<"exp1= "<<exp((double)(-(state[0][0]-min[0])/T))<<endl;
    cout<<"exp2= "<<exp((double)(-(state[1][0]-min[0])/T))<<endl;
    cout<<"exp3= "<<exp((double)(-(state[2][0]-min[0])/T))<<endl;
    cout<<"exp4= "<<exp((double)(-(state[3][0]-min[0])/T))<<endl;
    cout<<"Z= "<<Z<<endl;
    
    system("pause");
    //*/
    
    uniform_real_distribution<double> distribution_real(0,1); //вещественное равномерное распределение
    
    r = distribution_real(generator);
    ///cout<<"r= "<<r<<endl;
    
    //перебираем интервалы вероятностей
    double interval = P[0];
    for(int state_num=0; state_num < states_amount; state_num++)
    {
        if(r<=interval)
        {
            ///cout<<"r<="<<interval<<endl;
            
            array_1D[x] = state[state_num][2];
            array_1D[global_random_neighbor_num] = state[state_num][4];
            
            My += -(state[0][2]+state[0][4])+(array_1D[x] + array_1D[global_random_neighbor_num]);
            old_E += -state[0][0]+state[state_num][0];
            
            break;
        }
        else
        {
            interval+=P[state_num+1];
        }
    }
    
    ///cout<<"state: "<<array_1D[x]<<array_1D[global_random_neighbor_num]<<endl<<endl;
    ///system("pause");
    
    return old_E;
}



//переворот блока спинов
double Block_Spin_flip(int *array_1D, double old_E, float T, int **neighbors_array)
{
    
    /** перебор всех конфигураций границ и блока, сохранение E и M в массив --->>>> */
    
    int n = 3; //длина стороны блока
    
    int num_of_spins_in_block = n*n; //количество спинов в блоке
    ulong states_amount = 1<<num_of_spins_in_block;//количество конфигураций блока
    
    int num_of_spins_in_boundaries = n*4; //количество спинов на границе (длина стороны блока * количество соседей)
    ulong boundaries_amount = 1<<num_of_spins_in_boundaries; //количество вариантов границ
    
    //структура хранит энергию, намагниченность и конфигурации с учетом соседей
    struct EM_struct{
        int E;
        int M;
        vector <ulong> state_array;
        
        //конструктор
        EM_struct(int E, int M){EM_struct::E = E; EM_struct::M = M;}
    };
    
    vector < vector <EM_struct> > boundaries_and_states_array (boundaries_amount);//массив границ и состояний
    
    int block_array[states_amount][num_of_spins_in_block];
    int boundaries_array[num_of_spins_in_boundaries];
    
    int M, M_1, M_2, count=0;
    int E_block = 0; //энергия блока без учета границ
    int E_block_tot = 0; //энергия блока с учетом границ
    
    //перебор конфигураций блока побитовым сдвигом
    for(ulong state_num=0; state_num<states_amount; ++state_num)
    {
        ///cout << "State: " << state_num << " ( ";
        count=0;
        
        for(int spin_num_in_block=(num_of_spins_in_block-1); spin_num_in_block>=0; --spin_num_in_block)
        {
            if(state_num&(1<<spin_num_in_block))
                block_array[state_num][count]=1;
            else
                block_array[state_num][count]=-1;
            
            ///cout << block_array[state_num][count] << " "; //вывод конфигурации на экран
            count++;
        }
        ///cout<<")"<<endl;
    }
    ///system("pause");
    
    ulong num_of_unique_EM = 0;
    //перебор границ побитовым сдвигом
    for(ulong boundary_num=0; boundary_num<boundaries_amount; ++boundary_num)
    {
        ///cout << "\nBoundary: " << boundary_num << "( ";
        M_1=0;
        count=0;
        
        for(int spin_num_in_boundary=(num_of_spins_in_boundaries-1); spin_num_in_boundary>=0; --spin_num_in_boundary)
        {   
            if(boundary_num&(1<<spin_num_in_boundary)) 
                boundaries_array[count]=1;
            else 
                boundaries_array[count]=-1;
            
            ///cout << boundaries_array [count] << " "; //вывод границы на экран
            
            M_1+=boundaries_array[count];
            count++;
        }
        
        ///cout<<")"<<endl;
        
        //перебор состояний в рамках границ и подсчет энергии и намагниченности
        for(ulong state_num=0; state_num<states_amount; ++state_num)
        {
            ///cout << "State: " << state_num << " ( ";
            
            //подсчет энергии и намагниченности в блоке
            M_2=0;
            E_block = 0;
            for(int spin_num_in_block=0; spin_num_in_block<num_of_spins_in_block; ++spin_num_in_block)
            {
                ///cout << block_array[state_num][spin_num_in_block] << " "; //вывод конфигурации на экран
                
                //энергия блока без учета границ
                if(n==2)
                {
                    if(spin_num_in_block < num_of_spins_in_block-1)
                        E_block += -J*block_array[state_num][spin_num_in_block]*block_array[state_num][spin_num_in_block+1];
                    else if(spin_num_in_block == num_of_spins_in_block-1)
                        E_block += -J*block_array[state_num][spin_num_in_block]*block_array[state_num][0];
                } else
                    if(n==3)
                    {
                        if(spin_num_in_block < num_of_spins_in_block-2)
                            E_block += -J*block_array[state_num][spin_num_in_block]*block_array[state_num][spin_num_in_block+1];
                        else if(spin_num_in_block == num_of_spins_in_block-2)
                            E_block += -J*block_array[state_num][spin_num_in_block]*block_array[state_num][0];
                        else if(spin_num_in_block == num_of_spins_in_block-1)
                        {
                            E_block += -J*(
                                        block_array[state_num][spin_num_in_block]*block_array[state_num][1]
                                      + block_array[state_num][spin_num_in_block]*block_array[state_num][3]
                                      + block_array[state_num][spin_num_in_block]*block_array[state_num][5]
                                      + block_array[state_num][spin_num_in_block]*block_array[state_num][7]
                                          );
                        }
                    }
                
                //намагниченность блока
                M_2 += block_array[state_num][spin_num_in_block];
            }
            
            ///cout<<")"<<endl;
            ///system("pause");
            
            //энергия блока (4 спина) с учетом границ
            if(n==2){
                E_block_tot = E_block 
                        -J*(
                            boundaries_array[0]*block_array[state_num][0]
                          + boundaries_array[1]*block_array[state_num][1]
                          + boundaries_array[2]*block_array[state_num][1]
                          + boundaries_array[3]*block_array[state_num][2]
                          + boundaries_array[4]*block_array[state_num][2]
                          + boundaries_array[5]*block_array[state_num][3]
                          + boundaries_array[6]*block_array[state_num][3]
                          + boundaries_array[7]*block_array[state_num][0]
                           );
            } else
            
            //энергия блока (9 спинов) с учетом границ
            if(n==3){
                E_block_tot = E_block 
                        -J*(
                            boundaries_array[0]*block_array[state_num][0]
                          + boundaries_array[1]*block_array[state_num][1]
                          + boundaries_array[2]*block_array[state_num][2]
                          + boundaries_array[3]*block_array[state_num][2]
                          + boundaries_array[4]*block_array[state_num][3]
                          + boundaries_array[5]*block_array[state_num][4]
                          + boundaries_array[6]*block_array[state_num][4]
                          + boundaries_array[7]*block_array[state_num][5]
                          + boundaries_array[8]*block_array[state_num][6]
                          + boundaries_array[9]*block_array[state_num][6]
                          + boundaries_array[10]*block_array[state_num][7]
                          + boundaries_array[11]*block_array[state_num][0]
                           );
            } else cout<<"ERROR n!!!!!!!!!!!\n";
            
            M = M_1+M_2;
            //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
            //M и E блока посчитано, далее заполняем массив с учетом вырождений
            
            
            if(state_num == 0)
            {
                boundaries_and_states_array[boundary_num].push_back(EM_struct(E_block_tot,M));
                boundaries_and_states_array[boundary_num][0].state_array.push_back(0);
                
                //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
                //cout<<"E= "<<boundaries_and_states_array[boundary_num][boundaries_and_states_array[boundary_num].size()-1].E<<endl;
                //system("pause");
            }
            else
            for(ulong ii_4=0; ii_4<boundaries_and_states_array[boundary_num].size(); ++ii_4)
            {
                //cout<<"qq1"<<endl;
                //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
                
                if(boundaries_and_states_array[boundary_num][ii_4].E == E_block_tot && boundaries_and_states_array[boundary_num][ii_4].M == M)
                {
                    //cout<<"qq2"<<endl;
                    boundaries_and_states_array[boundary_num][ii_4].state_array.push_back(state_num);
                    //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
                    //cout << "E1= " <<boundaries_and_states_array[boundary_num][boundaries_and_states_array[boundary_num].size()-1].E<<endl;
                    //system("pause");
                    break;
                }
                else if(ii_4 == boundaries_and_states_array[boundary_num].size()-1)
                {
                    //cout<<"qq3"<<endl;
                    boundaries_and_states_array[boundary_num].push_back(EM_struct(E_block_tot,M));
                    boundaries_and_states_array[boundary_num][ii_4+1].state_array.push_back(state_num);
                    ii_4++;
                    //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
                    //cout<<"E2= "<<boundaries_and_states_array[boundary_num][boundaries_and_states_array[boundary_num].size()-1].E<<endl;
                    //system("pause");
                }
                //system("pause");
                ///cout << "M = " << states_array[ii][ii_4].M << endl;
                ///cout << "E = " << states_array[ii][ii_4].E << endl;
            }
            
            ///system("pause");
            //здесь заканчивается перебор состояний в рамках границ
        }
        
        ///cout << "Number of unique E&M = " << boundaries_and_states_array[boundary_num].size() << endl;
        num_of_unique_EM+=boundaries_and_states_array[boundary_num].size();
        
        /*
        //вывод E и M с учетом кратности вырождения
        for(int ii=0; ii<boundaries_and_states_array[boundary_num].size(); ++ii){
            cout<<"________\n";
            cout<<"E = "<<boundaries_and_states_array[boundary_num][ii].E;
            cout<<", M = "<<boundaries_and_states_array[boundary_num][ii].M<<endl;
            cout<<"States: { ";
            for(int jj=0; jj<boundaries_and_states_array[boundary_num][ii].state_array.size(); ++jj){
                cout<<boundaries_and_states_array[boundary_num][ii].state_array[jj]<<" ";
            }
            cout<<"}";
            cout<<endl;
        }
        //*/
        
        ///system("pause");
        //здесь заканчивается перебор границ
    }
    
    ///cout << "Number of boundaries = " << boundaries_and_states_array.size() << endl;
    //cout << "Number of states with unique E&M = " << boundaries_and_states_array[0][0].state_array.size() << endl;
    //cout << "M = " << boundaries_and_states_array[0][0].M << endl;
    //cout << "E = " << boundaries_and_states_array[0][0].E << endl;
    ///cout << "The number of all pairs EM = "<<boundaries_amount*states_amount<<endl;
    ///cout << "The number of unique pairs EM = "<<num_of_unique_EM<<endl;
    
    ///system("pause");
    
    
    /** 
     * на выходе получаем массив boundaries_and_states_array[][].{E,M,state_array[]}, 
     * который хранит все конфигурации границ, структуру уникальных значений E,M, и соответствующие им конфигурации блоков
     */
    /** <<<<----- здесь заканчивается перебор всех конфигураций границ и блока ----- */
    
    
    
    
    
    /** --------------- подставляем блок в систему --------------->>>> */
    
    vector <int> block_temp (n*n,0);
    vector <int> boundary_temp (n*4,0);
    
    ///добавить температуру
    ///float T = 0.0001;
    /// double C, hi;
    ///ofstream C_data("C.txt");
    ///ofstream hi_data("hi.txt");
    ///ofstream E_data("E.txt");
    
    
    //проход блока по системе
    for(ulong MCS = 0; MCS<10000*N*N; ++MCS)
    {
        ulong spin_num = rand()%(N*N);
        ///for(int neighbor_num=0; neighbor_num<4; ++neighbor_num) cout<<neighbors[spin_num][neighbor_num]<<",";
        
        if(n==2)
        {
            ///cout<<spin_num<<": ";
            // определение блока в системе
            ///cout<<"block: ";
            int next_spin=spin_num; //0
            block_temp[0]=array_1D[next_spin]; //0
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][0]; //1
            block_temp[1]=array_1D[next_spin]; //1
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][1]; //2
            block_temp[2]=array_1D[next_spin]; //2
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][2]; //3
            block_temp[3]=array_1D[next_spin]; //3
            ///cout<<next_spin;
            
            ///cout<<" ( ";
            ///for(ulong ii=0; ii<block_temp.size(); ++ii) cout<<block_temp[ii]<<" ";
            ///cout<<")\n";
            
            
            // определение границы в системе
            ///cout<<"boundary: ";
            next_spin=neighbors_array[spin_num][3]; //0
            boundary_temp[0] = array_1D[next_spin]; //0
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][0]; //1
            boundary_temp[1] = array_1D[next_spin]; //1
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[neighbors_array[next_spin][0]][1]; //2
            boundary_temp[2] = array_1D[next_spin]; //2
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][1]; //3
            boundary_temp[3] = array_1D[next_spin]; //3
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[neighbors_array[next_spin][1]][2]; //4
            boundary_temp[4] = array_1D[next_spin]; //4
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][2]; //5
            boundary_temp[5] = array_1D[next_spin]; //5
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[neighbors_array[next_spin][2]][3]; //6
            boundary_temp[6] = array_1D[next_spin]; //6
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][3]; //7
            boundary_temp[7] = array_1D[next_spin]; //7
            ///cout<<next_spin<<endl;
        } else 
            if(n==3)
            {
                ///cout<<spin_num<<": ";
                // определение блока в системе
                ///cout<<"block: ";
                int next_spin=spin_num;//0
                block_temp[0]=array_1D[next_spin];//0
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0];//1
                block_temp[1]=array_1D[next_spin];//1
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0];//2
                block_temp[2]=array_1D[next_spin];//2
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][1];//3
                block_temp[3]=array_1D[next_spin];//3
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][1];//4
                block_temp[4]=array_1D[next_spin];//4
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][2];//5
                block_temp[5]=array_1D[next_spin];//5
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][2];//6
                block_temp[6]=array_1D[next_spin];//6
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][3];//7
                block_temp[7]=array_1D[next_spin];//7
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0];//8
                block_temp[8]=array_1D[next_spin];//8
                ///cout<<next_spin;
                
                ///cout<<" ( ";
                
                ///for(ulong ii=0; ii<block_temp.size(); ++ii) cout<<block_temp[ii]<<" ";
                
                ///cout<<")\n";
                
                // определение границы в системе
                ///cout<<"boundary: ";
                next_spin=neighbors_array[spin_num][3]; //0
                boundary_temp[0] = array_1D[next_spin]; //0
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0]; //1
                boundary_temp[1] = array_1D[next_spin]; //1
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0]; //2
                boundary_temp[2] = array_1D[next_spin]; //2
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[neighbors_array[next_spin][0]][1]; //3
                boundary_temp[3] = array_1D[next_spin]; //3
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][1]; //4
                boundary_temp[4] = array_1D[next_spin]; //4
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][1]; //5
                boundary_temp[5] = array_1D[next_spin]; //5
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[neighbors_array[next_spin][1]][2]; //6
                boundary_temp[6] = array_1D[next_spin]; //6
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][2]; //7
                boundary_temp[7] = array_1D[next_spin]; //7
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][2]; //8
                boundary_temp[8] = array_1D[next_spin]; //8
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[neighbors_array[next_spin][2]][3]; //9
                boundary_temp[9] = array_1D[next_spin]; //9
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][3]; //10
                boundary_temp[10] = array_1D[next_spin]; //10
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][3]; //11
                boundary_temp[11] = array_1D[next_spin]; //11
                ///cout<<next_spin<<endl;
            }
            else cout<<"ERROR n2!!!!!!!!";
        ///cout<<endl;
        
        ulong boundary_dec=0; //конфигурация границы в десятичной системе
        int bit;
        for (ulong ii=0; ii<boundary_temp.size(); ++ii)
        {
            if (boundary_temp[ii]==-1) bit = 0; 
            else bit = 1;
            boundary_dec = (boundary_dec | bit);
            if (ii<boundary_temp.size()-1)
                boundary_dec = boundary_dec<<1;
        }
        ///cout<<"b_n = "<<dec<<endl;
        
        /*
        cout<<"boundary: ";
        for(ulong ii=0; ii<boundary_temp.size(); ++ii){
            cout<<boundary_temp[ii]<<" ";
        }
        cout<<endl<<endl;
        //*/
        
        
        ulong block_dec=0; //конфигурация блока в десятичной системе
        for (ulong ii=0; ii<block_temp.size(); ++ii)
        {
            if (block_temp[ii]==-1) bit = 0; 
            else bit = 1;
            block_dec = (block_dec | bit);
            if (ii<block_temp.size()-1)
                block_dec = block_dec<<1;
        }
        
        int E0=0, M0=0;
        int Emin = boundaries_and_states_array[boundary_dec][0].E;
        ulong num_of_unique_EM_in_boundary =boundaries_and_states_array[boundary_dec].size();
        
        for (ulong ii=0; ii<num_of_unique_EM_in_boundary; ++ii)
        {
            for(ulong jj=0; jj<boundaries_and_states_array[boundary_dec][ii].state_array.size(); ++jj)
            {
                if(boundaries_and_states_array[boundary_dec][ii].state_array[jj] == block_dec)
                {
                    E0 = boundaries_and_states_array[boundary_dec][ii].E;
                    M0 = boundaries_and_states_array[boundary_dec][ii].M;
                    //ii = num_of_unique_EM_in_boundary;
                    break;
                }
            }
            if(Emin > boundaries_and_states_array[boundary_dec][ii].E) 
                Emin = boundaries_and_states_array[boundary_dec][ii].E;
        }
        ///cout<<"E0 = "<<E0<<", M0 = "<<M0<<endl;
        ///cout << "Emin = " << Emin << endl;
        
        old_E-=E0;
        My-=M0;
        
        
        
        double Z = 0;//статсумма
        double Z_i[num_of_unique_EM_in_boundary];//массив статсумм по конфигурациям
        double P[num_of_unique_EM_in_boundary];//массив вероятностей энергий
        double r;//число от 0 до 1
        
        //считаем статсумму
        Z = 0;
        for(ulong set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; set_of_states++)
        {
            //вырождение * exp
            Z_i[set_of_states] = boundaries_and_states_array[boundary_dec][set_of_states].state_array.size()*
                    exp((double)(-(boundaries_and_states_array[boundary_dec][set_of_states].E-Emin)/T));
            Z += Z_i[set_of_states];
            ///cout<<exp((double)(-(boundaries_and_states_array[boundary_dec][set_of_states].E-Emin)/T))<<endl;
        }
        
        ///system("pause");
        
        ///cout << "Z = " << Z << endl;
        
        double E_aver_i = 0;
        //double sum=0;
        //считаем вероятности
        for(ulong set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; set_of_states++)
        {
            P[set_of_states]= Z_i[set_of_states]/Z;
            //sum+=P[set_of_states];
            ///cout << "P"<<set_of_states<<"= " << P[set_of_states] << endl;
            E_aver_i+= P[set_of_states]*(old_E+boundaries_and_states_array[boundary_dec][set_of_states].E);
            ///cout<<"P = "<<P[set_of_states]<<endl;
        }
        ///cout<<"E_aver_i = "<<E_aver_i<<endl;
        //cout<<"P_sum = " << sum << endl;
        E_aver += (E_aver_i - E_aver) / (MCS + 1);
        E_aver2 += (E_aver_i*E_aver_i - E_aver2) / (MCS + 1);
        ///cout<<"E_aver = "<<E_aver<<endl;
        ///system("pause");
        
        uniform_real_distribution<double> distribution_real(0,1); //вещественное равномерное распределение
        
        r = distribution_real(generator);
        ///cout<<"r= "<<r<<endl;
        
        int rand_state=0;//случайная конфигурация с одинаковыми параметрами
        
        //перебираем интервалы вероятностей
        double interval = P[0];
        for(ulong set_of_states=0; set_of_states < num_of_unique_EM_in_boundary; ++set_of_states)
        {
            if(r<=interval)
            {
                ///cout<<"r<="<<interval<<endl;
                
                My += boundaries_and_states_array[boundary_dec][set_of_states].M;
                old_E += boundaries_and_states_array[boundary_dec][set_of_states].E;
                
                //выбираем случайную конфигурацию с одинаковыми параметрами
                rand_state = rand() % boundaries_and_states_array[boundary_dec][set_of_states].state_array.size();
                
                //запоминаем конфигурацию
                rand_state = boundaries_and_states_array[boundary_dec][set_of_states].state_array[rand_state];
                
                //rand_state = 2; //проверка работоспособности
                
                //переводим в двоичную систему
                ///cout<<"rand_state = "<<rand_state<<endl;
                if(n==2)
                {
                    count=0;
                    for(bit=num_of_spins_in_block-1; bit>=0; --bit)
                    {
                        block_temp[count] = (1 & rand_state >> bit);
                        if(block_temp[count]==0) block_temp[count] = -1;
                        ///cout << block_temp[count]<< " ";
                        count++;
                    }
                    ///cout<<endl;
                    
                    // заменяем блок внутри границы
                    int next_spin=spin_num; //0
                    array_1D[next_spin]=block_temp[0]; //0
                    ///cout<<next_spin<<";";
                    next_spin=neighbors_array[next_spin][0]; //1
                    array_1D[next_spin]=block_temp[1]; //1
                    ///cout<<next_spin<<";";
                    next_spin=neighbors_array[next_spin][1]; //2
                    array_1D[next_spin]=block_temp[2]; //2
                    ///cout<<next_spin<<";";
                    next_spin=neighbors_array[next_spin][2]; //3
                    array_1D[next_spin]=block_temp[3]; //3
                    ///cout<<next_spin<<endl;
                } else 
                    if(n==3)
                    {
                        count=0;
                        for(bit=num_of_spins_in_block-1; bit>=0; --bit)
                        {
                            block_temp[count] = (1 & rand_state >> bit);
                            if(block_temp[count]==0) block_temp[count] = -1;
                            ///cout << block_temp[count]<< " ";
                            count++;
                        }
                        ///cout<<endl;
                        
                        // заменяем блок внутри границы
                        int next_spin=spin_num;//0
                        array_1D[next_spin]=block_temp[0];//0
                        ///cout<<next_spin<<";";
                        next_spin=neighbors_array[next_spin][0];//1
                        array_1D[next_spin]=block_temp[1];//1
                        ///cout<<next_spin<<";";
                        next_spin=neighbors_array[next_spin][0];//2
                        array_1D[next_spin]=block_temp[2];//2
                        ///cout<<next_spin<<";";
                        next_spin=neighbors_array[next_spin][1];//3
                        array_1D[next_spin]=block_temp[3];//3
                        ///cout<<next_spin<<";";
                        next_spin=neighbors_array[next_spin][1];//4
                        array_1D[next_spin]=block_temp[4];//4
                        ///cout<<next_spin<<";";
                        next_spin=neighbors_array[next_spin][2];//5
                        array_1D[next_spin]=block_temp[5];//5
                        ///cout<<next_spin<<";";
                        next_spin=neighbors_array[next_spin][2];//6
                        array_1D[next_spin]=block_temp[6];//6
                        ///cout<<next_spin<<";";
                        next_spin=neighbors_array[next_spin][3];//7
                        array_1D[next_spin]=block_temp[7];//7
                        ///cout<<next_spin<<";";
                        next_spin=neighbors_array[next_spin][0];//8
                        array_1D[next_spin]=block_temp[8];//8
                        ///cout<<next_spin;
                        ///cout<<endl;
                    } else cout<<"ERROR n3!!!!!!!!!!";
                break;
            }
            else
            {
                interval+=P[set_of_states+1];
            }
        }
       //здесь заканчивается проход блока по системе 
    }
    /** <<<<--------------------------------------------------- */
    
    
    ///cout << "\n" << (double)clock()/CLOCKS_PER_SEC << " sec.\n";
    ///system("pause");
    
    return old_E;
}
