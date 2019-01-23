/** Реализация блочного алгоритма **/

/** 
 * Список функций:
 * - Init_1D_Array - заполнение одномерного массива
 * - Full_Energy_Square_Ising_1D_array - функция подсчета энергии системы (4 соседа) для одномерного массива
 * - Full_Energy_Square_Ising_recalculation_1D_array - функция пересчета энергии системы (4 соседа) для одномерного массива
 * 
 * - Metropolis_Single_Spin_flip - алгоритм Метрополиса для одного спина
 * - Block_Spin_flip - переворот блока спинов для одномерного массива
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

#define N 100 // MAX 360
#define NumMC 20000*N*N //*20
#define HEATING 0
#define MONTE_CARLO 0
#define J 1 //('1'-ferro, '-1'-antiferro)

#define Tmin 0.0001
#define Tmax 4.01
float Tstep = 0.1;

double C1=0;

void Init_1D_Array(int *array_1D); //заполнение одномерного массива
double Full_Energy_Square_Ising_1D_array(int *array_1D); //функция подсчета энергии системы (4 соседа)
double Full_Energy_Square_Ising_recalculation_1D_array(int *array_1D, double old_E, int x); //функция пересчета энергии системы (4 соседа)

void Neighbors_searching(int *k, int *l, int i, int j); //определение соседей одного спина
void Neighbors_saving(int **neighbors_array); //функция формирование массива соседей каждого спина

void Metropolis_Single_Spin_flip(int *array_1D, double E_prev); //алгоритм Метрополиса для одного спина
void Block_Spin_flip(int *array_1D, double old_E, int **neighbors_array); //переворот блока спинов

using namespace std;

default_random_engine generator; //тип генератора случайных чисел

int My = 0;


int main()
{
    cout<<N*N<<" spins\n";
    time_t s_time = time(NULL);
    generator.seed(1); //стартовое число
    srand(1);
    
    //выделение места в памяти под динамический одномерный массив спинов:
    int *array_1D = new int [N*N];
    
    //выделение места в памяти под динамический двумерный массив соседей:
    int **neighbors_array = new int* [N*N];
    for (int i = 0; i < N*N; ++i)
        neighbors_array[i] = new int [4];
    
    //заполняем массив спинов
    Init_1D_Array(array_1D);
    
    //заполняем массив соседей
    Neighbors_saving(neighbors_array);
    
    //cout << "E_1D = " << Full_Energy_Square_Ising_1D_array(array_1D) << endl;
    //system("pause");
    
    
    //вычисляем начальную энергию
    double E_prev = Full_Energy_Square_Ising_1D_array(array_1D);
    cout << "E_start = " << E_prev << endl << endl;
    
    if(MONTE_CARLO)
    {
        cout<<"*** Starting Monte-Carlo Single Spin Flip ***\n\n";
        Metropolis_Single_Spin_flip(array_1D, E_prev);
    }
    
    //переворот блоков
    if(!MONTE_CARLO)
    {
        cout<<"*** Starting Block Spin Flip ***\n\n";
        Block_Spin_flip(array_1D, E_prev, neighbors_array);
    }
    
    //cout << "\nE_prev = " << E_prev << endl;
    //cout << "Full Energy = " << Full_Energy_Square_Ising_1D_array(array_1D) << endl;
    
    /** Очистка памяти */
    
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


//алгоритм Метрополиса для одного спина
void Metropolis_Single_Spin_flip(int *array_1D, double E_prev)
{
    double E_new; //хранит новую энергию
    
    int x; //номер случайного спина
    double r, ver; //вероятности для монте-карло
    
    uniform_int_distribution<int> distribution_int(0,N*N-1); //целочисленное равномерное распределение
    uniform_real_distribution<double> distribution_real(0,1); //вещественное равномерное распределение
    
    ofstream C_data("C.txt"); //файл теплоемкости
    ofstream hi_data("hi.txt"); //файл магнитной восприимчивости
    ofstream E_data("E.txt"); //файл энергии
    
    double E_aver = 0, E2_aver = 0; //средняя энергия и ее квадрат
    double Mpr_aver = 0, Mpr2_aver = 0; //средняя проекция намагниченности и ее квадрат

    for(float T = Tmin; T<Tmax; T+=Tstep)
    {
        cout << "T = " << T << endl << "----------"<< endl;
        
        E_aver = 0;
        E2_aver = 0;
        Mpr_aver = 0;
        Mpr2_aver = 0;
        
        /** ----- односпиновое Монте-Карло ----->>> **/
        //прогрев
        if(HEATING)
        {
            cout << "Heating...\n";
            
            for(ulong MCS = 0; MCS<NumMC; MCS++)
            {
                //Metropolis
                //--start
                
                x = distribution_int(generator);
                
                array_1D[x] *= (-1);
                E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_prev, x);
                
                if(E_new >= E_prev){
                    r = distribution_real(generator);
                    ver = exp(-((E_new-E_prev)/T));
                    if (r>ver){
                        array_1D[x] *= (-1);
                        //E_new = E_prev;
                        E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_new, x);
                    }
                }
                E_prev = E_new;
                
                //--finish
            }
        }
        
        //Монте-Карло
        cout << "Monte Carlo...\n";
        for(ulong MCS = 0; MCS<NumMC; MCS++)
        {
            //Metropolis
            //--start
            x = distribution_int(generator);
            
            ///cout << "E_new = "<<E_new<<endl;
            ///cout << "My = "<<My<<endl;
            ///cout<<"E_new(real) = "<<Full_Energy_Square_Ising_1D_array(array_1D)<<endl;
            ///cout << "My(real) = "<<My<<endl<<endl;
            
            //system("pause");
            
            array_1D[x] *= (-1);
            E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_prev, x);
            
            if(E_new >= E_prev)
            {
                r = distribution_real(generator);
                ver = exp(-(E_new-E_prev)/T);
                if (r>ver)
                {
                    array_1D[x] *= (-1);
                    //E_new = E_prev;
                    E_new = Full_Energy_Square_Ising_recalculation_1D_array(array_1D, E_new, x);
                }
            }
            E_prev = E_new;
            //--finish
            
            //средние
            E_aver += (E_new - E_aver) / (MCS + 1);
            E2_aver += (E_new*E_new - E2_aver) / (MCS + 1);
            
            Mpr_aver += (abs(My) - Mpr_aver) / (double)(MCS + 1);
            Mpr2_aver += ((My*My) - Mpr2_aver) / (double)(MCS + 1);
        }
        /** <<<----- односпиновое Монте-Карло ----- **/
        
        //cout << "Mpr_aver = " << Mpr_aver << endl << endl;
        
        cout << "E_aver = " << E_aver << endl;
        
        double C = ((E2_aver - E_aver*E_aver)/(T*T))/(N*N);
        cout << "C = " << C << endl;
        C_data << T << "\t" << C << endl;
        
        double hi = ((Mpr2_aver - Mpr_aver*Mpr_aver)/T) / (double)(N*N);
        cout << "hi = " << hi << endl << endl;
        hi_data << T << "\t" << hi << endl;
        
        E_data << T << "\t" << E_aver << endl;
        
        ///system("pause");
    }
    
    C_data.close();
    hi_data.close();
    E_data.close();
}


//переворот блока спинов
void Block_Spin_flip(int *array_1D, double old_E, int **neighbors_array)
{
    
    /** перебор всех конфигураций границ и блока, сохранение E и M в массив --->>>> */
    
    ulong n = 3; //длина стороны блока
    cout<<"Block size: "<<n<<"x"<<n<<endl<<endl;
    
    int num_of_spins_in_block = n*n; //количество спинов в блоке
    ulong states_amount = 1<<num_of_spins_in_block;//количество конфигураций блока
    
    int num_of_spins_in_boundaries = n*4; //количество спинов на границе (длина стороны блока * количество соседей)
    ulong boundaries_amount = 1<<num_of_spins_in_boundaries; //количество вариантов границ
    
    int Emin_array[boundaries_amount]; //минимальные значанияэнергий для каждой границы
    
    //структура хранит энергию, намагниченность и конфигурации с учетом соседей
    struct EM_state_struct{
        int E;
        int M;
        vector <ulong> state_array;
        
        //конструктор
        EM_state_struct(int E, int M){EM_state_struct::E = E; EM_state_struct::M = M;}
    };
    
    //структура хранит энергию, намагниченность и конфигурации с учетом соседей
    struct EM_struct{
        int E;
        int M;
    };
    
    vector < vector <EM_state_struct> > boundaries_EM_blocks (boundaries_amount);//массив границ, энергий, намагниченности и соответствующих состояний
    vector < vector <EM_struct> > boundaries_blocks_EM (boundaries_amount, vector <EM_struct> (states_amount));//массив границ, состояний и их энергий и намагниченности
    
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
                //ищем минимальную энергию для границы
                Emin_array[boundary_num]=E_block_tot;
                
                boundaries_EM_blocks[boundary_num].push_back(EM_state_struct(E_block_tot,M));
                boundaries_EM_blocks[boundary_num][0].state_array.push_back(0);
                
                //заполняем массив границ и соответствующие блоки
                boundaries_blocks_EM[boundary_num][state_num].E = E_block_tot;
                boundaries_blocks_EM[boundary_num][state_num].M = M;
                
                //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
                //cout<<"E= "<<boundaries_and_states_array[boundary_num][boundaries_and_states_array[boundary_num].size()-1].E<<endl;
                //system("pause");
            }
            else{
                //ищем минимальную энергию для границы
                if(Emin_array[boundary_num] > E_block_tot)
                    Emin_array[boundary_num] = E_block_tot;
                
                //заполняем массив границ и соответствующие блоки
                boundaries_blocks_EM[boundary_num][state_num].E = E_block_tot;
                boundaries_blocks_EM[boundary_num][state_num].M = M;
                
                for(ulong ii_4=0; ii_4<boundaries_EM_blocks[boundary_num].size(); ++ii_4)
                {
                    //cout<<"qq1"<<endl;
                    //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
                    
                    if(boundaries_EM_blocks[boundary_num][ii_4].E == E_block_tot && boundaries_EM_blocks[boundary_num][ii_4].M == M)
                    {
                        //cout<<"qq2"<<endl;
                        boundaries_EM_blocks[boundary_num][ii_4].state_array.push_back(state_num);
                        //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
                        //cout << "E1= " <<boundaries_and_states_array[boundary_num][boundaries_and_states_array[boundary_num].size()-1].E<<endl;
                        //system("pause");
                        break;
                    }
                    else if(ii_4 == boundaries_EM_blocks[boundary_num].size()-1)
                    {
                        //cout<<"qq3"<<endl;
                        boundaries_EM_blocks[boundary_num].push_back(EM_state_struct(E_block_tot,M));
                        boundaries_EM_blocks[boundary_num][ii_4+1].state_array.push_back(state_num);
                        ii_4++;
                        //cout<<"E_block_tot = "<<E_block_tot<<", M= "<<M<<endl;
                        //cout<<"E2= "<<boundaries_and_states_array[boundary_num][boundaries_and_states_array[boundary_num].size()-1].E<<endl;
                        //system("pause");
                    }
                    //system("pause");
                    ///cout << "M = " << states_array[ii][ii_4].M << endl;
                    ///cout << "E = " << states_array[ii][ii_4].E << endl;
                }
            }
            
            ///system("pause");
            //здесь заканчивается перебор состояний в рамках границ
        }
        
        ///cout << "Number of unique E&M = " << boundaries_and_states_array[boundary_num].size() << endl;
        num_of_unique_EM+=boundaries_EM_blocks[boundary_num].size();
        
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
    
    ulong block_coordinates_for_spins[N*N][num_of_spins_in_block]; //[номер спина][координаты блока]
    ulong boundary_coordinates_for_spins[N*N][num_of_spins_in_boundaries]; //[номер спина][координаты границ]
    
    for(int spin_num=0; spin_num<N*N; ++spin_num)
    {
        if(n==2)
        {
            ///cout<<spin_num<<": ";
            // определение координат блока в системе
            ///cout<<"block: ";
            int next_spin=spin_num; //0
            block_coordinates_for_spins[spin_num][0] = next_spin; //0
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][0]; //1
            block_coordinates_for_spins[spin_num][1] = next_spin; //1
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][1]; //2
            block_coordinates_for_spins[spin_num][2] = next_spin; //2
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][2]; //3
            block_coordinates_for_spins[spin_num][3] = next_spin; //3
            ///cout<<next_spin;
            
            ///cout<<" ( ";
            ///for(ulong ii=0; ii<block_temp.size(); ++ii) cout<<block_temp[ii]<<" ";
            ///cout<<")\n";
            
            
            // определение координат границы в системе
            ///cout<<"\nboundary: ";
            next_spin=neighbors_array[spin_num][3]; //0
            boundary_coordinates_for_spins[spin_num][0] = next_spin; //0
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][0]; //1
            boundary_coordinates_for_spins[spin_num][1] = next_spin; //1
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[neighbors_array[next_spin][0]][1]; //2
            boundary_coordinates_for_spins[spin_num][2] = next_spin; //2
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][1]; //3
            boundary_coordinates_for_spins[spin_num][3] = next_spin; //3
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[neighbors_array[next_spin][1]][2]; //4
            boundary_coordinates_for_spins[spin_num][4] = next_spin; //4
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][2]; //5
            boundary_coordinates_for_spins[spin_num][5] = next_spin; //5
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[neighbors_array[next_spin][2]][3]; //6
            boundary_coordinates_for_spins[spin_num][6] = next_spin; //6
            ///cout<<next_spin<<";";
            next_spin=neighbors_array[next_spin][3]; //7
            boundary_coordinates_for_spins[spin_num][7] = next_spin; //7
            ///cout<<next_spin<<endl;
        } else 
            if(n==3)
            {
                ///определить блок для каждого спина в начале программы
                ///cout<<spin_num<<": ";
                // определение координат блока в системе
                ///cout<<"block: ";
                int next_spin=spin_num;//0
                block_coordinates_for_spins[spin_num][0] = next_spin;//0
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0];//1
                block_coordinates_for_spins[spin_num][1] = next_spin;//1
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0];//2
                block_coordinates_for_spins[spin_num][2] = next_spin;//2
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][1];//3
                block_coordinates_for_spins[spin_num][3] = next_spin;//3
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][1];//4
                block_coordinates_for_spins[spin_num][4] = next_spin;//4
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][2];//5
                block_coordinates_for_spins[spin_num][5] = next_spin;//5
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][2];//6
                block_coordinates_for_spins[spin_num][6] = next_spin;//6
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][3];//7
                block_coordinates_for_spins[spin_num][7] = next_spin;//7
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0];//8
                block_coordinates_for_spins[spin_num][8] = next_spin;//8
                ///cout<<next_spin;
                
                ///cout<<" ( ";
                ///for(ulong ii=0; ii<block_temp.size(); ++ii) cout<<block_temp[ii]<<" ";
                ///cout<<")\n";
                
                // определение координат границы в системе
                ///cout<<"\nboundary: ";
                next_spin=neighbors_array[spin_num][3]; //0
                boundary_coordinates_for_spins[spin_num][0] = next_spin; //0
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0]; //1
                boundary_coordinates_for_spins[spin_num][1] = next_spin; //1
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][0]; //2
                boundary_coordinates_for_spins[spin_num][2] = next_spin; //2
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[neighbors_array[next_spin][0]][1]; //3
                boundary_coordinates_for_spins[spin_num][3] = next_spin; //3
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][1]; //4
                boundary_coordinates_for_spins[spin_num][4] = next_spin; //4
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][1]; //5
                boundary_coordinates_for_spins[spin_num][5] = next_spin; //5
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[neighbors_array[next_spin][1]][2]; //6
                boundary_coordinates_for_spins[spin_num][6] = next_spin; //6
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][2]; //7
                boundary_coordinates_for_spins[spin_num][7] = next_spin; //7
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][2]; //8
                boundary_coordinates_for_spins[spin_num][8] = next_spin; //8
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[neighbors_array[next_spin][2]][3]; //9
                boundary_coordinates_for_spins[spin_num][9] = next_spin; //9
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][3]; //10
                boundary_coordinates_for_spins[spin_num][10] = next_spin; //10
                ///cout<<next_spin<<";";
                next_spin=neighbors_array[next_spin][3]; //11
                boundary_coordinates_for_spins[spin_num][11] = next_spin; //11
                ///cout<<next_spin<<endl;
            }
            else 
                cout<<"ERROR n2!!!!!!!!";
        
        ///system("pause");
    }
    
    ofstream C_data("C.txt"); //файл теплоемкости
    ofstream hi_data("hi.txt"); //файл магнитной восприимчивости
    ofstream E_data("E.txt"); //файл энергии
    
    double E_aver = 0, E2_aver = 0; //средняя энергия и ее квадрат
    double Mpr_aver = 0, Mpr2_aver = 0; //средняя проекция намагниченности и ее квадрат
    
    ulong set_of_states, total_set_of_states=0;
    bool flag=0;
    
    uniform_real_distribution<double> distribution_real(0,1); //вещественное равномерное распределение
    
    int rand_state=0; //случайная конфигурация с одинаковыми параметрами EM
    
    double Z = 0; //статсумма
    double r; //число от 0 до 1
    double E_aver_i = 0;
    double E2_aver_i = 0;
    double M_aver_i = 0;
    double M2_aver_i = 0;
    
    //интервалы вероятностей
    double interval=0;
    
    ulong spin_num; //случайный спин
    
    //проход блока по системе
    for(float T = Tmin; T<Tmax; T+=Tstep) //цикл по температуре
    {
        if(T>2.1 && T<2.5) Tstep=0.01;
        else Tstep=0.1;
        cout << "T = " << T << endl << "----------"<<endl;
        
        E_aver = 0;
        E2_aver = 0;
        Mpr_aver = 0;
        Mpr2_aver = 0;
        
        for(ulong MCS = 0; MCS<NumMC; ++MCS)
        {
            spin_num = rand()%(N*N);
            /**
            if(MCS<N*N) spin_num = MCS;
            else spin_num = MCS-(int)(MCS/(N*N-1))*(N*N-1);
            //*/
            
            ///for(int neighbor_num=0; neighbor_num<4; ++neighbor_num) cout<<neighbors[spin_num][neighbor_num]<<",";
            
            ///cout<<endl;
            
            
            ///ПЕРЕДЕЛАЛ!!! ПРОВЕРИТЬ!!!
            /// ulong array_to_dec(int *array_1D) //Оформить как функцию
            ulong boundary_dec=0; //конфигурация границы в десятичной системе
            int bit;
            for (int ii=0; ii<num_of_spins_in_boundaries; ++ii)
            {
                if (array_1D[boundary_coordinates_for_spins[spin_num][ii]]==-1) 
                    bit = 0;
                else 
                    bit = 1;
                
                boundary_dec = (boundary_dec | bit);
                
                if (ii<num_of_spins_in_boundaries-1) 
                    boundary_dec = boundary_dec<<1;
            }
            ///cout<<"boundary_dec = "<<boundary_dec<<endl;
            
            ulong block_dec=0; //конфигурация блока в десятичной системе
            for (int ii=0; ii<num_of_spins_in_block; ++ii)
            {
                if (array_1D[block_coordinates_for_spins[spin_num][ii]]==-1) 
                    bit = 0;
                else 
                    bit = 1;
               
                block_dec = (block_dec | bit);
                
                if (ii<num_of_spins_in_block-1)
                    block_dec = block_dec<<1;
            }
            
            ulong num_of_unique_EM_in_boundary = boundaries_EM_blocks[boundary_dec].size(); //g
            double P[num_of_unique_EM_in_boundary]; //массив вероятностей энергий
            
            ///cout<<"E0 = "<<boundaries_blocks_EM[boundary_dec][block_dec].E<<", M0 = "<<boundaries_blocks_EM[boundary_dec][block_dec].M<<endl;
            ///cout << "Emin = " << Emin_array[boundary_dec] << endl;
            
            old_E-= boundaries_blocks_EM[boundary_dec][block_dec].E;
            My-= boundaries_blocks_EM[boundary_dec][block_dec].M;
            
            Z = 0; //статсумма
            E_aver_i = 0;
            E2_aver_i = 0;
            M_aver_i = 0;
            M2_aver_i = 0;
            
            //считаем статсумму
            for(set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; ++set_of_states)
            {
                //вырождение * exp
                P[set_of_states] = boundaries_EM_blocks[boundary_dec][set_of_states].state_array.size()*
                        exp((double)(-(boundaries_EM_blocks[boundary_dec][set_of_states].E-Emin_array[boundary_dec])/T));
                Z += P[set_of_states];
                ///cout<<exp((double)(-(boundaries_and_states_array[boundary_dec][set_of_states].E-Emin)/T))<<endl;
                /// 
                /// ПЕРЕНЕСТИ СЮДА РАСЧЕТ СРЕДНИХ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
            }
            
            ///system("pause");
            
            ///cout << "Z = " << Z << endl;
            
            r = distribution_real(generator);
            ///cout<<"r= "<<r<<endl;
            
            //перебираем интервалы вероятностей
            interval=0;
            flag=0;
            //double sum=0;
            
            //считаем вероятности
            for(set_of_states=0; set_of_states<num_of_unique_EM_in_boundary; ++set_of_states)
            {
                P[set_of_states] /= Z;
                //sum+=P[set_of_states];
                ///cout << "P"<<set_of_states<<"= " << P[set_of_states] << endl;
                E_aver_i += P[set_of_states]*(old_E+boundaries_EM_blocks[boundary_dec][set_of_states].E);
                M_aver_i += P[set_of_states]*(My+boundaries_EM_blocks[boundary_dec][set_of_states].M);
                
                E2_aver_i += P[set_of_states]*(old_E+boundaries_EM_blocks[boundary_dec][set_of_states].E)*(old_E+boundaries_EM_blocks[boundary_dec][set_of_states].E);
                M2_aver_i += P[set_of_states]*(My+boundaries_EM_blocks[boundary_dec][set_of_states].M)*(My+boundaries_EM_blocks[boundary_dec][set_of_states].M);
                
                interval+=P[set_of_states];
                if(r<=interval && flag==0)
                {
                    total_set_of_states = set_of_states;
                    flag = 1;
                }
            }
            ///cout<<"E_aver_i = "<<E_aver_i<<endl;
            //cout<<"P_sum = " << sum << endl;
            
            //средние
            E_aver += (E_aver_i - E_aver) / (double)(MCS + 1);
            //E2_aver += (E_aver_i*E_aver_i - E2_aver) / (double)(MCS + 1);
            E2_aver += (E2_aver_i - E2_aver) / (double)(MCS + 1);
            
            Mpr_aver += (abs(M_aver_i) - Mpr_aver) / (double)(MCS + 1);
            //Mpr2_aver += ((M_aver_i*M_aver_i) - Mpr2_aver) / (double)(MCS + 1);
            Mpr2_aver += (M2_aver_i - Mpr2_aver) / (double)(MCS + 1);
            
            ///cout<<"E_aver = "<<E_aver<<endl;
            ///system("pause");
            
            ///cout<<"r<="<<interval<<endl;
            
            My += boundaries_EM_blocks[boundary_dec][total_set_of_states].M;
            old_E += boundaries_EM_blocks[boundary_dec][total_set_of_states].E;
            
            //выбираем случайную конфигурацию с одинаковыми параметрами
            rand_state = rand() % boundaries_EM_blocks[boundary_dec][total_set_of_states].state_array.size();
            
            //запоминаем конфигурацию
            rand_state = boundaries_EM_blocks[boundary_dec][total_set_of_states].state_array[rand_state];
            
            //переводим в двоичную систему
            ///cout<<"rand_state = "<<rand_state<<endl;
            
            ///ПЕРЕДЕЛАЛ!!! ПРОВЕРИТЬ!!!
            count=0;
            for(bit=num_of_spins_in_block-1; bit>=0; --bit)
            {
                array_1D[block_coordinates_for_spins[spin_num][count]] = (1 & rand_state >> bit);
                
                if(array_1D[block_coordinates_for_spins[spin_num][count]]==0) 
                    array_1D[block_coordinates_for_spins[spin_num][count]] = -1;
                
                ///cout << block_temp[count]<< " ";
                count++;
            }
            ///cout<<endl;
            //здесь заканчивается проход блока по системе 
        }
        
        //cout << "Mpr_aver = " << Mpr_aver << endl << endl;
        
        cout << "E_aver = " << E_aver << endl;
        
        double C = ((E2_aver - E_aver*E_aver)/(T*T))/(N*N);
        cout << "C = " << C << endl;
        C_data << T << "\t" << C << endl;
        if(C1>C) Tstep /= 1.1; 
        else Tstep *= 1.1; 
        ///C1=C;
        
        double hi = ((Mpr2_aver - Mpr_aver*Mpr_aver)/T) / (double)(N*N);
        cout << "hi = " << hi << endl << endl;
        hi_data << T << "\t" << hi << endl;
        
        E_data << T << "\t" << E_aver << endl;
        
        ///system("pause");
        
        
        if(old_E < (-2*N*N) || old_E > (2*N*N)) cout<<"ERROR2 E = "<<old_E<<endl;
        
        //здесь заканчивается перебор температур 
    }
    
    C_data.close();
    hi_data.close();
    E_data.close();
    
    /** <<<<--------------------------------------------------- */
    
    ///cout << "\n" << (double)clock()/CLOCKS_PER_SEC << " sec.\n";
    ///system("pause");
}
