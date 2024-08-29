#include<iostream>
#include <stdio.h>
#include <math.h>
#include<stdlib.h>
#include<limits.h>
using namespace std;
#define Max_Row 10
#define Max_Column 10
int Num_Of_Constraints, Num_Of_Var, NOPTIMAL, Pivot_X, Pivot_Y, checker;

double Matrix[Max_Row][Max_Column];
void Taking_Input();
void Implementation();
void Calculate_Pivot_xy();
void Table_Transformantion();
void Optimization();


void Taking_Input()
{

    double R2;

    int i, j;

    printf("\n SIMPLEX METHOD FOR MAXIMISING\n\n");

    printf("\n NUMBER OF VARIABLES ? ");

    scanf("%d", &Num_Of_Var);

    printf("\n NUMBER OF CONSTRAINTS ? ");
    scanf("%d", &Num_Of_Constraints);

    printf("\n COEFFICIENTS OF OBJECTIVE FUNCTION:\n");

    for (j = 1; j <= Num_Of_Var; j++)
    {

        printf(" #%d ? ", j);
        scanf("%lf", &R2);

        Matrix[1][j + 1] = R2; // objective function 
    }

    printf(" constant d in objective function ? ");
    scanf("%lf", &R2);  // if there is a constant added in objective function 

    Matrix[1][1] = R2;

    for (i = 1; i <= Num_Of_Constraints; i++)
    {

        printf("\n CONSTRAINT #%d:\n", i);

        for (j = 1; j <= Num_Of_Var; j++)
        {

            printf(" #%d ? ", j);
            scanf("%lf", &R2);

            Matrix[i + 1][j + 1] = -R2;  // negative value in simplex table 
        }

        printf(" Right hand side ? ");
        scanf("%lf", &Matrix[i + 1][1]);
    }

    printf("\n\n RESULT OBTAINED :\n\n");

    for (j = 1; j <= Num_Of_Var; j++)
        Matrix[0][j + 1] = j; // numbering the interested columns 

    for (i = Num_Of_Var + 1; i <= Num_Of_Var + Num_Of_Constraints; i++)
        Matrix[i - Num_Of_Var + 1][0] = i; // numbering the interested rows
    // for(int i=0;i<10;i++)
    // {
    //     for(int j=0;j<10;j++)
    //     {
    //         cout<<Matrix[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
}



void Implementation()
{
    
top:
    
    Calculate_Pivot_xy();

    Table_Transformantion();

    Optimization();

    if (NOPTIMAL == 1)
        goto top;
}

void Calculate_Pivot_xy()
{

    double min, V, max;

    int i, j;

    max = INT_MIN;

    for (j = 2; j <= Num_Of_Var + 1; j++)
    {

        if (Matrix[1][j] > 0.0 && Matrix[1][j] > max)
        {

            max = Matrix[1][j];

            Pivot_Y = j; 
        }
    }

    min = INT_MAX;

    for (i = 2; i <= Num_Of_Constraints + 1; i++)
    {

        if (Matrix[i][Pivot_Y] >= 0.0)
            continue;

        V = fabs(Matrix[i][1] / Matrix[i][Pivot_Y]); //minimum positive ratio calculation 

        if (V < min)
        {

            min = V;

            Pivot_X = i;
        }

    }

    V = Matrix[0][Pivot_Y]; // interchanging incoming and outgoing numbering of rows and columns 
    Matrix[0][Pivot_Y] = Matrix[Pivot_X][0];
    Matrix[Pivot_X][0] = V;
    // for(int i=0;i<10;i++)
    // {
    //     for(int j=0;j<10;j++)
    //     {
    //         cout<<Matrix[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
    
}

void Table_Transformantion()
{

    int i, j;

    for (i = 1; i <= Num_Of_Constraints + 1; i++)
    {

        if (i == Pivot_X)
           { 
           continue;
           }

        for (j = 1; j <= Num_Of_Var + 1; j++)
        {

            if (j == Pivot_Y)
             
             { 
                continue;}

            Matrix[i][j] = Matrix[i][j] - (Matrix[Pivot_X][j] * Matrix[i][Pivot_Y] / Matrix[Pivot_X][Pivot_Y]);  // q=(pq -rs)/p =>  q=q-(rs/p)

        
        }

    
    }

    Matrix[Pivot_X][Pivot_Y] = 1.0 / Matrix[Pivot_X][Pivot_Y]; // p=1/p

    for (j = 1; j <= Num_Of_Var + 1; j++)
    {

        if (j == Pivot_Y)
            continue;

        Matrix[Pivot_X][j] = Matrix[Pivot_X][j]*(fabs(Matrix[Pivot_X][Pivot_Y]));  // row division containing pivot element

    // e100:;
    }

    for (i = 1; i <= Num_Of_Constraints + 1; i++)
    {

        if (i == Pivot_X)
            continue;

        Matrix[i][Pivot_Y] =  Matrix[i][Pivot_Y] * (Matrix[Pivot_X][Pivot_Y]);  // column division containing pivot element

    }
    // for(int i=0;i<10;i++)
    // {
    //     for(int j=0;j<10;j++)
    //     {
    //         cout<<Matrix[i][j]<<" ";
    //     }
    //     cout<<endl;
    // }
}

void Optimization()
{

    int i, j;

    for (i = 2; i <= Num_Of_Constraints + 1; i++)

        if (Matrix[i][1] < 0.0)
            checker = 1;

    NOPTIMAL = 0;

    if (checker == 1)
        return;

    for (j = 2; j <= Num_Of_Var + 1; j++)

        if (Matrix[1][j] > 0.0) // if any z-c is positive we need more optimization 
            NOPTIMAL = 1;
}

void Results()
{

    int i, j;
    if(checker==0){
    for (i = 1; i <= Num_Of_Var; i++){

        for (j = 2; j <= Num_Of_Constraints + 1; j++)
        {

            if (Matrix[j][0] != 1.0 * i)
                continue;
 
            printf(" VARIABLE #%d: %f\n", i, Matrix[j][1]);

        }
    }
    printf("\n OBJECTIVE FUNCTION VALUE:  %f\n", Matrix[1][1]);
    }
    else
    {
        printf("NO SOLUTION\n");
        return;
    }

}

int main()
{

    Taking_Input();

    Implementation();

    Results();

    return 0;
}
