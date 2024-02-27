#include <stdio.h>
#include <math.h>
int main()
{
    int t;

    printf("The nth Fibonacci number");
    scanf("%d", &t);

    int myArray[] = {0, 1};
    int temp;
    if (t > -1 && t < 2)
    {
        temp = myArray[t - 1];
    }
    else if (t >= 2)
    {
        for (int i = 0; i < t - 1; i++)
        {
            temp = myArray[0] + myArray[1];
            myArray[0] = myArray[1];
            myArray[1] = temp;
        }
    }
    printf("%d", temp);
    return 0;
}