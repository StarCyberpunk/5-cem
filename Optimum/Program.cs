using System;

namespace Optimum
{
    
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Ответ шаговый {0}",Extremum.Shagoviy(1,0.5,0.0001,x=>(x*x/2)+(8/(x*x))));
            Console.WriteLine("Золотое сечение {0}",Extremum.GoldSechenie(1, 5, 0.0001, x => (x * x / 2) + (8 / (x * x))));
            Console.WriteLine("Квад Апрокси {0}", Extremum.MetodKvadAproksim(-5,0,5, 0.0001, x => x * x+5*x));
            Vector xn = new Vector(2);
            xn[0] = 1;xn[1] = 1;
            Console.WriteLine("Градиент {0}", Extremum.Grad(xn,0.001, x => x[0]*x[0]+x[1]*x[0]+2*x[1]*x[1],0.1));
        }
    }
}
