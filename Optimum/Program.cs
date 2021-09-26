using System;

namespace Optimum
{
    
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("Ответ шаговый {0}",Extremum.Shagoviy(1,0.5,0.0001,x=>(x*x/2)+(8/(x*x))));
            Console.WriteLine("Золотое сечение {0}",Extremum.GoldSechenie(1, 5, 0.0001, x => (x * x / 2) + (8 / (x * x))));
            Console.WriteLine("Квад Апрокси {0}", Extremum.MetodKvadAproksim(1,2,5, 0.0001, x => (x * x)));
        }
    }
}
