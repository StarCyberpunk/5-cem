using System;
using System.Collections.Generic;

namespace Optimum
{
    
    class Program
    {
        static void Main(string[] args)
        {
            /*Console.WriteLine("Ответ шаговый {0}",Extremum.Shagoviy(1,0.5,0.0001,x=>(x*x/2)+(8/(x*x))));
            Console.WriteLine("Золотое сечение {0}",Extremum.GoldSechenie(1, 5, 0.0001, x => (x * x / 2) + (8 / (x * x))));
            Console.WriteLine("Квад Апрокси {0}", Extremum.MetodKvadAproksim(-5,0,5, 0.0001, x => x * x+5*x));
            Vector xn = new Vector(2);
            xn[0] = 1;xn[1] = 1;
            
            Console.WriteLine("Градиент {0}", Extremum.Grad(xn,0.001, x => x[0]*x[0]+x[1]*x[0]+2*x[1]*x[1],0.1));
*/
            Vector xn = new Vector(2);
            xn[0] = 1; xn[1] = 1;
            Console.WriteLine("Градиент Модификация {0}", Extremum.ModifiyGrad(xn, 0.001, x => x[0] * x[0] + x[1] * x[0] + 2* x[1] * x[1]+x[0]));
            Console.WriteLine("Градиент Модификация ++ {0}", Extremum.GradSoprIspra(xn, 0.001, x => x[0] * x[0] + x[1] * x[0] + 2 * x[1] * x[1] +x[0]));
            Vector[] xnn = new Vector[3];
            xnn[0] = xn;
           double[] xnnn = new double[2];
            xnnn[0] = 0;xnnn[1] = 1;
            xnn[1] = new Vector(xnnn);
            xnnn[0] = 0;xnnn[1] = 0;
            xnn[2] = new Vector(xnnn);


            Console.WriteLine("MSP {0}", Extremum.MSP(xn, 20, 0.1, 0.001, x => x[0] * x[0] + x[1] * x[0] + 2 * x[1] * x[1]+x[0]));
            Console.WriteLine("Метод деформированных многоугольноков {0}", Extremum.MDN(xnn, 0.001, x => x[0] * x[0] + x[1] * x[0] + 2 * x[1] * x[1] + x[0]));
            

        }
    }
}
