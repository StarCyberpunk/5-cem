using System;
using System.Collections.Generic;

namespace Optimum
{
    
    class Program
    {
        static void Main(string[] args)
        {
            #region First sem
            /*Console.WriteLine("Ответ шаговый {0}",Extremum.Shagoviy(1,0.5,0.0001,x=>(x*x/2)+(8/(x*x))));
            Console.WriteLine("Золотое сечение {0}",Extremum.GoldSechenie(1, 5, 0.0001, x => (x * x / 2) + (8 / (x * x))));
            Console.WriteLine("Квад Апрокси {0}", Extremum.MetodKvadAproksim(-5,0,5, 0.0001, x => x * x+5*x));
            Vector xn = new Vector(2);
            xn[0] = 1;xn[1] = 1;
            
            Console.WriteLine("Градиент {0}", Extremum.Grad(xn,0.001, x => x[0]*x[0]+x[1]*x[0]+2*x[1]*x[1],0.1));
*/
            /*Vector xn = new Vector(2);
            xn[0] = 0; xn[1] = 1;
            Console.WriteLine("Градиент Модификация {0}", Extremum.ModifiyGrad(xn, 0.001, x => x[0] * x[0] + x[1] * x[0] + 2* x[1] * x[1]+x[0]));
            Console.WriteLine("Градиент Модификация ++ {0}", Extremum.GradSoprIspra(xn, 0.001, x => x[0] * x[0] + x[1] * x[0] + 2 * x[1] * x[1] +x[0]));
            Vector[] xnn = new Vector[3];
            xnn[0] = xn;
           double[] xnnn = new double[2];
            xnnn[0] = 0;xnnn[1] = 0;
            xnn[1] = new Vector(xnnn);
            xnnn[0] = 1;xnnn[1] = 0;
            xnn[2] = new Vector(xnnn);

             
            Console.WriteLine("MSP {0}", Extremum.MSP(xn, 20, 0.1, 0.001, x => x[0] * x[0] + x[1] * x[0] + 2 * x[1] * x[1]+x[0]));
            Console.WriteLine("Метод деформированных многоугольноков {0}", Extremum.NelderaMida(xnn, 0.001, x => x[0] * x[0] + x[1] * x[0] + x[1]* x[1] -9* x[1] -6*x[0]));
            Vector ffff(Vector x)
            {
                Vector fc = new Vector(3);
                fc[0] = x[0] * x[0]+x[1]*x[0]+x[1]*x[1];
                fc[1] = x[0] * x[0]+x[1];
                fc[2] = x[0] * x[1]*x[1]*2;
                return fc;
            }

            xn[0] = 1;
            xn[1] = 3;
            
            Vector tipof = new Vector(3);
            tipof[0] = 3;tipof[1] = 2;tipof[2] = 1;
            Vector fv = new Vector(3);
            fv[0] = 3;fv[1] = 2;fv[2] = 10;
            Vector fn = new Vector(3);
            fn[0] = 1;fn[1] = -10;fn[2] = 1;
            Vector ai = new Vector(2);
            ai[0] = -2;ai[1] = 2;
            Vector bi = new Vector(2);
            bi[0] = -2;bi[1] = 2;

            Console.WriteLine("Метод ОЗУ {0}", Extremum.ozuMethod(xn, 0.0001, 0.1, fv, fn, tipof, ffff));*/
            #endregion
            Vector postav = new Vector(new double[] { 30, 40, 20 });
            Vector poluch = new Vector(new double[] { 20, 25, 30, 15 });
            Matrix stoimost = new Matrix(new double[,] { { 5, 3, 2, 4 }, { 6, 1, 7, 3 }, { 2, 5, 11, 8 } });
            Matrix res = Extremum.Init(postav, poluch, stoimost);
            Console.ReadKey();
        }
    }
}
