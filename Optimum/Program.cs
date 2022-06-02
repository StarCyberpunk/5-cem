
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



            Graph graph = new Graph();
            Vertex v1 = new Vertex("1");
            Vertex v2 = new Vertex("2");
            Vertex v3 = new Vertex("3");
            Vertex v4 = new Vertex("4");
            Vertex v5 = new Vertex("5");
           

            graph.AddVertex(v1);
            graph.AddVertex(v2);
            graph.AddVertex(v3);
            graph.AddVertex(v4);
            graph.AddVertex(v5);
          


            graph.AddEdge(v1, v2,12);
            graph.AddEdge(v1, v4,8);
            graph.AddEdge(v1, v5,9);
            graph.AddEdge(v1, v3,5);
           graph.AddEdge(v2, v1,12);
            graph.AddEdge(v2, v3,10);
            graph.AddEdge(v2, v4,6);
            graph.AddEdge(v2, v5,4);
            graph.AddEdge(v3, v1,5);
            graph.AddEdge(v3, v2,10);
            graph.AddEdge(v3, v4,8);
            graph.AddEdge(v3, v5,7);
            graph.AddEdge(v4, v1,8);
           graph.AddEdge(v4, v2,6);
            graph.AddEdge(v4, v3,8);
            graph.AddEdge(v4, v5,12);
           graph.AddEdge(v5, v1,9);
            graph.AddEdge(v5, v2,4);
            graph.AddEdge(v5, v3,7);
           graph.AddEdge(v5, v4,12);


            graph.ViewGraph();

            Extremum.Bliz2(graph, v2);


            Graph potok2=new Graph();
            Vertex ps = new Vertex("s");
            Vertex p2 = new Vertex("1");
            Vertex p3 = new Vertex("2");
            Vertex p4 = new Vertex("3");
            Vertex p5 = new Vertex("4");
            Vertex pt = new Vertex("t");

            potok2.AddVertex(ps);
            potok2.AddVertex(p2);
            potok2.AddVertex(p3);
            potok2.AddVertex(p4);
            potok2.AddVertex(p5);
            potok2.AddVertex(pt);


            potok2.AddEdge(ps, p2, 3);
            potok2.AddEdge(ps, p3, 5);
            potok2.AddEdge(p3, p2, 4);
            potok2.AddEdge(p2, p4, 5);
            potok2.AddEdge(p3, p5, 2);
            potok2.AddEdge(p4, p5, 6);
            potok2.AddEdge(p5, pt, 7);
            potok2.AddEdge(p4, pt, 5);



            /* List<int[]> pi = new List<int[]>();
             int[] p1 = new int[2] { 5,2};
             int[] p2 = new int[2] { 6, 3 };
             int[] p3 = new int[2] { 3, 4 };
             int[] p4 = new int[2] { 2, 1 };
             int[] p5 = new int[2] { 8, 5 };
             int[] p6 = new int[2] { 4, 2 };
             pi.Add(p1);
             pi.Add(p2);
             pi.Add(p3);
             pi.Add(p4);
             pi.Add(p5);
             pi.Add(p6);
             int vob = 10;
             Extremum.Bag(pi, vob);
            
            double[,] n12 = new double[3, 3] { { 1, -3, 2 },{ 0, 5, 4 },{ 2, 3, 2 } };
            Console.WriteLine("\nПервый пример");

            var matrix = new double[2, 4];

            matrix[0, 0] = 4;

            matrix[0, 1] = 2;

            matrix[0, 2] = 3;

            matrix[0, 3] = -1;

            matrix[1, 0] = -4;

            matrix[1, 1] = 0;

            matrix[1, 2] = -2;

            matrix[1, 3] = 2;

            Console.WriteLine($"Vсм = {Extremum.GameTheory(matrix)}");

            Console.WriteLine("\nВторой пример");

            var matrix2 = new double[2, 4];

            matrix2[0, 0] = 2;

            matrix2[0, 1] = 3;

            matrix2[0, 2] = 1;

            matrix2[1, 0] = 1;

            matrix2[1, 1] = 2;

            matrix2[1, 2] = 5;

            Console.WriteLine($"Vсм = {Extremum.GameTheory(matrix2)}");

            Console.ReadKey();
            Matrix m12 = new Matrix(n12);
            Extremum.TeoriaIgrFirst(m12);*/
            Console.ReadKey();
        }
    }
}
