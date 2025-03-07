﻿using System;
using System.Collections.Generic;
using System.Text;

namespace Optimum
{
    delegate double Fun(double x);
    delegate double Fun2(Vector x);
    delegate Vector Fun3(Vector x);
    class Extremum
    {
        #region первый семестр
        public static double Shagoviy(double xn, double h, double eps, Fun fun)
        {
            int k = 1;
            double fn = fun(xn);//вычисление функции начального значения
            double xs = xn + h;//прибавление приближения
            double fs = fun(xs);//вычисление функции нового х
            while (Math.Abs(h) > eps)
            {
                if (fs > fn)
                {
                    h = -h / 2;//шаг не в ту сторону

                }
                else
                {

                    h = h * 1.2;//шаг в правильную сторону

                }
                xn = xs;
                fn = fs;
                xs += h;
                fs = fun(xs);
                k++;
            }
            /*Console.WriteLine("k={0}", k);*/
            return xn;
        }
        public static double GoldSechenie(double a, double b, double eps, Fun fun)
        {
            int k = 2;
            double fi = (Math.Sqrt(5) - 1) / 2;
            double l = b - a;
            double v = a + l * (1 - fi);
            double w = a + l * fi;
            double fv = fun(v);
            double fw = fun(w);
            while (l > eps)
            {
                if (fv < fw)
                {
                    b = w;
                    w = v;
                    fw = fv;
                    l = b - a;
                    v = a + l * (1 - fi);
                    fv = fun(v);
                }
                else
                {
                    a = v;
                    v = w;
                    fv = fw;
                    l = b - a;
                    w = a + l * fi;
                    fw = fun(w);
                }
                k++;
            }
            /*Console.WriteLine("k={0}", k);*/
            return (a + b) / 2;
        }
        public static double MetodKvadAproksim(double x1, double x2, double x3, double eps, Fun fun)
        {
            double xpred = DopolenieMedoda(x1, x2, x3, fun);
            double xtek = 0;
            while (Math.Abs(xpred - xtek) < eps)
            {
                if ((x1 > x2) && (x1 > x3))
                {
                    x1 = xpred;
                    xtek = DopolenieMedoda(x1, x2, x3, fun);
                }
                if ((x2 > x1) && (x2 > x3))
                {
                    x2 = xpred;
                    xtek = DopolenieMedoda(x1, x2, x3, fun);
                }
                if ((x3 > x1) && (x3 > x2))
                {
                    x3 = xpred;
                    xtek = DopolenieMedoda(x1, x2, x3, fun);
                }

            }
            return xpred;
        }

        private static double DopolenieMedoda(double x1, double x2, double x3, Fun fun)
        {
            double f1 = fun(x1);
            double f2 = fun(x2);
            double f3 = fun(x3);
            double c = ((f3 - f1) * (x2 - x1) - (f2 - f1) * (x3 - x1)) / ((x3 * x3 - x1 * x1) * (x2 - x1) - (x2 * x2 - x1 * x1) * (x3 - x1));
            double b = ((f2 - f1) - c * (x2 * x2 - x1 * x1)) / (x2 - x1);
            double a = f1 - b * x1 - c * x1 * x1;
            double xpred = -b / 2 * c;
            double fm = xpred * xpred * c + b * xpred + a;
            return xpred;
        }
        public static Vector Grad(Vector xn, double eps, Fun2 func, double h) {
            Vector xnach = new Vector(xn);
            int k = 1;
            int n = xn.Size;
            double fn = func(xn);
            Vector xs = new Vector(xn);
            Vector gr;
            Vector dx;
            double delta = 0.5 * eps;
            double fs = fn + func(xs);
            do
            {
                gr = new Vector(n);
                for (int i = 0; i < n; i++)
                {

                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    gr[i] = (func(xg) - fn) / delta;
                }
                dx = -h * gr;
                xs = xn - h * gr;
                fs = func(xs);
                if (fs > fn)
                {
                    h = -h / 2;
                }
                else
                {
                    h = 1.2 * h;

                }
                xn = xs;
                fn = fs;
                k++;

            }
            while (Math.Abs(dx.NormaE()) > eps);
            Console.WriteLine(k);
            return xs;

        }
        public static Vector ModifiyGrad(Vector xn, double eps, Fun2 func)
        {
            Vector xnach = new Vector(xn);
            int k = 1;
            int n = xn.Size;
            double fn = func(xn);
            Vector xs = new Vector(xn);
            Vector gr;
            Vector dx;
            double delta = 0.5 * eps;
            double fs = fn + func(xs);
            double h = 0.1;

            do
            {
                //формирование вектора градиента
                gr = new Vector(n);
                for (int i = 0; i < n; i++)
                {
                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    gr[i] = (func(xg) - fn) / delta;

                }

                h = Shagoviy(h, 0.5, eps, he => func(xn - he * gr));//нахождение оптимального шага путем нахождения экстремума 1 переменной
                dx = -h * gr;
                xs = xn - h * gr;//вычистывание нового х
                fs = func(xs);//новое значение функции

                dx = -h * gr;
                xn = xs;
                fn = fs;
                k++;

            }
            while (Math.Abs(dx.NormaE()) > eps);
            Console.WriteLine(k);
            return xs;



        }

        public static Vector GradSoprIspra(Vector xn, double eps, Fun2 func)
        {
            Vector xnach = new Vector(xn);
            int k = 1;
            int n = xn.Size;
            double fn = func(xn);
            Vector xs = new Vector(xn);
            Vector gr = new Vector(n);
            Vector dx;
            double delta = 0.5 * eps;
            double fs = fn + func(xs);
            double h = 0.01;

            Vector grold = new Vector(n);//последнее направление
            double kof;

            for (int i = 0; i < n; i++)
            {
                Vector xg = xn.Copy();
                xg[i] = xg[i] + delta;
                gr[i] = (func(xg) - fn) / delta;

            }


            h = Shagoviy(h, 0.5, eps, he => func(xn - he * gr));//Вычисление опримального шага с помощью эксремума

            dx = -h * gr;
            xs = xn - h * gr;
            fs = func(xs);
            grold = gr.Copy();
            xn = xs;
            fn = fs;
            k++;
            //выше высчиталось первое направление для определение следующего
            do
            {

                for (int i = 0; i < n; i++)
                {
                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    gr[i] = (func(xg) - fn) / delta;

                }


                h = Shagoviy(h, 0.5, eps, he => func(xn - he * gr));
                //опредление оптимального шага

                kof = gr.NormaE() * gr.NormaE() / (grold.NormaE() * grold.NormaE());//высчитывание коэффицента для определения значимости предыдушего направления
                gr = gr + kof * grold;//высчитываение нового направление с учетом прошлого
                dx = -h * gr;
                xs = xn - h * gr;
                fs = func(xs);
                grold = gr.Copy();
                xn = xs;
                fn = fs;
                k++;

            }
            while (Math.Abs(dx.NormaE()) > eps);
            Console.WriteLine(k);
            return xs;


        }
        public static Vector MSP(Vector xn, int n, double h, double eps, Fun2 func)
        {
            Random rnd = new Random();
            int ni = xn.Size;
            Vector xs = new Vector(ni);
            Vector[] xpr = new Vector[3 * n];
            Vector fxp = new Vector(3 * n);
            double k = 0;

            do
            {
                xs = xn - EdRandR(2, rnd);
                double fn = func(xn);
                double fss = func(xs);


                for (int i = 0; i < 3 * n; i++) {
                    Vector rnvec = EdRandR(2, rnd);
                    Vector po = xn - h * rnvec;
                    xpr[i] = po;
                    double fs = func(xpr[i]);
                    fxp[i] = fs;
                    if (fss > fs) { fss = fs; xs = xpr[i]; }

                }

                if (fss < fn)
                {
                    xn = xs;
                    h = h * 1.2;
                }
                else
                {
                    h = h * 0.5;
                }
                k++;
            } while (h > eps);
            Console.WriteLine(k);
            return xn;

        }
        public static Vector NelderaMida(Vector[] xn, double eps, Fun2 func)
        {
            int n = xn.Length;
            double m = n + 1;


            Vector xotr = new Vector(xn[0].Size);
            double cps = 1.0 / (n - 1.0);

            double k = 0;
            double eeeps = 10;
            double[] fmas = new double[n];

            do
            {
                Vector xc = new Vector(xn[0].Size);//нахождение среднего 
                double kof = 1;
                for (int i = 0; i < n; i++)
                {
                    fmas[i] = func(xn[i]);
                }
                fmas = ShellSort(fmas);//Определение максимальных и минимальных значений
                //Высчитывание среднего значения
                for (int i = 0; i < n; i++)
                {
                    if (func(xn[i]) == fmas[n - 1]) { }
                    else
                    {
                        xc += xn[i] * cps;
                    }
                }
                //определение х отраженного. Максимальный ищется из функции FindXN
                xotr = xc + kof * (xc - FindXN(xn, fmas[n - 1], func));
                double fotr = func(xotr);
                //правильный шаг,растяжение с увлечением шага
                if (fotr < fmas[0])
                {
                    double newfotr = fotr;

                    eeeps = Math.Abs(func(FindXN(xn, fmas[n - 1], func)) - func(xotr));
                    if (fotr < fmas[1])
                    {
                        kof = kof * 2;
                        xotr = xc + kof * (xotr - xc);
                    }
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] == FindXN(xn, fmas[n - 1], func))
                        {
                            xn[i] = xotr.Copy();
                        }
                    }
                }
                //Расстяжение не полное
                else if (fotr > fmas[0])
                {
                    kof = kof / 2;
                    xotr = xc + kof * (FindXN(xn, fmas[n - 1], func) - xc);
                    eeeps = Math.Abs(-func(FindXN(xn, fmas[n - 1], func)) + func(xotr));
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] == FindXN(xn, fmas[n - 1], func))
                        {
                            xn[i] = xotr.Copy();
                        }
                    }

                }
                //Выполнение сжатия
                else
                {
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] != FindXN(xn, fmas[0], func))
                        {
                            for (int j = 0; j < xn[i].Size; j++)
                            {
                                xn[i][j] = xn[i][j] / 2;
                            }
                        }
                    }
                }
                k++;

            } while (eps < eeeps);// https://habr.com/ru/post/332092/
            Console.WriteLine(k);
            return FindXN(xn, fmas[0], func);




        }
        static Random rnd = new Random();
        public static Vector ozuMethod(Vector x, double eps, double h, Vector fn, Vector fv, Vector tipf, Fun3 func)
        {
            int k = 0;
            int n = x.Size;
            int m = 3 * n;
            int pmin = int.MinValue;
            Matrix xpTemp = new Matrix(m, n);
            Vector xp = new Vector(n);
            Vector fp = new Vector(m);
            double ft = fNorm(x, fn, fv, tipf, func);//функция текущая вычисляется по типу огранижений
            double fpmin = double.MaxValue;
            double temp, len;
            do
            {
                k++;
                for (int i = 0; i < m; i++)
                {
                    //Вычисление нового рандомного значения
                    for (int j = 0; j < n; j++)
                    {
                        temp = rnd.NextDouble() - 0.5;
                        xpTemp[i, j] = temp;
                    }
                    len = 0;
                    //Нахождение длинны,и временного значения xp
                    for (int j = 0; j < n; j++)
                        len += (xpTemp[i, j] - x[j] * (xpTemp[i, j] - x[j]));
                    len = Math.Sqrt(len);
                    for (int j = 0; j < n; j++)
                        xpTemp[i, j] = x[j] + h * (xpTemp[i, j] - x[j]) / len;
                }
                for (int i = 0; i < m; i++)
                {
                    //Присваниевение xp временные значения
                    for (int j = 0; j < n; j++)
                    {
                        xp[j] = xpTemp[i, j];
                    }
                    fp[i] = fNorm(xp, fn, fv, tipf, func);//вычисление новых значений функции
                    if (fp[i] < fpmin)
                    {
                        fpmin = fp[i];
                        pmin = i;
                    }

                }
                //правильный шаг
                if (fpmin < ft)
                {
                    for (int j = 0; j < n; j++)
                        x[j] = xpTemp[pmin, j];
                    ft = fpmin;
                    h = h * 1.2;

                }
                else h = h / 2.0;//неправльный шаг
            } while (h > eps);
            return x;
        }
        #endregion
        #region Второй семестр
        public static Matrix Init(Vector a, Vector b, Matrix c)
        {
            //a количество товара у поставщиков
            //b количество товара у потребителя
            //с стоимость перевозки ij

            if (!Equaleble(a, b, c))
            {
                Console.WriteLine("Error");
                return new Matrix(1, 1);
            }
            Matrix x = new Matrix(a.Size, b.Size);
            double s = 0;

            List<double[]> min_znach = FindMinInMatrix(c.Copy());
            foreach (double[] kord in min_znach)
            {
                if (AllEmpty(a)) break;
                int i = (int)kord[0];
                int j = (int)kord[1];
                if (a[i] == 0) continue;
                if (b[j] == 0) continue;
                if (a[i] > b[j])
                {
                    x[i, j] = b[j];
                    s += b[j] * c[i, j];
                    a[i] -= b[j];
                    b[j] = 0;
                }
                else if (a[i] < b[j])
                {
                    x[i, j] = a[i];
                    b[j] -= a[i];
                    s += a[i] * c[i, j];
                    a[i] = 0;
                }
                else
                {
                    x[i, j] = a[i];
                    s += b[j] * c[i, j];
                    a[i] = 0;
                    b[j] = 0;
                }
            }
            Matrix delta = new Matrix(a.Size, b.Size);
            List<double[]> notnull = FindNotNullMat(x);


            Console.WriteLine("Сумма: {0} ", s);
            return x;
        }
        
        public static void Bliz2(Graph gr,Vertex start)
        {
        
            Console.WriteLine("........................");
            double mi = 200;
            Edge tem = null;
            Vertex next = start;
            start.visited=true;
            List<Edge> edgesnn = new List<Edge>();
            int k = 0;
            do
            {
                mi = 200;
                if (gr.allvertexs.Count - 1 == k)
                {
                    foreach (Edge ee in gr.alledges)
                    {
                        if (ee.First == next && ee.End == start)
                        {
                            edgesnn.Add(ee);
                            break;
                        }
                    }
                    break;

                }
                foreach (Edge e in next.Edges)
                {
                    if (e.First == next && e.End.visited != true)
                    {
                        if (mi > e.Length)
                        {
                            mi = e.Length;
                            tem = e;

                        }
                    }
                    
                }
                /*}*/
                edgesnn.Add(tem);
                tem.End.prev = next;
                next = tem.End;
                
                next.visited = true;
                k++;
            } while (start != next);
            double summ = 0;  
            foreach(Edge e in edgesnn)
            {
               Console.Write(String.Format( e.ToString()+"-->" ));
                summ += e.Length;
            }
            Console.WriteLine("Сумма:");
            Console.Write(summ);
           
        }
        
        public static int Bag(List<int[]> pi, int vob)
        {
            int n = pi.Count;
            List<int> used = new List<int>();
            Matrix m = new Matrix(3, n);
            int sum = 0;
            int vi = 0;
            for (int i = 0; i < n; i++)
            {
                m[0, i] = pi[i][0];
                m[1, i] = pi[i][1];
                if (pi[i][1] == 0) throw new Exception();
                m[2, i] = (double)pi[i][0] / pi[i][1];
            }
            while (vob >= vi)
            {

                int i = FindMaxiesNext(m, used);
                if (i == -1) break;
                if (vi + (int)m[1, i] > vob) { used.Remove(i); break; }
                sum += (int)m[0, i];
                vi += (int)m[1, i];
            }

            foreach (int use in used)
            {

                for (int i = 0; i < n; i++)
                {

                    if (!HaveList(i, used))
                    {
                        int vinew = vi;
                        int sumnew = sum;
                        vinew -= (int)m[1, use];
                        sumnew -= (int)m[0, use];
                        vinew += (int)m[1, i];
                        sumnew += (int)m[0, i];
                        if (vob >= vinew && sumnew > sum) { vi = vinew; sum = sumnew; }
                    }

                }
            }

            return sum;

        }
        public static void TeoriaIgrFirst(Matrix A)
        {
            int n = A.GetCountRows();
            int m = A.GetCountColumns();
            double[] a = new double[2];
            double[] b = new double[2];
            Vector k = new Vector(2);
            Vector k1 = new Vector(3);
            double[] z = new double[2];
            for (int i = 0; i < n; i++)
            {

                k1[i] = FindMinIG(A.GetRow(i))[0];
            }
            a = FindMaxIG(k1);

            for (int i = 0; i < m; i++)
            {
                k1[i] = FindMaxIG(A.GetColumn(i))[0];
            }

            b = FindMinIG(k1);

            if (a[0] == b[0]) Console.WriteLine(String.Format("Седловая с координатами {0},{1}", a[1], b[1]));

        }

        public static double GameTheory(double[,] matrix)
        {
            if (matrix.GetLength(0) > 2)
            {
                throw new Exception("Рядов больше 2х");

                return 0;
            }

            double min = double.MaxValue;

            int iMin = 0;

            int jMin = 0;

            for (var i = 0; i < matrix.GetLength(1); i++)
            {
                var squareMatrix = new double[2, 2];

                squareMatrix[0, 0] = matrix[0, i];

                squareMatrix[1, 0] = matrix[1, i];

                for (var j = 0; j < matrix.GetLength(1); j++)
                {
                    if (i != j && i < j)
                    {
                        squareMatrix[0, 1] = matrix[0, j];

                        squareMatrix[1, 1] = matrix[1, j];

                        var del = (squareMatrix[0, 0] + squareMatrix[1, 1] - squareMatrix[0, 1] - squareMatrix[1, 0]);

                        if (del != 0)
                        {
                            var p1 = (squareMatrix[1, 1] - squareMatrix[1, 0]) / del;

                            var p2 = 1 - p1;

                            var q1 = (squareMatrix[1, 1] - squareMatrix[0, 1]) / del;

                            var q2 = 1 - q1;

                            if (p1 < 0 || p2 < 0 || q1 < 0 || q2 < 0)
                            {
                                continue;
                            }

                            var v = squareMatrix[0, 0] * p1 * q1 + squareMatrix[0, 1] * p1 * q2 + squareMatrix[1, 0] * p2 * q1 + squareMatrix[1, 1] * p2 * q2;

                            if (v < min)
                            {
                                min = v;

                                iMin = i;

                                jMin = j;
                            }
                        }
                    }
                }
            }

            Console.WriteLine($"Столбец 1: {iMin + 1}");

            Console.WriteLine($"Ряд 1: {jMin + 1}");

            return min;
        }

        public static void Potok(Graph gr, Vertex start)
        {
            double potok = 0;
            foreach (var v in gr.allvertexs)
            {

                v.prev = null;
                v.visited = false;
            }
            start.Weight = 0;
            start.visited = true;
            start.prev = null;
            Queue<Vertex> que = new Queue<Vertex>();
            que.Enqueue(start);//включать De выключать
            while (que.Count > 0)
            {

                Vertex u = que.Dequeue();
                double min = 0;
                foreach (Edge e in u.Edges)
                {
                    Vertex r = e.End;
                    if (r.visited == false)
                    {
                        if (min < e.Length)
                        {
                            min = e.Length;
                        }
                        r.visited = true;
                        r.Weight = u.Weight + 1;
                        r.prev = u;
                        que.Enqueue(r);
                    }
                }
                u.visited = true;
            }
        }

        public static void Potok2(Graph gr, Vertex startVertex)
        {

            foreach (Vertex vv in gr.allvertexs)
            {
                vv.visited = false;
                vv.prev = null;
            }

            startVertex.time = 0;

            foreach (Vertex vv in gr.allvertexs)
            {
                if (vv.visited == false)
                {
                    DFS_Visit2(startVertex);
                }
            }



        }
        private static void DFS_Visit2(Vertex u) //Для DFS
        {
            u.visited = true;
            u.discovered = u.time;
            u.time += 1; //Счетчик времени увеличиваем на 1

            foreach (Edge ee in u.Edges)
            {
                Vertex v = ee.End;
                if (v.visited == false)
                {
                    v.prev = u;
                    DFS_Visit2(v);
                }

            }
        }
        private static double[] FindMinIG(Vector v)
        {
            double[] res = new double[2];
            double min = Double.MaxValue;
            int j = -1;
            for (int i = 0; i < v.size; i++)
            {
                if (v[i] < min) { min = v[i]; j = i; }
            }
            res[0] = min;
            res[1] = j;
            return res;
        }
        private static double[] FindMaxIG(Vector v)
        {
            double[] res = new double[2];
            double max = Double.MinValue;
            int j = -1;
            for (int i = 0; i < v.size; i++)
            {
                if (v[i] > max) { max = v[i]; j = i; }
            }
            res[0] = max;
            res[1] = j;
            return res;
        }
        private static int FindMaxiesNext(Matrix m,List<int> list)
        {
            int max = 0;
            int add = -1;
            for (int i = 0; i < m.GetCountColumns(); i++)
            {
                if (max < m[2, i] && !HaveList(i,list)) {//эквалис лучше

                    max = (int)m[2, i];
                    add = i;
                }
            }
            if (add != -1) list.Add(add);
            return add;
        }
        private static bool HaveList(int i,List<int> lsit )
        {
            /*if (lsit.Count == 0) return true;*/
            foreach(int ii in lsit)
            {
                if (i == ii) return true;
            }
            return false;
        }

        private static List<double[]> FindMinInMatrix(Matrix m)
        { 
            int sizeMas = m.GetCountRows() * m.GetCountColumns();
            
            List<double[]> ress = new List<double[]>();
            double min = 100;
            int k = 0;
            while (k < sizeMas)
            {
                double[] res = new double[2];
                min = 100;
                for (int i = 0; i < m.GetCountRows(); i++)
                {
                    for (int j = 0; j < m.GetCountColumns(); j++)
                    {
                        if (min > m[i, j])
                        {
                            min = m[i, j];
                            res[0] = i;
                            res[1] = j;
                        }
                    }
                }
                
                ress.Add(res);
                 m[(int)res[0], (int)res[1]] = 100;
                k++;
                        
                    
                
            }
            return ress;
        }
        private static List<double[]> FindNotNullMat(Matrix m)
        {
            int sizeMas = m.GetCountRows() * m.GetCountColumns();

            List<double[]> ress = new List<double[]>();
            
            
              double[] res = new double[2];

            for (int i = 0; i < m.GetCountRows(); i++)
            {
                for (int j = 0; j < m.GetCountColumns(); j++)
                {
                    if (0 < m[i, j])
                    {
                        res[0] = i;
                        res[1] = j;
                        ress.Add(res);
                    }
                }
            }
            return ress;
        }
        private static bool AllEmpty(Vector a)
        {
            bool yeah=true;
            for(int i = 0; i < a.Size; i++)
            {
                if (a[i] != 0)
                {
                    yeah = false;
                }
            }
            return yeah;
        }
        private static bool Equaleble(Vector a,Vector b,Matrix c)
        {
            double suma = 0;
            double sumb = 0;
            bool res = true;
            for (int i=0; i < a.Size; i++)
            {
                suma += a[i];
            }
            for (int i = 0; i < b.Size; i++)
            {
                sumb += b[i];
            }
            if (suma != sumb) res = false;
            if((a.Size!=c.GetCountRows())&& (b.Size != c.GetCountColumns())){
                Console.WriteLine("Неправильно заполнена матрица");
                res = false;
            }
            return res;
        }

        #endregion
        #region Первый семестр Приватные методы
        private static double fNorm(Vector x,Vector fn,Vector fv,Vector tipf,Fun3 func)
        {
            Vector fx = func(x);
            int sizeCriteria = fx.Size;
            Vector g = new Vector(sizeCriteria);
            int maxg =0;
            double g1;
            double g2;
            for(int i = 0; i < sizeCriteria; i++)
            {
                int t = (int)tipf[i];
                if (t == 1)//niz
                {
                    if (fn[i] > 0)
                        g[i] = 2 - fx[i] / fn[i];
                    if (fn[i] == 0)
                        g[i] = -fx[i] / fn[i] + 1;
                    if (fn[i] < 0)
                        g[i] = fx[i] / fn[i];
                }
                if (t == 2)//vverx
                {

                    if (fn[i] > 0)
                        g[i] = fx[i] / fv[i];
                    if (fn[i] == 0)
                        g[i] = fx[i] / fv[i] + 1;
                    if (fn[i] < 0)
                        g[i] =2- fx[i] / fv[i];
                }
                if (t == 3)//2-x storon
                {

                    g1 = (fx[i] - fn[i]) / (fv[i] - fn[i]);
                    g2 = (fv[i] - fx[i]) / (fv[i] - fn[i]);
                    if (g1 > g2)
                        g[i] = g1;
                    else
                        g[i] = g2;
                }
                if (g[i] >= g[maxg])
                    maxg = i;
            }
            return g[maxg];


        }
        private static double[] ShellSort(double[] list) 
        {
            //расстояние между элементами, которые сравниваются
            var d = list.Length / 2;
            while (d >= 1)
            {
                for (var i = d; i < list.Length; i++)
                {
                    var j = i;
                    while ((j >= d) && (list[j - d].CompareTo(list[j]) > 0))
                    {
                        Swap(ref list[j], ref list[j - d]);
                        j = j - d;
                    }
                }

                d = d / 2;
            }
            return list;
        }
        private static void Swap(ref double a, ref double b)
        {
            double c = a; a = b; b = c;
        }
        private static Vector FindXN(Vector[] xn,double res,Fun2 func)
        {
            Vector n = new Vector(xn[0].Size);
            for (int i = 0; i < xn.Length; i++)
            {
                if (func(xn[i]) == res) {
                    n = xn[i];
                }
                
                
            }
            if (n == new Vector(xn[0].Size)) throw new Exception();
            return n;
        }
        private static Vector EdRand(int n,Random rnd)
        {
            double[] e = new double[n];
            for(int i = 0; i < n; i++)
            {
                e[i] = rnd.NextDouble()-0.5 ;
            }
            Vector ee = new Vector(e);
            double[] eee = new double[1];
            eee[0] = ee.NormaE();
            return new Vector(eee);
        }
        private static Vector EdRandR(int n, Random rnd)
        {
            double[] e = new double[n];
            for (int i = 0; i < n; i++)
            {
                e[i] = rnd.NextDouble() - 0.5;
            }
            Vector ee = new Vector(e);
            return ee;
        }
        #endregion

    }
}
