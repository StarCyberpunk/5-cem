using System;
using System.Collections.Generic;
using System.Text;

namespace Optimum
{
    delegate double Fun(double x);
    delegate double Fun2(Vector x);
    class Extremum
    {
        public static double Shagoviy(double xn, double h, double eps, Fun fun)
        {
            int k = 1;
            double fn = fun(xn);
            double xs = xn + h;
            double fs = fun(xs);
            while (Math.Abs(h) > eps)
            {
                if (fs > fn)
                {
                    h = -h / 2;

                }
                else
                {

                    h = h * 1.2;

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
                dx= -h * gr;
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
                gr = new Vector(n);
                for (int i = 0; i < n; i++)
                {
                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    gr[i] = (func(xg) - fn) / delta;
                   
                }
                    h = Shagoviy(h, 0.5, eps, he => func(xn - he * gr));
                    dx = -h * gr;
                    xs = xn - h * gr;
                    fs = func(xs);
                
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
            Vector gr=new Vector(n);
            Vector dx;
            double delta = 0.5 * eps;
            double fs = fn + func(xs);
            double h = 0.01;
            
            Vector grold=new Vector(n);
            double kof;

            for (int i = 0; i < n; i++)
            {
                Vector xg = xn.Copy();
                xg[i] = xg[i] + delta;
                gr[i] = (func(xg) - fn) / delta;
                
            }


            h = Shagoviy(h, 0.5, eps, he => func(xn - he * gr));

            dx = -h * gr;
            xs = xn - h * gr;
            fs = func(xs);
            grold = gr.Copy();
            xn = xs;
            fn = fs;
            k++;

            do
            {
                
                for (int i = 0; i < n; i++)
                {
                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    gr[i] = (func(xg) - fn) / delta;
                    
                }


                h = Shagoviy(h, 0.5, eps, he => func(xn - he * gr));
                //настройка коэфов
                
                kof = gr.NormaE()*gr.NormaE()/(grold.NormaE()*grold.NormaE());
                gr = gr + kof * grold;
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
        public static Vector MSP(Vector xn,int n,double h, double eps,Fun2 func)
        {
            Random rnd = new Random();
            int ni = xn.Size;
            Vector xs=new Vector(ni);
            Vector[] xpr = new Vector[3*n];
            double fpr = 0;
            
            do
            {
                xs = xn - EdRandR(2, rnd);
                double fn = func(xn);
                double fss = func(xs);
                
                for (int i = 0; i<3*n;i++) {
                    Vector rnvec = EdRandR(2, rnd);
                    Vector po= xn - h * rnvec;
                    xpr[i] = po;
                    double fs = func(xpr[i]);
                    if (fss > fs) { fss = fs; xs = xpr[i]; }
                    
                }
                if (fpr == fss) break;
                if (fss < fn)
                {
                    xn = xs;
                    h = h * 1.2;
                }
                else
                {
                    h = h * 0.5;
                }
            } while (h > eps);
            return xn;

        }
        public static Vector MDN(Vector[] xn, double eps, Fun2 func)
        {
            double n = xn.Length;
            double m = n + 1;
            Vector xmax = FindMax(xn, func);
            Vector xmin = FindMin(xn, func);
            Vector xc = new Vector(xn[0].Size);
            double cps = 1 / n;
            double kof = 2;
            do
            {
                for (int i = 0; i < n; i++)
                {
                    if (xn[i] == xmax) { }
                    else
                    {
                        xc += xn[i] * (cps - 1);
                    }
                }
                Vector xotr = xmax + kof * (xc - xmax);
                double fotr = func(xotr);
                double fmin = func(FindMin(xn, func));
                if (fotr < fmin)
                {
                    kof = kof * 2;
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] == xmax)
                        {
                            xn[i] = xotr.Copy();
                        }
                    }
                }
                else
                {
                    
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] == xmax)
                        {
                            xn[i] = xotr.Copy();
                        }
                    }
                    //переместить
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] != xmin)
                        {
                            for(int j = 0; j < xn[i].Size; j++)
                            {
                                xn[i][j] = xn[i][j] / 2;
                            }
                        }
                    }
                }
            } while (xmax[0]-xmin[0]<eps);//заглушка

            return FindMin(xn, func);




        }
       /* public static Vector MetodIskluchenia(int n,Fun2 func,Fun ogra)
        {
            Vector answer = new Vector(n);
            double delta = 0.5 * eps;
            double fn = func(xn);
            for (int i = 0; i < n; i++)
            {
                answer[i] = ogra();
            }
            for (int i = 0; i < n; i++)
            {
                Vector xg = xn.Copy();
                xg[i] = xg[i] + delta;
                gr[i] = (func(xg) - fn) / delta;

            }

            double z = func(answer);
        }*/
        private static Vector FindMax(Vector[] xn, Fun2 func)
        {
            Vector max = new Vector(xn[0].size);
            for (int i = 0; i < xn[0].size; i++)
            {
                max[i] = 0;//заглушка
            }
            for (int j = 0; j < xn.Length; j++)
            {
                double z = func(xn[j]);
                double zz = func(max);
                if (func(xn[j]) >= func(max))
                {
                    
                    max = xn[j];
                }
            }
            return max;
        }
        private static Vector FindMin(Vector[] xn, Fun2 func)
        {
            Vector min = new Vector(xn[0].size);
            for (int i = 0; i < xn[0].size; i++)
            {
                min[i] = 1; //заглушка
            }
            for (int j = 0; j < xn.Length; j++)
            {
                double z = func(xn[j]);
                double zz = func(min);
                if (func(xn[j]) <= func(min))
                {
                    min = xn[j];
                }
            }
            return min;
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

    } 
}
