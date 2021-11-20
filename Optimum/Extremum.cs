using System;
using System.Collections.Generic;
using System.Text;

namespace Optimum
{
    delegate double Fun(double x);
    delegate double Fun2(Vector x);
    delegate Vector Fun3(Vector x);
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
            Vector fxp = new Vector(3*n);
            double k = 0;

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
                Vector xc = new Vector(xn[0].Size);
                double kof = 1;
                for (int i = 0; i < n; i++)
                {
                    fmas[i] = func(xn[i]);
                }
                fmas = ShellSort(fmas);
                for (int i = 0; i < n; i++)
                {
                    if (func(xn[i]) == fmas[n - 1]) { }
                    else
                    {
                        xc += xn[i] * cps;
                    }
                }
                 xotr = xc + kof * (xc - FindXN(xn, fmas[n - 1], func));
                double fotr = func(xotr);

                if (fotr < fmas[0])
                {
                   double newfotr = fotr;
                    
                    eeeps =Math.Abs( func(FindXN(xn, fmas[n - 1], func)) - func(xotr));
                    if (fotr < fmas[1])
                    {
                        kof = kof * 2;
                        xotr = xc + kof * (xotr-xc );
                    }
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] == FindXN(xn, fmas[n - 1], func))
                        {
                            xn[i] = xotr.Copy();
                        }
                    }
                }
                else if (fotr > fmas[0])
                {
                    kof = kof / 2;
                    xotr = xc + kof * ( FindXN(xn, fmas[n - 1], func)-xc);
                    eeeps = Math.Abs(-func(FindXN(xn, fmas[n - 1], func)) + func(xotr));
                    for (int i = 0; i < n; i++)
                    {
                        if (xn[i] == FindXN(xn, fmas[n - 1], func))
                        {
                            xn[i] = xotr.Copy();
                        }
                    }

                }
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

            } while (eps< eeeps);// https://habr.com/ru/post/332092/
            Console.WriteLine(k);
            return FindXN(xn, fmas[0], func);




        }
        static Random rnd = new Random();
       public static Vector ozuMethod(Vector x,double eps,double h,Vector fn,Vector fv,Vector tipf,Fun3 func)
        {
            int k = 0; 
            int n = x.Size;
            int m = 3 * n;
            int pmin = int.MinValue;
            Matrix xpTemp = new Matrix(m, n);
            Vector xp = new Vector(m);
            Vector fp = new Vector(m);
            double ft = fNorm(xp, fn, fv, tipf, func);
            double fpmin = double.MaxValue;
            double temp,len;
            do
            {
                k++;
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        temp = rnd.NextDouble() - 0.5;
                        xpTemp[i, j] = temp;
                    }
                    len = 0;
                    for (int j = 0; j < n; j++)
                        len += (xpTemp[i, j] - x[j] * (xpTemp[i, j] - x[j]));
                    len = Math.Sqrt(len);
                    for (int j = 0; j < n; j++)
                        xpTemp[i, j] = x[j] + h * (xpTemp[i, j] - x[j]) / len;
                }
                for (int i = 0; i < m; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        xp[j] = xpTemp[i, j];
                    }
                    fp[i] = fNorm(xp,fn,fv,tipf,func);
                    if (fp[i] < fpmin)
                    {
                        fpmin = fp[i];
                        pmin = i;
                    }

                }
                if (fpmin < ft)
                {
                    for (int j = 0; j < n; j++)
                        x[j] = xpTemp[pmin, j];
                    ft = fpmin;
                    h = h * 1.2;

                }
                else h = h / 2.0;
            } while (h > eps);
            return x;
        }
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
                if (tipf[i] == 1)
                {
                    if (fn[i] > 0)
                        g[i] = 2 - fx[i] / fn[i];
                    if (fn[i] == 0)
                        g[i] = -fx[i] / fn[i] + 1;
                    if (fn[i] < 0)
                        g[i] = fx[i] / fn[i];
                }
                if (tipf[i] == 2)
                {

                    if (fn[i] > 0)
                        g[i] = fx[i] / fv[i];
                    if (fn[i] == 0)
                        g[i] = fx[i] / fv[i] + 1;
                    if (fn[i] < 0)
                        g[i] =2- fx[i] / fv[i];
                }
                if (tipf[i] == 3)
                {

                    g1 = (fx[i] - fn[i]) / (fv[i] - fn[i]);
                    g2 = (fv[i] - fx[i]) / (fv[i] - fn[i]);
                    if (g1 > g2)
                        g[i] = g1;
                    else
                        g[i] = g2;
                }
                if (g[i] > g[maxg])
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

    } 
}
