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
            double fs = fn + fun(xs);
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
            Console.WriteLine("k={0}", k);
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
                    /*Vector t = new Vector(1);
                    t[0] = xn[i];
                    Vector delta2 = new Vector(1);
                    delta2[0] = delta+xn[i];
                    gr[i+1] = (func(delta2)-func(t))/delta;*/
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
            double h = eps;
            Vector hopt = new Vector(xn);
            do
            {
                gr = new Vector(n);
                for (int i = 0; i < n; i++)
                {
                    Vector xg = xn.Copy();
                    xg[i] = xg[i] + delta;
                    gr[i] = (func(xg) - fn) / delta;
                    hopt[i] = h;
                }
                //Условие + границы
                while (fs > fn) {
                    h = GoldSechenie(h, 100, eps, he => func(xn-he * gr));
                }
                dx = -h * gr;
                xs = xn - h * gr;
                fs = func(xs);
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
            double h = eps;
            Vector hopt = new Vector(xn);
            Vector grold=new Vector(n);
            double kof;

            for (int i = 0; i < n; i++)
            {
                Vector xg = xn.Copy();
                xg[i] = xg[i] + delta;
                gr[i] = (func(xg) - fn) / delta;
                hopt[i] = h;
            }
            //Условие + границы
            while (fs > fn)
            {
                h = GoldSechenie(h, 100, eps, he => func(xn - he * gr));
            }
            dx = -h * gr;
            xs = xn - h * gr;
            fs = func(xs);
            grold = gr;
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
                    hopt[i] = h;
                }
                //Условие + границы
                while (fs > fn)
                {
                    h = GoldSechenie(h, 100, eps, he => func(xn - he * gr));
                }
                kof = (gr.NormaE() * gr.NormaE()) / (grold.NormaE() * grold.NormaE());
                gr = gr + kof * grold;
                dx = -h * gr;
                xs = xn - h * gr;
                fs = func(xs);
                grold = gr;
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
            do
            {
                
                Vector rnvec = EdRand(n, rnd);
                for(int i = 0; i < ni; i++)
                {
                    xs[i] = xn[i] - h * rnvec[0];
                }
                
                double fn = func(xn);
                double fs = func(xs);
                if (fs < fn)
                {
                    xn = xs;
                }
                else
                {
                    h = h / 2;
                }
            } while (h > eps);
            return xn;

        }
        public static Vector MDN(Vector xn,double eps,Fun2 func)
        {
           double n = xn.size;
           double m = xn.size - 1;
            double fn = func(xn);
            double sum = 0;
            double xs =0;
            for(int i=0; i < n; i++)
            {
                sum += xn[i];
            }
            xs = 1 / n * sum;


        }
        private static Vector FindMax(Vector xn,Fun2 func)
        {
            
            for
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

    } 
}
