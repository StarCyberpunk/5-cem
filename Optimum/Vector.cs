using System;
using System.Collections.Generic;
using System.Text;

namespace ConsoleTestVectorMatrix
{
    class Vector
    {
        /// Класс вектора
        /// </summary>

        public double[] vector;
        public int size = 0;

        public Vector(int size)
        {
            this.size = size;
            vector = new double[size];
        }

        public Vector(double[] m)
        {
            this.size = m.Length;
            vector = new double[size];
            for (int i = 0; i < size; i++)
                vector[i] = m[i];

        }

        public Vector(Vector m)
        {
            this.size = m.size;
            vector = new double[size];
            for (int i = 0; i < size; i++)
                vector[i] = m[i];

        }


        public double this[int i] // Индексатор
        {
            get
            {
                if (i < 0 && i >= size)
                {
                    Console.WriteLine("Индексы вышли за пределы матрицы ");
                    return 0;
                }
                else
                    return vector[i];
            }
            set
            {
                if (i < 0 && i >= size)
                {
                    Console.WriteLine("Индексы вышли за пределы матрицы ");
                }
                else
                    vector[i] = value;
            }
        }

        public override string ToString()
            => $"{{{string.Join(";", this.vector)}}}";

        public int Size { get { return size; } }
        /// Копирование вектора
        /// </summary>
        public Vector Copy()
        {
            Vector rez = new Vector(vector);
            return rez;
        }

        /// <summary>
        /// Вывод вектора
        /// </summary>
        public void View()
        {
            Console.Write("( ");
            for (int i = 0; i < this.size; i++)
                Console.Write("{0} ", this[i]);
            Console.WriteLine(")");
        }
        public void Clear()
        {
            for (int i = 0; i < this.size; i++) vector[i] = 0.0;
        }
        public double NormaE()
        {
            return Math.Sqrt(this * this);
        }
        /// Умножение вектора на число
        /// </summary>
        public Vector Multiplication(double x)
        {
            Vector rez = new Vector(size);
            for (int i = 0; i < size; i++)
                rez[i] = vector[i] * x;
            return rez;

        }
        /// Сложение векторов
        /// </summary>
        public Vector Addition(Vector a)
        {
            if (size == a.size)
            {
                Vector rez = new Vector(size);
                for (int i = 0; i < size; i++)
                    rez[i] = vector[i] + a[i];
                return rez;
            }

            return null;
        }
        public static Vector operator + (Vector a, Vector b)
        {
            if (a.Size == b.Size)
            {
                Vector c=new Vector(a.Size);
                for (int i = 0; i < a.Size; i++)
                    c[i] += a[i] + b[i];
                return c;
            }
            return null;
        }
        public static Vector operator - (Vector a, Vector b)
        {
            if (a.Size == b.Size)
            {
                Vector c = new Vector(a.Size);
                for (int i = 0; i < a.Size; i++)
                    c[i] += a[i] - b[i];
                return c;
            }
            return null;
        }
        public static Vector operator * (Vector a, double c)
        {          
                Vector r = new Vector(a.Size);
                for (int i = 0; i < a.Size; i++)
                    r[i] = a[i]*c;
                return r;          
        }
        public static Vector operator *(double c,Vector a)
        {
            Vector r = new Vector(a.Size);
            for (int i = 0; i < a.Size; i++)
                r[i] = a[i] * c;
            return r;
        }
        public static double operator * (Vector a, Vector b)
            {
            if(a.Size==b.Size)
            {
                double s = 0.0;
                for (int i = 0; i < a.Size; i++)
                    s += a[i] * b[i];
                return s;
            }
            return Double.NaN;
            }

        /// <summary>
        /// Вычитание векторов
        /// </summary>
        public Vector Subtraction(Vector a)
        {
            if (size == a.size)
            {
                Vector rez = new Vector(size);
                for (int i = 0; i < size; i++)
                    rez[i] = vector[i] - a[i];
                return rez;
            }

            return null;
        }

        /// <summary>
        /// Длина вектора
        /// </summary>
        /// <returns></returns>
        public double Len()
        {

            double x = 0;
            for (int i = 0; i < size; i++)
                x += Math.Pow(vector[i], 2);
            x = Math.Sqrt(x);

            return x;
        }

        /// <summary>
        /// Нормализация вектора
        /// </summary>
        public Vector Normalization()
        {
            Vector rez = new Vector(vector);
            double x = Len();
            if (x > 0)
            {
                for (int i = 0; i < size; i++)
                    rez[i] = rez[i] / x;
                return rez;
            }
            return null;
        }
        private Vector ForEach(Action2 action, Vector vector)
        {
            for (int index = 0; index < this.Size; ++index)
                action(ref this.vector[index], ref vector.vector[index]);

            return this;
        }

        public double ScalarProduct(Vector vector)
        {
            if (vector.Size != this.Size)
                throw new InvalidOperationException("Size of both vectors should be equals");

            var result = 0.0;
            this.ForEach((ref double f, ref double l) => result += f * l, vector);
            return result;
        }
        private delegate void Action2(ref double v, ref double l);
        private delegate void Action(ref double v);
    }
    }
