using System;
using System.Collections.Generic;
using System.Text;

namespace Optimum
{
    class Matrix
    {
        protected int rows, columns;
        protected double[,] data;
        public Matrix(int r, int c)
        {
            this.rows = r; this.columns = c;
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++) data[i, j] = 0;
        }
        public Matrix(double[,] mm)
        {
            this.rows = mm.GetLength(0); this.columns = mm.GetLength(1);
            data = new double[rows, columns];
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    data[i, j] = mm[i, j];
        }
        public int GetCountRows()
        {
            return rows;
        }
        public int GetCountColumns()
        {
            return columns;
        }
        public double this[int i, int j]
        {
            get
            {
                if (i < 0 || j < 0 || i >= rows || j >= columns)
                {
                    Console.WriteLine(" Индексы вышли за пределы матрицы ");
                    return Double.NaN;
                }
                else
                    return data[i, j];
            }
            set
            {
                if (i < 0 && j < 0 && i >= rows && j >= columns)
                {
                    //Console.WriteLine(" Индексы вышли за пределы матрицы ");
                }
                else
                    data[i, j] = value;
            }
        }
        public Vector GetRow(int r)
        {
            if (r >= 0 && r < rows)
            {
                Vector row = new Vector(columns);
                for (int j = 0; j < columns; j++) row[j] = data[r, j];
                return row;
            }
            return null;
        }
        public Vector GetColumn(int c)
        {
            if (c >= 0 && c < columns)
            {
                Vector column = new Vector(rows);
                for (int i = 0; i < rows; i++) column[i] = data[i, c];
                return column;
            }
            return null;
        }
        public void SetRow(int r, Vector rr) 
        {
            if (r >= 0 && r < rows)
            {
                if (columns == rr.Size)
                    for (int j = 0; j < columns; j++) data[r, j] = rr[j];
            }
        }
        public void SetColumn(int c, Vector cc)  
        {
            if (c >= 0 && c < columns)
            {
                if (rows == cc.Size)
                    for (int i = 0; i < rows; i++) data[i, c] = cc[i];
            }
        }
        public Matrix Copy()
        {
            Matrix c = new Matrix(columns, rows);
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < columns; j++)
                    c.data[i, j] = data[i, j];
            return c;

        }
        //нормы
        public double Norma1()
        {
            double s = 0;
            for (int i = 0; i < rows; i++)
                for (int j = 0; j < rows; j++)
                    s += data[i, j] * data[i, j];
            return Math.Sqrt(s);
        }
        public double Norma2()
        {
            double max = 0, s = 0;
            for (int i = 0; i < rows; i++)
            {
                s = 0;
                for (int j = 0; j < rows; j++)
                    s += Math.Abs(data[i, j]);
                if (s > max) max = s;
            }
            return max;
        }
        //умножение матрицы на вектор
        public static Vector operator *(Matrix a, Vector b)
        {
            if (a.columns != b.Size) return null;
            Vector r = new Vector(a.rows);
            for (int i = 0; i < a.rows; i++)
            {
                r[i] = a.GetRow(i) * b;
            }
            return r;
        }

        //печать
        public void Print()
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    Console.Write($"{data[i, j]} \t");             //t - горизонтальная табуляция
                }
                Console.WriteLine();
            }
            Console.WriteLine("\n");
        }
        public override string ToString()
        {
            string s = "{";
            for (int i = 0; i < rows - 1; i++)
            {
                Vector v = this.GetRow(i);
                s += v.ToString() + ";";
            }
            Vector vc = this.GetRow(rows - 1);
            s += vc.ToString() + "}";
            return s;
        }
        //Транспонированние
        public Matrix Transpose()
        {
            Matrix transposeMatrix = new Matrix(rows, columns);
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    transposeMatrix.data[i, j] = data[j, i];
                }
            }
            return transposeMatrix;
        }
        //Умножение матрицы на чило
        public static Matrix MultByNum(Matrix m, double n)   //Умножаем матрицу на число
        {
            Matrix result = new Matrix(m.rows, m.columns);
            for (int i = 0; i < m.rows; i++)
            {
                for (int j = 0; j < m.columns; j++)
                {
                    result[i, j] = m[i, j] * n;
                }
            }
            return result;
        }
        public static Matrix operator *(Matrix m, double n)
        {
            return Matrix.MultByNum(m, n);
        }

        public static Matrix operator *(double n, Matrix m)
        {
            return Matrix.MultByNum(m, n);
        }
        //Умножение матриц
        public static Matrix operator *(Matrix m1, Matrix m2)
        {
            if (m1.columns != m2.rows)
            {
                throw new Exception("Количество столбцов первой матрицы не равно количеству строк второй");
            }
            Matrix result = new Matrix(m1.rows, m2.columns);
            for (int i = 0; i < m1.rows; i++)
            {
                for (int j = 0; j < m2.columns; j++)
                {
                    for (int k = 0; k < m1.columns; k++)
                    {
                        result[i, j] += m1[i, k] * m2[k, j];
                    }
                }
            }
            return result;
        }
        //Сложение матриц
        public static Matrix operator +(Matrix m1, Matrix m2)     //Сложение матриц
        {
            if (m1.rows != m2.rows || m1.columns != m2.columns)
            {
                throw new Exception("Матрицы не совпадают по размерности");
            }
            Matrix result = new Matrix(m1.rows, m1.columns);

            for (int i = 0; i < m1.rows; i++)
            {
                for (int j = 0; j < m2.columns; j++)
                {
                    result[i, j] = m1[i, j] + m2[i, j];
                }
            }
            return result;
        }
        //Вычитание матриц
        public static Matrix operator -(Matrix m1, Matrix m2)     //Вычитание матриц
        {
            if (m1.rows != m2.rows || m1.columns != m2.columns)
            {
                throw new Exception("Матрицы не совпадают по размерности");
            }
            Matrix result = new Matrix(m1.rows, m2.columns);

            for (int i = 0; i < m1.rows; i++)
            {
                for (int j = 0; j < m2.columns; j++)
                {
                    result[i, j] = m1[i, j] - m2[i, j];
                }
            }
            return result;
        }
        //Нижняя треугольная матрица
        public static Vector Lower(Matrix a, Vector b)
        {
            if (a.rows != b.Size)
                throw new Exception("Матрицы и вектор не совпадают по размерности");
            for (int i = 0; i < a.rows; i++)
                if (a[i, i] == 0)
                    throw new Exception("Матрицы null diaganal");

            Vector x = new Vector(b.Size);
            x[0] = b[0] / a[0, 0];
            for (int i = 1; i < a.rows; i++)
            {
                double c = 0;
                for (int j = 0; j < i; j++) { c += x[j] * a[i, j]; }
                x[i] = (b[i] - c) / (a[i, i]);
            }
            return x;

        }

        //Верхняя треугольная матрица
        public static Vector Upper(Matrix a, Vector b)
        {
            if (a.rows != b.Size)
                throw new Exception("Матрицы и вектор не совпадают по размерности");
            for (int i = 0; i < a.rows; i++)
                if (a[i, i] == 0)
                    throw new Exception("Матрицы null diaganal");

            Vector x = new Vector(b.Size);
            x[b.Size - 1] = b[b.Size - 1] / a[a.rows - 1, a.columns - 1];
            for (int i = a.rows - 2; i >= 0; i--)
            {
                double c = 0;
                for (int j = b.Size - 1; j > i; j--) { c += a[i, j] * x[j]; }
                x[i] = (b[i] - c) / (a[i, i]);
            }
            return x;

        }


        //Метод Гаусса
        public static Vector Gauss(Matrix aa, Vector bb)
        {
            Matrix a = aa.Copy();
            Vector b = bb.Copy();
            Vector x = new Vector(b.Size);

            if (a.columns != b.Size) return null;
            for (int j = 0; j < a.rows - 1; j++)
            {
                double maxc = 0;
                int jmax = j;
                for (int i = j; i < a.rows; i++)
                {
                    if (Math.Abs(a[i, j]) > maxc) { maxc = Math.Abs(a[i, j]); jmax = i; }
                }

                if (maxc == 0) return null;//вырождена 
                //change rows
                Vector ri = a.GetRow(j);
                Vector rjmax = a.GetRow(jmax);
                for (int k = 0; k < a.columns; k++)
                {
                    a[j, k] = rjmax[k]; a[jmax, k] = ri[k];
                }
                double c = b[j]; b[j] = b[jmax]; b[jmax] = c;
                for (int i = j + 1; i < a.rows; i++)
                {
                    double cc = a[i, j] / a[j, j];
                    for (int k = 0; k < a.columns; k++)
                    {
                        a[i, k] = a[i, k] - cc * a[j, k];
                    }
                    b[i] = b[i] - cc * b[j];
                }
            }
            x = Upper(a, b);

            //  xx = SLU_DOWN(a, b);

            return x;
        }
        //Гивенс
        public static Vector Givens(Matrix aa, Vector bb)
        {
            Matrix a = aa.Copy();
            Vector b = bb.Copy();
            Vector x = new Vector(b.Size);
            if (a.columns != b.Size) return null;
            if (a.columns != b.Size) return null;
            double z = 0;
            double cos = 0;
            double sin = 0;
            double[,] edmat = new double[a.rows, a.columns];
            for (int k = 0; k < a.rows; k++)
                for (int l = 0; l < a.columns; l++)
                    if (k == l)
                        edmat[k, l] = 1;
                    else
                        edmat[k, l] = 0;
            Matrix aaa = new Matrix(edmat);
            for (int i = 0; i < a.rows - 1; i++)
                for (int j = i + 1; j < a.rows; j++)
                {
                    Matrix ed = new Matrix(edmat);
                    /* if (Math.Abs(a[i, i]) > Math.Abs(a[j, i]))
                        z = a[i, i];
                     else
                        z = a[j, i];
                     if (z == 0)
                     {
                        cos = 1;
                        sin = 0;
                     }
                     if (z == a[i,i])
                     {
                         cos = 1 / Math.Sqrt(1 + a[j, i] * a[j, i] / (a[i, i] * a[i, i]);
                         sin = -cos * a[j, i] / a[i, i];
                     }
                     else
                     {
                         sin = 
                     }*/
                    cos = a[i, i] / Math.Sqrt(a[i, i] * a[i, i] + a[j, i] * a[j, i]);
                    sin = -a[j, i] / Math.Sqrt(a[i, i] * a[i, i] + a[j, i] * a[j, i]);

                    ed[i, i] = cos;
                    ed[i, j] = sin;
                    ed[j, i] = -sin;
                    ed[j, j] = cos;
                    aaa = ed.Transpose() * aaa;
                    a = ed.Transpose() * a;
                    b = ed.Transpose() * b;
                }
            Console.WriteLine(aaa);
            Console.WriteLine(a);
            x = Upper(a, b);
            return x;
        }

        //Метод последовательного приближения
        public static Vector SuccessiveApproximation(Matrix aa, Vector bb, double eps, int kmax)
        {

            Matrix a = aa.Copy();

            Vector b = bb.Copy();
            Vector x = bb.Copy();
            Vector xs = bb.Copy();
            Vector r = bb.Copy();
            if (a.columns != b.Size || a.rows != a.columns) return null;

            for (int i = 0; i < a.rows; i++)
            {
                if (a[i, i] != 0)
                {
                    for (int j = 0; j < a.rows; j++)
                    {
                        if (i != j)
                            a[i, j] = -a[i, j] / a[i, i];

                    }
                    b[i] = b[i] / a[i, i];
                }
                a[i, i] = 0;
            }

            x = b;
            int k = 0;
            do
            {
                k++;
                xs = a * x + b;
                r = xs - x;
                x = xs.Copy();
                if (k > kmax) return null;
            }
            while (r.NormaE() > eps);

            return x;
        }
        //Метод прогонки
        public static Vector Tridiagonal(Vector d, Vector m, Vector u, Vector xx)
        {
            Vector a = d.Copy();
            Vector b = m.Copy();
            Vector c = u.Copy();
            Vector x = xx.Copy();
            Vector otv = new Vector(b.size);
            Vector alpha = new Vector(b.size);
            Vector beta = new Vector(b.size);
            Vector y = new Vector(b.size);

            if (a.size + 1 != b.size || b.size != c.size + 1 || a.size != c.size) return null;

            y[0] = b[0];
            alpha[0] = -c[0] / y[0];
            beta[0] = x[0] / y[0];

            for (int i = 1; i < b.size - 1; i++)
            {
               
                y[i] = b[i] + ((a[i-1]) * alpha[i - 1]);
                alpha[i] = -c[i] / y[i];
                beta[i] = (x[i] - a[i-1] * beta[i - 1]) / y[i];
            }

            y[b.size-1] = b[b.size-1] + a[b.size-2] * alpha[b.size - 2];
            beta[b.size-1] = (x[b.size-1] - a[b.size-2] * beta[b.size - 2]) / y[b.size-1];

            otv[b.size - 1] = beta[b.size - 1];

            for (int i = b.size - 2; i > -1; i--)
            {
                otv[i] = alpha[i] * otv[i + 1] + beta[i];
            }

            return otv;
        }
         
    }
}
